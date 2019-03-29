package conspos

import (
	"bytes"
	"fmt"
	"os"
	"strconv"
	"strings"

	fa "github.com/kentwait/gofasta"
)

// ConsistentAlignmentPositions returns the list of positions in the alignment that are considered consistent given by the alignment pattern per site across all given alignments.
func ConsistentAlignmentPositions(gapChar string, matrices ...[][]int) []bool {
	// Transpose matrices and combine as string
	patternSetMap := make(map[string]int)
	var templatePattern []string
	var patternBuffer bytes.Buffer
	/*
		Loops over all ungapped position matrices passed to the function.
		Assumes that all matrices have the same shape.
		For two matrices given, each matrix will look like this:
		a := [[0, 1, 2, 3,-1, 4, 5, 6,-1, 7,...]
			  [0, 1, 2, 3, 4,-1, 5, 6,-1, 7,...]
			  [0, 1, 2, 3,-1,-1, 4, 5, 6, 7,...]]
		b := [[0, 1, 2, 3,-1, 4, 5, 6,-1, 7,...]
			  [0, 1, 2, 3,-1, 4, 5, 6,-1, 7,...]
			  [0, 1, 2, 3,-1,-1, 4, 5, 6, 7,...]]

		This loop iterates over each matrix, and transposes the matrix such that:
		at := [[ 0, 0, 0],
			   [ 1, 1, 1],
			   [ 2, 2, 2],
			   [ 3, 3, 3],
			   [-1, 4,-1],
			   [ 5, 5, 4],
			   [ 6, 6, 5],
			   [-1,-1, 6],
			   [ 7, 7, 7],
				   ...   ]
		bt := [[ 0, 0, 0],
			   [ 1, 1, 1],
			   [ 2, 2, 2],
			   [ 3, 3, 3],
			   [-1,-1,-1],
			   [ 4, 4,-1],
			   [ 5, 5, 4],
			   [ 6, 6, 5],
			   [-1,-1, 6],
			   [ 7, 7, 7],
				   ...   ]

		The goal is to compare whether the pattern in the rows are consistent across all matrices.
		Comparing rows means comparing slices. In Golang, slices cannot be immediately compared, and must use a loop to do an itemwise comparison.
		As a workaround, each row is encoded as a string whose values are separated by a comma "," like this:
		at := ["0,0,0",
			   "1,1,1",
			   "2,2,2",
			   "3,3,3",
			   "-1,4,-1",
			   "5,5,4",
			   "6,6,5",
			   "-1,-1,6",
			   "7,7,7",
				   ...   ]
		bt := ["0,0,0",
			   "1,1,1",
			   "2,2,2",
			   "3,3,3",
			   "-1,-1,-1",
			   "4,4,-1",
			   "5,5,4",
			   "6,6,5",
			   "-1,-1,6",
			   "7,7,7",
				   ...   ]

		Through this method, each row becomes a string and multiple comparisons becomes straightforward.

		"0,0,0" == "0,0,0" = true
		"1,1,1" == "1,1,1" = true
		"1,1,1" == "1,1,1" = true
		"2,2,2" == "2,2,2" = true
		"3,3,3" == "3,3,3" = true
		"-1,4,-1" == "-1,-1,-1" = false
		...

	*/
	for k, matrix := range matrices {
		for j := 0; j < len(matrix[0]); j++ {
			// Gets the pattern for each column in the matrix
			for i := 0; i < len(matrix); i++ {
				patternBuffer.WriteString(strconv.Itoa(matrix[i][j]) + ",")
			}

			// If this is the first matrix, store its pattern in a slice.
			// This pattern will be used as the template later.
			if k == 0 {
				templatePattern = append(templatePattern, patternBuffer.String())
			}
			// Adds the pattern to a mapping.
			// Increments to record the number of times the pattern has been observed.
			patternSetMap[patternBuffer.String()]++
			// Clears the buffer for the next column
			patternBuffer.Reset()
		}
	}

	pos := make([]bool, len(matrices[0][0]))
	for j, pattern := range templatePattern {
		// For the current column, get how many times the template pattern was observed.
		// If the number of times is equal to the number of matrices, then this pattern was observed over all the matrices and is consistent.
		// However, if the number of times is less than the number of matrices, then some matrices showed a different pattern.
		// This deficit makes the current column inconsistent.
		if patternSetMap[pattern] == len(matrices) {
			pos[j] = true
		} else {
			pos[j] = false
		}
	}
	return pos
}

// ConsistentCodonAlignmentPositions returns the list of positions in the codon alignment that are considered consistent given by the alignment pattern per site across all given alignments.
func ConsistentCodonAlignmentPositions(gapChar string, matrices ...[][]int) []bool {
	// ConsistentAlignmentPositions constructs a boolean slice for the codon matrices
	var codonPos []bool
	for _, pos := range ConsistentAlignmentPositions(gapChar, matrices...) {
		// Added 3 times to because each codon has 3 nucleotide sites
		codonPos = append(codonPos, pos, pos, pos)
	}
	return codonPos
}

// ConsistentAlnPipeline aligns using global, local, and affine-local alignment strategies to determine positions that have a consistent alignment pattern over the three different strategies.
func ConsistentAlnPipeline(inputPath, gapChar, markerID, consistentMarker, inconsistentMarker string, iterations int, toUpper, toLower, saveTempAlns bool) bytes.Buffer {
	// TODO: Allow this to be a parameter instead of being hard-coded
	const mafftCmd = "mafft"

	os.Stderr.WriteString(fmt.Sprintf("%s: ", inputPath))

	/* Align using 3 different strategies implemented in MAFFT
	   - global alignment (GinsiAlign)
	   - local alignment (LinsiAlign)
	   - affine-gap local alignment (EinsiAlign)

	   These calls run sequentially with MAFFT saturating all cores.
	   MAFFT outputs results to stdout and these functions capture stdout to return a string.
	*/
	// TODO: Propagate ExecMafft error into *Align
	// TODO: *Align should probably output a buffer instead of a string
	ginsiString := CharAlign(mafftCmd, inputPath, "ginsi", iterations)
	linsiString := CharAlign(mafftCmd, inputPath, "linsi", iterations)
	einsiString := CharAlign(mafftCmd, inputPath, "einsi", iterations)

	// Check if alignments are not empty.
	// If empty, print error message to stderr and exit with code 1
	if len(ginsiString) == 0 {
		EmptyAlnError("G-INSI", inputPath)
	}
	if len(linsiString) == 0 {
		EmptyAlnError("L-INSI", inputPath)
	}
	if len(einsiString) == 0 {
		EmptyAlnError("E-INSI", inputPath)
	}

	// Each result (string) is parsed to create character alignments
	ginsiAln := fa.FastaToAlignment(strings.NewReader(ginsiString), false)
	linsiAln := fa.FastaToAlignment(strings.NewReader(linsiString), false)
	einsiAln := fa.FastaToAlignment(strings.NewReader(einsiString), false)
	os.Stderr.WriteString(".")

	// Writes temp alignments if necessary
	// TODO: ToFasta is not necessary, just write the ginsiString, linsiString, einsiString
	if saveTempAlns == true {
		einsiAln.ToFastaFile(inputPath + ".einsi.aln")
		ginsiAln.ToFastaFile(inputPath + ".ginsi.aln")
		linsiAln.ToFastaFile(inputPath + ".linsi.aln")
	}

	// consistentPos is a boolean slice indicating per position whether it is consistent or not.
	consistentPos := ConsistentAlignmentPositions(
		gapChar,
		einsiAln.UngappedPositionMatrix(gapChar),
		ginsiAln.UngappedPositionMatrix(gapChar),
		linsiAln.UngappedPositionMatrix(gapChar),
	)
	os.Stderr.WriteString(".")

	if toUpper == true {
		einsiAln.ToUpper()
	} else if toLower == true {
		einsiAln.ToLower()
	}

	os.Stderr.WriteString(" Done.\n")

	// TODO: Add aiblity to select what alignment is outputted
	return MarkedAlignmentToBuffer(einsiAln, consistentPos, markerID, consistentMarker, inconsistentMarker)
}

// ConsistentCodonAlnPipeline aligns codon sequences using global, local, and  affine-local alignment strategies to determine positions that have a consistent alignment pattern over the three different strategies.
func ConsistentCodonAlnPipeline(inputPath, gapChar, markerID, consistentMarker, inconsistentMarker string, iterations int, toUpper, toLower, saveTempAlns bool) bytes.Buffer {
	// TODO: Allow this to be a parameter instead of being hard-coded
	const mafftCmd = "mafft"

	os.Stderr.WriteString(fmt.Sprintf("%s: ", inputPath))

	// Create an Alignment of CodonSequence to generate translated protein sequence from nucleotide sequence
	c := fa.FastaFileToCodonAlignment(inputPath)

	// Read protein sequences from Alignment of CodonSequences and create a Fasta string in buffer
	protReader := strings.NewReader(c.ToFasta())

	// Pass buff to each of the three alignment strategies.
	// These will align protein sequences in MAFFT.
	// Based on the protein alignment, the original codon alignment is adjusted using the AlignCodonsUsingProtAlignment function.
	ginsiString := CodonAlignStdin(mafftCmd, protReader, "ginsi", iterations, c)
	linsiString := CodonAlignStdin(mafftCmd, protReader, "linsi", iterations, c)
	einsiString := CodonAlignStdin(mafftCmd, protReader, "einsi", iterations, c)

	// Check if string alignment is not empty
	if len(ginsiString) == 0 {
		EmptyAlnError("G-INSI", inputPath)
	}
	if len(linsiString) == 0 {
		EmptyAlnError("L-INSI", inputPath)
	}
	if len(einsiString) == 0 {
		EmptyAlnError("E-INSI", inputPath)
	}

	// The FASTA outputs are parsed to create codon alignments.
	ginsiAln := fa.FastaToAlignment(strings.NewReader(ginsiString), true)
	linsiAln := fa.FastaToAlignment(strings.NewReader(linsiString), true)
	einsiAln := fa.FastaToAlignment(strings.NewReader(einsiString), true)
	os.Stderr.WriteString(".")

	// TODO: ToFasta conversion is unnecessary.
	// *insiString is already in FASTA format
	if saveTempAlns == true {
		einsiAln.ToFastaFile(inputPath + ".einsi.aln")
		ginsiAln.ToFastaFile(inputPath + ".ginsi.aln")
		linsiAln.ToFastaFile(inputPath + ".linsi.aln")
	}

	// consistentPos is a boolean slice indicating per position whether it is consistent or not. Length of consistentPos is the length of the codon alignment as single characters.
	consistentPos := ConsistentCodonAlignmentPositions(
		gapChar,
		einsiAln.UngappedPositionMatrix(gapChar),
		ginsiAln.UngappedPositionMatrix(gapChar),
		linsiAln.UngappedPositionMatrix(gapChar),
	)

	if toUpper == true {
		einsiAln.ToUpper()
	} else if toLower == true {
		einsiAln.ToLower()
	}

	os.Stderr.WriteString(" Done.\n")

	return MarkedAlignmentToBuffer(einsiAln, consistentPos, markerID, consistentMarker, inconsistentMarker)
}
