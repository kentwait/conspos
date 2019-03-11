package main

import (
	"bytes"
	"fmt"
	"os"
	"strconv"

	aln "github.com/kentwait/conspos/alignment"
)

// ConsistentAlignmentPositions returns the list of positions in the alignment
// that are considered consistent given by the alignment pattern per site
// across all given alignments.
func ConsistentAlignmentPositions(gapChar string, matrices ...[][]int) []bool {
	// Transpose matrices and combine as string
	patternSetMap := make(map[string]int)
	var templatePattern []string
	var patternBuffer bytes.Buffer
	for k, matrix := range matrices {
		for j := 0; j < len(matrix[0]); j++ {
			for i := 0; i < len(matrix); i++ {
				patternBuffer.WriteString(strconv.Itoa(matrix[i][j]) + ",")
			}

			if k == 0 {
				templatePattern = append(templatePattern, patternBuffer.String())
			}
			patternSetMap[patternBuffer.String()]++
			patternBuffer.Reset()
		}
	}

	pos := make([]bool, len(matrices[0][0]))
	for j, pattern := range templatePattern {
		// For the current column, compare the column pattern across the
		// matrices. If a difference in corresponding values are detected,
		// then the pattern is deemed inconsistent and the current column
		// is therefore also inconsistent.
		if patternSetMap[pattern] == len(matrices) {
			pos[j] = true
		} else {
			pos[j] = false
		}
	}
	return pos
}

func ConsistentCodonAlignmentPositions(gapChar string, matrices ...[][]int) []bool {
	var codonPos []bool
	for _, protPos := range ConsistentAlignmentPositions(gapChar, matrices...) {
		for i := 0; i < 3; i++ {
			codonPos = append(codonPos, protPos)
		}
	}
	return codonPos
}

// ConsistentAlnPipeline aligns using global, local, and affine-local alignment
// strategies to determine positions that have a consistent alignment pattern over
// the three different strategies.
func ConsistentAlnPipeline(inputPath, gapChar, markerID, consistentMarker, inconsistentMarker string, iterations int, toUpper, toLower, saveTempAlns bool) bytes.Buffer {

	const mafftCmd = "mafft"

	os.Stderr.WriteString(fmt.Sprintf("%s: ", inputPath))

	ginsiString := GinsiAlign(mafftCmd, inputPath, iterations)
	linsiString := LinsiAlign(mafftCmd, inputPath, iterations)
	einsiString := EinsiAlign(mafftCmd, inputPath, iterations)

	// Check if string alignment is not empty
	if len(ginsiString) < 1 {
		EmptyAlnError("G-INSI", inputPath)
		os.Exit(1)
	}
	if len(linsiString) < 1 {
		EmptyAlnError("L-INSI", inputPath)
		os.Exit(1)
	}
	if len(einsiString) < 1 {
		EmptyAlnError("E-INSI", inputPath)
		os.Exit(1)
	}

	ginsiAln := aln.StringToCharAlignment(ginsiString)
	linsiAln := aln.StringToCharAlignment(linsiString)
	einsiAln := aln.StringToCharAlignment(einsiString)
	os.Stderr.WriteString(".")

	if saveTempAlns == true {
		einsiAln.ToFasta(inputPath + ".einsi.aln")
		ginsiAln.ToFasta(inputPath + ".ginsi.aln")
		linsiAln.ToFasta(inputPath + ".linsi.aln")
	}

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

	return MarkedAlignmentToBuffer(einsiAln, consistentPos, markerID, consistentMarker, inconsistentMarker)
}

// ConsistentCodonAlnPipeline aligns codon sequences using global, local, and
// affine-local alignment strategies to determine positions that have a
// consistent alignment pattern over the three different strategies.
func ConsistentCodonAlnPipeline(inputPath, gapChar, markerID, consistentMarker, inconsistentMarker string, iterations int, toUpper, toLower, saveTempAlns bool) bytes.Buffer {

	const mafftCmd = "mafft"

	os.Stderr.WriteString(fmt.Sprintf("%s: ", inputPath))

	// Create CodonSequences to generate translated protein sequence from nucleotide sequence
	c := FastaToCodonAlignment(inputPath)

	// Read protein sequences from Alignment of CodonSequences and create a Fasta string in buffer
	buff := aln.ProtAlignmentToBuffer(c)

	ginsiString := GinsiCodonAlign(mafftCmd, buff, iterations, c)
	linsiString := LinsiCodonAlign(mafftCmd, buff, iterations, c)
	einsiString := EinsiCodonAlign(mafftCmd, buff, iterations, c)

	// Check if string alignment is not empty
	if len(ginsiString) < 1 {
		EmptyAlnError("G-INSI", inputPath)
		os.Exit(1)
	}
	if len(linsiString) < 1 {
		EmptyAlnError("L-INSI", inputPath)
		os.Exit(1)
	}
	if len(einsiString) < 1 {
		EmptyAlnError("E-INSI", inputPath)
		os.Exit(1)
	}

	// Translate FASTA to protein then align
	ginsiAln := aln.StringToCodonAlignment(ginsiString)
	linsiAln := aln.StringToCodonAlignment(linsiString)
	einsiAln := aln.StringToCodonAlignment(einsiString)
	os.Stderr.WriteString(".")

	if saveTempAlns == true {
		einsiAln.ToFasta(inputPath + ".einsi.aln")
		ginsiAln.ToFasta(inputPath + ".ginsi.aln")
		linsiAln.ToFasta(inputPath + ".linsi.aln")
	}

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
