package main

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"os"
	"reflect"
	"strconv"
	"strings"
)

/* Single character
 */

// StringToCharSequences loads a string generated from properly formatted
// FASTA file into a SequenceAlignment struct.
func StringToCharSequences(s string) (sequences SequenceAlignment) {
	lines := strings.Split(s, "\n")

	var id, title string
	var seqBuffer bytes.Buffer
	var splitted []string

	for _, line := range lines {
		if strings.HasPrefix(line, ">") {
			if seqBuffer.Len() > 0 {
				sequences = append(sequences, &CharSequence{id, title, seqBuffer.String()})
				seqBuffer.Reset()
			}
			splitted = strings.SplitN(line[1:], " ", 2)
			id = splitted[0]
			if len(splitted) == 2 {
				title = splitted[1]
			}
		} else if strings.HasPrefix(line, "\n") {
			continue
		} else if strings.HasPrefix(line, "#") {
			continue
		} else {
			seqBuffer.WriteString(line)
		}
	}
	if seqBuffer.Len() > 0 {
		sequences = append(sequences, &CharSequence{id, title, seqBuffer.String()})
	}
	return
}

/* Protein utilities
 */

// ProtSequencesToString converts sequences in the sequene alignment into a FASTA
// formatted string.
func ProtSequencesToString(a SequenceAlignment) string {
	buffer := SequencesToBuffer(a)
	return buffer.String()
}

// ProtSequencesToBuffer converts sequences in the sequene alignment into a buffered
// stream which can then be converted to bytes or a string.
func ProtSequencesToBuffer(a SequenceAlignment) bytes.Buffer {
	var buffer bytes.Buffer
	// Append each Sequence in SequenceAlignment
	for _, s := range a {
		if len(s.Title()) > 0 {
			buffer.WriteString(fmt.Sprintf(">%s %s\n", s.ID(), s.Title()))
		} else {
			buffer.WriteString(fmt.Sprintf(">%s\n", s.ID()))
		}
		v := reflect.ValueOf(s).Elem().FieldByName("prot").String()
		buffer.WriteString(v + "\n")
	}
	return buffer
}

/* Codon utilities
 */

// StringToCodonSequences loads a string generated from properly formatted
// FASTA file into a SequenceAlignment struct.
func StringToCodonSequences(s string) (sequences SequenceAlignment) {
	lines := strings.Split(s, "\n")

	var id, title string
	var seqBuffer bytes.Buffer
	var splitted []string

	for _, line := range lines {
		if strings.HasPrefix(line, ">") {
			if seqBuffer.Len() > 0 {
				sequences = append(sequences, NewCodonSequence(id, title, seqBuffer.String()))
				seqBuffer.Reset()
			}
			splitted = strings.SplitN(line[1:], " ", 2)
			id = splitted[0]
			if len(splitted) == 2 {
				title = splitted[1]
			}
		} else if strings.HasPrefix(line, "\n") {
			continue
		} else if strings.HasPrefix(line, "#") {
			continue
		} else {
			seqBuffer.WriteString(line)
		}
	}
	if seqBuffer.Len() > 0 {
		sequences = append(sequences, NewCodonSequence(id, title, seqBuffer.String()))
	}
	return
}

// Translate naively converts nucleotides into amino acids without regard
// for the reading frame.
func Translate(s string) *bytes.Buffer {
	var buff bytes.Buffer
	var trans string
	for i := 0; i < len(s); i += 3 {
		trans = geneticCode[string(s[i:i+3])]
		if trans == "" {
			buff.WriteString("X")
		} else {
			buff.WriteString(trans)
		}
	}
	return &buff
}

// FastaToCodonSequence reads a FASTA file and converts it into a sequence
// alignment of codon sequences.
func FastaToCodonSequence(path string) SequenceAlignment {
	b, err := ioutil.ReadFile(path)
	if err != nil {
		panic(err)
	}
	return StringToCodonSequences(string(b))
}

// AlignCodonsToBuffer takes an unaligned set of codon sequences and an
// aligned set of protein sequences to create a new sequence alignment
// by using the aligned protein sequences as a guide. Outputs a buffer
// containing a the new codon alignment as a FASTA-formatted string.
func AlignCodonsToBuffer(c, p SequenceAlignment) bytes.Buffer {
	gapRune := []rune("-")[0]
	var newFasta bytes.Buffer
	for i := range p {
		newSeq := bytes.Buffer{}

		if len(c[i].Title()) > 0 {
			newFasta.WriteString(fmt.Sprintf(">%s %s\n", c[i].ID(), p[i].Title()))
		} else {
			newFasta.WriteString(fmt.Sprintf(">%s\n", c[i].ID()))
		}

		ucnt := 0
		for _, char := range p[i].Sequence() {
			seq := reflect.ValueOf(c[i]).Elem().FieldByName("seq").String()
			cStart := ucnt * 3
			cEnd := (ucnt + 1) * 3
			if char == gapRune {
				newSeq.WriteString("---")
			} else {
				newSeq.WriteString(seq[cStart:cEnd])
				ucnt++
			}
		}
		newSeq.WriteString("\n")
		newFasta.Write(newSeq.Bytes())

		newSeq.Reset()
	}
	return newFasta
}

// AlignCodonsToString takes an unaligned set of codon sequences and an
// aligned set of protein sequences to create a new sequence alignment
// by using the aligned protein sequences as a guide. Outputs a FASTA-
// formatted string.
func AlignCodonsToString(c, p SequenceAlignment) string {
	b := AlignCodonsToBuffer(c, p)
	return b.String()
}

/* IO utils
 */

// SequencesToBuffer converts sequences in the sequene alignment into a buffered
// stream which can then be converted to bytes or a string.
func SequencesToBuffer(a SequenceAlignment) bytes.Buffer {
	var buffer bytes.Buffer
	// Append each Sequence in SequenceAlignment
	for _, s := range a {
		if len(s.Title()) > 0 {
			buffer.WriteString(fmt.Sprintf(">%s %s\n", s.ID(), s.Title()))
		} else {
			buffer.WriteString(fmt.Sprintf(">%s\n", s.ID()))
		}
		buffer.WriteString(s.Sequence() + "\n")
	}
	return buffer
}

// SequencesToString converts sequences in the sequene alignment into a FASTA
// formatted string.
func SequencesToString(a SequenceAlignment) string {
	buffer := SequencesToBuffer(a)
	return buffer.String()
}

// BufferedMarkedAlignment writes a marked multiple sequence alignment
// in the FASTA format to the buffer.
func BufferedMarkedAlignment(template SequenceAlignment, consistentPos []bool, markerID, consistentMarker, inconsistentMarker string) bytes.Buffer {
	var buffer bytes.Buffer

	// Append marker sequence
	buffer.WriteString(fmt.Sprintf(">%s\n", markerID))
	for _, t := range consistentPos {
		if t == true {
			buffer.WriteString(consistentMarker)
		} else {
			buffer.WriteString(inconsistentMarker)
		}
	}
	buffer.WriteString("\n")

	// Append each Sequence in SequenceAlignment
	for _, s := range template {
		if len(s.Title()) > 0 {
			buffer.WriteString(fmt.Sprintf(">%s %s\n", s.ID(), s.Title()))
		} else {
			buffer.WriteString(fmt.Sprintf(">%s\n", s.ID()))
		}
		buffer.WriteString(s.Sequence() + "\n")
	}
	return buffer
}

// WriteBufferToFile writes the contents of a buffer to file.
func WriteBufferToFile(path string, b bytes.Buffer) {
	f, err := os.Create(path)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	b.WriteTo(f)
	f.Sync()
}

/* Pipeline
 */

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

	ginsiAln := StringToCharSequences(ginsiString)
	linsiAln := StringToCharSequences(linsiString)
	einsiAln := StringToCharSequences(einsiString)
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

	return BufferedMarkedAlignment(einsiAln, consistentPos, markerID, consistentMarker, inconsistentMarker)
}

// ConsistentCodonAlnPipeline aligns codon sequences using global, local, and
// affine-local alignment strategies to determine positions that have a
// consistent alignment pattern over the three different strategies.
func ConsistentCodonAlnPipeline(inputPath, gapChar, markerID, consistentMarker, inconsistentMarker string, iterations int, toUpper, toLower, saveTempAlns bool) bytes.Buffer {

	const mafftCmd = "mafft"

	os.Stderr.WriteString(fmt.Sprintf("%s: ", inputPath))

	// Create CodonSequences to generate translated protein sequence from nucleotide sequence
	c := FastaToCodonSequence(inputPath)

	// Read protein sequences from SequenceAlignment of CodonSequences and create a Fasta string in buffer
	buff := ProtSequencesToBuffer(c)

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
	ginsiAln := StringToCodonSequences(ginsiString)
	linsiAln := StringToCodonSequences(linsiString)
	einsiAln := StringToCodonSequences(einsiString)
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

	if toUpper == true {
		einsiAln.ToUpper()
	} else if toLower == true {
		einsiAln.ToLower()
	}

	os.Stderr.WriteString(" Done.\n")

	return BufferedMarkedAlignment(einsiAln, consistentPos, markerID, consistentMarker, inconsistentMarker)
}
