package main

import (
	"bytes"
	"fmt"
	"reflect"

	fa "github.com/kentwait/gofasta"
)

// MarkedAlignmentToBuffer writes a marked multiple sequence alignment
// in the FASTA format to the buffer.
func MarkedAlignmentToBuffer(template fa.Alignment, consistentPos []bool, markerID, consistentMarker, inconsistentMarker string) bytes.Buffer {
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

	// Append each Sequence in Alignment
	for _, s := range template {
		if len(s.Description()) > 0 {
			buffer.WriteString(fmt.Sprintf(">%s %s\n", s.ID(), s.Description()))
		} else {
			buffer.WriteString(fmt.Sprintf(">%s\n", s.ID()))
		}
		buffer.WriteString(s.Sequence() + "\n")
	}
	return buffer
}

// AlignCodonsUsingProtAlignment takes an unaligned set of codon sequences and an
// aligned set of protein sequences to create a new sequence alignment
// by using the aligned protein sequences as a guide. Outputs a buffer
// containing a the new codon alignment as a FASTA-formatted string.
func AlignCodonsUsingProtAlignment(c, p fa.Alignment) bytes.Buffer {
	gapRune := []rune("-")[0]
	var newFasta bytes.Buffer
	for i := range p {
		newSeq := bytes.Buffer{}

		if len(c[i].Description()) > 0 {
			newFasta.WriteString(fmt.Sprintf(">%s %s\n", c[i].ID(), p[i].Description()))
		} else {
			newFasta.WriteString(fmt.Sprintf(">%s\n", c[i].ID()))
		}

		ucnt := 0
		for _, char := range p[i].Sequence() {
			seq := reflect.ValueOf(c[i]).Elem().FieldByName("sequence").String()
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
