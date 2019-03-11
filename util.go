package main

import (
	"bytes"
	"fmt"

	aln "github.com/kentwait/conspos/alignment"
)

// MarkedAlignmentToBuffer writes a marked multiple sequence alignment
// in the FASTA format to the buffer.
func MarkedAlignmentToBuffer(template aln.Alignment, consistentPos []bool, markerID, consistentMarker, inconsistentMarker string) bytes.Buffer {
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
		if len(s.Title()) > 0 {
			buffer.WriteString(fmt.Sprintf(">%s %s\n", s.ID(), s.Title()))
		} else {
			buffer.WriteString(fmt.Sprintf(">%s\n", s.ID()))
		}
		buffer.WriteString(s.Sequence() + "\n")
	}
	return buffer
}
