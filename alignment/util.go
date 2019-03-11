package alignment

import (
	"bytes"
	"fmt"
)

// SequenceAlignmentToString converts sequences in the sequene alignment into a FASTA
// formatted string.
func SequenceAlignmentToString(a Alignment) string {
	buffer := SequenceAlignmentToBuffer(a)
	return buffer.String()
}

// SequenceAlignmentToBuffer converts sequences in the sequene alignment into a buffered
// stream which can then be converted to bytes or a string.
func SequenceAlignmentToBuffer(a Alignment) bytes.Buffer {
	var buffer bytes.Buffer
	// Append each Sequence in Alignment
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
