package alignment

import (
	"bytes"
	"fmt"
)

// ToString converts sequences in the sequene alignment into a FASTA
// formatted string.
func ToString(a Alignment) string {
	buffer := ToBuffer(a)
	return buffer.String()
}

// ToBuffer converts sequences in the sequene alignment into a buffered
// stream which can then be converted to bytes or a string.
func ToBuffer(a Alignment) bytes.Buffer {
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
