package alignment

import (
	"bytes"
	"fmt"
	"reflect"
	"strings"

	"github.com/kentwait/conspos/sequence"
)

// StringToCharAlignment loads a string generated from properly formatted
// FASTA file into a Alignment struct.
func StringToCharAlignment(s string) (sequences Alignment) {
	lines := strings.Split(s, "\n")

	var id, title string
	var seqBuffer bytes.Buffer
	var splitted []string

	for _, line := range lines {
		if strings.HasPrefix(line, ">") {
			if seqBuffer.Len() > 0 {
				sequences = append(sequences, sequence.NewCharSequence(id, title, seqBuffer.String()))
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
		sequences = append(sequences, sequence.NewCharSequence(id, title, seqBuffer.String()))
	}
	return
}

// ProtAlignmentToBuffer converts sequences in the sequene alignment into a buffered
// stream which can then be converted to bytes or a string.
func ProtAlignmentToBuffer(a Alignment) bytes.Buffer {
	var buffer bytes.Buffer
	// Append each Sequence in Alignment
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
