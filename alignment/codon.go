package alignment

import (
	"bytes"
	"fmt"
	"reflect"
	"strings"

	"github.com/kentwait/conspos/sequence"
)

// AlignCodonsUsingProtAlignment takes an unaligned set of codon sequences and an
// aligned set of protein sequences to create a new sequence alignment
// by using the aligned protein sequences as a guide. Outputs a buffer
// containing a the new codon alignment as a FASTA-formatted string.
func AlignCodonsUsingProtAlignment(c, p Alignment) bytes.Buffer {
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

// StringToCodonAlignment loads a string generated from properly formatted
// FASTA file into a Alignment struct.
func StringToCodonAlignment(s string) (sequences Alignment) {
	lines := strings.Split(s, "\n")

	var id, title string
	var seqBuffer bytes.Buffer
	var splitted []string

	for _, line := range lines {
		if strings.HasPrefix(line, ">") {
			if seqBuffer.Len() > 0 {
				sequences = append(sequences, sequence.NewCodonSequence(id, title, seqBuffer.String()))
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
		sequences = append(sequences, sequence.NewCodonSequence(id, title, seqBuffer.String()))
	}
	return
}
