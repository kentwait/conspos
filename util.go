package conspos

import (
	"bytes"
	"fmt"
	"os"

	fa "github.com/kentwait/gofasta"
)

// Exists returns whether the given file or directory Exists or not,
// and accompanying errors.
func Exists(path string) (bool, error) {
	_, err := os.Stat(path)
	if err == nil {
		return true, nil
	}
	if os.IsNotExist(err) {
		return false, nil
	}
	return true, err
}

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

// OffsetAlignCodons aligns codons by adding gaps to offset positions based on its aligned translated protein sequence.
func OffsetAlignCodons(c []*fa.CodonSequence, p []*fa.CharSequence) (sequences []*fa.CodonSequence) {
	gapRune := []rune("-")[0]
	for i := range p {
		newSeq := bytes.Buffer{}
		ucnt := 0
		for _, char := range p[i].Sequence() {
			seq := c[i].Sequence()
			cStart := ucnt * 3
			cEnd := (ucnt + 1) * 3
			if char == gapRune {
				newSeq.WriteString("---")
			} else {
				newSeq.WriteString(seq[cStart:cEnd])
				ucnt++
			}
		}
		sequence := fa.NewCodonSequence(c[i].ID(), c[i].Description(), newSeq.String())
		newSeq.Reset()
		sequences = append(sequences, sequence)
	}
	return
}
