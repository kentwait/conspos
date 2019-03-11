package sequence

import "bytes"

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
