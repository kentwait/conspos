package main

import (
	"io/ioutil"

	aln "github.com/kentwait/conspos/alignment"
)

// FastaToCodonAlignment reads a FASTA file and converts it into a sequence
// alignment of codon sequences.
func FastaToCodonAlignment(path string) aln.Alignment {
	b, err := ioutil.ReadFile(path)
	if err != nil {
		panic(err)
	}
	return StringToCodonAlignment(string(b))
}
