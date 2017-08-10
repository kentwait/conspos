package main

import (
	"strings"
	"testing"
)

func TestTranslate(t *testing.T) {
	var codons = []string{
		"TTT", "TTC", "TTA", "TTG",
		"TCT", "TCC", "TCA", "TCG",
		"TAT", "TAC", "TAA", "TAG",
		"TGT", "TGC", "TGA", "TGG",
		"CTT", "CTC", "CTA", "CTG",
		"CCT", "CCC", "CCA", "CCG",
		"CAT", "CAC", "CAA", "CAG",
		"CGT", "CGC", "CGA", "CGG",
		"ATT", "ATC", "ATA", "ATG",
		"ACT", "ACC", "ACA", "ACG",
		"AAT", "AAC", "AAA", "AAG",
		"AGT", "AGC", "AGA", "AGG",
		"GTT", "GTC", "GTA", "GTG",
		"GCT", "GCC", "GCA", "GCG",
		"GAT", "GAC", "GAA", "GAG",
		"GGT", "GGC", "GGA", "GGG",
		"---", "A-A",
	}
	var aa = []string{
		"F", "F", "L", "L",
		"L", "L", "L", "L",
		"I", "I", "I", "M",
		"V", "V", "V", "V",
		"S", "S", "S", "S",
		"P", "P", "P", "P",
		"T", "T", "T", "T",
		"A", "A", "A", "A",
		"Y", "Y", "*", "*",
		"H", "H", "Q", "Q",
		"N", "N", "K", "K",
		"D", "D", "E", "E",
		"C", "C", "*", "W",
		"R", "R", "R", "R",
		"S", "S", "R", "R",
		"G", "G", "G", "G",
		"-", "X",
	}
	s := strings.Join(codons, "")
	exp := strings.Join(aa, "")
	res := Translate(s).String()

	for i := range exp {
		if res[i] != exp[i] {
			t.Errorf("Translate(%s): expected %s, actual %s", s, string(exp[i]), string(res[i]))
		}
	}
}
