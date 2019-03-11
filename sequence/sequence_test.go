package sequence

import (
	"testing"
)

func TestCharSequence_UngappedCoords(t *testing.T) {
	seq := "TTT---TTCTTATTG"
	s := CharSequence{"test", "", seq}
	exp := []int{0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14}

	res := s.UngappedCoords("-")

	for i, expValue := range exp {
		if expValue != res[i] {
			t.Errorf("UngappedCoords(\"-\"): expected (%d) %d, actual %d",
				i, expValue, res[i],
			)
		}
	}
}

func TestCharSequence_UngappedPositionSlice(t *testing.T) {
	seq := "TTT---TTCTTATTG"
	s := CharSequence{"test", "", seq}
	exp := []int{0, 1, 2, -1, -1, -1, 3, 4, 5, 6, 7, 8, 9, 10, 11}

	res := s.UngappedPositionSlice("-")

	for i, expValue := range exp {
		if expValue != res[i] {
			t.Errorf("UngappedCoords(\"-\"): expected (%d) %d, actual %d",
				i, expValue, res[i],
			)
		}
	}
}

func TestCodonSequence_SetSequence_seq(t *testing.T) {
	s := CodonSequence{CharSequence{"test", "", ""}, "", []string{}}
	seq := "TTT---TTCTTATTG"
	s.SetSequence(seq)

	if s.seq != seq {
		t.Errorf("SetSequence(\"%s\"): expected %s, actual %s", seq, seq, s.seq)
	}
}

func TestCodonSequence_SetSequence_prot(t *testing.T) {
	s := CodonSequence{CharSequence{"test", "", ""}, "", []string{}}
	seq := "TTTTTCTTATTGTCTTCCTCATCGTATTACTAATAGTGTTGCTGATGGCTTCTCCTACTGCCTCCCCCACCGCATCACCAACAGCGTCGCCGACGGATTATCATAATGACTACCACAACGAATAACAAAAAGAGTAGCAGAAGGGTTGTCGTAGTGGCTGCCGCAGCGGATGACGAAGAGGGTGGCGGAGGG---NNN"
	s.SetSequence(seq)
	exp := "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG-X"

	if s.prot != exp {
		t.Errorf("SetSequence(\"%s\"): expected %s, actual %s", seq, exp, s.prot)
	}
}

func TestCodonSequence_SetSequence_codon(t *testing.T) {
	s := CodonSequence{CharSequence{"test", "", ""}, "", []string{}}
	seq := "TTTTTCTTATTGTCTTCCTCATCGTATTACTAATAGTGTTGCTGATGGCTTCTCCTACTGCCTCCCCCACCGCATCACCAACAGCGTCGCCGACGGATTATCATAATGACTACCACAACGAATAACAAAAAGAGTAGCAGAAGGGTTGTCGTAGTGGCTGCCGCAGCGGATGACGAAGAGGGTGGCGGAGGG---NNN"
	s.SetSequence(seq)
	exp := []string{
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
		"---", "NNN",
	}

	for i, expValue := range exp {
		if s.codons[i] != expValue {
			t.Errorf("SetSequence(\"%s\"): expected codon (%d) %s, actual %s", seq, i, expValue, s.codons[i])
		}
	}
}

func TestCodonSequence_SetCodons_seq(t *testing.T) {
	s := CodonSequence{CharSequence{"test", "", ""}, "", []string{}}
	codons := []string{
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
		"---", "NNN",
	}
	exp := "TTTTTCTTATTGTCTTCCTCATCGTATTACTAATAGTGTTGCTGATGGCTTCTCCTACTGCCTCCCCCACCGCATCACCAACAGCGTCGCCGACGGATTATCATAATGACTACCACAACGAATAACAAAAAGAGTAGCAGAAGGGTTGTCGTAGTGGCTGCCGCAGCGGATGACGAAGAGGGTGGCGGAGGG---NNN"
	s.SetCodons(codons)

	if s.seq != exp {
		t.Errorf("SetCodons(\"%v\"): expected %s, actual %s", codons, exp, s.seq)
	}
}

func TestCodonSequence_SetCodons_prot(t *testing.T) {
	s := CodonSequence{CharSequence{"test", "", ""}, "", []string{}}
	codons := []string{
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
		"---", "NNN",
	}
	s.SetCodons(codons)
	exp := "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG-X"

	if s.prot != exp {
		t.Errorf("SetCodons(\"%v\"): expected %s, actual %s", codons, exp, s.prot)
	}
}

func TestCodonSequence_SetCodons_codon(t *testing.T) {
	s := CodonSequence{CharSequence{"test", "", ""}, "", []string{}}
	codons := []string{
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
		"---", "NNN",
	}
	s.SetCodons(codons)

	for i, expValue := range codons {
		if s.codons[i] != expValue {
			t.Errorf("SetCodons(\"%v\"): expected codon (%d) %s, actual %s", codons, i, expValue, s.codons[i])
		}
	}
}

func TestCodonSequence_UngappedCoords(t *testing.T) {
	seq := "TTT---TTCTTATTG"
	s := NewCodonSequence("test", "", seq)
	gapChar := "---"
	exp := []int{0, 2, 3, 4}

	res := s.UngappedCoords(gapChar)

	for i, expValue := range exp {
		if expValue != res[i] {
			t.Errorf("UngappedCoords(\"%s\"): expected (%d) %d, actual %d",
				gapChar, i, expValue, res[i],
			)
		}
	}
}

func TestCodonSequence_UngappedPositionSlice(t *testing.T) {
	seq := "TTT---TTCTTATTG"
	s := NewCodonSequence("test", "", seq)
	gapChar := "---"
	exp := []int{0, -1, 1, 2, 3}

	res := s.UngappedPositionSlice(gapChar)

	for i, expValue := range exp {
		if expValue != res[i] {
			t.Errorf("UngappedCoords(\"%s\"): expected (%d) %d, actual %d",
				gapChar, i, expValue, res[i],
			)
		}
	}
}
