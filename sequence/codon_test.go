package sequence

import (
	"strings"
	"testing"
)

func TestCodonSequence_Properties(t *testing.T) {
	id := "a"
	title := "test"
	seq := "ATGGCGTGG"
	prot := "MAW"
	codons := []string{"ATG", "GCG", "TGG"}
	actual := CodonSequence{CharSequence{id, title, seq}, prot, codons}

	if id != actual.ID() {
		t.Errorf("ID: expected %#v, actual %#v", id, actual.ID())
	}
	if title != actual.Title() {
		t.Errorf("Title: expected %#v, actual %#v", title, actual.Title())
	}
	if seq != actual.Sequence() {
		t.Errorf("Sequence: expected %#v, actual %#v", seq, actual.Sequence())
	}
	if prot != actual.Prot() {
		t.Errorf("Prot: expected %#v, actual %#v", prot, actual.Prot())
	}
	if strings.Join(codons, "") != strings.Join(actual.Codons(), "") {
		t.Errorf("Codons: expected %#v, actual %#v", codons, actual.Codons())
	}
}

func TestCodonSequence_ToUpper(t *testing.T) {
	actual := CodonSequence{CharSequence{"a", "test", "atggcgtgg"}, "maw", []string{"atg", "gcg", "tgg"}}
	actual.ToUpper()
	if exp := "ATGGCGTGG"; exp != actual.seq {
		t.Errorf("ToUpper: expected %#v, actual %#v", exp, actual.seq)
	}
	if exp := "MAW"; exp != actual.prot {
		t.Errorf("ToUpper: expected %#v, actual %#v", exp, actual.prot)
	}
	if exp := "ATGGCGTGG"; exp != strings.Join(actual.codons, "") {
		t.Errorf("ToUpper: expected %#v, actual %#v", exp, actual.codons)
	}
}
func TestCodonSequence_ToLower(t *testing.T) {
	actual := CodonSequence{CharSequence{"a", "test", "ATGGCGTGG"}, "MAW", []string{"ATG", "GCG", "TGG"}}
	actual.ToLower()
	if exp := "atggcgtgg"; exp != actual.seq {
		t.Errorf("ToUpper: expected %#v, actual %#v", exp, actual.seq)
	}
	if exp := "maw"; exp != actual.prot {
		t.Errorf("ToUpper: expected %#v, actual %#v", exp, actual.prot)
	}
	if exp := "atggcgtgg"; exp != strings.Join(actual.codons, "") {
		t.Errorf("ToUpper: expected %#v, actual %#v", exp, actual.codons)
	}
}

func TestCodonSequence_Char(t *testing.T) {
	seq := CharSequence{id: "a", title: "test", seq: "GTGGCGTAG"}
	exp := "A"
	actual := seq.Char(7)
	if exp != actual {
		t.Errorf("Char: expected %#v, actual %#v", exp, actual)
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
