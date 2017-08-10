package main

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
	seq := "TTTTTCTTATTGTCTTCCTCATCGTATTACTAATAGTGTTGCTGATGGCTTCTCCTACTGCCTCCCCCACCGCATCACCAACAGCGTCGCCGACGGATTATCATAATGACTACCACAACGAATAACAAAAAGAGTAGCAGAAGGGTTGTCGTAGTGGCTGCCGCAGCGGATGACGAAGAGGGTGGCGGAGGG---"
	s.SetSequence(seq)
	exp := "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG-"

	if s.prot != exp {
		t.Errorf("SetSequence(\"%s\"): expected %s, actual %s", seq, exp, s.prot)
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
