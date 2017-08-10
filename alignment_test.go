package main

import (
	"testing"
)

func TestSequenceAlignment_UngappedCoords(t *testing.T) {
	seq1 := "TTT---TTCTTATTG"
	seq2 := "TTT---TTCTTTTTG"
	seq3 := "TTTTTCTTC---TTG"
	a := SequenceAlignment{
		&CharSequence{"test", "", seq1},
		&CharSequence{"test", "", seq2},
		&CharSequence{"test", "", seq3},
	}
	expR := []int{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}
	expC := []int{0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14}

	r, c := a.UngappedCoords("-")

	for i, expValue := range expR {
		if r[i] != expValue {
			t.Errorf("UngappedCoords(\"-\"): expected row value at (%d) %d, actual %d",
				i, expValue, r[i],
			)
		}
	}
	for i, expValue := range expC {
		if c[i] != expValue {
			t.Errorf("UngappedCoords(\"-\"): expected column value at (%d) %d, actual %d",
				i, expValue, c[i],
			)
		}
	}
}

func TestSequenceAlignment_UngappedPositionMatrix(t *testing.T) {
	seq1 := "TTT---TTCTTATTG"
	seq2 := "TTT---TTCTTTTTG"
	seq3 := "TTTTTCTTC---TTG"
	a := SequenceAlignment{
		&CharSequence{"test", "", seq1},
		&CharSequence{"test", "", seq2},
		&CharSequence{"test", "", seq3},
	}
	exp := [][]int{
		[]int{0, 1, 2, -1, -1, -1, 3, 4, 5, 6, 7, 8, 9, 10, 11},
		[]int{0, 1, 2, -1, -1, -1, 3, 4, 5, 6, 7, 8, 9, 10, 11},
		[]int{0, 1, 2, 3, 4, 5, 6, 7, 8, -1, -1, -1, 9, 10, 11},
	}

	m := a.UngappedPositionMatrix("-")

	for i, expRow := range exp {
		for j, expValue := range expRow {
			if m[i][j] != expValue {
				t.Errorf("UngappedPositionMatrix(\"-\"): expected value at (%d,%d) %d, actual %d",
					i, j, expValue, m[i][j],
				)
			}
		}
	}
}

func TestSequenceAlignment_ToFastaString(t *testing.T) {
	seq1 := "TTT---TTCTTATTG"
	seq2 := "TTT---TTCTTTTTG"
	seq3 := "TTTTTCTTC---TTG"
	a := SequenceAlignment{
		&CharSequence{"test", "", seq1},
		&CharSequence{"test", "", seq2},
		&CharSequence{"test", "", seq3},
	}
	exp := ">test\nTTT---TTCTTATTG\n>test\nTTT---TTCTTTTTG\n>test\nTTTTTCTTC---TTG\n"

	res := a.ToFastaString()

	if res != exp {
		t.Errorf("ToFastaString(): expected:\n%s\n actual:\n%s", exp, res)
	}
}
