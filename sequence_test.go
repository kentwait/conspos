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
