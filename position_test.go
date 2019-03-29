package conspos

import (
	"testing"
)

func TestScorePositions(t *testing.T) {
	m1 := [][]int{
		[]int{0, 1, -1, 2, -1, 3, -1, 4},
		[]int{0, -1, 1, 2, -1, 3, -1, 4},
		[]int{0, 1, -1, 2, -1, 3, 4, -1},
	}
	m2 := [][]int{
		[]int{0, 1, -1, 2, -1, 3, -1, 4},
		[]int{0, -1, 1, 2, -1, 3, 4, -1},
		[]int{0, 1, -1, 2, -1, 3, 4, -1},
	}
	m3 := [][]int{
		[]int{0, 1, -1, 2, -1, 3, -1, 4},
		[]int{0, -1, 1, 2, -1, 3, 4, -1},
		[]int{0, -1, 1, 2, -1, 3, -1, 4},
	}
	exp := []uint{
		3, 2, 2, 3, 3, 3, 1, 1,
	}
	res := ScorePositions(m1, m2, m3)
	// Test if scores are equal
	for i := range res {
		if res[i] != exp[i] {
			t.Errorf("Expected and actual scores are not equal: %v != %v", res[i], exp[i])
		}
	}
}

func TestMarkConsistent(t *testing.T) {
	c := "C"
	i := "I"
	n := "N"
	scoreSlice := []uint{
		3, 2, 2, 3, 3, 3, 1, 1,
	}
	exp := []string{
		c, i, i, c, c, c, n, n,
	}
	res := MarkConsistent(scoreSlice, 3, c, i, n)
	// Test if scores are equal
	for i := range res {
		if res[i] != exp[i] {
			t.Errorf("Expected and actual characters are not the same: %v != %v", res[i], exp[i])
		}
	}
}

func TestMarkOnlyConsistent(t *testing.T) {
	c := "C"
	n := "N"
	scoreSlice := []uint{
		3, 2, 2, 3, 3, 3, 1, 1,
	}
	exp := []string{
		c, n, n, c, c, c, n, n,
	}
	res := MarkOnlyConsistent(scoreSlice, 3, c, n)
	// Test if scores are equal
	for i := range res {
		if res[i] != exp[i] {
			t.Errorf("Expected and actual characters are not the same: %v != %v", res[i], exp[i])
		}
	}
}

func TestMarkConsistentCodon(t *testing.T) {
	c := "C"
	i := "I"
	n := "N"
	scoreSlice := []uint{
		3, 2, 2, 3, 3, 3, 1, 1,
	}
	exp := []string{
		c, c, c, i, i, i, i, i, i, c, c, c, c, c, c, c, c, c, n, n, n, n, n, n,
	}
	res := MarkConsistentCodon(scoreSlice, 3, c, i, n)
	// Test if scores are equal
	for i := range res {
		if res[i] != exp[i] {
			t.Errorf("Expected and actual characters are not the same: %v != %v", res[i], exp[i])
		}
	}
}

func TestMarkOnlyConsistentCodon(t *testing.T) {
	c := "C"
	n := "N"
	scoreSlice := []uint{
		3, 2, 2, 3, 3, 3, 1, 1,
	}
	exp := []string{
		c, c, c, n, n, n, n, n, n, c, c, c, c, c, c, c, c, c, n, n, n, n, n, n,
	}
	res := MarkOnlyConsistentCodon(scoreSlice, 3, c, n)
	// Test if scores are equal
	for i := range res {
		if res[i] != exp[i] {
			t.Errorf("Expected and actual characters are not the same: %v != %v", res[i], exp[i])
		}
	}
}
