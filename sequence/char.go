package sequence

import (
	"sort"
	"strings"
)

// CharSequence struct for storing single-character biological sequences such
// as nucleotides and single-letter amino acids. However, any sequence that
// whose element can be represented as a single string character can be stored
// in CharSequence.
type CharSequence struct {
	id    string
	title string
	seq   string
}

// ID returns the id field of CharSequence.
func (s *CharSequence) ID() string {
	return s.id
}

// Title returns the title field of CharSequence.
func (s *CharSequence) Title() string {
	return s.title
}

// Sequence returns the seq field of CharSequence.
func (s *CharSequence) Sequence() string {
	return s.seq
}

// Char returns a single character from the seq field of CharSequence.
func (s *CharSequence) Char(i int) string {
	return string([]rune(s.seq)[i])
}

// SetSequence assigns a string to the seq field of CharSequence.
func (s *CharSequence) SetSequence(seq string) {
	s.seq = seq
}

// UngappedCoords returns the positions in the sequence where the character
// does not match the gap character.
func (s *CharSequence) UngappedCoords(gapChar string) (colCoords []int) {
	set := make(map[int]struct{})
	// Assumes gapChar contains only a "single character"
	// Convert single character string to rune slice, taking the first item
	gapRune := []rune(gapChar)[0]
	// Range over rune slice, j counts by Unicode code points, s is the rune representation of the character
	for j, s := range []rune(s.seq) {
		// If sequence rune is not a gap character rune, add to rune position to set, 0-indexed
		if s != gapRune {
			set[j] = struct{}{} // Uses empty anonymous struct
		}
	}
	// Range over set of positions
	// Since this is a map, order is scrambled
	for key := range set {
		colCoords = append(colCoords, key)
	}
	sort.Ints(colCoords)
	return
}

// UngappedPositionSlice returns a slice that counts only over characters
// that does not match the gap character in the sequence.
// If a character matches the gap character, -1 is inserted instead of the
// ungapped count.
func (s *CharSequence) UngappedPositionSlice(gapChar string) (arr []int) {
	// Assumes gapChar contains only a "single character"
	// Convert single character string to rune slice, taking the first item
	gapRune := []rune(gapChar)[0]
	cnt := 0
	for _, s := range []rune(s.seq) {
		// If sequence rune is not a gap character rune, append current count value to array and increment
		if s != gapRune {
			arr = append(arr, cnt)
			cnt++
			// If it is equal to the gap character rune, then append a -1.
			// Do not increment.
		} else {
			arr = append(arr, -1)
		}
	}
	return
}

// ToUpper changes the case of the sequence to all uppercase letters.
func (s *CharSequence) ToUpper() {
	s.seq = strings.ToUpper(s.seq)
}

// ToLower changes the case of the sequence to all lowercase letters.
func (s *CharSequence) ToLower() {
	s.seq = strings.ToLower(s.seq)
}
