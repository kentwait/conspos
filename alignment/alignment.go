package alignment

import (
	"os"

	"github.com/kentwait/conspos/sequence"
)

// Alignment is a slice of Sequence pointers.
type Alignment []sequence.Sequence

// UngappedCoords returns the row and column positions in the sequence alignment
// where the character does not match the gap character.
func (a Alignment) UngappedCoords(gapChar string) (rowCoords, colCoords []int) {
	var currColCoords []int
	for i, s := range a {
		currColCoords = s.UngappedCoords(gapChar)
		for c := 0; c < len(currColCoords); c++ {
			rowCoords = append(rowCoords, i)
		}
		colCoords = append(colCoords, currColCoords...)
	}
	return
}

// UngappedPositionMatrix returns a matrix that counts only over characters
// that does not match the gap character for each sequence in the alignment.
// If a character in a sequence matches the gap character, -1 is inserted
// instead of the ungapped count.
func (a Alignment) UngappedPositionMatrix(gapChar string) (m [][]int) {
	for _, s := range a {
		m = append(m, s.UngappedPositionSlice(gapChar))
	}
	return
}

// ToUpper changes the case of all sequences to all uppercase letters.
func (a Alignment) ToUpper() {
	for _, s := range a {
		s.ToUpper()
	}
}

// ToLower changes the case of of all sequences to all lowercase letters.
func (a Alignment) ToLower() {
	for _, s := range a {
		s.ToLower()
	}
}

// ToFastaString returns the FASTA-formatted string of the sequence alignment.
func (a Alignment) ToFastaString() string {
	return ToString(a)
}

// ToFasta saves the sequence alignment to a FASTA file.
func (a Alignment) ToFasta(path string) {
	buff := ToBuffer(a)
	f, err := os.Create(path)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	buff.WriteTo(f)
	f.Sync()
}
