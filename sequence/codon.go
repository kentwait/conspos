package sequence

import (
	"fmt"
	"sort"
	"strings"
)

// CodonSequence is a struct for specifically designed for triplet nucleotide
// codon sequences. It embeds the CharSequence struct which also gives it
// id, title and seq fields. Additionally, CodonSequence has a prot field which
// stores a string and a codon string field which stores a slice of strings.
// The seq, prot and codons fields follow a positional correspondence.
// The first item in the codons slice translates to the first character
// in the prot string. The first item in the codons slice is equal to
// the first three characters of the seq string. This codon-seq correspondence
// should be consistent across the entire sequence.
type CodonSequence struct {
	CharSequence
	prot   string
	codons []string // TODO: Change to *string to avoid duplication
}

// NewCodonSequence is a constructor that creates a new CodonSequence where
// prot and codons field values are automatically computed from the provided
// nucleotide sequence.
func NewCodonSequence(id, title, seq string) *CodonSequence {
	if len(seq)%3 != 0 {
		panic(fmt.Sprintf("Given seq's length (%d) not divisible by 3", len(seq)))
	}
	s := new(CodonSequence)
	s.id = id
	s.title = title
	s.SetSequence(seq)
	return s
}

// ID returns the id field of CodonSequence.
func (s *CodonSequence) ID() string {
	return s.id
}

// Title returns the title field of CodonSequence.
func (s *CodonSequence) Title() string {
	return s.title
}

// Sequence returns the seq field of CodonSequence. The seq field contains
// a nucleotide sequence stored as a string.
func (s *CodonSequence) Sequence() string {
	return s.seq
}

// Codons returns the codon field of CodonSequence. The codon field
// contains a nucleotide sequence delimited by codon. This is stored
// as a slice of 3-character strings.
func (s *CodonSequence) Codons() []string {
	return s.codons
}

// Prot returns the prot field of CodonSequence. The prot field
// contains the translated amino acid sequence based on the seq
// field using the standard genetic code. The amino acid sequence
// is encoded as single-character amino acids and stored as a
// string.
func (s *CodonSequence) Prot() string {
	return s.prot
}

// Char returns a single nucleotide from the seq field of CodonSequence.
func (s *CodonSequence) Char(i int) string {
	return string([]rune(s.seq)[i])
}

// ProtChar returns a single amino acid from the prot field of CodonSequence.
func (s *CodonSequence) ProtChar(i int) string {
	return string(s.prot[i])
}

// Codon returns a single codon 3 nucleotides long from the codons field of
// CodonSequence.
func (s *CodonSequence) Codon(i int) string {
	return string(s.codons[i])
}

/* The following two methods are setters for sequence fields in CodonSequence.
   Note that there is not method to set a protein sequence in the prot field.
   Because of the relationships between seq, prot, and codons, it is impossible
   to compute the values of seq and codons from the protein sequence alone.
   Although a protein sequence can be set literally, this is not recommended as
   there is no way to ensure that the relationships between seq, prot, and
   codons are maintained.
*/

// SetSequence assigns a nucleotide sequence to the seq field of CodonSequence.
// It also automatically fills the codons and prot fields by splitting the
// nucleotide sequence into triplets and translating each codon into its
// corresponding amino acid using the standard genetic code respectively.
func (s *CodonSequence) SetSequence(seq string) {
	// Converts sequence to rune slice to deal with unicode chars
	seqRune := []rune(seq)
	if len(seqRune)%3 != 0 {
		panic(fmt.Sprintf("Length of given seq \"%s\" is not divisible by 3", seq))
	}
	// Overwrite value of .seq
	s.seq = seq
	// Overwrites value of .codons
	var codons []string
	for i := 0; i < len(seqRune); i += 3 {
		codons = append(codons, string(seqRune[i:i+3]))
	}
	s.codons = codons
	// Overwrites the value of .prot
	s.prot = Translate(seq).String()
}

// SetCodons assigns a nucleotide sequence delimited by codon to the codons
// field of CodonSequence. It also automatically fills the seq and prot
// fields by joining the codons into a single continuous string and
// translating each codon into its corresponding amino acid using the
// standard genetic code respectively.
func (s *CodonSequence) SetCodons(seq []string) {
	// Overwrites value of .codons
	s.codons = seq
	// Overwrite value of .seq
	s.seq = strings.Join(seq, "")
	// Overwrites the value of .prot
	s.prot = Translate(s.seq).String()
}

// UngappedCoords returns the positions in the sequence where the character
// does not match the gap character.
func (s *CodonSequence) UngappedCoords(gapChar string) (colCoords []int) {
	// Counts length of the rune slice instead of byte length of the string
	if len([]rune(gapChar))%3 != 0 {
		panic(fmt.Sprintf("Length of given gapChar \"%s\" is not equal to 3", gapChar))
	}
	// Range over the slice of codons
	// If the codon does not match the gapChar string, then it is not a gap
	// Adds its position to the set map.
	set := make(map[int]struct{})
	for j, codon := range s.codons {
		if codon != gapChar {
			set[j] = struct{}{}
		}
	}
	// Range of set. Since this is a map, order is scrambled.
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
func (s *CodonSequence) UngappedPositionSlice(gapChar string) (arr []int) {
	// Counts length of the rune slice instead of byte length of the string
	if len([]rune(gapChar))%3 != 0 {
		panic(fmt.Sprintf("Length of given gapChar \"%s\" is not equal to 3", gapChar))
	}
	cnt := 0
	// Range over the slice of codons
	// If the codon does not match the gapChar string, then it is not a gap
	// Adds the current ungapped count to the array.
	for _, codon := range s.codons {
		if codon != gapChar {
			arr = append(arr, cnt)
			cnt++
			// If it is a gap, adds -1 to the array instead
		} else {
			arr = append(arr, -1)
		}
	}
	return
}

// ToUpper changes the case of the sequence to all uppercase letters.
func (s *CodonSequence) ToUpper() {
	s.seq = strings.ToUpper(s.seq)
	s.prot = strings.ToUpper(s.prot)
	for i := range s.codons {
		s.codons[i] = strings.ToUpper(s.codons[i])
	}
}

// ToLower changes the case of the sequence to all lowercase letters.
func (s *CodonSequence) ToLower() {
	s.seq = strings.ToLower(s.seq)
	s.prot = strings.ToLower(s.prot)
	for i := range s.codons {
		s.codons[i] = strings.ToLower(s.codons[i])
	}
}
