package sequence

// Sequence is an interface for single character sequences stored as a string
// and multi-character sequences stored as a slice.
type Sequence interface {
	Meta
	Getter
	Setter
}

// Meta defines methods to retrieve sequence metadata.
type Meta interface {
	ID() string
	Title() string
}

// Getter contains methods to retrieve information about sequence data.
type Getter interface {
	Sequence() string
	Char(int) string
	UngappedCoords(string) []int
	UngappedPositionSlice(string) []int
}

// Setter contains methods to set/modify sequence data.
type Setter interface {
	SetSequence(string)
	ToUpper()
	ToLower()
}
