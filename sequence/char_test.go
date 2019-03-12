package sequence

import "testing"

func TestNewCharSequence(t *testing.T) {
	id := "a"
	title := "test"
	seq := "ATGGCGTAG"
	exp := CharSequence{id: id, title: title, seq: seq}
	actual := *NewCharSequence(id, title, seq)

	if exp.id != actual.id && exp.title != actual.title && exp.seq != actual.seq {
		t.Errorf("NewCharSequence: expected %#v, actual %#v", exp, actual)
	}
}

func TestCharSequence_Properties(t *testing.T) {
	id := "a"
	title := "test"
	seq := "ATGGCGTAG"
	actual := CharSequence{id: id, title: title, seq: seq}

	if id != actual.ID() {
		t.Errorf("ID: expected %#v, actual %#v", id, actual.ID())
	}
	if title != actual.Title() {
		t.Errorf("Title: expected %#v, actual %#v", title, actual.Title())
	}
	if seq != actual.Sequence() {
		t.Errorf("Sequence: expected %#v, actual %#v", seq, actual.Sequence())
	}
}

func TestCharSequence_SetID(t *testing.T) {
	seq := CharSequence{id: "a", title: "test", seq: "ATGGCGTAG"}
	exp := "x"
	seq.SetID(exp)
	if exp != seq.id {
		t.Errorf("SetID: expected %#v, actual %#v", exp, seq.id)
	}
}

func TestCharSequence_SetTitle(t *testing.T) {
	seq := CharSequence{id: "a", title: "test", seq: "ATGGCGTAG"}
	exp := "test again"
	seq.SetTitle(exp)
	if exp != seq.title {
		t.Errorf("SetTitle: expected %#v, actual %#v", exp, seq.title)
	}
}

func TestCharSequence_SetSequence(t *testing.T) {
	seq := CharSequence{id: "a", title: "test", seq: "ATGGCGTAG"}
	exp := "CCCCCCCCC"
	seq.SetSequence(exp)
	if exp != seq.seq {
		t.Errorf("SetSequence: expected %#v, actual %#v", exp, seq.seq)
	}
}

func TestCharSequence_ToUpper(t *testing.T) {
	seq := CharSequence{id: "a", title: "test", seq: "atggcgtag"}
	exp := "ATGGCGTAG"
	seq.ToUpper()
	if exp != seq.seq {
		t.Errorf("ToUpper: expected %#v, actual %#v", exp, seq.seq)
	}
}

func TestCharSequence_ToLower(t *testing.T) {
	seq := CharSequence{id: "a", title: "test", seq: "ATGGCGTAG"}
	exp := "atggcgtag"
	seq.ToLower()
	if exp != seq.seq {
		t.Errorf("ToLower: expected %#v, actual %#v", exp, seq.seq)
	}
}

func TestCharSequence_Char(t *testing.T) {
	seq := CharSequence{id: "a", title: "test", seq: "GTGGCGTAG"}
	exp := "A"
	actual := seq.Char(7)
	if exp != actual {
		t.Errorf("Char: expected %#v, actual %#v", exp, actual)
	}
}

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
