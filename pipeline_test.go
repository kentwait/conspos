package conspos

import (
	"strings"
	"testing"

	fa "github.com/kentwait/gofasta"
)

func TestScoringPipeline(t *testing.T) {
	aln1 := []*fa.CharSequence{
		fa.NewCharSequence("test1", "desc1", "ATGAA-TTTT"),
		fa.NewCharSequence("test2", "", "ATGAATT-TT"),
		fa.NewCharSequence("test3", "desc3", "ATGAA-TTTT"),
	}
	aln2 := []*fa.CharSequence{
		fa.NewCharSequence("test1", "desc1", "ATGAA-TTTT"),
		fa.NewCharSequence("test2", "", "ATGAAT-TTT"),
		fa.NewCharSequence("test3", "desc3", "ATGAA-TTTT"),
	}
	aln3 := []*fa.CharSequence{
		fa.NewCharSequence("test1", "desc1", "ATGAAT-TTT"),
		fa.NewCharSequence("test2", "", "ATGAAT-TTT"),
		fa.NewCharSequence("test3", "desc3", "ATGAAT-TTT"),
	}
	c := "C"
	i := "I"
	n := "N"
	exp := []string{
		c, c, c, c, c, i, n, n, c, c,
	}
	res := ScoringPipeline("-", c, i, n, aln1, aln2, aln3)
	if exp, res := strings.Join(exp, ""), strings.Join(res, ""); exp != res {
		t.Errorf("expected and actual score marks are not the same: %s != %s", exp, res)
	}
}

func TestScoringPipeline_Binary(t *testing.T) {
	aln1 := []*fa.CharSequence{
		fa.NewCharSequence("test1", "desc1", "ATGAA-TTTT"),
		fa.NewCharSequence("test2", "", "ATGAATT-TT"),
		fa.NewCharSequence("test3", "desc3", "ATGAA-TTTT"),
	}
	aln2 := []*fa.CharSequence{
		fa.NewCharSequence("test1", "desc1", "ATGAA-TTTT"),
		fa.NewCharSequence("test2", "", "ATGAAT-TTT"),
		fa.NewCharSequence("test3", "desc3", "ATGAA-TTTT"),
	}
	aln3 := []*fa.CharSequence{
		fa.NewCharSequence("test1", "desc1", "ATGAAT-TTT"),
		fa.NewCharSequence("test2", "", "ATGAAT-TTT"),
		fa.NewCharSequence("test3", "desc3", "ATGAAT-TTT"),
	}
	c := "C"
	i := ""
	n := "N"
	exp := []string{
		c, c, c, c, c, n, n, n, c, c,
	}
	res := ScoringPipeline("-", c, i, n, aln1, aln2, aln3)
	if exp, res := strings.Join(exp, ""), strings.Join(res, ""); exp != res {
		t.Errorf("expected and actual score marks are not the same: %s != %s", exp, res)
	}
}
