package conspos

import (
	"io"
	"os"
	"strings"

	fa "github.com/kentwait/gofasta"
)

// SingleCharMafftPipeline aligns sequences using global, local and affine-local strategies in MAFFT.
func SingleCharMafftPipeline(inputFile io.Reader, mafftCmd string, saveAlns string, iterations, verbosity int) ([]*fa.CharSequence, []*fa.CharSequence, []*fa.CharSequence) {
	/* Align using 3 different strategies implemented in MAFFT
	   - global alignment (GinsiAlign)
	   - local alignment (LinsiAlign)
	   - affine-gap local alignment (EinsiAlign)

	   These calls run sequentially with MAFFT saturating all cores.
	   MAFFT outputs results to stdout and these functions capture stdout to return a string.
	*/
	var err error

	// MAFFT global alignment
	ginsiString, err := CharAlignStdio(mafftCmd, inputFile, "ginsi", iterations, verbosity)
	if err != nil {
		// Write error to stderr and exit
		os.Stderr.WriteString("[ERROR] Could not perform global alignment (ginsi)\n")
		os.Stderr.WriteString(err.Error() + "\n")
		os.Exit(1)
	}
	if len(ginsiString) == 0 {
		// Print error message to stderr and exit with code 1
		os.Stderr.WriteString("[ERROR] Empty global alignment (ginsi)\n")
		os.Exit(1)
	}

	// MAFFT local alignment
	linsiString, err := CharAlignStdio(mafftCmd, inputFile, "linsi", iterations, verbosity)
	if err != nil {
		// Write error to stderr and exit
		os.Stderr.WriteString("[ERROR] Could not perform local alignment (linsi)\n")
		os.Stderr.WriteString(err.Error() + "\n")
		os.Exit(1)
	}
	if len(linsiString) == 0 {
		// Print error message to stderr and exit with code 1
		os.Stderr.WriteString("[ERROR] Empty local alignment (linsi)\n")
		os.Exit(1)
	}

	// MAFFT affine gap alignment
	einsiString, err := CharAlignStdio(mafftCmd, inputFile, "einsi", iterations, verbosity)
	if err != nil {
		// Write error to stderr and exit
		os.Stderr.WriteString("[ERROR] Could not perform affine-gap penalty alignment (einsi)\n")
		os.Stderr.WriteString(err.Error() + "\n")
		os.Exit(1)
	}
	if len(einsiString) == 0 {
		// Print error message to stderr and exit with code 1
		os.Stderr.WriteString("[ERROR] Empty affine-gap penalty alignment (einsi)\n")
		os.Exit(1)
	}

	// Each result (string) is parsed to create character alignments
	ginsiAln := fa.FastaToCharSequences(strings.NewReader(ginsiString))
	linsiAln := fa.FastaToCharSequences(strings.NewReader(linsiString))
	einsiAln := fa.FastaToCharSequences(strings.NewReader(einsiString))

	if len(saveAlns) > 0 {
		// Save alignments to file
	}

	return ginsiAln, linsiAln, einsiAln
}

// CodonMafftPipeline aligns translated protein sequences using global, local and affine-local strategies in MAFFT.
func CodonMafftPipeline(inputFile io.Reader, mafftCmd string, saveProtAlns, saveCodonAlns string, iterations, verbosity int) ([]*fa.CodonSequence, []*fa.CodonSequence, []*fa.CodonSequence) {
	// TODO: Whats happening here?
	// Create an Alignment of CodonSequence to generate translated protein sequence from nucleotide sequence
	codonSequences := fa.FastaToCodonSequences(inputFile)

	// Translate codon sequences into proteins sequences
	// Read protein sequences from Alignment of CodonSequences and create a Fasta string in buffer
	var protSequences fa.Alignment
	for _, s := range codonSequences {
		protSequences = append(protSequences, s.ToProtSequence())
	}
	protReader := strings.NewReader(protSequences.ToFasta())

	// Pass buff to each of the three alignment strategies.
	// These will align protein sequences in MAFFT.
	// Based on the protein alignment, the original codon alignment is adjusted using the AlignCodonsUsingProtAlignment function.
	protGinsiAln, protLinsiAln, protEinsiAln := SingleCharMafftPipeline(protReader, mafftCmd, saveProtAlns, iterations, verbosity)

	// Create codon alignments based on protein alignments
	codonGinsiAln := OffsetAlignCodons(codonSequences, protGinsiAln)
	codonLinsiAln := OffsetAlignCodons(codonSequences, protLinsiAln)
	codonEinsiAln := OffsetAlignCodons(codonSequences, protEinsiAln)

	if len(saveCodonAlns) > 0 {
		// Save codon alignments to file
	}

	return codonGinsiAln, codonLinsiAln, codonEinsiAln
}

// ScoringPipeline compares multiple alignments and marks consistent alignment patterns.
func ScoringPipeline(gapChar, consChar, interChar, inconsChar string, alns ...fa.Alignment) (markSlice []string) {
	var ungappedMatrices [][][]int
	for _, sequences := range alns {
		var aln fa.Alignment
		for _, sequence := range sequences {
			seq := fa.Sequence(sequence)
			aln = append(aln, seq)
		}
		ungappedMatrix := aln.UngappedPositionMatrix(gapChar)
		ungappedMatrices = append(ungappedMatrices, ungappedMatrix)
	}
	scoreSlice := ScorePositions(ungappedMatrices...)

	if len(interChar) > 0 {
		markSlice = MarkConsistent(scoreSlice, len(alns), consChar, interChar, inconsChar)
	} else {
		markSlice = MarkOnlyConsistent(scoreSlice, len(alns), consChar, inconsChar)
	}
	return
}
