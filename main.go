package main

import (
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"runtime"
	"strconv"
	"strings"
)

// Translate naively converts nucleotides into amino acids without regard
// for the reading frame.
func Translate(s string) *bytes.Buffer {
	var buff bytes.Buffer
	var trans string
	for i := 0; i < len(s); i += 3 {
		trans = geneticCode[string(s[i:i+3])]
		if trans == "" {
			buff.WriteString("X")
		} else {
			buff.WriteString(trans)
		}
	}
	return &buff
}

// ExecMafft calls the MAFFT program with the given arguments.
func ExecMafft(mafftCmd string, args []string) string {
	absPath, lookErr := exec.LookPath(mafftCmd)
	if lookErr != nil {
		panic(lookErr)
	}

	threads := runtime.NumCPU() - 1
	args = append([]string{"--thread", strconv.Itoa(threads)}, args...)

	cmd := exec.Command(absPath, args...)
	stdout, err := cmd.Output()
	if err != nil {
		fmt.Println(fmt.Sprint(err) + ": " + string(stdout))
		return ""
	}
	return string(stdout)
}

// EinsiAlign calls MAFFT to align sequences by local alignment with
// affine-gap scoring.
func EinsiAlign(mafftCmd, fastaPath string, iterations int) (stdout string) {
	args := []string{
		"--quiet",
		"--genafpair",
		"--maxiterate",
		strconv.Itoa(iterations),
		fastaPath,
	}
	stdout = ExecMafft(mafftCmd, args)
	return
}

// LinsiAlign calls MAFFT to align sequences by local alignment.
func LinsiAlign(mafftCmd, fastaPath string, iterations int) (stdout string) {
	args := []string{
		"--quiet",
		"--localpair",
		"--maxiterate",
		strconv.Itoa(iterations),
		fastaPath,
	}
	stdout = ExecMafft(mafftCmd, args)
	return
}

// GinsiAlign calls MAFFT to align sequences by global alignment.
func GinsiAlign(mafftCmd, fastaPath string, iterations int) (stdout string) {
	args := []string{
		"--quiet",
		"--globalpair",
		"--maxiterate",
		strconv.Itoa(iterations),
		fastaPath,
	}
	stdout = ExecMafft(mafftCmd, args)
	return
}

// ExecMafftStdin calls the MAFFT program with the given arguments and using
// standard input as input.
func ExecMafftStdin(mafftCmd string, buff bytes.Buffer, args []string) string {
	absPath, lookErr := exec.LookPath(mafftCmd)
	if lookErr != nil {
		panic(lookErr)
	}

	args = append(args, "-")
	cmd := exec.Command(absPath, args...)

	os.Stdin.Write(buff.Bytes())
	cmd.Stdin = os.Stdin
	stdout, err := cmd.Output()
	if err != nil {
		fmt.Println(fmt.Sprint(err) + ": " + string(stdout))
		return ""
	}
	return string(stdout)
}

// FastaToCodonSequence reads a FASTA file and converts it into a sequence
// alignment of codon sequences.
func FastaToCodonSequence(path string) SequenceAlignment {
	b, err := ioutil.ReadFile(path)
	if err != nil {
		panic(err)
	}
	return StringToCodonSequences(string(b))
}

// AlignCodonsToBuffer takes an unaligned set of codon sequences and an
// aligned set of protein sequences to create a new sequence alignment
// by using the aligned protein sequences as a guide. Outputs a buffer
// containing a the new codon alignment as a FASTA-formatted string.
func AlignCodonsToBuffer(c, p SequenceAlignment) bytes.Buffer {
	gapRune := []rune("-")[0]
	var newFasta bytes.Buffer
	for i := range p {
		newSeq := bytes.Buffer{}

		if len(c[i].Title()) > 0 {
			newFasta.WriteString(fmt.Sprintf(">%s %s\n", c[i].ID(), p[i].Title()))
		} else {
			newFasta.WriteString(fmt.Sprintf(">%s\n", c[i].ID()))
		}

		ucnt := 0
		for _, char := range p[i].Sequence() {
			seq := reflect.ValueOf(c[i]).Elem().FieldByName("seq").String()
			cStart := ucnt * 3
			cEnd := (ucnt + 1) * 3
			if char == gapRune {
				newSeq.WriteString("---")
			} else {
				newSeq.WriteString(seq[cStart:cEnd])
				ucnt++
			}
		}
		newSeq.WriteString("\n")
		newFasta.Write(newSeq.Bytes())

		newSeq.Reset()
	}
	return newFasta
}

// AlignCodonsToString takes an unaligned set of codon sequences and an
// aligned set of protein sequences to create a new sequence alignment
// by using the aligned protein sequences as a guide. Outputs a FASTA-
// formatted string.
func AlignCodonsToString(c, p SequenceAlignment) string {
	b := AlignCodonsToBuffer(c, p)
	return b.String()
}

// EinsiCodonAlign calls MAFFT to align sequences by local alignment with
// affine-gap scoring.
func EinsiCodonAlign(mafftCmd, fastaPath string, iterations int) string {
	args := []string{
		"--quiet",
		"--genafpair",
		"--maxiterate",
		strconv.Itoa(iterations),
	}
	// Create CodonSequences to generate translated protein sequence from nucleotide sequence
	c := FastaToCodonSequence(fastaPath)

	// Read protein sequences from SequenceAlignment of CodonSequences and create a Fasta string in buffer
	buff := ProtSequencesToBuffer(c)
	// Pass Fasta string as stdin to mafft then capture stdout string
	stdout := ExecMafftStdin(mafftCmd, buff, args)

	// Create CharSequences from protein alignment
	p := StringToCharSequences(stdout)

	// Use protein alignment to offset codons and match alignment. Output as Fasta string
	newStdout := AlignCodonsToString(c, p)

	return newStdout
}

// LinsiCodonAlign calls MAFFT to align sequences by local alignment.
func LinsiCodonAlign(mafftCmd, fastaPath string, iterations int) (stdout string) {
	args := []string{
		"--quiet",
		"--localpair",
		"--maxiterate",
		strconv.Itoa(iterations),
		fastaPath,
	}
	stdout = ExecMafft(mafftCmd, args)
	return
}

// GinsiCodonAlign calls MAFFT to align sequences by global alignment.
func GinsiCodonAlign(mafftCmd, fastaPath string, iterations int) (stdout string) {
	args := []string{
		"--quiet",
		"--globalpair",
		"--maxiterate",
		strconv.Itoa(iterations),
		fastaPath,
	}
	stdout = ExecMafft(mafftCmd, args)
	return
}

// StringToCharSequences loads a string generated from properly formatted
// FASTA file into a SequenceAlignment struct.
func StringToCharSequences(s string) (sequences SequenceAlignment) {
	lines := strings.Split(s, "\n")

	var id, title string
	var seqBuffer bytes.Buffer
	var splitted []string

	for _, line := range lines {
		if strings.HasPrefix(line, ">") {
			if seqBuffer.Len() > 0 {
				sequences = append(sequences, &CharSequence{id, title, seqBuffer.String()})
				seqBuffer.Reset()
			}
			splitted = strings.SplitN(line[1:], " ", 2)
			id = splitted[0]
			if len(splitted) == 2 {
				title = splitted[1]
			}
		} else if strings.HasPrefix(line, "\n") {
			continue
		} else if strings.HasPrefix(line, "#") {
			continue
		} else {
			seqBuffer.WriteString(line)
		}
	}
	if seqBuffer.Len() > 0 {
		sequences = append(sequences, &CharSequence{id, title, seqBuffer.String()})
	}
	return
}

// StringToCodonSequences loads a string generated from properly formatted
// FASTA file into a SequenceAlignment struct.
func StringToCodonSequences(s string) (sequences SequenceAlignment) {
	lines := strings.Split(s, "\n")

	var id, title string
	var seqBuffer bytes.Buffer
	var splitted []string

	for _, line := range lines {
		if strings.HasPrefix(line, ">") {
			if seqBuffer.Len() > 0 {
				sequences = append(sequences, NewCodonSequence(id, title, seqBuffer.String()))
				seqBuffer.Reset()
			}
			splitted = strings.SplitN(line[1:], " ", 2)
			id = splitted[0]
			if len(splitted) == 2 {
				title = splitted[1]
			}
		} else if strings.HasPrefix(line, "\n") {
			continue
		} else if strings.HasPrefix(line, "#") {
			continue
		} else {
			seqBuffer.WriteString(line)
		}
	}
	if seqBuffer.Len() > 0 {
		sequences = append(sequences, NewCodonSequence(id, title, seqBuffer.String()))
	}
	return
}

// SequencesToBuffer converts sequences in the sequene alignment into a buffered
// stream which can then be converted to bytes or a string.
func SequencesToBuffer(a SequenceAlignment) bytes.Buffer {
	var buffer bytes.Buffer
	// Append each Sequence in SequenceAlignment
	for _, s := range a {
		if len(s.Title()) > 0 {
			buffer.WriteString(fmt.Sprintf(">%s %s\n", s.ID(), s.Title()))
		} else {
			buffer.WriteString(fmt.Sprintf(">%s\n", s.ID()))
		}
		buffer.WriteString(s.Sequence() + "\n")
	}
	return buffer
}

// SequencesToString converts sequences in the sequene alignment into a FASTA
// formatted string.
func ProtSequencesToString(a SequenceAlignment) string {
	buffer := SequencesToBuffer(a)
	return buffer.String()
}

// SequencesToBuffer converts sequences in the sequene alignment into a buffered
// stream which can then be converted to bytes or a string.
func ProtSequencesToBuffer(a SequenceAlignment) bytes.Buffer {
	var buffer bytes.Buffer
	// Append each Sequence in SequenceAlignment
	for _, s := range a {
		if len(s.Title()) > 0 {
			buffer.WriteString(fmt.Sprintf(">%s %s\n", s.ID(), s.Title()))
		} else {
			buffer.WriteString(fmt.Sprintf(">%s\n", s.ID()))
		}
		v := reflect.ValueOf(s).Elem().FieldByName("prot").String()
		buffer.WriteString(v + "\n")
	}
	return buffer
}

// SequencesToString converts sequences in the sequene alignment into a FASTA
// formatted string.
func SequencesToString(a SequenceAlignment) string {
	buffer := SequencesToBuffer(a)
	return buffer.String()
}

// ConsistentAlignmentPositions returns the list of positions in the alignment
// that are considered consistent given by the alignment pattern per site
// across all given alignments.
func ConsistentAlignmentPositions(gapChar string, matrices ...[][]int) []bool {
	// Transpose matrices and combine as string
	patternSetMap := make(map[string]int)
	var templatePattern []string
	var patternBuffer bytes.Buffer
	for k, matrix := range matrices {
		for j := 0; j < len(matrix[0]); j++ {
			for i := 0; i < len(matrix); i++ {
				patternBuffer.WriteString(strconv.Itoa(matrix[i][j]) + ",")
			}

			if k == 0 {
				templatePattern = append(templatePattern, patternBuffer.String())
			}
			patternSetMap[patternBuffer.String()]++
			patternBuffer.Reset()
		}
	}

	pos := make([]bool, len(matrices[0][0]))
	for j, pattern := range templatePattern {
		// For the current column, compare the column pattern across the
		// matrices. If a difference in corresponding values are detected,
		// then the pattern is deemed inconsistent and the current column
		// is therefore also inconsistent.
		if patternSetMap[pattern] == len(matrices) {
			pos[j] = true
		} else {
			pos[j] = false
		}
	}
	return pos
}

// BufferedMarkedAlignment writes a marked multiple sequence alignment
// in the FASTA format to the buffer.
func BufferedMarkedAlignment(template SequenceAlignment, consistentPos []bool, markerID, consistentMarker, inconsistentMarker string) bytes.Buffer {
	var buffer bytes.Buffer

	// Append marker sequence
	buffer.WriteString(fmt.Sprintf(">%s\n", markerID))
	for _, t := range consistentPos {
		if t == true {
			buffer.WriteString(consistentMarker)
		} else {
			buffer.WriteString(inconsistentMarker)
		}
	}
	buffer.WriteString("\n")

	// Append each Sequence in SequenceAlignment
	for _, s := range template {
		if len(s.Title()) > 0 {
			buffer.WriteString(fmt.Sprintf(">%s %s\n", s.ID(), s.Title()))
		} else {
			buffer.WriteString(fmt.Sprintf(">%s\n", s.ID()))
		}
		buffer.WriteString(s.Sequence() + "\n")
	}
	return buffer
}

// ConsistentAlnPipeline aligns using global, local, and affine-local alignment
// strategies to determine positions that have a consistent alignment pattern over
// the three different strategies.
func ConsistentAlnPipeline(inputPath, gapChar, markerID, consistentMarker, inconsistentMarker string, iterations int, toUpper, toLower, saveTempAlns bool) bytes.Buffer {

	const mafftCmd = "mafft"

	ginsiAln := StringToCharSequences(GinsiAlign(mafftCmd, inputPath, iterations))
	linsiAln := StringToCharSequences(LinsiAlign(mafftCmd, inputPath, iterations))
	einsiAln := StringToCharSequences(EinsiAlign(mafftCmd, inputPath, iterations))

	if saveTempAlns == true {
		einsiAln.ToFasta(inputPath + ".einsi.aln")
		ginsiAln.ToFasta(inputPath + ".ginsi.aln")
		linsiAln.ToFasta(inputPath + ".linsi.aln")
	}

	consistentPos := ConsistentAlignmentPositions(
		gapChar,
		einsiAln.UngappedPositionMatrix(gapChar),
		ginsiAln.UngappedPositionMatrix(gapChar),
		linsiAln.UngappedPositionMatrix(gapChar),
	)

	if toUpper == true {
		einsiAln.ToUpper()
	} else if toLower == true {
		einsiAln.ToLower()
	}

	return BufferedMarkedAlignment(einsiAln, consistentPos, markerID, consistentMarker, inconsistentMarker)
}

func ConsistentCodonAlnPipeline(inputPath, gapChar, markerID, consistentMarker, inconsistentMarker string, iterations int, toUpper, toLower, saveTempAlns bool) bytes.Buffer {

	const mafftCmd = "mafft"

	// Translate FASTA to protein then align
	ginsiAln := StringToCodonSequences(GinsiCodonAlign(mafftCmd, inputPath, iterations))
	linsiAln := StringToCodonSequences(LinsiCodonAlign(mafftCmd, inputPath, iterations))
	einsiAln := StringToCodonSequences(EinsiCodonAlign(mafftCmd, inputPath, iterations))

	if saveTempAlns == true {
		einsiAln.ToFasta(inputPath + ".einsi.aln")
		ginsiAln.ToFasta(inputPath + ".ginsi.aln")
		linsiAln.ToFasta(inputPath + ".linsi.aln")
	}

	consistentPos := ConsistentAlignmentPositions(
		gapChar,
		einsiAln.UngappedPositionMatrix(gapChar),
		ginsiAln.UngappedPositionMatrix(gapChar),
		linsiAln.UngappedPositionMatrix(gapChar),
	)

	if toUpper == true {
		einsiAln.ToUpper()
	} else if toLower == true {
		einsiAln.ToLower()
	}

	return BufferedMarkedAlignment(einsiAln, consistentPos, markerID, consistentMarker, inconsistentMarker)
}

// WriteBufferToFile writes the contents of a buffer to file.
func WriteBufferToFile(path string, b bytes.Buffer) {
	f, err := os.Create(path)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	b.WriteTo(f)
	f.Sync()
}

// exists returns whether the given file or directory exists or not,
// and accompanying errors.
func exists(path string) (bool, error) {
	_, err := os.Stat(path)
	if err == nil {
		return true, nil
	}
	if os.IsNotExist(err) {
		return false, nil
	}
	return true, err
}

func main() {
	toUpper := false
	toLower := false

	// isCodonPtr := flag.Bool("codon", 0, "Create a codon-based alignment.")
	maxIterPtr := flag.Int("maxiterate", 1, "Maximum number of iterative refinement that MAFFT will perform.")
	gapCharPtr := flag.String("gapchar", "-", "Character in the alignment used to represent a gap.")
	markerIDPtr := flag.String("marker_id", "marker", "Name of marker sequence.")
	cMarkerPtr := flag.String("consistent_marker", "C", "Character to indicate a site is consistent across all alignment strategies.")
	icMarkerPtr := flag.String("inconsistent_marker", "N", "Character to indicate a site is inconsistent in at least one alignment strategy.")
	changeCasePtr := flag.String("change_case", "upper", "Change the case of the sequences. {upper|lower|no}")
	saveTempAlnPtr := flag.Bool("save_temp_alignments", false, "Save G-INSI, L-INSI, E-INSI alignments generated by MAFFT.")
	isBatchPtr := flag.String("batch", "", "Run in batch mode which reads files found in the specified folder.")
	inSuffixPtr := flag.String("input_suffix", ".fa", "Only files ending with this suffix will be processed. Used in conjunction with -batch.")
	outSuffixPtr := flag.String("output_suffix", ".aln", "Suffix to be appended to the end of the filename of resulting alignments. Used in conjunction with -batch.")
	outDirPtr := flag.String("outdir", "", "Output directory where alignments will be saved. Used in conjunction with -batch.")

	flag.Parse()

	if _, lookErr := exec.LookPath("mafft"); lookErr != nil {
		os.Stderr.WriteString("Error: \"mafft\" is not found in $PATH. Please make sure that mafft is installed and is accessible through the command \"mafft\".\n")
		os.Exit(1)
	}

	if len(*isBatchPtr) < 1 {
		// Single file mode
		args := flag.Args()

		if len(args) < 1 {
			os.Stderr.WriteString("Error: Missing path to FASTA file.\n")
			os.Exit(1)
		} else if len(args) > 1 {
			os.Stderr.WriteString("Error: More than 1 positional argument passed.\n")
			os.Exit(1)
		}

		if doesExist, _ := exists(args[0]); doesExist == false {
			os.Stderr.WriteString("Error: file does not exist.\n")
			os.Exit(1)
		}

		switch {
		case *changeCasePtr == "lower":
			toLower = true
		case *changeCasePtr == "upper":
			toUpper = true
		}

		buffer := ConsistentAlnPipeline(args[0], *gapCharPtr, *markerIDPtr, *cMarkerPtr, *icMarkerPtr, *maxIterPtr, toUpper, toLower, *saveTempAlnPtr)

		fmt.Print(buffer.String())

	} else {
		// Batch mode
		if doesExist, _ := exists(*isBatchPtr); doesExist == false {
			os.Stderr.WriteString("Error: Specified directory containing FASTA files does not exist.\n")
			os.Exit(1)
		}

		// Check if outdir flag used
		// Check if folder exists
		if len(*outDirPtr) < 1 {
			os.Stderr.WriteString("Error: Missing output directory.\nUse -outdir to specify an output directory where alignments will be saved.\n")
			os.Exit(1)
		}

		if doesExist, _ := exists(*outDirPtr); doesExist == false {
			os.Stderr.WriteString("Error: Specified output directory does not exist.\n")
			os.Exit(1)
		}

		// Read all fasta files in directory matching suffix
		files, err := filepath.Glob(*isBatchPtr + "/*" + *inSuffixPtr)
		if err != nil {
			panic(err)
		}
		var outputPath string
		for _, f := range files {
			fmt.Println(f)
			buffer := ConsistentAlnPipeline(f, *gapCharPtr, *markerIDPtr, *cMarkerPtr, *icMarkerPtr, *maxIterPtr, toUpper, toLower, *saveTempAlnPtr)
			outputPath = *outDirPtr + "/" + filepath.Base(f) + *outSuffixPtr
			WriteBufferToFile(outputPath, buffer)
		}
	}
}
