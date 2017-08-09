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
	"sort"
	"strconv"
	"strings"
)

var bases = [4]string{"T", "C", "A", "G"}
var codons = [64]string{
	"TTT", "TTC", "TTA", "TTG",
	"TCT", "TCC", "TCA", "TCG",
	"TAT", "TAC", "TAA", "TAG",
	"TGT", "TGC", "TGA", "TGG",
	"CTT", "CTC", "CTA", "CTG",
	"CCT", "CCC", "CCA", "CCG",
	"CAT", "CAC", "CAA", "CAG",
	"CGT", "CGC", "CGA", "CGG",
	"ATT", "ATC", "ATA", "ATG",
	"ACT", "ACC", "ACA", "ACG",
	"AAT", "AAC", "AAA", "AAG",
	"AGT", "AGC", "AGA", "AGG",
	"GTT", "GTC", "GTA", "GTG",
	"GCT", "GCC", "GCA", "GCG",
	"GAT", "GAC", "GAA", "GAG",
	"GGT", "GGC", "GGA", "GGG",
}
var stopCodons = [3]string{"TGA", "TAG", "TAA"}
var aminoAcids = [20]string{
	"A",
	"R",
	"N",
	"D",
	"C",
	"Q",
	"E",
	"G",
	"H",
	"I",
	"L",
	"K",
	"M",
	"F",
	"P",
	"S",
	"T",
	"W",
	"Y",
	"V",
}
var geneticCode = map[string]string{
	"TTT": "F",
	"TTC": "F",
	"TTA": "L",
	"TTG": "L",
	"TCT": "L",
	"TCC": "L",
	"TCA": "L",
	"TCG": "L",
	"TAT": "I",
	"TAC": "I",
	"TAA": "I",
	"TAG": "M",
	"TGT": "V",
	"TGC": "V",
	"TGA": "V",
	"TGG": "V",
	"CTT": "S",
	"CTC": "S",
	"CTA": "S",
	"CTG": "S",
	"CCT": "P",
	"CCC": "P",
	"CCA": "P",
	"CCG": "P",
	"CAT": "T",
	"CAC": "T",
	"CAA": "T",
	"CAG": "T",
	"CGT": "A",
	"CGC": "A",
	"CGA": "A",
	"CGG": "A",
	"ATT": "Y",
	"ATC": "Y",
	"ATA": "*",
	"ATG": "*",
	"ACT": "H",
	"ACC": "H",
	"ACA": "Q",
	"ACG": "Q",
	"AAT": "N",
	"AAC": "N",
	"AAA": "K",
	"AAG": "K",
	"AGT": "D",
	"AGC": "D",
	"AGA": "E",
	"AGG": "E",
	"GTT": "C",
	"GTC": "C",
	"GTA": "*",
	"GTG": "W",
	"GCT": "R",
	"GCC": "R",
	"GCA": "R",
	"GCG": "R",
	"GAT": "S",
	"GAC": "S",
	"GAA": "R",
	"GAG": "R",
	"GGT": "G",
	"GGC": "G",
	"GGA": "G",
	"GGG": "G",
	"---": "-",
}

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

// ExecMafft calls the MAFFT program with the given arguments
func ExecMafft(mafftCmd string, args []string) string {
	absPath, lookErr := exec.LookPath(mafftCmd)
	if lookErr != nil {
		panic(lookErr)
	}

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

func FastaToCodonSequence(path string) SequenceAlignment {
	b, err := ioutil.ReadFile(path)
	if err != nil {
		panic(err)
	}
	return StringToCodonSequences(string(b))
}

func AlignCodonsByProt(codonSeq, protAln SequenceAlignment) SequenceAlignment {
	gapRune := []rune("-")[0]
	var newAln SequenceAlignment
	for i := range protAln {
		newCodons := []string{}
		newSeq := bytes.Buffer{}
		for j, c := range protAln[i].Sequence() {
			codons := reflect.ValueOf(codonSeq[i]).Elem().FieldByName("codons").String()
			cStart := j * 3
			cEnd := (j + 1) * 3
			if c == gapRune {
				newCodons = append(newCodons, "---")
				newSeq.WriteString("---")
			} else {
				newCodons = append(newCodons, codons[cStart:cEnd])
				newSeq.WriteString(codons[cStart:cEnd])
			}
		}
		newAln = append(newAln, &CodonSequence{CharSequence{codonSeq[i].ID(), codonSeq[i].Title(), newSeq.String()}, protAln[i].Sequence(), newCodons})
	}
	return newAln
}

// EinsiCodonAlign calls MAFFT to align sequences by local alignment with
// affine-gap scoring.
func EinsiCodonAlign(mafftCmd, fastaPath string, iterations int) (stdout string) {
	args := []string{
		"--quiet",
		"--genafpair",
		"--maxiterate",
		strconv.Itoa(iterations),
	}
	s := FastaToCodonSequence(fastaPath)
	buff := ProtSequencesToBuffer(s)
	stdout = ExecMafftStdin(mafftCmd, buff, args)

	return
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

// Sequence is an interface for single character sequences stored as a string
// and multi-character sequences stored as a slice.
type Sequence interface {
	ID() string
	Title() string
	Sequence() string
	Char(int) string
	SetSequence(string)
	ToUpper()
	ToLower()
	UngappedCoords(string) []int
	UngappedPositionSlice(string) []int
}

// CharSequence is a struct for nucleotide and single-letter protein sequences.
type CharSequence struct {
	id    string
	title string
	seq   string
}

func (s *CharSequence) ID() string {
	return s.id
}

func (s *CharSequence) Title() string {
	return s.title
}

func (s *CharSequence) Sequence() string {
	return s.seq
}

func (s *CharSequence) Char(i int) string {
	return string(s.seq[i])
}

func (s *CharSequence) SetSequence(seq string) {
	s.seq = seq
}

// UngappedCoords returns the positions in the sequence where the character
// does not match the gap character.
func (s *CharSequence) UngappedCoords(gapChar string) (colCoords []int) {
	set := make(map[int]struct{})
	for j := 0; j < len(s.seq); j++ {
		if string(s.seq[j]) != gapChar {
			set[j] = struct{}{}
		}
	}
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
	cnt := 0
	for j := 0; j < len(s.seq); j++ {
		if string(s.seq[j]) != gapChar {
			arr = append(arr, cnt)
			cnt++
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

// CodonSequence is a struct for triplet nucleotide codon sequences
type CodonSequence struct {
	CharSequence
	prot   string
	codons []string
}

func NewCodonSequence(id, title, seq string) *CodonSequence {
	if len(seq)%3 == 0 {
		panic("seq length not divisible by 3")
	}
	s := new(CodonSequence)
	s.id = id
	s.title = title
	s.SetSequence(seq)
	return s
}

func (s *CodonSequence) ID() string {
	return s.id
}

func (s *CodonSequence) Title() string {
	return s.title
}

func (s *CodonSequence) Sequence() string {
	return s.seq
}

func (s *CodonSequence) Codons() []string {
	return s.codons
}

func (s *CodonSequence) Prot() string {
	return s.prot
}

func (s *CodonSequence) Char(i int) string {
	return string(s.seq[i])
}

func (s *CodonSequence) ProtChar(i int) string {
	return string(s.prot[i])
}

func (s *CodonSequence) Codon(i int) string {
	return string(s.codons[i])
}

func (s *CodonSequence) SetSequence(seq string) {
	for i := 0; i < len(seq); i += 3 {
		s.codons = append(s.codons, string(seq[i:i+3]))
	}
	s.prot = (*Translate(seq)).String()
}

func (s *CodonSequence) SetCodons(seq []string) {
	s.codons = seq
}

func (s *CodonSequence) SetProt(seq string) {
	s.prot = seq
}

// UngappedCoords returns the positions in the sequence where the character
// does not match the gap character.
func (s *CodonSequence) UngappedCoords(gapChar string) (colCoords []int) {
	set := make(map[int]struct{})
	for j := 0; j < len(s.seq); j++ {
		if string(s.seq[j]) != gapChar {
			set[j] = struct{}{}
		}
	}
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
	cnt := 0
	for j := 0; j < len(s.seq); j++ {
		if string(s.seq[j]) != gapChar {
			arr = append(arr, cnt)
			cnt++
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
	for i := 0; i < len(s.seq); i++ {
		s.codons[i] = strings.ToUpper(s.codons[i])
	}
}

// ToLower changes the case of the sequence to all lowercase letters.
func (s *CodonSequence) ToLower() {
	s.seq = strings.ToLower(s.seq)
	s.prot = strings.ToLower(s.prot)
	for i := 0; i < len(s.seq); i++ {
		s.codons[i] = strings.ToLower(s.codons[i])
	}
}

// SequenceAlignment is a slice of Sequence pointers.
type SequenceAlignment []Sequence

// UngappedCoords returns the row and column positions in the sequence alignment
// where the character does not match the gap character.
func (a SequenceAlignment) UngappedCoords(gapChar string) (rowCoords, colCoords []int) {
	var currColCoords []int
	for i, s := range a {
		// s := reflect.ValueOf(sPtr)
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
func (a SequenceAlignment) UngappedPositionMatrix(gapChar string) (m [][]int) {
	for _, s := range a {
		m = append(m, s.UngappedPositionSlice(gapChar))
	}
	return
}

// ToUpper changes the case of all sequences to all uppercase letters.
func (a SequenceAlignment) ToUpper() {
	for _, s := range a {
		s.ToUpper()
	}
}

// ToLower changes the case of of all sequences to all lowercase letters.
func (a SequenceAlignment) ToLower() {
	for _, s := range a {
		s.ToLower()
	}
}

// ToFastaString returns the FASTA-formatted string of the sequence alignment.
func (a SequenceAlignment) ToFastaString() string {
	return SequencesToString(a)
}

// ToFasta saves the sequence alignment to a FASTA file.
func (a SequenceAlignment) ToFasta(path string) {
	WriteBufferToFile(path, SequencesToBuffer(a))
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
			// fmt.Println(pattern)
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

// ConsistentAlignmentPipeline aligns using global, local, and affine-local alignment
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
	maxIterPtr := flag.Int("maxiterate", 0, "Maximum number of iterative refinement that MAFFT will perform.")
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
