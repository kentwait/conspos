package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"os/exec"
	"strings"

	. "github.com/kentwait/conspos"
	fa "github.com/kentwait/gofasta"
)

func main() {
	// Required flags
	modePtr := flag.String("mode", "char", "Perform conspos for single-character sites (char) or codons (codon). {char|codon}")

	// Optional conspos flags
	markerIDPtr := flag.String("marker_id", "marker", "Name of marker sequence.")
	markerCharsPtr := flag.String("markers", "CIN", "Characters to indicate a site is fully consistent across all alignment strategies (3/3), partially consistent (2/3), or completely inconsistent (1/3). If only 2 characters are included, sites are classified as either consistent (3/3) or inconsistent (not 3/3).")
	gapCharPtr := flag.String("gap_char", "-", "Character in the alignment used to represent a gap.")
	changeCasePtr := flag.String("change_case", "upper", "Change the case of the sequences. {upper|lower|no}")
	verbosityPtr := flag.Int("verbose", 0, "Amount of notifications to show. (0-5)")

	// Optional MAFFT-related flags
	maxIterPtr := flag.Int("max_iterate", 1, "Maximum number of iterative refinement that MAFFT will perform.")
	tempPathPtr := flag.String("temp", "", "Save path for G-INSI, L-INSI, E-INSI alignments generated by MAFFT. Appends .ginsi.aln/.linsi.aln/.einsi.aln suffix to the given path.")
	mafftPathPtr := flag.String("mafft_path", "mafft", "Path to MAFFT executable. If MAFFT is registered in $PATH, you can use \"mafft\".")

	flag.Parse()

	// Checks if values of arguments are valid.

	// Validates supplied path for MAFFT executable.
	// Raises an error and exits if the path does not exist.
	if _, lookErr := exec.LookPath(*mafftPathPtr); lookErr != nil {
		os.Stderr.WriteString("[InputError]\tInvalid MAFFT path.\n" + "            \tMake sure that the MAFFT executable is installed and is accessible at the path specified in -mafft_path.\n")
		os.Exit(1)
	}

	// Checks whether there is at least one positional argument present.
	// Raises an error and exits if no positional arguments are present, or when more than one is given.
	args := flag.Args()
	if len(args) == 0 {
		os.Stderr.WriteString("[InputError]\tMissing path to FASTA file.\n")
		os.Exit(1)
	} else if len(args) > 1 {
		os.Stderr.WriteString("[InputError]\tMore than 1 positional argument passed.\n")
		os.Exit(1)
	}
	// Given that there is only one positional argument supplied, checks whether a file exists at that path.
	// This does not check whether the file is a FASTA file though.
	if doesExist, _ := Exists(args[0]); doesExist == false {
		os.Stderr.WriteString(fmt.Sprintf("[FileError]\tFile does not exist: %s\n", args[0]))
		os.Exit(1)
	}

	// Converts case change choices to boolean variables.
	var toUpper, toLower bool
	switch *changeCasePtr {
	case "lower":
		toLower = true
	case "upper":
		toUpper = true
	case "no":
	default:
		os.Stderr.WriteString(fmt.Sprintf("[InputError]\tInvalid -change_case value: %s\n", *changeCasePtr) + "            \t{upper|lower|no}.\n")
		os.Exit(1)
	}

	// Convert marker chars string to individual characters
	var markerChars []string
	for _, runeChar := range *markerCharsPtr {
		markerChars = append(markerChars, string(runeChar))
	}
	c, i, n := markerChars[0], markerChars[1], markerChars[2]

	switch *modePtr {
	case "char":
		file, err := os.Open(args[0])
		if err != nil {
			log.Fatal(err)
		}
		defer file.Close()

		// Mafft pipeline
		ginsiAln, linsiAln, einsiAln := SingleCharMafftPipeline(file, *mafftPathPtr, *tempPathPtr, *maxIterPtr, *verbosityPtr)
		// Scoring pipeline
		markerSequence := ScoringPipeline(*gapCharPtr, c, i, n, ginsiAln, linsiAln, einsiAln)
		// Final alignment
		var finalAln []*fa.CharSequence
		finalAln = append(finalAln, fa.NewCharSequence(*markerIDPtr, "", strings.Join(markerSequence, "")))
		finalAln = append(finalAln, einsiAln...)

	case "codon":
		file, err := os.Open(args[0])
		if err != nil {
			log.Fatal(err)
		}
		defer file.Close()

		// Mafft pipeline
		ginsiAln, linsiAln, einsiAln := CodonMafftPipeline(file, *mafftPathPtr, *tempPathPtr, *tempPathPtr, *maxIterPtr, *verbosityPtr)
		// Scoring pipeline
		markerSequence := ScoringPipeline(*gapCharPtr, c, i, n, ginsiAln, linsiAln, einsiAln)
		// Final alignment
		var finalAln []*fa.CharSequence
		finalAln = append(finalAln, fa.NewCharSequence(*markerIDPtr, "", strings.Join(markerSequence, "")))
		finalAln = append(finalAln, einsiAln...)

	case "":
		os.Stderr.WriteString("[InputError]\tSelect an alignment mode {char|codon} using -mode.\n")
		os.Exit(1)
	default:
		os.Stderr.WriteString(fmt.Sprintf("[InputError]\tInvalid -mode value: %s\n", *modePtr) + "            \t{char|codon}.\n")
		os.Exit(1)
	}

}