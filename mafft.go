package main

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"runtime"
	"strconv"
)

// MafftError writes to stderr that something has gone wrong with MAFFT.
func MafftError(err error) {
	os.Stderr.WriteString("Error: MAFFT did not exit properly (" + fmt.Sprint(err) + ")\nCheck if input sequences are valid.\n")
}

// MAFFT functions for nucleotide and protein alignments

// ExecMafft calls the MAFFT program with the given arguments.
func ExecMafft(mafftCmd string, args []string) (string, error) {
	absPath, lookErr := exec.LookPath(mafftCmd)
	if lookErr != nil {
		panic(lookErr)
	}

	threads := runtime.NumCPU() - 1
	args = append([]string{"--thread", strconv.Itoa(threads)}, args...)

	// Output MAFFT call
	// os.Stderr.WriteString(absPath + " " + strings.Join(args, " ") + "\n")

	cmd := exec.Command(absPath, args...)

	stdout, err := cmd.Output()
	if err != nil {
		MafftError(err)
		os.Exit(1)
	}

	return string(stdout), nil
}

// EinsiAlign calls MAFFT to align sequences by local alignment with
// affine-gap scoring.
func EinsiAlign(mafftCmd, fastaPath string, iterations int) string {
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		"--genafpair",
		"--quiet",
		fastaPath,
	}...)
	os.Stderr.WriteString("E")
	stdout, _ := ExecMafft(mafftCmd, args)
	return stdout
}

// LinsiAlign calls MAFFT to align sequences by local alignment.
func LinsiAlign(mafftCmd, fastaPath string, iterations int) string {
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		"--localpair",
		"--quiet",
		fastaPath,
	}...)
	os.Stderr.WriteString("L")
	stdout, _ := ExecMafft(mafftCmd, args)
	return stdout
}

// GinsiAlign calls MAFFT to align sequences by global alignment.
func GinsiAlign(mafftCmd, fastaPath string, iterations int) string {
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		"--globalpair",
		"--quiet",
		fastaPath,
	}...)
	os.Stderr.WriteString("G")
	stdout, _ := ExecMafft(mafftCmd, args)
	return stdout
}

// ExecMafftStdin calls the MAFFT program with the given arguments and using
// standard input as input.
func ExecMafftStdin(mafftCmd string, buff bytes.Buffer, args []string) (string, error) {
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
		MafftError(err)
		os.Exit(1)
	}
	return string(stdout), nil
}

// MAFFT functions for codon alignemnts

// EinsiCodonAlign calls MAFFT to align sequences by local alignment with
// affine-gap scoring.
func EinsiCodonAlign(mafftCmd string, buffer bytes.Buffer, iterations int, c SequenceAlignment) string {
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		"--genafpair",
		"--quiet",
	}...)

	// Pass Fasta string as stdin to mafft then capture stdout string
	stdout, _ := ExecMafftStdin(mafftCmd, buffer, args)

	// Create CharSequences from protein alignment
	p := StringToCodonSequences(stdout)

	// Use protein alignment to offset codons and match alignment. Output as Fasta string
	newStdout := AlignCodonsToString(c, p)

	return newStdout
}

// LinsiCodonAlign calls MAFFT to align sequences by local alignment.
func LinsiCodonAlign(mafftCmd string, buffer bytes.Buffer, iterations int, c SequenceAlignment) string {
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		"--localpair",
		"--quiet",
	}...)

	// Pass Fasta string as stdin to mafft then capture stdout string
	stdout, _ := ExecMafftStdin(mafftCmd, buffer, args)

	// Create CharSequences from protein alignment
	p := StringToCharSequences(stdout)

	// Use protein alignment to offset codons and match alignment. Output as Fasta string
	newStdout := AlignCodonsToString(c, p)

	return newStdout
}

// GinsiCodonAlign calls MAFFT to align sequences by global alignment.
func GinsiCodonAlign(mafftCmd string, buffer bytes.Buffer, iterations int, c SequenceAlignment) string {
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		"--globalpair",
		"--quiet",
	}...)

	// Pass Fasta string as stdin to mafft then capture stdout string
	stdout, _ := ExecMafftStdin(mafftCmd, buffer, args)

	// Create CharSequences from protein alignment
	p := StringToCharSequences(stdout)

	// Use protein alignment to offset codons and match alignment. Output as Fasta string
	newStdout := AlignCodonsToString(c, p)

	return newStdout
}
