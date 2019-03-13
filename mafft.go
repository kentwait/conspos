package main

import (
	"io"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"

	fa "github.com/kentwait/gofasta"
)

// MAFFT functions for nucleotide and protein alignments

// ExecMafft calls the MAFFT program with the given arguments.
// Returns stdout as a string and nil if no errors are encountered.
// If an error occurs, returns a nil string and the error encountered.
func ExecMafft(mafftCmd string, args []string) (string, error) {
	// Check if mafftCmd is in $PATH and returns its absolute path.
	// However, if mafftCmd contains slashes, exec.LookPath assumes this is the absolute/relative path to the program.
	// Because of this, exec.LookPath will not look in $PATH and directly try to run mafftCmd using the given path.
	// Panics if LookPath returns an error
	absPath, lookErr := exec.LookPath(mafftCmd)
	if lookErr != nil {
		// TODO: Handle this without panicking, Return error instead
		panic(lookErr)
	}

	// Sets the number of threads MAFFT will use by getting the number of CPUs minus 1
	// TODO: Allow the user to set this manually and only use all CPUs when set to -1
	threads := runtime.NumCPU() - 1
	args = append([]string{"--thread", strconv.Itoa(threads)}, args...)

	// TODO: Add debug/verbose state which outputs the MAFFT call to stderr
	// Output MAFFT call
	// os.Stderr.WriteString(absPath + " " + strings.Join(args, " ") + "\n")

	// Builds the command to execute the MAFFT program from the path, and for the given set of arguments
	// This does not execute the program and the arguments yet.
	// It only creates the struct containing the necessary information to execute the program.
	cmd := exec.Command(absPath, args...)
	// Calling .Output() executes the command and captures stdout, discards sterr.
	stdout, err := cmd.Output()
	// Check if the program returned an error
	if err != nil {
		// TODO: Do not exit, return error instead.
		MafftError(err)
		os.Exit(1)
	}

	// No errors encountered, returns stdout as a string and nil error.
	return string(stdout), nil
}

// ExecMafftStdin calls the MAFFT program with the given arguments and using standard input as input.
// Returns stdout as a string and nil if no errors are encountered.
// If an error occurs, returns a nil string and the error encountered.
func ExecMafftStdin(mafftCmd string, stdin io.Reader, args []string) (string, error) {
	// TODO: Combine with ExecMafft?
	absPath, lookErr := exec.LookPath(mafftCmd)
	if lookErr != nil {
		panic(lookErr)
	}
	args = append(args, "-")

	cmd := exec.Command(absPath, args...)
	cmd.Stdin = stdin

	stdout, err := cmd.Output()
	if err != nil {
		MafftError(err)
		os.Exit(1)
	}
	return string(stdout), nil
}

// CharAlign calls MAFFT to align sequences depending on the specified alignment method.
func CharAlign(mafftCmd, fastaPath string, method string, iterations int) string {
	var methodFlag, indicatorChar string
	if method == "einsi" {
		methodFlag = "--genafpair"
		indicatorChar = "E"
	} else if method == "linsi" {
		methodFlag = "--localpair"
		indicatorChar = "L"
	} else if method == "ginsi" {
		methodFlag = "--globalpair"
		indicatorChar = "G"
	}
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		methodFlag,
		"--quiet",
		fastaPath,
	}...)
	// TODO: Add verbosity level to silence output
	os.Stderr.WriteString(indicatorChar)
	stdout, _ := ExecMafft(mafftCmd, args)
	return stdout
}

// CharAlignStdin calls MAFFT to align sequences depending on the specified alignment method coming from standard input.
func CharAlignStdin(mafftCmd string, r io.Reader, method string, iterations int) string {
	var methodFlag, indicatorChar string
	if method == "einsi" {
		methodFlag = "--genafpair"
		indicatorChar = "E"
	} else if method == "linsi" {
		methodFlag = "--localpair"
		indicatorChar = "L"
	} else if method == "ginsi" {
		methodFlag = "--globalpair"
		indicatorChar = "G"
	}
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		methodFlag,
		"--quiet",
	}...)
	// TODO: Add verbosity level to silence output
	os.Stderr.WriteString(indicatorChar)
	stdout, _ := ExecMafftStdin(mafftCmd, r, args)
	return stdout
}

// MAFFT functions for codon alignemnts

// CodonAlign calls MAFFT to align codon sequences depending on the specified alignment method.
func CodonAlign(mafftCmd, method string, fastaPath string, iterations int, c fa.Alignment) string {
	var methodFlag, indicatorChar string
	if method == "einsi" {
		methodFlag = "--genafpair"
		indicatorChar = "E"
	} else if method == "linsi" {
		methodFlag = "--localpair"
		indicatorChar = "L"
	} else if method == "ginsi" {
		methodFlag = "--globalpair"
		indicatorChar = "G"
	}
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		methodFlag,
		"--quiet",
		fastaPath,
	}...)
	// TODO: Add verbosity level to silence output
	os.Stderr.WriteString(indicatorChar)
	stdout, _ := ExecMafft(mafftCmd, args)
	// Create CharSequences from protein alignment
	p := fa.FastaToAlignment(strings.NewReader(stdout), false)

	// Use protein alignment to offset codons and match alignment. Output as Fasta string
	buff := AlignCodonsUsingProtAlignment(c, p)
	newStdout := buff.String()
	os.Stderr.WriteString("C")

	return newStdout
}

// CodonAlignStdin calls MAFFT to align codon sequences depending on the specified alignment method coming from standard input.
func CodonAlignStdin(mafftCmd string, r io.Reader, method string, iterations int, c fa.Alignment) string {
	var methodFlag, indicatorChar string
	if method == "einsi" {
		methodFlag = "--genafpair"
		indicatorChar = "E"
	} else if method == "linsi" {
		methodFlag = "--localpair"
		indicatorChar = "L"
	} else if method == "ginsi" {
		methodFlag = "--globalpair"
		indicatorChar = "G"
	}
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		methodFlag,
		"--quiet",
	}...)
	// TODO: Add verbosity level to silence output
	os.Stderr.WriteString(indicatorChar)
	stdout, _ := ExecMafftStdin(mafftCmd, r, args)

	// Create CharSequences from protein alignment
	p := fa.FastaToAlignment(strings.NewReader(stdout), false)

	// Use protein alignment to offset codons and match alignment. Output as Fasta string
	buff := AlignCodonsUsingProtAlignment(c, p)
	newStdout := buff.String()
	os.Stderr.WriteString("C")

	return newStdout
}
