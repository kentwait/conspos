package main

import (
	"bytes"
	"os"
	"os/exec"
	"runtime"
	"strconv"

	aln "github.com/kentwait/conspos/alignment"
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

// EinsiAlign calls MAFFT to align sequences by local alignment with affine-gap scoring.
func EinsiAlign(mafftCmd, fastaPath string, iterations int) string {
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		"--genafpair",
		"--quiet",
		fastaPath,
	}...)
	// TODO: Add verbosity level to silence output
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
	// TODO: Add verbosity level to silence output
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
	// TODO: Add verbosity level to silence output
	os.Stderr.WriteString("G")
	stdout, _ := ExecMafft(mafftCmd, args)
	return stdout
}

// ExecMafftStdin calls the MAFFT program with the given arguments and using standard input as input.
// Returns stdout as a string and nil if no errors are encountered.
// If an error occurs, returns a nil string and the error encountered.
func ExecMafftStdin(mafftCmd string, buff bytes.Buffer, args []string) (string, error) {
	// TODO: Combine with ExecMafft?
	absPath, lookErr := exec.LookPath(mafftCmd)
	if lookErr != nil {
		panic(lookErr)
	}
	args = append(args, "-")

	cmd := exec.Command(absPath, args...)
	cmd.Stdin = &buff

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
func EinsiCodonAlign(mafftCmd string, buffer bytes.Buffer, iterations int, c aln.Alignment) string {
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		"--genafpair",
		"--quiet",
	}...)

	// TODO: Add verbosity level to silence output
	os.Stderr.WriteString("E")
	// Pass Fasta string as stdin to mafft then capture stdout string
	stdout, _ := ExecMafftStdin(mafftCmd, buffer, args)

	// Create CharSequences from protein alignment
	p := StringToCharSequences(stdout)

	// Use protein alignment to offset codons and match alignment. Output as Fasta string
	newStdout := AlignCodonsToString(c, p)
	os.Stderr.WriteString("C")

	return newStdout
}

// LinsiCodonAlign calls MAFFT to align sequences by local alignment.
func LinsiCodonAlign(mafftCmd string, buffer bytes.Buffer, iterations int, c aln.Alignment) string {
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		"--localpair",
		"--quiet",
	}...)

	// TODO: Add verbosity level to silence output
	os.Stderr.WriteString("L")
	// Pass Fasta string as stdin to mafft then capture stdout string
	stdout, _ := ExecMafftStdin(mafftCmd, buffer, args)

	// Create CharSequences from protein alignment
	p := StringToCharSequences(stdout)

	// Use protein alignment to offset codons and match alignment. Output as Fasta string
	newStdout := AlignCodonsToString(c, p)
	os.Stderr.WriteString("C")

	return newStdout
}

// GinsiCodonAlign calls MAFFT to align sequences by global alignment.
func GinsiCodonAlign(mafftCmd string, buffer bytes.Buffer, iterations int, c aln.Alignment) string {
	var args []string
	args = append(args, []string{
		"--maxiterate",
		strconv.Itoa(iterations),
		"--globalpair",
		"--quiet",
	}...)

	// TODO: Add verbosity level to silence output
	os.Stderr.WriteString("G")
	// Pass Fasta string as stdin to mafft then capture stdout string
	stdout, _ := ExecMafftStdin(mafftCmd, buffer, args)

	// Create CharSequences from protein alignment
	p := StringToCharSequences(stdout)

	// Use protein alignment to offset codons and match alignment. Output as Fasta string
	newStdout := AlignCodonsToString(c, p)
	os.Stderr.WriteString("C")

	return newStdout
}
