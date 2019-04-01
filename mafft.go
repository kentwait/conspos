package conspos

import (
	"io"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"
)

// ExecMafft calls the MAFFT program with the given arguments.
// Returns stdout as a string and nil if no errors are encountered.
// If an error occurs, returns a nil string and the error encountered.
func ExecMafft(mafftCmd string, args []string, verbosity int) (string, error) {
	// Check if mafftCmd is in $PATH and returns its absolute path.
	// However, if mafftCmd contains slashes, exec.LookPath assumes this is the absolute/relative path to the program.
	// Because of this, exec.LookPath will not look in $PATH and directly try to run mafftCmd using the given path.
	// Panics if LookPath returns an error
	absPath, lookErr := exec.LookPath(mafftCmd)
	if lookErr != nil {
		return "", lookErr
	}

	// Sets the number of threads MAFFT will use by getting the number of CPUs minus 1
	// TODO: Allow the user to set this manually and only use all CPUs when set to -1
	threads := runtime.NumCPU() - 1
	args = append([]string{"--thread", strconv.Itoa(threads)}, args...)

	// TODO: Add debug/verbose state which outputs the MAFFT call to stderr
	// Output MAFFT call
	if verbosity > 0 {
		os.Stderr.WriteString(absPath + " " + strings.Join(args, " ") + "\n")
	}

	// Builds the command to execute the MAFFT program from the path, and for the given set of arguments
	// This does not execute the program and the arguments yet.
	// It only creates the struct containing the necessary information to execute the program.
	cmd := exec.Command(absPath, args...)
	// Calling .Output() executes the command and captures stdout, discards sterr.
	stdout, err := cmd.Output()
	// Check if the program returned an error
	if err != nil {
		return string(stdout), err
	}

	// No errors encountered, returns stdout as a string and nil error.
	return string(stdout), nil
}

// CharAlign calls MAFFT to align sequences depending on the specified alignment method.
func CharAlign(mafftCmd, fastaPath string, method string, iterations, verbosity int) (string, error) {
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
	return ExecMafft(mafftCmd, args, verbosity)
}

// ExecMafftStdio calls the MAFFT program with the given arguments and using standard input as input.
// Returns stdout as a string and nil if no errors are encountered.
// If an error occurs, returns a nil string and the error encountered.
func ExecMafftStdio(mafftCmd string, stdin io.Reader, args []string, verbosity int) (string, error) {
	// TODO: Combine with ExecMafft?
	absPath, lookErr := exec.LookPath(mafftCmd)
	if lookErr != nil {
		return "", lookErr
	}
	args = append(args, "-")

	if verbosity > 0 {
		os.Stderr.WriteString(absPath + " " + strings.Join(args, " ") + "\n")
	}

	cmd := exec.Command(absPath, args...)
	cmd.Stdin = stdin

	stdout, err := cmd.Output()
	if err != nil {
		return string(stdout), err
	}
	return string(stdout), nil
}

// CharAlignStdio calls MAFFT to align sequences depending on the specified alignment method coming from standard input.
func CharAlignStdio(mafftCmd string, r io.Reader, method string, iterations, verbosity int) (string, error) {
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
	return ExecMafftStdio(mafftCmd, r, args, verbosity)
}
