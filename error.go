package main

import (
	"fmt"
	"os"
)

// EmptyAlnError writes to stderr that the resulting alignment from MAFFT is
// empty.
func EmptyAlnError(alnType, inputPath string) {
	msg := fmt.Sprintf("Error: %s alignment from %s is empty.\nMAFFT may have encountered an error. Check if input sequences are valid.\n", alnType, inputPath)
	os.Stderr.WriteString(msg)
}
