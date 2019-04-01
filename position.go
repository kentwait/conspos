package conspos

import (
	"bytes"
	"strconv"
)

// ScorePositions takes a list of position matrices and returns the the consistency scores using the first matrix for positional reference.
func ScorePositions(matrices ...[][]int) []int {
	patternMap := make(map[string]int)
	var templatePattern []string
	var patternBuffer bytes.Buffer
	/*
		Loops over all ungapped position matrices passed to the function.
		Assumes that all matrices have the same shape.
		For two matrices given, each matrix will look like this:
		a := [[0, 1, 2, 3,-1, 4, 5, 6,-1, 7,...]
			  [0, 1, 2, 3, 4,-1, 5, 6,-1, 7,...]
			  [0, 1, 2, 3,-1,-1, 4, 5, 6, 7,...]]
		b := [[0, 1, 2, 3,-1, 4, 5, 6,-1, 7,...]
			  [0, 1, 2, 3,-1, 4, 5, 6,-1, 7,...]
			  [0, 1, 2, 3,-1,-1, 4, 5, 6, 7,...]]

		This loop iterates over each matrix, and transposes the matrix such that:
		at := [[ 0, 0, 0],
			   [ 1, 1, 1],
			   [ 2, 2, 2],
			   [ 3, 3, 3],
			   [-1, 4,-1],
			   [ 5, 5, 4],
			   [ 6, 6, 5],
			   [-1,-1, 6],
			   [ 7, 7, 7],
				   ...   ]
		bt := [[ 0, 0, 0],
			   [ 1, 1, 1],
			   [ 2, 2, 2],
			   [ 3, 3, 3],
			   [-1,-1,-1],
			   [ 4, 4,-1],
			   [ 5, 5, 4],
			   [ 6, 6, 5],
			   [-1,-1, 6],
			   [ 7, 7, 7],
				   ...   ]

		The goal is to compare whether the pattern in the rows are consistent across all matrices.
		Comparing rows means comparing slices. In Golang, slices cannot be immediately compared, and must use a loop to do an itemwise comparison.
		As a workaround, each row is encoded as a string whose values are separated by a comma "," like this:
		at := ["0,0,0",
			   "1,1,1",
			   "2,2,2",
			   "3,3,3",
			   "-1,4,-1",
			   "5,5,4",
			   "6,6,5",
			   "-1,-1,6",
			   "7,7,7",
				   ...   ]
		bt := ["0,0,0",
			   "1,1,1",
			   "2,2,2",
			   "3,3,3",
			   "-1,-1,-1",
			   "4,4,-1",
			   "5,5,4",
			   "6,6,5",
			   "-1,-1,6",
			   "7,7,7",
				   ...   ]

		Through this method, each row becomes a string and multiple comparisons becomes straightforward.

		"0,0,0" == "0,0,0" = true
		"1,1,1" == "1,1,1" = true
		"1,1,1" == "1,1,1" = true
		"2,2,2" == "2,2,2" = true
		"3,3,3" == "3,3,3" = true
		"-1,4,-1" == "-1,-1,-1" = false
		...

	*/

	// i is rows, j is columns, k is matrix number
	for k, matrix := range matrices {
		for j := 0; j < len(matrix[0]); j++ {
			// Gets the pattern for each column in the matrix
			for i := 0; i < len(matrix); i++ {
				patternBuffer.WriteString(strconv.Itoa(matrix[i][j]) + ",")
			}

			// If this is the first matrix, store the pattern as a key
			// If it is NOT the first matrix, add only is the pattern exists
			patternString := patternBuffer.String()
			if k == 0 {
				patternMap[patternString]++
				templatePattern = append(templatePattern, patternString)
			} else if _, ok := patternMap[patternString]; ok {
				patternMap[patternString]++
			}
			patternBuffer.Reset()
		}
	}

	scoreSlice := make([]int, len(matrices[0][0]))
	for j, pattern := range templatePattern {
		// For the current column, get how many times the template pattern was observed.
		scoreSlice[j] = patternMap[pattern]
	}
	return scoreSlice
}

// MarkConsistent returns a slice of 1-character strings denoting whether the position is consistent, always inconsistent, or sometimes inconsistent.
func MarkConsistent(scoreSlice []int, consScore int, consChar, interChar, inconsChar string) (markSlice []string) {
	for _, score := range scoreSlice {
		switch {
		case score == consScore:
			markSlice = append(markSlice, consChar)
		case score == 1:
			markSlice = append(markSlice, inconsChar)
		default:
			markSlice = append(markSlice, interChar)
		}
	}
	return
}

// MarkOnlyConsistent returns a slice of 1-character strings denoting whether the position is consistent or inconsistent.
func MarkOnlyConsistent(scoreSlice []int, consScore int, consChar, inconsChar string) (markSlice []string) {
	for _, score := range scoreSlice {
		switch {
		case score == consScore:
			markSlice = append(markSlice, consChar)
		default:
			markSlice = append(markSlice, inconsChar)
		}
	}
	return
}
