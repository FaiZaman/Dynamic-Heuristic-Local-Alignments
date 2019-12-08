# Dynamic-Heuristic-Local-Alignments

Use dynamic programming and heuristic procedures to find DNA sequence alignments. There are three main files:

# dynprog.py

Uses dynamic programming (Smith-Waterman algorithm) in quadratic space and time to find the best local alignment of two sequences

# dynproglin.py

Uses dynamic programming (Smith-Waterman algorithm) in linear space and time (using a variation of Hirschberg's algorithm) to find the best local alignment of two sequences

# heuralign.py

Uses the FASTA algorithm with heuristic procedures to find an alignment as good as possible while trading off time

# To run the programs

Open one of the files and change the variables at the bottom (alphabet, substitution matrix, sequence 1, sequence 2) and then run the program. It will output the best score found and the indices of the sequences for the alignment
