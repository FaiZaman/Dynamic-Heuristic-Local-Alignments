import numpy as np
import sys
from heapq import nlargest
np.set_printoptions(threshold=sys.maxsize)

def heuralign(alphabet, substitution_matrix, seq1, seq2):

	# defining parameters
	ktup = 2  # length of matches
	cutoff_score = -3  # cutoff score when scoring diagonals
	width = 2	 # width of band for banded DP
	
	# get the index table and seeds
	index_table = get_index_table(ktup, seq1)
	diagonal_seeds = get_seeds(ktup, index_table, seq2)

	# score the diagonals
	diagonal_score = score_diagonals(alphabet, substitution_matrix, seq1, seq2, ktup, cutoff_score, diagonal_seeds)

	# get the best 3 diagonals and run banded DP on them
	best_diagonals = nlargest(3, diagonal_score, key=diagonal_score.get)
	banded_DP(alphabet, substitution_matrix, seq1, seq2, best_diagonals, width)


def get_index_table(ktup, seq1):

	index_table = {}

	# go through seq1 and store ktup strings in index table
	for letter_index in range(0, len(seq1) - ktup + 1):      
		match = seq1[letter_index:letter_index + ktup]

		if match in index_table:
			index_table[match].append(letter_index)
		else:
			index_table[match] = [letter_index]

	return index_table


def get_seeds(ktup, index_table, seq2):

	# go through seq2 and get seeds by matching to index table
	seeds = {}

	for letter_index in range(0, len(seq2) - ktup + 1):
		match = seq2[letter_index:letter_index + ktup]

		if match in index_table:
			for seq1_position in index_table[match]:
				diagonal = seq1_position - letter_index
				if diagonal in seeds:
					seeds[diagonal].append((seq1_position, letter_index))
				else:
					seeds[diagonal] = [(seq1_position, letter_index)]

	return seeds


def score_diagonals(alphabet, substitution_matrix, seq1, seq2, ktup, cutoff_score, diagonal_seeds):

	# initialise diagonal score dictionary
	diagonal_score = {}

	# initalise diagonal entries
	for diagonal in diagonal_seeds:
		diagonal_score[diagonal] = 0

	# loop through diagonals and find best score
	for diagonal in diagonal_seeds:
		for (seed_i, seed_j) in diagonal_seeds[diagonal]:
			updated = True

			# get the score for matching the current seed
			current_score = 0
			for k in range(0, ktup):
				seq1_letter = seq1[seed_i + k]
				seq2_letter = seq2[seed_j + k]
				current_score += substitution_matrix[alphabet.index(seq1_letter)][alphabet.index(seq2_letter)]

			# initialise max score and current and best endpoints of alignments
			max_score = current_score

			seq1_current_start_index = seed_i
			seq2_current_start_index = seed_j
			seq1_current_end_index = seed_i + ktup
			seq2_current_end_index = seed_j + ktup

			seq1_best_start_index = seed_i
			seq2_best_start_index = seed_j
			seq1_best_end_index = seed_i + ktup
			seq2_best_end_index = seed_j + ktup      
			
			while updated:

				updated = False
				
				while current_score > cutoff_score:	# extend left

					seq1_current_start_index -= 1
					seq2_current_start_index -= 1
					if seq1_current_start_index < 0 or seq2_current_start_index < 0:
						break
					
					seq1_letter = seq1[seq1_current_start_index]
					seq2_letter = seq2[seq2_current_start_index]
					# get score for a match
					current_score += substitution_matrix[alphabet.index(seq1_letter)][alphabet.index(seq2_letter)]

					if current_score > max_score:
						updated = True
						max_score = current_score
						seq1_best_start_index = seq1_current_start_index
						seq2_best_start_index = seq2_current_start_index

				# reset to best score and indices
				seq1_current_start_index = seq1_best_start_index
				seq2_current_start_index = seq1_best_start_index

				while current_score > cutoff_score:	# extend right
		
					seq1_current_end_index += 1
					seq2_current_end_index += 1
					if seq1_current_end_index > len(seq1) - 1 or seq2_current_end_index > len(seq2) - 1:
						break
					
					seq1_letter = seq1[seq1_current_end_index]
					seq2_letter = seq2[seq2_current_end_index]
					# get score for a match
					current_score += substitution_matrix[alphabet.index(seq1_letter)][alphabet.index(seq2_letter)]

					if current_score > max_score:
						updated = True
						max_score = current_score
						seq1_best_end_index = seq1_current_end_index
						seq2_best_end_index = seq2_current_end_index

				seq1_current_end_index = seq1_best_end_index
				seq2_current_end_index = seq2_best_end_index
				
		# if seeds absorbed then remove them from the diagonal dictionary
		for (seed_k, seed_l) in diagonal_seeds[diagonal]:
			index = diagonal_seeds[diagonal].index((seed_k, seed_l))
			if seed_k != seed_i and seed_l != seed_j:  # not current seeds
				if seq1_best_start_index < seed_k < seq1_best_end_index:
					del diagonal_seeds[diagonal][index]

		diagonal_score[diagonal] = max(diagonal_score[diagonal], max_score)
	return diagonal_score

def banded_DP(alphabet, substitution_matrix, seq1, seq2, best_diagonals, width):

	best_scores = []
	for diagonal in best_diagonals:
		
		# initialise max score data and got band values
		max_score, max_score_row, max_score_column = 0, -1, -1
		diagonal *= -1
		upper_diagonal = diagonal + width
		lower_diagonal = diagonal - width
		
		# initialise matrices using np.empty therefore being subquadratic time
		scoring_matrix = np.empty((len(seq2) + 1, len(seq1) + 1))
		backtracking_matrix = np.empty((len(seq2) + 1, len(seq1) + 1))
		
		# initialise first row and column
		scoring_matrix[0] = 0
		scoring_matrix[:,0] = 0

		# run local alignment on the cells in the band
		for row in range(1, len(seq2) + 1):
			for column in range(max(row - upper_diagonal, 1), min(row - lower_diagonal, len(seq1)) + 1):

				# get the score and where it comes from
				score_data = calculate_score_data(row, column, alphabet, substitution_matrix, scoring_matrix, seq1, seq2)
				score, score_origin = score_data[0], score_data[1]

				# replace max score data if greater
				if score > max_score:
					max_score = score
					max_score_row = row
					max_score_column = column

				scoring_matrix[row][column] = score
				backtracking_matrix[row][column] = score_origin

		best_scores.append(max_score)

	print(scoring_matrix)


def calculate_score_data(row, column, alphabet, substitution_matrix, scoring_matrix, seq1, seq2):

	# calculate and return the best score and its origin for the current scoring matrix cell
    seq1_letter = seq1[column - 1]
    seq2_letter = seq2[row - 1]

    match_score = substitution_matrix[alphabet.index(seq1_letter)][alphabet.index(seq2_letter)]

    diagonal_score = scoring_matrix[row - 1][column - 1] + match_score
    left_score = scoring_matrix[row][column - 1] + substitution_matrix[alphabet.index(seq1_letter)][-1]
    up_score = scoring_matrix[row - 1][column] + substitution_matrix[alphabet.index(seq2_letter)][-1]
    
	# check if the max score is out of bounds
    score = max(diagonal_score, up_score, left_score, 0)
    score_origin = 0

    # 8 = DIAGONAL, 2 = UP, 4 = LEFT
    if score == diagonal_score:
        score_origin = 8
    elif score == up_score:
        score_origin = 2
    else:
        score_origin = 4

    return (score, score_origin)


def get_indices(backtracking_matrix, row, column):
    
    seq1_indices = []
    seq2_indices = []
    seq1_alignment = ""
    seq2_alignment = ""

    # iterate through backtracking matrix starting with cell which has the max score
    # iterate while collecting indices for the best alignment for both sequences
    while row > 0 and column > 0:
        score_origin = backtracking_matrix[row][column]

        if score_origin == 8:
            seq1_alignment += seq1[column - 1]
            seq2_alignment += seq2[row - 1]
            row = row - 1
            column = column - 1
            seq1_indices.append(column)
            seq2_indices.append(row)
        elif score_origin == 2:
            seq1_alignment += '-'
            seq2_alignment += seq2[row - 1]
            row = row - 1
        else:
            seq1_alignment += seq1[column - 1]
            seq2_alignment += '-'
            column = column - 1
    
    seq1_indices.sort()
    seq2_indices.sort()
    seq1_alignment = seq1_alignment[::-1]
    seq2_alignment = seq2_alignment[::-1]
    displayAlignment([seq1_alignment, seq2_alignment])
    return (seq1_indices, seq2_indices)


def displayAlignment(alignment):
    string1 = alignment[0]
    string2 = alignment[1]
    string3 = ''
    for i in range(min(len(string1), len(string2))):
        if string1[i] == string2[i]:
            string3 = string3 + "|"
        else:
            string3 = string3 + " "
    print('String1: ' + string1)
    print('         ' + string3)
    print('String2: ' + string2 + '\n\n')

alphabet = "ABCD"
substitution_matrix = [[1, -5, -5, -5, -1],
					   [-5, 1, -5, -5, -1],
					   [-5, -5, 5, -5, -4],
					   [-5, -5, -5, 6, -4],
					   [-1, -1, -4, -4, -9]]
seq1 = "DDCDDCCCDCAA"
seq2 = "DDCDDCCCDCBC"

alignments = heuralign(alphabet, substitution_matrix, seq1, seq2)