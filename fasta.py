import numpy as np

def heuralign(alphabet, substitution_matrix, seq1, seq2):

	# defining parameters
	ktup = 2  # length of matches
	cutoff_score = 0  # cutoff score when scoring diagonals
	
	# get the index table and seeds
	index_table = get_index_table(ktup, seq1)
	diagonal_seeds = get_seeds(ktup, index_table, seq2)
	
	# score the diagonals 
	diagonal_score = score_diagonals(ktup, cutoff_score, diagonal_seeds)
	print(diagonal_score)


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


def score_diagonals(ktup, cutoff_score, diagonal_seeds):

	# initialise diagonal score dictionary
	diagonal_score = {}

	# initalise diagonal entries
	for diagonal in diagonal_seeds:
		diagonal_score[diagonal] = 0

	# loop through diagonals and find best score
	for diagonal in diagonal_seeds:
		for (seed_i, seed_j) in diagonal_seeds[diagonal]:
			seeds = diagonal_seeds[diagonal]
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
			
			extending_left = True

			while updated:
				updated = False

				# extend left/right until score cutoff score and keep track of max score found
				if extending_left:
					seq1_current_start_index -= 1
					seq2_current_start_index -= 1
				else:
					seq1_current_end_index += 1
					seq2_current_end_index += 1

				if seq1_current_start_index < 0 or seq2_current_start_index < 0 or seq1_current_end_index > len(seq1) or seq2_current_end_index > len(seq2):
					break

				# get score for a match
				seq1_letter = seq1[seed_i]
				seq2_letter = seq2[seed_j]
				current_score += substitution_matrix[alphabet.index(seq1_letter)][alphabet.index(seq2_letter)]

				if current_score > max_score:
					updated = True
					max_score = current_score
					if extending_left:
						seq1_best_start_index = seq1_current_start_index
						seq2_best_start_index = seq2_current_start_index
					else:
						seq1_best_end_index = seq1_current_end_index
						seq2_best_end_index = seq2_current_end_index
				elif current_score < cutoff_score:  # backtrack to indices for best score seen so far
					if extending_left:
						seq1_current_start_index = seq1_best_start_index
						seq2_current_start_index = seq2_best_start_index
					else:
						seq1_current_end_index = seq1_best_end_index
						seq2_current_end_index = seq2_best_end_index
						extending_left = not(extending_left)
				
		# if seeds absorbed then remove them from the diagonal dictionary
		for seed in range(0, len(diagonal_seeds[diagonal]) - 1):
			seeds = diagonal_seeds[diagonal]
			if seeds[seed][0] != seed_i and seeds[seed][1] != seed_j:  # not current seeds
				if seq1_best_start_index < seeds[seed][0] < seq1_best_end_index:
					del diagonal_seeds[diagonal][seed]

		diagonal_score[diagonal] = max(diagonal_score[diagonal], max_score)
	return diagonal_score

	

alphabet = "ABCD"
substitution_matrix = [[1, -5, -5, -5, -1],
				[-5, 1, -5, -5, -1],
				[-5, -5, 5, -5, -4],
				[-5, -5, -5, 6, -4],
				[-1, -1, -4, -4, -9]]
seq1 = "DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD"
seq2 = "DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD"

alignments = heuralign(alphabet, substitution_matrix, seq1, seq2)