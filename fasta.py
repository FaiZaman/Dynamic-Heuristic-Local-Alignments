import numpy as np

def heuralign(alphabet, substitution_matrix, seq1, seq2):

  # defining parameters
  ktup = 2  # length of matches
  cutoff_score = 0  # cutoff score when scoring diagonals
  
  # get the index table and seeds
  index_table = get_index_table(ktup, seq1)
  diagonal_seeds = get_seeds(ktup, index_table, seq2)
  
  # score the diagonals 
  score_diagonals(ktup, cutoff_score, diagonal_seeds)


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
      updated = True

      # get the score for matching the current seed
      current_score = 0
      for k in range(0, ktup):
        seq1_letter = seq1[seed_i + k]
        seq2_letter = seq2[seed_j + k]
        current_score += substitution_matrix[alphabet.index(seq1_letter)][alphabet.index(seq2_letter)]

      max_score = current_score
      max_score_index_1 = seed_i
      max_score_index_2 = seed_j
      extending_left = True

      while updated:
        updated = False

        # extend left/right until score cutoff score and keep track of max score found
        if extending_left:
          seed_i -= 1
          seed_j -= 1
        else:
          seed_i += 1
          seed_j += 1
        if seed_i < 0 or seed_j < 0 or seed_i > len(seq1) or seed_j > len(seq2):
          break
        seq1_letter = seq1[seed_i]
        seq2_letter = seq2[seed_j]
        current_score += substitution_matrix[alphabet.index(seq1_letter)][alphabet.index(seq2_letter)]

        if current_score > max_score:
          updated = True
          max_score = current_score
          max_score_index_1 = seed_i
          max_score_index_2 = seed_j
        elif current_score < cutoff_score:
          extending_left = not(extending_left)
          seed_i = max_score_index_1
          seed_j = max_score_index_2
        print(max_score, max_score_index_1, max_score_index_2)

  

alphabet = "ACGT"
substitution_matrix = [[2, -1, -1, -1, -2], [-1, 2, -1, -1, -2], [-1, -1, 2, -1, -2], [-1, -1, -1, 2, -2], [-1, -1, -1, -1, -2]]
seq1 = "TATGCTA"
seq2 = "AGTACGCA"

alignments = heuralign(alphabet, substitution_matrix, seq1, seq2)