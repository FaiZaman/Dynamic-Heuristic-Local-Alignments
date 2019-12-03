import numpy as np

def heuralign(alphabet, substitution_matrix, seq1, seq2):

  # defining parameters
  ktup = 2  # length of matches
  cutoff = -3  # cutoff when scoring diagonals
  
  # get the index table and seeds
  index_table = get_index_table(ktup, seq1)
  diagonal_seeds = get_seeds(ktup, index_table, seq2)
  
  # score the diagonals 
  score_diagonals(ktup, cutoff, diagonal_seeds)


def get_index_table(ktup, seq1):

  index_table = {}

  # go through seq1 and store ktup strings in index table
  for letter_index in range(0, len(seq1) - ktup + 1):      
    match = get_match(ktup, letter_index, seq1)

    if match in index_table:
      index_table[match].append(letter_index)
    else:
      index_table[match] = [letter_index]

  return index_table


def get_seeds(ktup, index_table, seq2):

  # go through seq2 and get seeds by matching to index table
  seeds = {}
  for letter_index in range(0, len(seq2) - ktup + 1):
      match = get_match(ktup, letter_index, seq2)

      if match in index_table:
        for seq1_position in index_table[match]:
          diagonal = seq1_position - letter_index
          if diagonal in seeds:
            seeds[diagonal].append((seq1_position, letter_index))
          else:
            seeds[diagonal] = [(seq1_position, letter_index)]

  return seeds


def score_diagonals(ktup, cutoff, diagonal_seeds):

  # initialise max score data
  max_score = 0
  max_score_index_1 = -1
  max_score_index_2 = -1

  # loop through diagonals and find best score
  for diagonal in diagonal_seeds:
    for (seed_i, seed_j) in diagonal_seeds[diagonal]:
      updated = True

      # get the score for matching the current seed
      seed_score = 0
      for k in range(0, ktup):
        seq1_letter = seq1[seed_i + k]
        seq2_letter = seq2[seed_j + k]
        seed_score += substitution_matrix[alphabet.index(seq1_letter)][alphabet.index(seq2_letter)]
        
      while updated:
        updated = False

        # extend left until score cutoff and keep track of max score found

        

    



def get_match(ktup, letter_index, seq):

  # get string of length ktup
  match = seq[letter_index:letter_index + ktup]
  return match

alphabet = "ACGT"
substitution_matrix = [[2, -1, -1, -1, -2], [-1, 2, -1, -1, -2], [-1, -1, 2, -1, -2], [-1, -1, -1, 2, -2], [-1, -1, -1, -1, -2]]
seq1 = "TATGCTA"
seq2 = "AGTACGCA"

alignments = heuralign(alphabet, substitution_matrix, seq1, seq2)