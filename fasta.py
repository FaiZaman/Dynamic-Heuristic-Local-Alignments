import numpy as np

def heuralign(alphabet, substitution_matrix, seq1, seq2):

  # defining parameters
  ktup = 2
  
  index_table = get_index_table(ktup, seq1)
  seeds = get_seeds(ktup, index_table, seq2)
  


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
  seeds = []
  for letter_index in range(0, len(seq2) - ktup + 1):
      match = get_match(ktup, letter_index, seq2)

      if match in index_table:
        for seq1_position in index_table[match]:
          seeds.append((seq1_position, letter_index))

  return seeds


def get_match(ktup, letter_index, seq):

  # get string of length ktup
  match = seq[letter_index:letter_index + ktup]
  return match

alphabet = "ABC"
substitution_matrix = [[2, -1, -1, -1, -2], [-1, 2, -1, -1, -2], [-1, -1, 2, -1, -2], [-1, -1, -1, 2, -2], [-1, -1, -1, -1, -2]]
seq1 = "TATGCTA"
seq2 = "AGTACGCA"

alignments = heuralign(alphabet, substitution_matrix, seq1, seq2)