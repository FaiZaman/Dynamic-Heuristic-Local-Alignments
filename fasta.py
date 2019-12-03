import numpy as np

def heuralign(alphabet, substitution_matrix, seq1, seq2):

  # defining parameters
  ktup = 2
  
  index_table = get_index_table(ktup, seq1)
  seeds = get_seeds(ktup, index_table, seq2)
  print(index_table)
  print(seeds)


def get_index_table(ktup, seq1):

  index_table = {}

  # go through seq1 and store ktup strings in index table
  for letter_index in range(0, len(seq1) - ktup + 1):      
    match = seq1[letter_index]
    for i in range(1, ktup):
      match += seq1[letter_index + i]

    if match in index_table:
      index_table[match].append(letter_index)
    else:
      index_table[match] = [letter_index]

  return index_table


def get_seeds(ktup, index_table, seq2):

  # go through seq2 and get seeds by matching to index table
  seeds = []
  for letter_index in range(0, len(seq2) - ktup + 1):
      match = seq2[letter_index]
      for i in range(1, ktup):
        match += seq2[letter_index + i]
      if match in index_table:
        for seq1_position in index_table[match]:
          seeds.append((seq1_position, letter_index))
  return seeds



alphabet = "ABC"
substitution_matrix = [[2, -1, -1, -1, -2], [-1, 2, -1, -1, -2], [-1, -1, 2, -1, -2], [-1, -1, -1, 2, -2], [-1, -1, -1, -1, -2]]
seq1 = "TATGCTA"
seq2 = "AGTACGCA"

alignments = heuralign(alphabet, substitution_matrix, seq1, seq2)