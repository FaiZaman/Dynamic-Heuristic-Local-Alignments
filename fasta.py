import numpy as np

def heuralign(alphabet, substitution_matrix, seq1, seq2):

  # defining parameters
  ktup = 2
  index_table = {}

  # go through seq1 and store ktup strings in index table
  for letter in range(0, len(seq1) - ktup + 1):      
    match = seq1[letter]
    for i in range(1, ktup):
      match += seq1[letter + i]

    if match in index_table:
      index_table[match].append(letter)
    else:
      index_table[match] = [letter]

  print(index_table)





alphabet = "ABC"
substitution_matrix = [[2, -1, -1, -1, -2], [-1, 2, -1, -1, -2], [-1, -1, 2, -1, -2], [-1, -1, -1, 2, -2], [-1, -1, -1, -1, -2]]
seq1 = "TATGCTA"
seq2 = "AGTACGCA"

alignments = heuralign(alphabet, substitution_matrix, seq1, seq2)