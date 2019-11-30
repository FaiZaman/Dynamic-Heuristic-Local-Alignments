import numpy as np

def dynprog(alphabet, substitution_matrix, seq1, seq2):
    
    # gathering values
    p = len(alphabet)
    seq1_len = len(seq1)
    seq2_len = len(seq2)

    # initalise scoring matrix, backtracking matrix, and max score data
    scoring_matrix = np.zeros((seq2_len + 1, seq1_len + 1))
    max_score, max_score_row, max_score_column = 0, -1, -1
    backtracking_matrix = np.zeros((seq2_len + 1, seq1_len + 1))
    
    # iterate through scoring matrix to fill it in while filling backtracking matrix
    for row in range(1, seq2_len + 1):
        for column in range(1, seq1_len + 1):
            score_data = calculate_score_data(row, column, substitution_matrix, scoring_matrix)
            score, score_origin = score_data[0], score_data[1]
            if score > max_score:
                max_score = score
                max_score_row = row
                max_score_column = column

            scoring_matrix[row][column] = score
            backtracking_matrix[row][column] = score_origin
    
    indices = get_indices(backtracking_matrix, max_score_row, max_score_column)
    return (max_score, indices[0], indices[1])


def calculate_score_data(row, column, substitution_matrix, scoring_matrix):

    # calculate and return the best score and its origin for the current scoring matrix cell
    seq1letter = seq1[column - 1]
    seq2letter = seq2[row - 1]

    match_score = substitution_matrix[alphabet.index(seq1letter)][alphabet.index(seq2letter)]

    diagonal_score = scoring_matrix[row - 1][column - 1] + match_score
    left_score = scoring_matrix[row][column - 1] + substitution_matrix[alphabet.index(seq2letter)][-1]
    up_score = scoring_matrix[row - 1][column] + substitution_matrix[alphabet.index(seq1letter)][-1]
    
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
    print(row, column)

    # iterate through backtracking matrix starting with cell which has the max score
    # iterate while collecting indices for the best alignment for both sequences
    while row > 1 and column > 1:
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
        print(row, column)
    
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


alphabet = "ABC"
substitution_matrix = [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]]
seq1 = "AABBAACA"
seq2 = "CBACCCBA"

a = dynprog(alphabet, substitution_matrix, seq1, seq2)
print("Score:   ", a[0])
print("Indices: ", a[1], a[2])