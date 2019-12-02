import numpy as np

def dynproglin(alphabet, substitution_matrix, seq1, seq2):

    # gathering values
    p = len(alphabet)
    seq1_len = len(seq1)
    seq2_len = len(seq2)

   # initalise scoring matrix, backtracking matrix, and max score data
    seq1_alignment = ""
    seq2_alignment = ""

    max_score_data = [0, 0, 0]

    if seq1_len == 0:
        for i in range(0, seq2_len):
            seq1_alignment = seq1_alignment + '-'
            seq2_alignment = seq2_alignment + seq2[i]
    elif seq2_len == 0:
        for j in range(0, seq1_len):
            seq1_alignment = seq1_alignment + seq1[j]
            seq2_alignment = seq2_alignment + '-'
    elif seq1_len == 1 or seq2_len == 1:
        needlemen_alignments = needleman_wunsch(alphabet, substitution_matrix, seq1, seq2)
        seq1_alignment = needlemen_alignments[0]
        seq2_alignment = needlemen_alignments[1]
    else:
        seq2_mid = int(seq2_len/2)
        scoreL = NWScore(alphabet, substitution_matrix, seq1, seq2[0:seq2_mid], max_score_data)
        seq2_rev = seq2[seq2_mid:len(seq2)]
        seq2_rev = seq2_rev[::-1]
        seq1_rev = seq1[::-1]

        scoreR = NWScore(alphabet, substitution_matrix, seq1_rev, seq2_rev, max_score_data)
        scoreR = np.flip(scoreR)
        scoreSum = np.add(scoreL, scoreR)
        seq1_mid = np.argmax(scoreSum)

        alignment1 = dynproglin(alphabet, substitution_matrix, seq1[0:seq1_mid], seq2[0:seq2_mid])
        alignment2 = dynproglin(alphabet, substitution_matrix, seq1[seq1_mid:len(seq1)], seq2[seq2_mid:len(seq2)])

        seq1_alignment = alignment1[0] + alignment2[0]
        seq2_alignment = alignment1[1] + alignment2[1]

    return [seq1_alignment, seq2_alignment, max_score_data]


def NWScore(alphabet, substitution_matrix, seq1, seq2, max_score_data):

    scoring_matrix = np.zeros((2, len(seq1) + 1))
    
    for row in range(1, len(seq2) + 1):
        for column in range(1, len(seq1) + 1):
            score = calculate_score_data(row, column, substitution_matrix, scoring_matrix, seq1, seq2)
            scoring_matrix[1][column] = score[0]
            if score[0] > max_score_data[0]:
                max_score_data[0] = score[0]
                max_score_data[1] = row
                max_score_data[2] = column
        scoring_matrix[0,:] = scoring_matrix[1,:]

    last_line = scoring_matrix[1]
    return max_score_data


def calculate_score_data(row, column, substitution_matrix, scoring_matrix, seq1, seq2):

    # calculate and return the best score and its origin for the current scoring matrix cell
    seq1letter = seq1[column - 1]
    seq2letter = seq2[row - 1]

    match_score = substitution_matrix[alphabet.index(seq1letter)][alphabet.index(seq2letter)]

    if (len(seq1) > 1 and len(seq2) > 1):
        diagonal_score = scoring_matrix[0][column - 1] + match_score
        left_score = scoring_matrix[1][column - 1] + substitution_matrix[alphabet.index(seq1letter)][-1]
        up_score = scoring_matrix[0][column] + substitution_matrix[alphabet.index(seq2letter)][-1]
    else:
        diagonal_score = scoring_matrix[row - 1][column - 1] + match_score
        left_score = scoring_matrix[row][column - 1] + substitution_matrix[alphabet.index(seq1letter)][-1]
        up_score = scoring_matrix[row - 1][column] + substitution_matrix[alphabet.index(seq2letter)][-1]

    score = max(diagonal_score, up_score, left_score, 0)

    # 8 = DIAGONAL, 2 = UP, 4 = LEFT
    if score == diagonal_score:
        score_origin = 8
    elif score == up_score:
        score_origin = 2
    else:
        score_origin = 4

    return (score, score_origin)
  
def needleman_wunsch(alphabet, substitution_matrix, seq1, seq2):

    scoring_matrix = np.zeros((len(seq2) + 1, len(seq1) + 1))
    backtracking_matrix = np.zeros((len(seq2) + 1, len(seq1) + 1))

    for row in range(1, len(seq2) + 1):
        for column in range(1, len(seq1) + 1):
            score = calculate_score_data(row, column, substitution_matrix, scoring_matrix, seq1, seq2)
            scoring_matrix[row][column] = score[0]
            backtracking_matrix[row][column] = score[1]

    alignments = backtrack(len(seq2), len(seq1), backtracking_matrix, seq1, seq2)
    return alignments


def backtrack(row, column, backtracking_matrix, seq1, seq2):
    seq1_alignment = ""
    seq2_alignment = ""

    while row > 1 and column > 1:
        score_origin = backtracking_matrix[row][column]
        if score_origin == 8:
            seq1_alignment += seq1[column - 1]
            seq2_alignment += seq2[row - 1]
            row = row - 1
            column = column - 1
        elif score_origin == 2:
            seq1_alignment += '-'
            seq2_alignment += seq2[row - 1]
            row = row - 1
        elif score_origin == 4:
            seq1_alignment += seq1[column - 1]
            seq2_alignment += '-'
            column = column - 1
        else:
            break

    seq1_alignment = seq1_alignment[::-1]
    seq2_alignment = seq2_alignment[::-1]
    return (seq1_alignment, seq2_alignment)

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

alphabet = "AGCT"
substitution_matrix = [[2, -1, -1, -1, -2], [-1, 2, -1, -1, -2], [-1, -1, 2, -1, -2], [-1, -1, -1, 2, -2], [-1, -1, -1, -1, -2]]
#substitution_matrix = [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]]
#substitution_matrix = [[1, -1, -1, -1, 0], [-1, 1, -1, -1, 0], [-1, -1, 1, -1, 0], [-1, -1, -1, 1, 0], [-1, -1, -1, -1, 0]]
seq1 = "TATGC"
seq2 = "AGTACGCA"

alignments = NWScore(alphabet, substitution_matrix, seq1, seq2, [0, 0, 0])
print(alignments)
seq1_end_index = alignments[2]
seq2_end_index = alignments[1]

seq1_rev = seq1[::-1]
seq2_rev = seq2[::-1]
new_alignments = NWScore(alphabet, substitution_matrix, seq1_rev, seq2_rev, [0, 0, 0])
print(new_alignments)

start_score_row = new_alignments[1]
start_score_column = new_alignments[2]
seq1_start_index = len(seq1) - start_score_column
seq2_start_index = len(seq2) - start_score_row
print(seq1_start_index, seq2_start_index)
print(seq1_end_index, seq2_end_index)
print(seq1[seq1_start_index:seq1_end_index], seq2[seq2_start_index:seq2_end_index])

