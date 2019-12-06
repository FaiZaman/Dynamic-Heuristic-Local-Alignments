import numpy as np

def dynproglin(alphabet, substitution_matrix, seq1, seq2):
    
    # get score and end indices of local alignment
    end_score_data = NWLocalScore(alphabet, substitution_matrix, seq1, seq2)
    score = end_score_data[0]
    seq1_end_index = end_score_data[2]
    seq2_end_index = end_score_data[1]

    # reverse sequences and run again
    seq1_rev = seq1[::-1]
    seq2_rev = seq2[::-1]
    start_score_data = NWLocalScore(alphabet, substitution_matrix, seq1_rev, seq2_rev)

    # get start indices and build local sequences
    seq1_start_index = len(seq1) - start_score_data[2]
    seq2_start_index = len(seq2) - start_score_data[1]
    seq1_local = seq1[seq1_start_index:seq1_end_index]
    seq2_local = seq2[seq2_start_index:seq2_end_index]

    # run Hirschberg global alignment on local sequences
    alignments = Hirschberg(alphabet, substitution_matrix, seq1_local, seq2_local)
    seq1_alignment = alignments[0]
    seq2_alignment = alignments[1]

    (seq1_indices, seq2_indices) = get_indices(seq1_alignment, seq2_alignment, seq1_start_index, seq2_start_index)
    return int(score), seq1_indices, seq2_indices


def Hirschberg(alphabet, substitution_matrix, seq1, seq2):

    # gathering values
    p = len(alphabet)
    seq1_len = len(seq1)
    seq2_len = len(seq2)

   # initalise scoring matrix, backtracking matrix, and max score data
    seq1_alignment = ""
    seq2_alignment = ""

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
        scoreL = NWScore(alphabet, substitution_matrix, seq1, seq2[0:seq2_mid])
        seq2_rev = seq2[seq2_mid:len(seq2)]
        seq2_rev = seq2_rev[::-1]
        seq1_rev = seq1[::-1]

        scoreR = NWScore(alphabet, substitution_matrix, seq1_rev, seq2_rev)
        scoreR = np.flip(scoreR)
        scoreSum = np.add(scoreL, scoreR)
        seq1_mid = np.argmax(scoreSum)

        alignment1 = Hirschberg(alphabet, substitution_matrix, seq1[0:seq1_mid], seq2[0:seq2_mid])
        alignment2 = Hirschberg(alphabet, substitution_matrix, seq1[seq1_mid:len(seq1)], seq2[seq2_mid:len(seq2)])

        seq1_alignment = alignment1[0] + alignment2[0]
        seq2_alignment = alignment1[1] + alignment2[1]

    return [seq1_alignment, seq2_alignment]


def NWScore(alphabet, substitution_matrix, seq1, seq2):

    p = len(alphabet)
    scoring_matrix = np.zeros((2, len(seq1) + 1))
    
    for column in range(1, len(seq1) + 1):
        letter = seq1[column - 1]
        scoring_matrix[0][column] = scoring_matrix[0][column - 1] + substitution_matrix[alphabet.index(letter)][p]
        
    for row in range(1, len(seq2) + 1):
        letter = seq2[row - 1]
        scoring_matrix[1][0] = scoring_matrix[0][0] + substitution_matrix[alphabet.index(letter)][p]
        for column in range(1, len(seq1) + 1):
            score = calculate_score_data(row, column, substitution_matrix, scoring_matrix, seq1, seq2)
            scoring_matrix[1][column] = score[0]
        scoring_matrix[0,:] = scoring_matrix[1,:]

    last_line = scoring_matrix[1]
    return last_line


def NWLocalScore(alphabet, substitution_matrix, seq1, seq2):

    scoring_matrix = np.zeros((2, len(seq1) + 1))
    max_score_data = [0, 0, 0]  # max score, row, and column
    
    for row in range(1, len(seq2) + 1):
        for column in range(1, len(seq1) + 1):
            score = calculate_score_data(row, column, substitution_matrix, scoring_matrix, seq1, seq2)
            scoring_matrix[1][column] = score[0]
            if score[0] > max_score_data[0]:
                max_score_data[0] = score[0]
                max_score_data[1] = row
                max_score_data[2] = column
        scoring_matrix[0,:] = scoring_matrix[1,:]

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

    score = max(diagonal_score, up_score, left_score)

    # 8 = DIAGONAL, 2 = UP, 4 = LEFT
    if score == diagonal_score:
        score_origin = 8
    elif score == up_score:
        score_origin = 2
    else:
        score_origin = 4

    return (score, score_origin)
  

def needleman_wunsch(alphabet, substitution_matrix, seq1, seq2):

    p = len(alphabet)
    scoring_matrix = np.zeros((len(seq2) + 1, len(seq1) + 1))
    backtracking_matrix = np.zeros((len(seq2) + 1, len(seq1) + 1))

    for i in range(1, len(seq2) + 1):
        letter = seq2[i - 1]
        scoring_matrix[i][0] = substitution_matrix[alphabet.index(letter)][p] + scoring_matrix[i - 1][0]
        backtracking_matrix[i][0] = 2

    for j in range(1, len(seq1) + 1):
        letter = seq1[j - 1]
        scoring_matrix[0][j] = substitution_matrix[alphabet.index(letter)][p] + scoring_matrix[0][j - 1]
        backtracking_matrix[0][j] = 4

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

    while row != 0 or column != 0:
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
        else:
            seq1_alignment += seq1[column - 1]
            seq2_alignment += '-'
            column = column - 1

    seq1_alignment = seq1_alignment[::-1]
    seq2_alignment = seq2_alignment[::-1]
    return (seq1_alignment, seq2_alignment)

def get_indices(seq1_alignment, seq2_alignment, seq1_start_index, seq2_start_index):
    
    seq1_indices = []
    seq2_indices = []
    seq1_gaps = 0
    seq2_gaps = 0

    for letter in range(0, len(seq1_alignment)):
        if seq1_alignment[letter] == "-":
            seq1_gaps += 1
        elif seq2_alignment[letter] == "-":
            seq2_gaps += 1
        else:
            seq1_indices.append(letter + seq1_start_index - seq1_gaps)
            seq2_indices.append(letter + seq2_start_index - seq2_gaps)

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
substitution_matrix = [

[ 1,-5,-5,-5,-1],

[-5, 1,-5,-5,-1],

[-5,-5, 5,-5,-4],

[-5,-5,-5, 6,-4],

[-1,-1,-4,-4,-9]]

seq1 = "AAAAACCDDCCDDAAAAACC"
seq2 = "CCAAADDAAAACCAAADDCCAAAA"

alignments = dynproglin(alphabet, substitution_matrix, seq1, seq2)
print("Score:   ", alignments[0])
print("Indices: ", alignments[1], alignments[2])
