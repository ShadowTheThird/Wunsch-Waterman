# by Michal Walach
import argparse

# argument and flag 
modifiers = argparse.ArgumentParser(description="Set values for match, mismatch and gap accordingly")
modifiers.add_argument("-m", "--match", type = int, default = 3, help = "Specify the value added on match")
modifiers.add_argument("-s", "--missmatch", type = int, default = -3, help = "Specify the value added on missmatch")
modifiers.add_argument("-g", "--gap", type = int, default = -2, help = "Specify the value added on gap")
args = modifiers.parse_args()
match = args.match
missmatch = args.missmatch
gap = args.gap
# check
# print(f"match: {args.match}")
# print(f"missmatch: {args.missmatch}")
# print(f"gap: {args.gap}")

nucleotides = ['A','C','G','T']
while 1:
    valid = 1
    with open("Smith-Watermann\sequence.fasta", "r") as file:
        sequences = file.readlines()
    sequences = [sequence.strip() for sequence in sequences]
    sequences = [sequence.upper() for sequence in sequences]
    for sequence in sequences:
        for letter in sequence:
            if letter not in nucleotides:
                print("sequence contains invalid symbols")
                valid = 0
                break
    if valid:
        print(f"loaded sequences are:\n>SW1\t{sequences[0]}\n>SW2\t{sequences[1]}")
        break
matrix = [[0] * (len(sequences[1]) + 1) for _ in range(len(sequences[0]) + 1)]
# modification = missmatch
# for row in range(1,len(matrix)):
#     for cell in range(1,len(matrix[row])):
#         if sequences[0][row-1] == sequences[1][cell-1]:
#             modification = match
#         else:
#             modification = missmatch
#         if matrix[row-1][cell] + gap > matrix[row][cell-1] + gap:
#             if matrix[row-1][cell] + gap > matrix[row-1][cell-1] + modification:
#                 matrix[row][cell] = matrix[row-1][cell] + gap
#             else:
#                 matrix[row][cell] = matrix[row-1][cell-1] + modification
#         else:
#             if matrix[row][cell-1] + gap > matrix[row-1][cell-1] + modification:
#                 matrix[row][cell] = matrix[row][cell-1] + gap
#             else:
#                 matrix[row][cell] = matrix[row-1][cell-1] + modification
highest_score = 0
for row in range(1,len(matrix)):
    for cell in range(1,len(matrix[row])):
        if sequences[0][row-1] == sequences[1][cell-1]:
            if max(matrix[row-1][cell]+gap,matrix[row][cell-1]+gap) > matrix[row-1][cell-1]+match:
                matrix[row][cell] = max(matrix[row-1][cell]+gap, matrix[row][cell-1]+gap)
            else:
                matrix[row][cell] = matrix[row-1][cell-1]+match
        else:
            if max(matrix[row-1][cell]+gap,matrix[row][cell-1]+gap) > matrix[row-1][cell-1]+missmatch:
                matrix[row][cell] = max(matrix[row-1][cell]+gap, matrix[row][cell-1]+gap, 0)
            else:
                matrix[row][cell] = max(matrix[row-1][cell-1]+missmatch, 0)
        if matrix[row][cell] > highest_score:
            highest_score = matrix[row][cell]
            highest_score_position = [row, cell]
line = "X\t"
for letter in sequences[1]:
    line += "\t" + letter
print(line)
line = ""
for position in range(0,len(sequences[0])):
    line += "\t" + str(matrix[0][position])
print(line)
for position in range(0,len(sequences[0])):
    line = sequences[0][position]
    for cell in matrix[position+1]:
        line += "\t" + str(cell)
    print(line)
print(f"the match score is {highest_score}")
