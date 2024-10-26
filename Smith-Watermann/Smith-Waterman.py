# by Michal Walach
import argparse

class InnerLoopBreak(Exception):
    pass
def format_type(text):
    text = text.lower()
    if text not in ["fasta","text"]:
        raise argparse.ArgumentError(f"invalid format '{text}'. Must be 'fasta' or 'text'.")
    return text.lower()

# argument and flag 
modifiers = argparse.ArgumentParser(description="Set values for match, mismatch and gap accordingly")
modifiers.add_argument("-m", "--match", type = int, default = 3, help = "Specify the value added on match")
modifiers.add_argument("-s", "--missmatch", type = int, default = -3, help = "Specify the value added on missmatch")
modifiers.add_argument("-g", "--gap", type = int, default = -2, help = "Specify the value added on gap")
modifiers.add_argument("-f", "--format", type = format_type, default = "text", help = "Defines whether the program should check if loaded sequences consist of only nucletid symbols")
args = modifiers.parse_args()
match = args.match
missmatch = args.missmatch
gap = args.gap
format = args.format
# check
# print(f"match: {args.match}")
# print(f"missmatch: {args.missmatch}")
# print(f"gap: {args.gap}")

nucleotides = ['A','C','G','T']
valid = 1
with open("Smith-Watermann\sequences.fasta", "r") as file:
    file_text = file.readlines()
file_text = [sequence.strip().upper() for sequence in file_text]
file_text = [sequence.split() for sequence in file_text]
current_sequence = -1
sequences = ["",""]
try:
    for line in file_text:
        for segment in line:
            if segment[0] == '>':
                current_sequence += 1
                if current_sequence > 1:
                    raise InnerLoopBreak
            else:
                sequences[current_sequence] += segment
except:
    pass
if format == "fasta":
    for sequence in sequences:
        for letter in sequence:
            if letter not in nucleotides:
                print("sequence contains invalid symbols")
                valid = 0
                break
    if not valid:
        print("loaded sequences are invalid!\nexiting program")
        exit()
print(f"loaded sequences are:\n>seq1\t{sequences[0]}\n>seq2\t{sequences[1]}")
matrix = [[0] * (len(sequences[1]) + 1) for _ in range(len(sequences[0]) + 1)]
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
# for letter in sequences[1]:
#     line += "\t" + letter
# print(line)
# line = ""
# for position in range(0,len(sequences[0])):
#     line += "\t" + str(matrix[0][position])
# print(line)
# for position in range(0,len(sequences[0])):
#     line = sequences[0][position]
#     for cell in matrix[position+1]:
#         line += "\t" + str(cell)
#     print(line)
current_score = highest_score; current_position = highest_score_position
best_match = ["",""]
while current_score > 0:
    if matrix[current_position[0]-1][current_position[1]-1] == current_score - match and sequences[0][current_position[0]-1] == sequences[1][current_position[1]-1] or matrix[current_position[0]-1][current_position[1]-1] == current_score - missmatch:
        current_position[0] -= 1; current_position[1] -= 1
        best_match[0] += sequences[0][current_position[0]]; best_match[1] += sequences[1][current_position[1]]
    else:
        if matrix[current_position[0]-1][current_position[1]] == current_score - gap:
            current_position[0] -= 1
            best_match[0] += sequences[0][current_position[0]]; best_match[1] += "-"
        else:
            current_position[1] -= 1
            best_match[0] += "-"; best_match[1] += sequences[1][current_position[1]]
    current_score = matrix[current_position[0]][current_position[1]]
print(f"the match score is {highest_score}:\n{best_match[0][::-1]}\n{best_match[1][::-1]}")
