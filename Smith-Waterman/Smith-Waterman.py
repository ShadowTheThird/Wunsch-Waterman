# by Michal Walach
import argparse

# setting variables for estethic output
RED =   "\033[91m"
RESET = "\033[0m"
BOLD = "\033[1m"
DYELLOW = "\033[38;5;214m"
BLUE = "\033[38;5;4m"
GREEN = "\033[38;5;2m"
GREY = "\033[38;5;8m"

# case insensitive flag value collector for data format
def format_type(text):
    text = text.lower()
    if text not in ["fasta","text"]:
        raise argparse.ArgumentError(f"invalid format '{text}'. Must be 'fasta' or 'text'.")
    return text.lower()

# acceptable flags and their default values
modifiers = argparse.ArgumentParser(description="Set values for match, mismatch and gap accordingly")
modifiers.add_argument("-m", "--match", type = int, default = 3, help = "Specify the value added on match")
modifiers.add_argument("-s", "--missmatch", type = int, default = -3, help = "Specify the value added on missmatch")
modifiers.add_argument("-g", "--gap", type = int, default = -2, help = "Specify the value added on gap")
modifiers.add_argument("-f", "--format", type = format_type, default = "text", help = "Defines whether the program should check if loaded sequences consist of only nucletid symbols")
modifiers.add_argument("-o", "--output", type = str, default = "Smith-Waterman\output.txt", help = "Sets the name of the output file and its directory")
modifiers.add_argument("-i", "--input", type = str, default = "Smith-Waterman\sequences.fasta", help = "Sets the name of the input file and its directory")
args = modifiers.parse_args()
match = args.match
missmatch = args.missmatch
gap = args.gap
format = args.format
output_file = args.output
input_file = args.input

try:
    # data collector from selected file
    with open(input_file, "r") as file:
        file_text = file.readlines()
    # formatting text to acquire usable data
    file_text = [sequence.strip().upper() for sequence in file_text]
    file_text = [sequence.split() for sequence in file_text]
    # spliting aquired data into sequences
    current_sequence = -1
    sequences = ["",""]
    names = ["",""]
    for line in file_text:
        for segment in line:
            if segment[0] == '>':
                current_sequence += 1
                if current_sequence == len(sequences):
                    sequences.append("")
                    names.append("")
                names[current_sequence] += segment[1:]
            else:
                sequences[current_sequence] += segment
    if current_sequence == -1:
        raise ValueError("No sequences detected. Ensure all sequences are preceded by '>' with a non-empty name.")
    for name in names:
        if name == "":
            raise ValueError("One or more names are empty. Ensure all sequences are preceded by '>' with a non-empty name.")
    for sequence in sequences:
        if len(sequence) < 5:
            raise ValueError("One or more sequences are too short. Ensure all sequences are at least 5 symbols long.")

    # check for fasta format if sequences are composed entirely of ACGT letters
    if format == "fasta":
        nucleotides = ['A','C','G','T']
        valid = 1
        for sequence in sequences:
            for letter in sequence:
                if letter not in nucleotides:
                    print("sequence contains invalid symbols")
                    valid = 0
                    break
        # exit case if sequences had invalid symbols
        if not valid:
            raise ValueError("Sequence contains invalid symbols. Only A, C, G, and T are allowed for fasta format.")

except FileNotFoundError:
    print(f"{RED}{BOLD}ERROR:{RESET}\t{DYELLOW}The file '{input_file}' was not found. Please check the file path.{RESET}")
    exit()
except ValueError as ve:
    print(f"{RED}{BOLD}ERROR:{RESET}\t{DYELLOW}{ve}{RESET}")
    exit()
except Exception as e:
    print(f"{RED}{BOLD}ERROR:{RESET}\t{DYELLOW}an uncaught error occured: {e}{RESET}")
    exit()

# selection of sequences to be compared
with open(output_file, 'w') as ofile:
    for seq1 in range(0,len(sequences)):
        for seq2 in range(seq1+1, len(sequences)):
            comparables = [sequences[seq1], sequences[seq2]]
            print(f"loaded sequences are:\n{DYELLOW}>{names[seq1]}\t{BLUE}{comparables[0]}\n{DYELLOW}>{names[seq2]}\t{BLUE}{comparables[1]}{RESET}")
            ofile.write(f"loaded sequences are:\n>{names[seq1]}\t{comparables[0]}\n>{names[seq2]}\t{comparables[1]}\n")
            # creation of matrix
            matrix = [[0] * (len(comparables[1]) + 1) for _ in range(len(comparables[0]) + 1)]
            highest_score = 0
            for row in range(1,len(matrix)):
                for cell in range(1,len(matrix[row])):
                    if comparables[0][row-1] == comparables[1][cell-1]:
                        if max(matrix[row-1][cell]+gap,matrix[row][cell-1]+gap) > matrix[row-1][cell-1]+match:
                            matrix[row][cell] = max(matrix[row-1][cell]+gap, matrix[row][cell-1]+gap)
                        else:
                            matrix[row][cell] = matrix[row-1][cell-1]+match
                    else:
                        if max(matrix[row-1][cell]+gap,matrix[row][cell-1]+gap) > matrix[row-1][cell-1]+missmatch:
                            matrix[row][cell] = max(matrix[row-1][cell]+gap, matrix[row][cell-1]+gap, 0)
                        else:
                            matrix[row][cell] = max(matrix[row-1][cell-1]+missmatch, 0)
                    # saving highest score for code optimization
                    if matrix[row][cell] > highest_score:
                        highest_score = matrix[row][cell]
                        highest_score_position = [row, cell]
            # searching for the match within the generated matrix using the highest recorded score as starting point
            if highest_score == 0:
                print(f"no match was found for selected sequences\n{GREY}-----------------------------------------------------------------------{RESET}")
                continue
            current_score = highest_score; current_position = highest_score_position
            best_match = ["",""]
            while current_score > 0:
                if matrix[current_position[0]-1][current_position[1]-1] == current_score - match and comparables[0][current_position[0]-1] == comparables[1][current_position[1]-1] or matrix[current_position[0]-1][current_position[1]-1] == current_score - missmatch:
                    current_position[0] -= 1; current_position[1] -= 1
                    best_match[0] += comparables[0][current_position[0]]; best_match[1] += comparables[1][current_position[1]]
                else:
                    if matrix[current_position[0]-1][current_position[1]] == current_score - gap:
                        current_position[0] -= 1
                        best_match[0] += comparables[0][current_position[0]]; best_match[1] += "-"
                    else:
                        current_position[1] -= 1
                        best_match[0] += "-"; best_match[1] += comparables[1][current_position[1]]
                current_score = matrix[current_position[0]][current_position[1]]
            # printing the found sequences in reverse order as they have been searched for and recorded from back to front
            print(f"the match score is {highest_score}:\n{GREEN}{best_match[0][::-1]}\n{best_match[1][::-1]}\n{GREY}-----------------------------------------------------------------------{RESET}")
            ofile.write(f"the match score is {highest_score}:\n{best_match[0][::-1]}\n{best_match[1][::-1]}\n-----------------------------------------------------------------------\n")