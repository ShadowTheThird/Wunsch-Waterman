# by Michal Walach
nucleotides = ['A','C','G','T']
while 1:
    valid = 1
    with open("Smith-Watermann\sequence.fasta", "r") as file:
        sequences = file.readlines()
    sequences = [sequence.strip() for sequence in sequences]
    sequences = [str.upper(sequence) for sequence in sequences]
    for letter in sequences:
        if letter not in nucleotides:
            print("sequence contains invalid symbols")
            valid = 0
            break
    if valid:
        print("loaded sequence is " + sequences)
        break
