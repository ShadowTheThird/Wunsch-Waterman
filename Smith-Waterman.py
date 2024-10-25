# by Michal Walach
nucleotides = ['A','C','G','T']
valid = 1
while 1:
    sequence1 = input("Input 1st sequence: ")
    print("you have input a sequence of: " + sequence1)
    for letter in sequence1:
        if letter not in nucleotides:
            print("sequence contains invalid symbols")
            valid = 0
            break
    if valid == 0:
        print("loaded sequence is " + sequence1)
        break
    