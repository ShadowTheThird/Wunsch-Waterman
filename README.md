# Smith-Waterman-like algorithm
## The actual algorithm
This is a python program following the idea of Smith-Waterman algorithm for finding the local sequence alignment. 
It builds a score matrix based on the alignment of corresponding positions within the two sequences. 
Unlike global alignment algorithm it focuses only on high scoring regions and does not allow negative values within the matrix.
After building the matrix through backtracking from the region with the highest score it is possible to traceback the segment with highest similarity based on the provided gap and missmatch penalties.
## My implementation
My program allows for discovering the region of closest alignment between any two of the provided sequences and is heavily inspired by Smith-Waterman algorithm.
It takes the sequences provided and compares each pair until it finds a region of highest similarity where it begins to backtrack through the matrix to rebuild the most similar segment in the way it was aligned.
After such segment is found program outputs the discovered alignment together with its score and goes on to compare the next pair.
On top of outputting the results onto the terminal it also writes the same outputs to the file.
## Data Preperation
### Input formatting
My program is format sensitive and requires specific data formatting in order to work correctly.
Primarly a readable file with all of the sequences is required.
Such file has to have all the sequences preceded by its name written in '>{name}' format specifically.
If no '>' sign is provided the program will ignore the sequence or add it to a different sequence alltogether.
If '>' sign is not followed by the sequence name, the program will exit informing the user that such data formatting is not permitted.
The name must be input with no white signs before or within the name.
After the name is detected program will assume that whatever follows is the sequence for comparison, but it will ignore or white signs, whether that be a space of a new line symbol.
There are no requirements regarding the name outside of it being at least 1 sign long, whereas the sequence needs to be at least 5 symbols long. I honestly made that requirement just to play around with exceptions.
### Input/Output files
Once you have a dataset formatted to my programs requirements you will need to save it in a readable file.
There are no requirements regarding the files name or extension, but if you want the program to work without the need to supply it with any extra arguments on launch you would need to save your data into a 'sequences.fasta' file and locate the file in the same directory as the Smith-Waterman.py file.
The output file will be saved in the same directory as the Smith-Waterman.py and will be named 'output.txt' unless specified otherwise when launching the program.
## Making use of the program
If all data has been prepared accordingly you may run the program.
It was designed for use with Linux terminal or visual studio code python debugger, but does also run in cmd (althought the output is not very readable).
When launching the program you can specify several optional arguments: match/missmatch/gap penalties, input/output file location and name, as well as format.<br>
The penalties can be set with -m/--match -s/--missmatch and -g/--gap followed by an intiger number.
The default values are match: 3 missmatch: -3 gap:  -2.<br>
Input and output files can be set with -i/--input and -o/--output followed by the location and name of the file.<br>
Format is a setting that classifies the data as either text or fasta data.
If the fasta format is selected the program will run check if all sequences consists of only ACGT letters to avoid data errors.
## Use example
sequences.boi
```
>alfa
AGCTTAGGCTACGTAGCT
>beta
TCGGATCGTAGGCTTCAA
>omega
GCTAGGCTATAGCTGGTCGA
```
terminal command
```
python3 Smith-Waterman.py -i input.boi -o output.boi -m 5 -g -4
```
output.boi
```
loaded sequences are:
>ALFA	AGCTTAGGCTACGTAGCT
>BETA	TCGGATCGTAGGCTTCAA
the match score is 41:
TAGGCTACGTA-GCT
TCGGAT-CGTAGGCT
-----------------------------------------------------------------------
loaded sequences are:
>ALFA	AGCTTAGGCTACGTAGCT
>OMEGA	GCTAGGCTATAGCTGGTCGA
the match score is 58:
GCTTAGGCTACGTAGCT
GC-TAGGCTA--TAGCT
-----------------------------------------------------------------------
loaded sequences are:
>BETA	TCGGATCGTAGGCTTCAA
>OMEGA	GCTAGGCTATAGCTGGTCGA
the match score is 33:
CGTAGGCT-TCA
C-TAGGCTAT-A
-----------------------------------------------------------------------
```
