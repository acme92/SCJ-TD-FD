# SCJD

SCJD takes two genomes as input - a non-duplicated genome S and one with duplicates D, and returns the optimal distance between S and D. 


### Requirements

PhySca is composed of a set of Python scripts. It requires the following to be available:

* python (3)


### Usage

Currently, the genes are stored in a text file called genome.txt. We can change the input method later as and when required.
Running the following command gives the distance between the two input genomes.

python d_SCJDup.py genome.txt

### Files

If we plan to make the modular:
1. The file [create_pr.py] creates the S' and D' genomes. The input for this is from the file genome.txt
2. The genomes S' and D' (list of lists) then act as the input for [sym_diff.py]. 
