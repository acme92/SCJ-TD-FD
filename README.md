# SCJD

SCJD takes two genomes as input - a trivial genome A and non-trivial D, and returns the optimal distance between A and D. 


### Requirements

PhySca is composed of a set of Python scripts. It requires the following to be available:

* python (3)


### Usage

The input files should have the ancestor genome on line 1 and descendant on line 2. All the input files should be in a folder. Currently, the foldername is provided as the input on the command line. We can change the input method later as and when required.
Running the following command gives the distance between the two input genomes.

python SCJTDFD.py <foldername>

### Files

If we plan to make the modular:
1. The file [create_pr.py] creates the S' and D' genomes. The input for this is from the file genome.txt
2. The genomes S' and D' (list of lists) then act as the input for [sym_diff.py]. 
