# SCJD

SCJD takes two genomes as input - a non-duplicated genome A and one with duplicates D, and returns the optimal distance between A and D. 


### Requirements

SCJTDFD requires the following to be available:

* python (3)


### Usage

Currently, the genes are stored in a text file. We can change the input method later as and when required.
Running the following command gives the distance between the two input genomes.

python SCJTDFD.py genomefilename

### Files

If we plan to make the modular:
1. The genomes A'' and D'' (list of lists) act as the input for [sym_diff.py]. 
