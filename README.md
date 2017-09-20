## dSCJTDFD: the directed Single Cut or Join distance with Tandem Duplications and Floating Duplications
Pedro Feijao, Aniket Mane, Cedric Chauve  
School of Computing Science and Department of Mathematics, Simon Fraser University  
Contact: cedric.chauve@sfu.ca

### Overview

SCJTDFD is a program that implements genome rearrangement algorithms in the directed Single Cut or Join with Tandem Duplications and Floating Duplications (SCJTDFD) model. This model and the related algorithms are described in the paper "A tractable variant of the Single Cut or Join distance with duplicated genes" published in the proceedings of [RECOMB-CG 2017](http://www.crg.eu/en/event/15th-recomb-comparative-genomics-satellite-workshop), Lecture Notes in Boinformatics (LNBI), volume 10562, Springer (DOI: 10.1007/978-3-319-67979-2).

It is implements three algorithms:
* Computing the directed SCJTDFD distance between an ancestral genome A and a descendant genome D.
* Computing a parsimonious directed SCJTDFD scenario between an ancestral genome A and a descendant genome D.
* Computing a parsimonious directed median genome A of k descendant genomes D1, ..., Dk.

### Requirements

The code is composed of a set of Python scripts. It requires the following to be available:

* python (3)

### Usage

Use the following command to run the script:

python DSCJ.py -s/-m/-d <inputfile>

Use -d for computing the dSCJTDFD distance between the two input genomes
    -s for computing a parsimonious dSCJTDFD scenario between the two input genomes
    -m for computing a parsimonious dSCJTDFD median of k genomes
    
The input file should have the following format.
1. Each genome should start with the genome name followed by the genome itself in the next non-empty and non-commented line.
2. Each line that is read contains exactly one chromosome.
3. Linear chromosomes should end with '|' while circular chromosomes with ')'. 
4. Any gene in reverse orientation should be preceded by a '-' sign.
5. Lines containing the genome name should not end with '|' or ')'.
6. Any line that is a comment should start with a '#'. 
7. For -d and -s, the input file should contain exactly two genomes, both having the same set of gene families.
8. For -m, the input file should have at least two genomes.

Example:
Genome 1:
a -b c |
-d e f g )

Genome 2:
a -b c d )
e -f g |

### Files

Three main files:
The main python script DSCJ_smd.py.
A script to format the input as required.
The output file.
