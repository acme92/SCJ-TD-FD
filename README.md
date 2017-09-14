### SCJTDFD

The DSCJ_smd.py is able to do the following:
1. Given a trivial ancestor A and its nontrivial descendant D, it can compute the SCJTDFD distance between the two.
2. Given a trivial ancestor A and its nontrivial descendant D, it can compute an instance of a scenario having the optimal SCJTDFD distance between the two.
3. Given a set of genomes, it can compute the median genome.


### Requirements

The code is composed of a set of Python scripts. It requires the following to be available:

* python (3)


### Usage

Use the following command to run the script:

python DSCJ_smd.py -s/-m/-d <inputfile>

Use -s for computing the optimal scenario
    -m for computing the median
    -d for computing the distance between two genomes

The input file should have the following format.
1. Each genome should start with the genome name followed by the genome itself from the next non-empty and non-commented line.
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
Another script to run the main script recursively for all files in the folder.
The output file.
