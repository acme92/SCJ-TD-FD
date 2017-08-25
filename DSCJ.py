__author__ = "Aniket Mane"
__email__ = "amane@sfu.ca"

"""
DSCJ.py implements the algorithm mentioned in the paper "A tractable variant of the Single Cut or Join distance with duplicated genes".
The following code consists of three main functions:
1. 	Implementation: python DSCJ.py -d <inputfile>
	Given a genome A without duplicate genes and genome D having duplicate genes, 
	compute the SCJTDFD distance between the two.
2.	Implementation: python DSCJ.py -s <inputfile>
	Given a genome A without duplicate genes and genome D having duplicate genes, 
	compute one scenario that can transform A to D under the SCJTDFD model.
3.	Implementation: python DSCJ.py -m <inputfile>
	Given a set of genomes having duplicate genes,
	compute their median genome M having no duplicate genes.
"""

from sys import argv
import random
import networkx as nx
import matplotlib.pyplot as plt



#Using python 3 interpreter

#Function definitions
#---------------------------------------------------------------------------
#A genome is identified by its name and a list of chromosomes. Data type: List of lists.
#Return genome name.
def get_genome_name(genome):
	return genome[0]
#Return list of chromosomes in the genome.
def get_genome_chr_list(genome):
	return genome[1]

#A chromosome has a type (linear or circular) and an ordered set of genes. Data type: List.
#Return chromosome name.
def get_chr_name(chromosome):
	return chromosome[0]
#Return chromosome type.
def get_chr_type(chromosome):
	return chromosome[1]
#Return list of genes in the chromosome.	
def get_chr_gene_list(chromosome):
	return chromosome[2]

#Gene in reverse direction. If gene is a string '-g', returns string 'g'. 
def reverse(gene):
	return gene[1:] if gene[0] == '-' else str('-' + gene)

#Check if genome A is trivial. Maintains a list of unique gene families. If same gene family repeats, the genome in non-trivial.
def check_if_trivial(chromosome, seen, is_trivial):
	for gene in chromosome:
		if gene not in seen and reverse(gene) not in seen:
			seen.add(gene)
		else:
			is_trivial = False
		return (seen, is_trivial)	

#Forms a list of all gene families in input genome.
def get_gene_list(chr_list):
	gene_list = []
	for chromosome in chr_list:
		for gene in chromosome[2]:
			if gene[0] == '-':
				if reverse(gene) not in gene_list:
					gene_list.append(reverse(gene))
			else:
				if gene not in gene_list:
					gene_list.append(gene)
	return gene_list

#Output genome as a list of lists of genes.
def create_genome(string_list):
	genome_list = [[]*2 for i in range(2)]
	chr_list = [[]*3 for i in range(2)]
	seen = set()															#Set of genes in A to check for trivialness
	is_trivial = True
	gene_count = [0,0]														#Maintain a count of number of genes in A and D, respectively
	 
	i = -1
	for line in string_list:
		if line[-1] not in {')','|'}:
			i += 1
			genome_list[i].append(line)
			genome_list[i].append([])
		elif line[-1] == '|':
			line = line.split(' ')
			line = [x for x in line if x != '|']
			if i == 0:														#Check trivialness only for A.
				seen, trivial = check_if_trivial(line[1:], seen, is_trivial)
			genome_list[i][1].append(line[0])
			chr_list[i].append([line[0], 'L', line[1:]])
			gene_count[i] += len(line[1:])
		elif line[-1] == ')':
			line = line.split(' ')
			line[-1] = line[1]
			if i == 0:														#Check trivialness only for A.
				seen, trivial = check_if_trivial(line[1:-1], seen, is_trivial)
			genome_list[i][1].append(line[0])
			chr_list[i].append([line[0], 'C', line[1:]])
			gene_count[i] += len(line[1:-1])
	
	if trivial == False:													#If gene repeats, A is nontrivial. Terminate program.
		print("Error message: Ancestor genome must be trivial.")
		quit()
	if set(get_gene_list(chr_list[0])) != set(get_gene_list(chr_list[1])):		#If different set of gene families, terminate program.
		print("Error message: Ancestor and descendant genomes have different sets of gene families.")
		quit()		
	return (genome_list, chr_list, gene_count)	

#Counts TD from Arrays and reduces genome
def reduceGenome(chr_list):
	TD_from_arrays = 0														#Maintain a count of TD from arrays.
	for chromosome in chr_list:										
		for gene_idx in range(len(chromosome[2])):
			while gene_idx < len(chromosome[2]) - 1:
				if chromosome[2][gene_idx] == chromosome[2][gene_idx + 1]:		#Remove tandem arrays, if any.
					del chromosome[2][gene_idx]	
					TD_from_arrays += 1
				else:
					gene_idx += 1
	return (chr_list, TD_from_arrays)

#Forms adjacency list of input genome
def get_adj_list(chr_list):
	adj_list = []															#List of adjacencies where 
	for chromosome in chr_list:												#every adjacency is of the format:
		for gene_idx in range(len(chromosome[2]) - 1):							#[(g1,'h'/'t'),(g2,'h'/'t')]
			if chromosome[2][gene_idx][0] == '-':
				left = (chromosome[2][gene_idx][1:], 't')
			else:
				left = (chromosome[2][gene_idx], 'h')
			if chromosome[2][gene_idx + 1][0] == '-':
				right = (chromosome[2][gene_idx + 1][1:], 'h')
			else:
				right = (chromosome[2][gene_idx + 1], 't')
			adj_list.append([left, right])
	return adj_list


#Main functions
#---------------------------------------------------------------------------
#Finds distance between the two given genomes in the input file
def distance(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#']	#Read line only if it is nonempty and not a comment.

	genome_list, chr_list, gene_count = create_genome(string_list)

	A = chr_list[0]			
	D = chr_list[1]
	
	D, TD_from_arrays = reduceGenome(D)
	
	gene_count[1] -= TD_from_arrays											#Number of genes in D after removing tandem arrays.
	n_duplicates = gene_count[1] - gene_count[0]							#Number of genes in D - number of genes in A

	A_adj = get_adj_list(A)
	D_adj = get_adj_list(D)

	preserved_adj = [adj for adj in A_adj if adj in D_adj or list(reversed(adj)) in D_adj]		#Intersection of adjacency sets, A and D
	n_cuts = len(A_adj) - len(preserved_adj)								#Adjacencies seen in A but NOT preserved in D
	n_joins = len(D_adj) - len(preserved_adj)								#Adjacencies seen in D but NOT preserved from A

	d_DSCJ = n_cuts + n_joins + 2*n_duplicates + TD_from_arrays				#d_DSCJ(A,D) = |A-D| + |D-A| + 2*n_d + TDA.

	print(d_DSCJ)
	print(n_cuts)
	print(n_joins)
	print(n_duplicates)
	print(TD_from_arrays)


#Input
#---------------------------------------------------------------------------
if len(argv) < 3:															#Takes file with genomes as argument in command line		
	print('Usage: python DSCJ_smd.py -d/-s/-m <genome_file>')
	exit(1)

if argv[1] == '-d':
	distance(argv[2])
elif argv[1] == '-s':
	scenario(argv[2])
elif argv[1] == '-m':
	median(argv[2])
else:
	print('Incorrect usage')
	print('Usage: python DSCJ_smd.py -d/-s/-m <genome_file>')