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
#A genome is identified by its name and a list of chromosomes. Data type: List of lists of strings.
#Set genome data (use better words)
def make_genome_list(genome_list, genome_name, i): #if i = 0 => A else => D
	genome_list[i].append(genome_name)		
	genome_list[i].append([])
	return genome_list

def list_genome_chr(genome_list, chr_name, i):
	genome_list[i][1].append(chr_name)
	return genome_list

def list_chr_data(chr_list, chromosome, chr_type, i):
	chr_list[i].append([chromosome[0], chr_type, chromosome[1:]])
	return chr_list

#Return genome name.
def get_genome_name(genome):
	return genome[0]
#Return list of chromosomes in the genome.
def get_genome_chr_list(genome):
	return genome[1]

#A chromosome has a type (linear or circular) and an ordered set of genes. Data type: List of strings.
#Return chromosome name.
def get_chr_name(chromosome):
	return chromosome[0]
#Return chromosome type.
def get_chr_type(chromosome):
	return chromosome[1]
#Return list of genes in the chromosome.	
def get_chr_gene_list(chromosome):
	return chromosome[2]
#Return gene at specific index in the chromosome.
def get_gene_by_posn(chromosome, j):
	if j > len(chromosome[2]):
		print("Gene position out of bounds.")
	else:
		return chromosome[2][j]

#A gene is characterised by its family, a unique name within its family and its orientation. Data type: String.
# Return the family of a gene g
def get_gene_posn(g):
    return(g[0])
# Test the orientation of a gene
def get_gene_orientation(g):
    return(g[1])
def get_left_neighbor(g):
	return(g[2])
def get_right_neighbor(g):
	return(g[3])
def get_left_neighbor_posn(g):
	return(g[4])
def get_right_neighbor_posn(g):
	return(g[5])

#Gene in reverse direction. If gene is a string '-g', returns string 'g'. 
def reverse(gene):
	if gene:
		return gene[1:] if gene[0] == '-' else str('-' + gene)
	else:
		return None

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
			genome_list = make_genome_list(genome_list, line, i)
		elif line[-1] == '|':
			line = line.split(' ')
			line = [x for x in line if x != '|']
			genome_list = list_genome_chr(genome_list, line[0], i)
			chr_list = list_chr_data(chr_list, line, 'L', i)
			gene_count[i] += len(line[1:])
			if i == 0:														#Check trivialness only for A.
				seen, trivial = check_if_trivial(line[1:], seen, is_trivial)			
		elif line[-1] == ')':
			line = line.split(' ')
			line[-1] = line[1]
			genome_list = list_genome_chr(genome_list, line[0], i)
			chr_list = list_chr_data(chr_list, line, 'C', i)
			gene_count[i] += len(line[1:-1])
			if i == 0:														#Check trivialness only for A.
				seen, trivial = check_if_trivial(line[1:-1], seen, is_trivial)
				
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


#---------------------------------------------------------------------------
def set_gene_family(gene):
	if gene[0] == '-':
		return gene[1:]
	else:
		return gene

def set_gene_orientation(gene):
	if gene[0] == '-':
		return 'BACKWARD'
	else:
		return 'FORWARD'

def set_left_neighbor(chromosome, j):
	if j == 0:
		if get_chr_type(chromosome) == 'L':
			return None
		elif get_chr_type(chromosome) == 'C':
			return get_gene_by_posn(chromosome, -1)
	else:
		return get_gene_by_posn(chromosome, j-1)

def set_right_neighbor(chromosome, j):
	if j == len(chromosome[2]) - 1:
		if get_chr_type(chromosome) == 'L':
			return None
		elif get_chr_type(chromosome) == 'C':
			return get_gene_by_posn(chromosome, 0)
	else:
		return get_gene_by_posn(chromosome, j+1)

def left_neighbor_posn(chromosome, posn):
	i, j = posn[0], posn[1]
	if j == 0:
		if get_chr_type(chromosome) == 'L':
			return None
		elif get_chr_type(chromosome) == 'C':
			return (i, len(chromosome[2]) - 1)
	else:
		return (i, j-1)

def right_neighbor_posn(chromosome, posn):
	i, j = posn[0], posn[1]
	if j == len(chromosome[2]) - 1:
		if get_chr_type(chromosome) == 'L':
			return None
		elif get_chr_type(chromosome) == 'C':
			return (i, 0)
	else:
		return (i, j+1)

#Updates A_dict by relabeling gene matched with gene in A
def updateA(gene, A_dict):
	LN, RN = get_left_neighbor(A_dict[gene]), get_right_neighbor(A_dict[gene])					#Left and right neighbors in A
	if LN:																	#If Left neighbor exists, update corresponding entry.
		if LN[0] == '-':													#(If LN exists, its right neighbor will be the gene itself, 
			A_dict[LN[1:]][3] = A_dict[LN[1:]][3]+str('copy')+str(1)		#so needs to be relabeled) 
		else: 
			A_dict[LN][3] = A_dict[LN][3]+str('copy')+str(1) 
	if RN:																	#If Right neighbor exists, update corresponding entry.
		if RN[0] == '-': 													#(If RN exists, its left neighbor will be the gene itself,	
			A_dict[RN[1:]][2] = A_dict[RN[1:]][2]+str('copy')+str(1)		#so needs to be relabeled)
		else: 
			A_dict[RN][2] = A_dict[RN][2]+str('copy')+str(1)
	A_dict[gene+str('copy')+str(1)] = A_dict[gene]
	del A_dict[gene]
	return A_dict

#Updates Idx_dict and D_dict by relabeling gene matched with gene in A
def updateD(gene, Idx_dict, A_dict, D_dict, posn, i):
	Idx_dict[gene+str('copy')+str(i)] = [posn]								#Update Idx_dict by introducing entry gcopy'i'
	Idx_dict[gene].remove(posn)											#and removing corresponding index from Idx_dict[g]

	LNIdx, RNIdx = get_left_neighbor_posn(D_dict[posn]), get_right_neighbor_posn(D_dict[posn])			#Left and right neighbors in D
	if LNIdx: D_dict[LNIdx][3] = D_dict[LNIdx][3]+str('copy')+str(i) 		#If LN exists, relabel its right neighbor
	if RNIdx: D_dict[RNIdx][2] = D_dict[RNIdx][2]+str('copy')+str(i) 		#If RN exists, relabel its left neighbor
		
	return Idx_dict, D_dict

#Updates list of Floating Duplicates
def updateFD(FD, i, gene, Idx_dict, D_dict):
	for posn in Idx_dict[gene]:
		FD.append(gene+str('copy')+str(i))									
		Idx_dict[gene+str('copy')+str(i)] = [posn]
		LNIdx, RNIdx = get_left_neighbor_posn(D_dict[posn]), get_right_neighbor_posn(D_dict[posn])			#Left and right neighbors in D
		if LNIdx: D_dict[LNIdx][3] = D_dict[LNIdx][3]+str('copy')+str(i)		#If LN exists, relabel its right neighbor
		if RNIdx: D_dict[RNIdx][2] = D_dict[RNIdx][2]+str('copy')+str(i)		#If RN exists, relabel its left neighbor
		i += 1
	return(FD, Idx_dict, D_dict)


#---------------------------------------------------------------------------
#Finds scenario with optimal distance to obtain second genome from first genome in the input file

#Pseudocode:

#Remove Tandem arrays if any. Complexity O(n_D)
#Create a dictionary for A. Complexity O(n_A) 
#Create an index dictionary for D: Key=gene name, Val=list of positions of gene in D
#Create a dictionary for D. Created along with index dictionary. Complexity O(n_D)
#Shuffle positions in A for randomness.
#
#For every gene g in A: Complexity O(n_A)
#	If g is from nontrivial family:
#		Check for context strongly conserved. Complexity O(copy number of g) 
#			If instance found: 
#           	Update all the dictionaries.
#		If instance found: 
#			Update list of FDs. Complexity O(copy number of g)
#		
#		If not strongly conserved:
#			Check for conservation of adjacencies on either side. Complexity O(copy number of g)
#			Check if context weakly conserved: Complexity O(copy number of g)
#				If instance found:
# 					Update all the dictionaries.
#				Else check only if left adj conserved:  
#					If so, update all the dictionaries accordingly.
#				Else check only if right adj conserved:  
#					If so, update all the dictionaries accordingly.
#				Else: 
#					Match with the first copy of g in D. Update all the dictionaries.
#
#			If weakly conserved:
#				Update list of FDs for remaining copies. Complexity O(copy number of g)	
#			Else not conserved:
#				Update list of FDs for remaining copies. Complexity O(copy number of g)
#
#Create an adjacency set of D using the dictionary. Complexity O(n_D)
#Create an adjacency set of A using the dictionary. Also add FDs. Complexity O(n_D)
#
#Find the distance: |A-D| + |D-A| + n_d + TDA

def scenario(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#']	#Read line only if it is nonempty and not a comment.

	genome_list, chr_list, gene_count = create_genome(string_list)

	A = chr_list[0]			
	D = chr_list[1]
	
	D, TD_from_arrays = reduceGenome(D)
	
	A_dict = {}

	for i in range(len(A)):		
		if get_chr_type(A[i]) == 'C':
			A[i][2] = A[i][2][:-1]
		chromosome = get_chr_gene_list(A[i])
		for j in range(len(chromosome)):
			A_dict[set_gene_family(chromosome[j])]= [	(i,j),
														set_gene_orientation(chromosome[j]),
														set_left_neighbor(A[i],j),
														set_right_neighbor(A[i],j)
													]
			
	#print(A_dict)


	Idx_dict = {}	#Dictionary for indices. Key = Gene family, Value = List of positions of g in D. 
	D_dict = {}		#Dictionary for D. Key = Index. Value = (Sign, Left neighbor, Left neighbor index, Right neighbor, Right neighbor index)
	
	for i in range(len(D)):
		if get_chr_type(D[i]) == 'C':
			D[i][2] = D[i][2][:-1]
		chromosome = get_chr_gene_list(D[i])
		for j in range(len(chromosome)):
			try:
				Idx_dict[set_gene_family(chromosome[j])].append((i,j))
			except KeyError:
				Idx_dict[set_gene_family(chromosome[j])] = [(i,j)]
			D_dict[(i,j)] = [	(i,j), 
								set_gene_orientation(chromosome[j]),
								set_left_neighbor(D[i],j),
								set_right_neighbor(D[i],j),
								left_neighbor_posn(D[i],(i,j)),
								right_neighbor_posn(D[i],(i,j))		
							]

	#print(D_dict)

	coords = []										#List of co-ordinates of A, shuffled for randomness							
	for i in range(len(A)):												
		for j in range(len(A[i][2])):
			coords.append([i,j])		
	random.shuffle(coords)

	FD = []											#List of FD created
	TD = []											#List of TD created

	for i, j in coords:
		gene = get_gene_by_posn(A[i],j)
		gene = set_gene_family(gene)

		if len(Idx_dict[gene]) > 1:			#Enter only if copy number of gene > 1
			#CASE 1: Context strongly conserved
			strong_context = 0
			weak_context = 0 
			for posn in Idx_dict[gene]:
		 		if strong_context == 0:
		 			if get_gene_orientation(A_dict[gene]) == get_gene_orientation(D_dict[posn]):
		 				if (get_left_neighbor(A_dict[gene]) and get_left_neighbor(A_dict[gene]) == get_left_neighbor(D_dict[posn]) and
		 						get_right_neighbor(A_dict[gene]) and get_right_neighbor(A_dict[gene]) == get_right_neighbor(D_dict[posn])):
		 					strong_context = 1
		 					A_dict = updateA(gene, A_dict)									#Updating A_dict
		 					Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, posn, 1)
		 			else:
		 				if (get_left_neighbor(A_dict[gene]) and get_left_neighbor(A_dict[gene]) == reverse(get_right_neighbor(D_dict[posn])) and
		 						get_right_neighbor(A_dict[gene]) and get_right_neighbor(A_dict[gene])==reverse(get_left_neighbor(D_dict[posn]))):
		 					strong_context = 1
		 					A_dict = updateA(gene, A_dict)
		 					Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, posn, 1)

			if strong_context == 1:
				FD, Idx_dict, D_dict = updateFD(FD, 2, gene, Idx_dict, D_dict)

			if strong_context == 0:
				left_adj_at, right_adj_at = None, None
				for posn in Idx_dict[gene]:
					if not left_adj_at:
						if get_gene_orientation(A_dict[gene]) == get_gene_orientation(D_dict[posn]):
							if (get_left_neighbor(A_dict[gene]) and get_left_neighbor(A_dict[gene]) == get_left_neighbor(D_dict[posn])):
								left_adj_at = posn
						else:
							if (get_left_neighbor(A_dict[gene]) and get_left_neighbor(A_dict[gene]) == reverse(get_right_neighbor(D_dict[posn]))):
								left_adj_at = posn
					if not right_adj_at:
						if get_gene_orientation(A_dict[gene]) == get_gene_orientation(D_dict[posn]):
							if (get_right_neighbor(A_dict[gene]) and get_right_neighbor(A_dict[gene]) == get_right_neighbor(D_dict[posn])):
								right_adj_at = posn
						else:
							if (get_right_neighbor(A_dict[gene]) and get_right_neighbor(A_dict[gene])==reverse(get_left_neighbor(D_dict[posn]))):
								right_adj_at = posn

				#Case 2: Weak context
				if left_adj_at and right_adj_at:
					weak_context = 1
					TD.append(gene)

					Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, left_adj_at, 1)
					Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, right_adj_at, 2)

					LN, RN = get_left_neighbor(A_dict[gene]), get_right_neighbor(A_dict[gene])
					Posn, Orientation = get_gene_posn(A_dict[gene]), get_gene_orientation(A_dict[gene])

					A_dict[gene+str('copy')+str(1)] = [Posn, Orientation, LN, RN]
					A_dict[gene+str('copy')+str(2)] = [Posn, Orientation, LN, RN]

					if LN[0] == '-':
						A_dict[LN[1:]][3] = A_dict[LN[1:]][3]+str('copy')+str(1)
						A_dict[gene+str('copy')+str(2)][2] = A_dict[LN[1:]][3]
					else:
						A_dict[LN][3] = A_dict[LN][3]+str('copy')+str(1)
						A_dict[gene+str('copy')+str(2)][2] = A_dict[LN][3]
					if RN[0] == '-':
						A_dict[RN[1:]][2] = A_dict[RN[1:]][2]+str('copy')+str(2)
						A_dict[gene+str('copy')+str(1)][3] = A_dict[RN[1:]][2]
					else:
						A_dict[RN][2] = A_dict[RN][2]+str('copy')+str(2)
						A_dict[gene+str('copy')+str(1)][3] = A_dict[RN][2]
					del A_dict[gene]

				#Case 3:
				elif left_adj_at and not right_adj_at:
					Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, left_adj_at, 1)		#Updating D by relabeling left match (to original gene in A)						
					A_dict = updateA(gene, A_dict)													#Updating A_dict

				elif right_adj_at and not left_adj_at:
					Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, right_adj_at, 1)		#Updating D by relabeling right match (to original gene in A)						
					A_dict = updateA(gene, A_dict)													#Updating A_dict

				#CASE 4: Context not conserved. No adjacency conserved.			
				else:
					Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, posn, 1)	#Updating D by relabeling last copy (to original gene in A)						
					A_dict = updateA(gene, A_dict)											#Updating A_dict

				if left_adj_at and right_adj_at:									#If weakly conserved context, remaining genes matched with FDs
					FD, Idx_dict, D_dict = updateFD(FD, 3, gene, Idx_dict, D_dict)									
				else:																#If context not conserved, remaining genes matched with FDs
					FD, Idx_dict, D_dict = updateFD(FD, 2, gene, Idx_dict, D_dict)								

	#All Cases covered. Create adjacency lists from dictionaries

	#Logic:
	#Iterate through sorted dictionary.	
	#If right neighbor (RN) exists, adjacency is [left neighbor of RN, RN] (head or tail chosen according to orientation)
	D_adj = []													#Form adjacency set for relabeled genome D'							
	for x in sorted((k,v) for (k,v) in D_dict.items()):
		left, right = None, None 								#All adjacencies of the format: [(g1,'h'/'t'),(g2,'h'/'t')]
		if x[1][5]:
			RNIdx = x[1][5]
			if D_dict[RNIdx][2][0] == '-':					#Left extremity of adjacency
				left = (D_dict[RNIdx][2][1:], 't')
			else:
				left = (D_dict[RNIdx][2], 'h')
			if D_dict[x[0]][3][0] == '-': 					#Right extremity of adjacency
				right = (D_dict[x[0]][3][1:], 'h')
			else:
				right = (D_dict[x[0]][3], 't')
			D_adj.append([left, right])

	#print(D_adj)
	A_adj = []													#Form adjacency set for relabeled genome A'
	for x in sorted((v[0],k) for (k,v) in A_dict.items()):
		left, right = None, None 								#All adjacencies of the format: [(g1,'h'/'t'),(g2,'h'/'t')]	
		if A_dict[x[1]][3]:
			RN = A_dict[x[1]][3]
			if A_dict[x[1]][1] == 'BACKWARD':					#Left extremity of adjacency
				left = (x[1], 't')
			else:
				left = (x[1], 'h')
			if RN[0] == '-':									#Right extremity of adjacency
				right = (RN[1:], 'h')
			else:
				right = (RN, 't')
			A_adj.append([left, right])			
	for x in FD:												#Append (g_h g_t) adjacencies for each FD
		A_adj.append([(x, 'h'),(x, 't')])

	#print(A_adj)	

	preserved_adj = [adj for adj in A_adj if adj in D_adj or list(reversed(adj)) in D_adj]		#Intersection of adjacency sets, A' and D'
	n_cuts = len(A_adj) - len(preserved_adj)					#Adjacencies seen in A' but NOT preserved in D'
	n_joins = len(D_adj) - len(preserved_adj)					#Adjacencies seen in D' but NOT preserved from A'
	n_duplicates = len(FD) + len(TD)

	distance = n_cuts + n_joins + n_duplicates + TD_from_arrays	#d_DSCJ(A,D) = |A'-D'| + |D'-A'| + n_d + TDA.


	print(distance)
	print(n_cuts, n_joins)			
	print(len(FD), len(TD), TD_from_arrays)




#Adjacency weight function
def wtAdj(adj, adj_list):
	weight = 0
	for genome in adj_list:
		if adj in genome or adj[::-1] in genome:	
			weight += 1
	weight = 2*weight - len(adj_list)
	return weight

#Create MWM graph
def createGraph(adj_list, total_gene_list, total_adj_list):
	G = nx.Graph()
	for g in total_gene_list:
		G.add_node((g,'t'))
		G.add_node((g,'h'))
	edge_list = []
	for adj in total_adj_list:
		edge_list.append((adj[0],adj[1],wtAdj(adj, adj_list)))
	G.add_weighted_edges_from(edge_list)
	return G

#Max weight matching
def MWMedges(G, total_adj_list):
	kept_adj = []
	disc_adj = []
	adj_kept_mwm = {}
	adj_info = {}
	for adj in total_adj_list:
		adj_kept_mwm[tuple(adj)] = False
	M = nx.max_weight_matching(G)
	M_list = sorted(list(M.keys()))
	for m1 in M_list:
		m2 = M[m1]
		adj_kept_mwm[(m1,m2)] = True
	for adj in total_adj_list:
		if adj_kept_mwm[tuple(adj)] == True:
			kept_adj.append(tuple(adj))
		else:
			disc_adj.append(tuple(adj))
	return((kept_adj, disc_adj))


#---------------------------------------------------------------------------
#Finds the median of all given genomes in the input file
'''
def median(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#']
	genome_list = []
	i = -1
	for line in string_list:
		if line[-1] not in {')','|'}:
			i += 1
			genome_list = make_genome_list(genome_list, line, i)
		elif line[-1] == '|':
			line = line.split(' ')
			line = [x for x in line if x != '|']
			genome_list = list_genome_chr(genome_list, line[0], i)
		else:
			line = line.split(' ')
			line[-1] = line[0]
			genome_list = list_genome_chr(genome_list, line[0], i)

	TD_from_arrays = [0] * len(genomes)
	i = 0
	for genome in genome_list:
		for chromosome in genome:				
			for gene_idx in range(len(chromosome[2])):
				while gene_idx < len(chromosome[2]) - 1:
					if chromosome[2][gene_idx] == chromosome[2][gene_idx + 1]:
						del chromosome[2][gene_idx]
						TD_from_arrays[i] += 1
						print(TD_from_arrays[i])
					else:
						gene_idx += 1
		i += 1

	adj_list = []
	total_gene_list = []	
	for genome in genome_list:
		adj_list.append(get_adj_list(genome))
		total_gene_list = list(set(total_gene_list + get_gene_list(genome)))

	total_adj_list = []
	for genome in adj_list:
		for adj in genome:
			if adj in total_adj_list or adj[::-1] in total_adj_list:
				continue
			else:
				total_adj_list.append(adj)	

	G = createGraph(adj_list, total_gene_list, total_adj_list)
	M = nx.max_weight_matching(G)
	print(M)
	print(MWMedges(G, total_adj_list))
	
	pos=nx.spring_layout(G)
	nx.draw(G,pos,with_labels=True)
	nx.draw_networkx_edge_labels(G,pos)
	plt.axis('off')
	plt.show()
'''





					
						
					




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