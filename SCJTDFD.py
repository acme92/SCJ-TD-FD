from sys import argv
import random
from functools import reduce

#Using python 3 interpreter

#Function definitions
#---------------------------------------------------------------------------
#Reverses an individual sign, a gene or a gene segment
def negate(a):	
	assert(type(a)==str or type(a)==list)	
	if type(a) == str:
		if len(a)==1:	
			#negates an individual sign (+ or -)				
			return '-' if a == '+' else '+'
		else:
			#reverses a gene (eg. +a to -a, -b to +b)
			return negate(a[0]) + a[1:]
	elif type(a) == list:
		#reverses a gene segment (eg. [+a,-b,+c] to [-c,+b,-a])
		a_negative = ["" for i in range(len(a))]
		for idx in range(len(a)):
			a_negative[idx] = negate(a[idx])
		return a_negative[::-1]

#Locates chromosome and gene index in target genome
def locateGene(inputgene, genome):
	for chr_idx, chromosome in enumerate(genome):
		for gene_idx, gene in enumerate (chromosome):
			if gene == inputgene or gene == negate(inputgene):
				return(chr_idx, gene_idx)

#Obtain triplets to determine if context is strongly conserved
def getTriplets(chromosome, gene_idx):
	triplet1 = []
	triplet2 = []
	if gene_idx > 0 and gene_idx < len(chromosome)-1:
		triplet1 = [chromosome[gene_idx-1], chromosome[gene_idx], chromosome[gene_idx+1]]
		triplet2 = negate(triplet1)		#Triplet1 in opposite direction
	return triplet1, triplet2

#Converts list to string
def charListToString(charList):
	return reduce(lambda x, y: x + y, charList, '')

#Outputs adjacencies of a gene in a given chromosome
def getAdj(chromosome, gene_idx):
	adj1 = []
	adj2 = []
	if gene_idx > 0:
		adj1 = [chromosome[gene_idx-1], chromosome[gene_idx]]
	if gene_idx < len(chromosome)-1:
		adj2 = [chromosome[gene_idx], chromosome[gene_idx+1]]
	return adj1, adj2

#Compares adjacencies
def compareAdj(chromosome1, chromosome2):
	assert(type(chromosome1)==list and type(chromosome2==list))

	if charListToString(chromosome1) == charListToString(chromosome2):
		return True
	elif charListToString(negate(chromosome1)) == charListToString(chromosome2):
		return True
	else:
		return False

#Functions used while calculating symmetric difference
#Pairs up adjacent genes in the list of genes that forms the chromosome
def pairUp(a):
    set_a = set()
    for chrom in a:
        for i in range(len(chrom) - 1):
            set_a.add(chrom[i] + chrom[i+1])
    return set_a

#Adds the adjacencies in both directions to create a look-up dictionary
def duplicateOrder(A):
    A_double = []
    for i in A:
        A_double.append(i)
        # complement
        count = 0
        index = 0
        for idx, val in enumerate(i):
            if val == '+' or val == '-':
                count += 1
                if count == 2:
                    index = idx                    
        A_double.append(
            negate(i[index]) + i[index+1:] +
            negate(i[0]) + i[1:index]
        )
    return A_double

#Elements of set A that aren't in the dictionary of set B
def diff(A, B):
    B_dict = {i:None for i in duplicateOrder(B)}
    return [i for i in A if i not in B_dict]

#Symmetric difference
def sym_diff(A, B):
   	return len(diff(A,B)) + len(diff(B,A))


#Input and preprocessing
#---------------------------------------------------------------------------

if len(argv) < 2:															#Takes file with genomes as argument in command line
	print('Usage: python d_SCjDup.py <genome_file>')
	exit(1)

genome_file = argv[1]
string = open(genome_file, "r").read()		

genomes = string.split("\n")												#Genomes A and D split by "\n" and obtains as a list of chromosomes
genomes = [i.strip().split("|") for i in genomes if len(i)]					#Each chromosome in a genome split by "|"
A = [j.strip().split(" ") for j in genomes[0]]								#Each chromosome is obtained as a list of genes
D = [j.strip().split(" ") for j in genomes[1]]								#Currently all chromosomes treated linear

TD = 0																		#No. of tandem duplications
FD = 0																		#No. of floating duplications
TD_from_arrays = 0 															#This is NOT included in number of TDs

D_pr = [ [i for i in chromosome] for chromosome in D]						#Initializes D' to D and removes tandem arrays in D'
for chr_idx in range(len(D_pr)):
	gene_idx = 0
	while gene_idx < len(D_pr[chr_idx]) - 1:
		if D_pr[chr_idx][gene_idx] == D_pr[chr_idx][gene_idx + 1]:
			del D_pr[chr_idx][gene_idx]
			TD_from_arrays += 1
		else:
			gene_idx = gene_idx + 1

#Main program
#---------------------------------------------------------------------------

A_2 = [ [i for i in chromosome] for chromosome in A]						#Initializes A_2 and D_2 to A and D_pr respectively
D_2 = [ [i for i in chromosome] for chromosome in D_pr]

coords = []																	#Constructs a list of co-ordinates of genes in A
for chr_idx in range(len(A)):												#Shuffles the list for randomness
	for gene_idx in range(len(A[chr_idx])):
		coords.append([chr_idx,gene_idx])		
random.shuffle(coords)
print(coords)

for A_chr_idx, A_gene_idx in coords:										#Iterates randomly through genes in A
	gene = A[A_chr_idx][A_gene_idx]
																			#Find location of gene in A_2 (chromosome, gene)
	A2_chr_idx, A2_gene_idx = locateGene(gene, A_2)							#since A_2 keeps changing and gene may not have the same location as in A																			
	
	search = [None]*2														#Assigns elements - gene and its reverse to be searched in D_2
	search[0] = gene 														#For e.g.: If gene is '-a', search includes '-a' and '+a'
	search[1] = negate(gene[0]) + gene[1:]

	found_indices = []														#Finds occurrences of the gene in D_2
	for D_chr_idx in range(len(D_2)):										#Lists their indices - (chromosome, gene) in found_indices	
		D_chromosome = D_2[D_chr_idx]
		for D_gene_idx in range(len(D_chromosome)):
			D_gene = D_chromosome[D_gene_idx]
			if search[0] == D_gene or search[1] == D_gene:
				found_indices.append( (D_chr_idx, D_gene_idx) )
				
	count = len(found_indices)												#Finds copy number of gene
	
	if count > 1:															#Entered only if duplicates exist

		#Step 1: Check for strongly conserved context-----------------------
		
		A_idx_count = {}													#Keep a count of number of occurrences of the gene in A_2 and D_2
		D_idx_count = {}													#Required for incremental labeling of duplicates

		triplet1, triplet2 = getTriplets(A_2[A2_chr_idx], A2_gene_idx)		#Triplet -> List of [left neighbour, gene, right neighbour]
																			#Basically, to check for the context of gene in a genome

		is_triplet_matched = False
		unmatched_indices = []												#List of occurrences of gene in D_2 that are not yet matched

		if len(triplet1) and len(triplet2):									#If triplet found in A_2 (gene is NOT a telomere in A_2)
			for D_chr_idx, D_gene_idx in found_indices:
				D_triplet, _ = getTriplets(D_2[D_chr_idx], D_gene_idx)		#Gets triplet for each occurrence of the gene in D_2

				if len(D_triplet):											#If triplet found in D_2 (gene is NOT a telomere in D_2)
					triplet1 = charListToString(triplet1)					
					triplet2 = charListToString(triplet2)
					D_triplet = charListToString(D_triplet)
					A_gene = charListToString(gene)
																			#If a triplet from A_2 has already been matched to one in D_2	
					if is_triplet_matched == False:							#then it shouldn't look for more instances of a triplet
																			
						if triplet1 == D_triplet:							#Checks if triplet is matched to one in D_2
							is_triplet_matched = True
							A_idx_count[A_gene[1:]] = 1															#If match found, increments label
							A_2[A2_chr_idx][A2_gene_idx] = A_gene + str('copy') + str(A_idx_count[A_gene[1:]])	#Updates labeling of corresponding					
							D_2[D_chr_idx][D_gene_idx] += str('copy') + str(A_idx_count[A_gene[1:]])			#genes in A_2 and D_2

						elif triplet2 == D_triplet:							#Checks if reverse of the triplet is matched to one in D_2
							#print('Found triplet: ', triplet2)
							is_triplet_matched = True
							A_idx_count[A_gene[1:]] = 1															#If match found, increments label
							A_2[A2_chr_idx][A2_gene_idx] = A_gene + str('copy') + str(A_idx_count[A_gene[1:]])	#Updates labeling of corresponding						
							D_2[D_chr_idx][D_gene_idx] += str('copy') + str(A_idx_count[A_gene[1:]])			#genes in A_2 and D_2

						else:													#Triplets exist for gene in both A_2 and D_2 but don't match
							unmatched_indices.append((D_chr_idx, D_gene_idx))	#Implies that context NOT strongly conserved and corresponding genes NOT matched

					else:														#Strongly conserved context already found and no more matches required
						unmatched_indices.append((D_chr_idx, D_gene_idx))		#Gene in current iteration can NOT be matched to gene from A_2

				else:															#Triplet doesn't exist for gene in D_2 as gene is a telomere
					unmatched_indices.append((D_chr_idx, D_gene_idx))			#Gene in current iteration can NOT be matched to gene from A_2

		else:																	#Triplets don't exist for gene in A_2 as gene is a telomere				
			unmatched_indices = [i for i in found_indices]						#No gene in D_2 can be matched as a strongly conserved context
																				#This case will handled later as only one adjacency conserved technically
 
		if is_triplet_matched == True:														#If triplets are matched, context strongly conserved
			for unmatched_D_chr_idx, unmatched_D_gene_idx in unmatched_indices:				#All remaining copies in D_2 will be matched to a floating duplicate
				A_idx_count[A_gene[1:]] += 1												#Increments label, creates a floating duplicate and updates label in D_2
				A_2.append( [A_gene + str('copy') + str(A_idx_count[A_gene[1:]]), A_gene + str('copy') + str(A_idx_count[A_gene[1:]])] )
				D_2[unmatched_D_chr_idx][unmatched_D_gene_idx] += str('copy') + str(A_idx_count[A_gene[1:]])
				FD += 1																		

		#------------------------End of Step 1------------------------------	

		#Step 2: Check for weakly conserved context or adjacency conserved on one side

		if is_triplet_matched == False:										#Enters this block only if context is NOT strongly conserved
			adj_found_at = [None, None]										#Index 0: left adjacency, index 1: right adjacency

			for D_chr_idx, D_gene_idx in unmatched_indices:					
				gene_adj_in_D2 = getAdj(D_2[D_chr_idx], D_gene_idx)			#Gets adjacencies of gene in D_2
				gene_adj_in_A2 = getAdj(A_2[A2_chr_idx], A2_gene_idx)		#As well as in A_2

				for A_adj_side in range(len(gene_adj_in_A2)):				#gene_adj_in_*2[0] => left adj, gene_adj_in_*2[0] => right adj
					A_adj = gene_adj_in_A2[A_adj_side]

					if len(A_adj):											#If adjacency exists in A_2
						for D_adj_side in range(len(gene_adj_in_D2)):
							D_adj = gene_adj_in_D2[D_adj_side]

							if len(D_adj):									#If adjacency exists in D_2
								if compareAdj(A_adj, D_adj):				#Compares the two adjacencies
									if adj_found_at[A_adj_side] == None:	#Enters indices to adj_found_at ONLY if adjacency has NOT been found already to avoid 
										adj_found_at[A_adj_side] = (D_chr_idx, D_gene_idx, D_adj_side)				#multiple genes in D_2 matching the gene in A_2			

			A_gene = charListToString(gene)
			A_idx_count[A_gene[1:]] = 0
			left_adj_found = False
														
			if adj_found_at[0] != None:										#Some adjacency in D_2 matched with left adjacency in A_2
				left_adj_found = True
				A_idx_count[A_gene[1:]] += 1 																	#Increments label
				A_2[A2_chr_idx][A2_gene_idx] = A_gene + str('copy') + str(A_idx_count[A_gene[1:]])				#Updates labeling of corresponding
				D_2[adj_found_at[0][0]][adj_found_at[0][1]] +=  str('copy') + str(A_idx_count[A_gene[1:]])		#genes in A_2 and D_2
				unmatched_indices.remove((adj_found_at[0][0], adj_found_at[0][1]))								#Since gene in D_2 has been matched,	
																												#it is removed from unmatched_indices

			if adj_found_at[1] != None:										#Some adjacency in D_2 matched with right adjacency in A_2
				A_idx_count[A_gene[1:]] += 1								#Increments label	
				if left_adj_found == False:									#Only right adjacency matched
					A_2[A2_chr_idx][A2_gene_idx] = A_gene + str('copy') + str(A_idx_count[A_gene[1:]])			#Updates labeling of corresponding gene in A_2
					unmatched_indices.remove((adj_found_at[1][0], adj_found_at[1][1]))							#Since gene in D_2 has been matched
																												#It is removed from unmatched_indices

				else:																							#Weakly conserved context found	
					A_2[A2_chr_idx].insert(A2_gene_idx+1, A_gene + str('copy') + str(A_idx_count[A_gene[1:]]))	#Add and relabel tandem duplicate
					TD += 1																						
					unmatched_indices.remove((adj_found_at[1][0], adj_found_at[1][1]))							#Since gene in D_2 has been matched				
																												#it is removed from unmatched_indices
				D_2[adj_found_at[1][0]][adj_found_at[1][1]] += str('copy') + str(A_idx_count[A_gene[1:]])		#Updates labeling of corresponding gene in D_2

			no_adj_conserved = False

			if adj_found_at[0] == None and adj_found_at[1] == None:												#No adjacency conserved
				no_adj_conserved = True
				unmatched_D_chr_idx, unmatched_D_gene_idx = unmatched_indices[0]								#Chooses to match with first unmatched gene in D_2
				A_idx_count[A_gene[1:]] += 1																	#Increments label
				A_2[A2_chr_idx][A2_gene_idx] = A_gene + str('copy') + str(A_idx_count[A_gene[1:]])				#Updates labeling of corresponding
				D_2[unmatched_D_chr_idx][unmatched_D_gene_idx] += str('copy') + str(A_idx_count[A_gene[1:]])	#genes in A_2 and D_2

				for unmatched_D_chr_idx, unmatched_D_gene_idx in unmatched_indices[1:]:		#All remaining copies in D_2 will be matched to a floating duplicate					
					A_idx_count[A_gene[1:]] += 1											#Increments label, creates a floating duplicate and updates label in D_2
					A_2.append( [A_gene + str('copy') + str(A_idx_count[A_gene[1:]]), A_gene + str('copy') + str(A_idx_count[A_gene[1:]])] )
					D_2[unmatched_D_chr_idx][unmatched_D_gene_idx] += str('copy') + str(A_idx_count[A_gene[1:]])
					FD += 1

			if no_adj_conserved == False:													#At least one adjacency conserved AND 
				for unmatched_D_chr_idx, unmatched_D_gene_idx in unmatched_indices:			#at least one gene in D_2 matched with original copy in A_2
					A_idx_count[A_gene[1:]] += 1											#Creates FD for all unmatched copies in this case
					A_2.append( [A_gene + str('copy') + str(A_idx_count[A_gene[1:]]), A_gene + str('copy') + str(A_idx_count[A_gene[1:]])] )
					D_2[unmatched_D_chr_idx][unmatched_D_gene_idx] += str('copy') + str(A_idx_count[A_gene[1:]])
					FD += 1

		#------------------------End of Step 2------------------------------	

#print("\n\n\n\n")
#print("Final A_2:", A_2)
#print("Final D_2:", D_2)


#Distance computation
#---------------------------------------------------------------------------

A_adj_set = pairUp(A_2)														#Forms adjacency sets of A_2 and D_2 respectively														
D_adj_set = pairUp(D_2)

print("Number of FD:", FD)
print("Number of TD:", TD)
print("Number of TD from arrays:", TD_from_arrays)
print("Number of cuts:", len(diff(A_adj_set, D_adj_set)))
print("Number of joins:", len(diff(D_adj_set, A_adj_set)))					
print("Total distance:", sym_diff(A_adj_set, D_adj_set))					#Computes the symmetric difference between the two adj sets