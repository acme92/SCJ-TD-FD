from functools import reduce
from sys import argv
import random

#Function definitions
#-----------------------------------------------------------------------------------------

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
				return([chr_idx, gene_idx])

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

#Outputs adjacency of a gene in a given chromosome
def getAdj(chromosome, gene_idx):
	'''
	get adjacency of a gene  
	index 0: left adjacency, index 1: right adjacency
	eg: for chromosome +a+b+c, 	the adjacencies of gene b (idx: 1) are +a+b (or -b-a) and +b+c (or -c-b)
								the adjacency of gene a (idx: 0) is +a+b (or -b-a)
								the adjacency of gene c (idx: 2) is +b+c (or -c-b) 
	''' 
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
#pairs up adjacent genes in the list of genes that forms the chromosome
def pairUp(a):
    set_a = set()
    for chrom in a:
        for i in range(len(chrom) - 1):
            set_a.add(chrom[i] + chrom[i+1])
    return set_a

#adds the adjacencies in opposite direction to create a look-up dictionary
def duplicateOrder(A):
    A_double = []
    for i in A:
        A_double.append(i)
        # complement
        count = 0
        index = 0
        for idx, val in enumerate(i):
            #print(idx, val)
            if val == '+' or val == '-':
                #print(idx, val)
                count += 1
                if count == 2:
                    index = idx         
        #print("Opp of adj:", i, ":", negate(i[index]) + i[index+1:] + negate(i[0]) + i[1:index])            
        A_double.append(
            negate(i[index]) + i[index+1:] +
            negate(i[0]) + i[1:index]
        )
    #print(A_double)
    return A_double

#A-B
def diff(A, B):
    B_dict = {i:None for i in duplicateOrder(B)}
    return [i for i in A if i not in B_dict]

def sym_diff(A, B):
    return len(diff(A, B)) + len(diff(B, A))


#Input and preprocessing
#-----------------------------------------------------------------------------------------

#Takes file with genomes as argument in command line
if len(argv) < 2:
	print('Usage: python d_SCjDup.py <genome_file>')
	exit(1)

genome_file = argv[1]

#Genomes A and D split by "\n". Each chromosome in a genome split by "|"
#Each chromosome is obtained as a list of genes. Currently all chromosomes treated linear
string = open(genome_file, "r").read()
genomes = string.split("\n")								
genomes = [i.strip().split("|") for i in genomes if len(i)]

#Obtains genomes A and D as a list of chromosomes (list of list of genes)
A = [j.strip().split(" ") for j in genomes[0]]
D = [j.strip().split(" ") for j in genomes[1]]
#print("A:", A)
#print("D:", D)

#Initializes D' to D and removes tandem arrays in D_pr
D_pr = [ [i for i in chromosome] for chromosome in D]
for chr_idx in range(len(D_pr)):
	gene_idx = 0
	while gene_idx < len(D_pr[chr_idx]) - 1:
		if D_pr[chr_idx][gene_idx] == D_pr[chr_idx][gene_idx + 1]:
			del D_pr[chr_idx][gene_idx]
		else:
			gene_idx = gene_idx + 1
#print("D':", D_pr)


#Main program
#-----------------------------------------------------------------------------------------

#Initializes A_2 and D_2 to A and D_pr respectively
A_2 = [ [i for i in chromosome] for chromosome in A]
D_2 = [ [i for i in chromosome] for chromosome in D_pr]
print("A_2:", A_2)
print("D_2:", D_2)

#Constructs a list of co-ordinates of genes in A and shuffles the list for randomness
coords = []
for chr_idx in range(len(A)):
	for gene_idx in range(len(A[chr_idx])):
		coords.append([chr_idx,gene_idx])
#print("Before shuffle:", coords)		
random.shuffle(coords)
#print("After shuffle:", coords)

#Iterates randomly through genes in A
for chr_idx, gene_idx in coords:
	gene = A[chr_idx][gene_idx]

	#Locate chromosome in A_2
	chromosome = A_2[chr_idx]
	#print("Chromosome searched: ", chromosome)

	#Find location of gene in A_2 (chromosome, gene)
	A_2_idx = locateGene(gene, A_2)

	#Assigns elements - gene and its reverse to be searched in D_2
	search = [None]*2
	search[0] = gene 							
	search[1] = negate(gene[0]) + gene[1:]

	#Finds occurrences of the gene in D_2 and lists their indices (chromosome, gene) in D_2
	found_indices = []
	for D_chr_idx in range(len(D_2)):	
			D_chromosome = D_2[D_chr_idx]
			for D_gene_idx in range(len(D_chromosome)):
				D_gene = D_chromosome[D_gene_idx]
				if search[0] == D_gene or search[1] == D_gene:
					found_indices.append( (D_chr_idx, D_gene_idx) )
	
	#Finds copy number of gene				
	count = len(found_indices)
	#print("Copy number of", gene[1:], "is", count)

	#If duplicates exist
	if count > 1:
		
		#STEP 1: If context of gene is strongly conserved

		triplet1, triplet2 = getTriplets(chromosome, A_2_idx[1])		#Triplet -> List of [left neighbour, gene, right neighbour]
		#print("Triplets for", gene, "are:", triplet1, triplet2)

		#Keep a count of number of occurences of the gene in A_2 and D_2
		A_idx_count = {}
		D_idx_count = {}

		# where the left and right adjacencies were found
		# index 0: left adjacency, index 1: right adjacency
		is_triplet_matched = False
		unmatched_triplet_indices = []
		if len(triplet1) and len(triplet2):

			for D_chr_idx, D_gene_idx in found_indices:
				D_chromosome = D_2[D_chr_idx]
				D_triplet, _ = getTriplets(D_chromosome, D_gene_idx)
				#print ('D_triplet for', gene, ":", D_triplet)

				if len(D_triplet):
					#converts lists to string in order to compare
					triplet1 = charListToString(triplet1)
					triplet2 = charListToString(triplet2)
					D_triplet = charListToString(D_triplet)
					A_gene = charListToString(gene)
					#print(A_gene)

					if triplet1 == D_triplet:
						#print('Found triplet: ', triplet1)
						is_triplet_matched = True
						
						A_idx_count[A_gene[1:]] = 1
						#print(A_gene + str(A_idx_count[A_gene[1:]]))
						A_2[A_2_idx[0]][A_2_idx[1]] = A_gene + str(A_idx_count[A_gene[1:]])
						#print ("A_2 is:", A_2)						
						D_2[D_chr_idx][D_gene_idx] = D_2[D_chr_idx][D_gene_idx][0] + A_gene[1:] + str(A_idx_count[A_gene[1:]])
						#print ("D_2 is:", D_2)

					elif triplet2 == D_triplet:
						#print('Found triplet: ', triplet2)
						is_triplet_matched = True	
						
						A_idx_count[A_gene[1:]] = 1
						#print(A_gene + str(A_idx_count[A_gene[1:]]))
						A_2[A_2_idx[0]][A_2_idx[1]] = A_gene + str(A_idx_count[A_gene[1:]])
						#print ("A_2 is:", A_2)						
						D_2[D_chr_idx][D_gene_idx] = D_2[D_chr_idx][D_gene_idx][0] + A_gene[1:] + str(A_idx_count[A_gene[1:]])
						#print ("D_2 is:", D_2)

					else:
						unmatched_triplet_indices.append( (D_chr_idx, D_gene_idx) )

				else:
					unmatched_triplet_indices.append( (D_chr_idx, D_gene_idx) )

		else:				
			unmatched_triplet_indices = [i for i in found_indices]	

		#print("Unmatched indices for", gene[1:], ":", unmatched_triplet_indices)	
		
		if is_triplet_matched == True:
			for unmatched_D_chr_idx, unmatched_D_gene_idx in unmatched_triplet_indices:
				A_idx_count[A_gene[1:]] += 1
				#print("Index count for", gene, ":", A_idx_count)
				A_2.append( [A_gene + str(A_idx_count[A_gene[1:]]), A_gene + str(A_idx_count[A_gene[1:]])] )
				D_2[unmatched_D_chr_idx][unmatched_D_gene_idx] += str(A_idx_count[A_gene[1:]])
				#print("Added floating duplication for:", gene, A_2)	


		#STEP 2: If context of gene is weakly conserved or only one adjacency is conserved

		if not is_triplet_matched:
			# where the left and right adjacencies were found
			# index 0: left adjacency, index 1: right adjacency
			adj_found_at = [None, None]

			for D_chr_idx, D_gene_idx in unmatched_triplet_indices:
				D_chr = D_2[D_chr_idx]

				D_adjacencies = getAdj(D_chr, D_gene_idx)
				A_adjacencies = getAdj(chromosome, gene_idx)
				#print("D_adj for", gene, "are:", D_adjacencies)
				#print("A_adj for", gene, "are:", A_adjacencies)

				for A_adj_idx in range(len(A_adjacencies)):
					A_adj = A_adjacencies[A_adj_idx]

					if len(A_adj):
						for D_adj_idx in range(len(D_adjacencies)):
							D_adj = D_adjacencies[D_adj_idx]

							if len(D_adj):
								if compareAdj(A_adj, D_adj):
									if adj_found_at[A_adj_idx] == None:
										adj_found_at[A_adj_idx] = (D_chr_idx, D_gene_idx, D_adj_idx)
									break

			#print('Adjacencies of ', chromosome[gene_idx])
			#print('at ', adj_found_at)

			A_gene = charListToString(gene)
			A_idx_count[A_gene[1:]] = 0

			left_adj_found = False			
			if adj_found_at[0] != None:
				left_adj_found = True
				A_idx_count[A_gene[1:]] += 1 
				#print("Index count for", gene, ":", A_idx_count)

				# left tandem (according to A)

				#A_idx_count[A_gene[1:]] += 1
				A_2[A_2_idx[0]][A_2_idx[1]] = A_gene + str(A_idx_count[A_gene[1:]])
				#print("A_2 is:", A_2)
				D_2[adj_found_at[0][0]][adj_found_at[0][1]] = D_2[adj_found_at[0][0]][adj_found_at[0][1]][0] + A_gene[1:] + str(A_idx_count[A_gene[1:]])
				#print("D_2 is:", D_2)
				unmatched_triplet_indices.remove((adj_found_at[0][0], adj_found_at[0][1]))
				#print("After removal:", unmatched_triplet_indices, "for:", gene)

			if adj_found_at[1] != None:
				A_idx_count[A_gene[1:]] += 1
				#print("Index count for", gene, ":", A_idx_count)

				# right tandem (according to A)
				if left_adj_found == False:
					#A_idx_count[A_gene[1:]] += 1
					A_2[A_2_idx[0]][A_2_idx[1]] =  A_gene + str(A_idx_count[A_gene[1:]])
					#print("A_2 is:", A_2)
					unmatched_triplet_indices.remove((adj_found_at[1][0], adj_found_at[1][1]))
					#print("After removal:", unmatched_triplet_indices, "for:", gene)
				else:
					A_2[A_2_idx[0]].insert(A_2_idx[1]+1, A_gene + str(A_idx_count[A_gene[1:]]))
					#print("Added tandem:", A_2)
					unmatched_triplet_indices.remove((adj_found_at[1][0], adj_found_at[1][1]))
					#print("After removal:", unmatched_triplet_indices, "for:", gene)					
				D_2[adj_found_at[1][0]][adj_found_at[1][1]] = D_2[adj_found_at[1][0]][adj_found_at[1][1]][0] + A_gene[1:] + str(A_idx_count[A_gene[1:]])
				#print("D_2 is:", D_2)				

			if adj_found_at[0] == None and adj_found_at[1] == None:

				unmatched_D_chr_idx, unmatched_D_gene_idx = unmatched_triplet_indices[0]
				A_idx_count[A_gene[1:]] += 1
				#print("Index count for", gene, ":", A_idx_count)

				A_2[A_2_idx[0]][A_2_idx[1]] = A_gene + str(A_idx_count[A_gene[1:]])
				#print("A_2 is:", A_2)
				D_2[unmatched_D_chr_idx][unmatched_D_gene_idx] += str(A_idx_count[A_gene[1:]])
				#print("D_2 is:", D_2)

				for unmatched_D_chr_idx, unmatched_D_gene_idx in unmatched_triplet_indices[1:]:
					A_idx_count[A_gene[1:]] += 1
					#print("Index count for", gene, ":", A_idx_count[A_gene[1:]])

					A_2.append( [A_gene + str(A_idx_count[A_gene[1:]]), A_gene + str(A_idx_count[A_gene[1:]])] )
					#print("A_2 is:", A_2)
					D_2[unmatched_D_chr_idx][unmatched_D_gene_idx] += str(A_idx_count[A_gene[1:]])
					#print("D_2 is:", D_2)
			
			for unmatched_D_chr_idx, unmatched_D_gene_idx in unmatched_triplet_indices:
				A_idx_count[A_gene[1:]] += 1
				A_2.append( [A_gene + str(A_idx_count[A_gene[1:]]), A_gene + str(A_idx_count[A_gene[1:]])] )
				D_2[unmatched_D_chr_idx][unmatched_D_gene_idx] += str(A_idx_count[A_gene[1:]])
				#print("Added floating duplication for:", gene, A_2)
			

#A_2 = [ [ gene for tandem in chromosome for gene in tandem ] for chromosome in A_2  ]
print("\n\n\n\n")
print("Final A_2:", A_2)
print("Final D_2:", D_2)

#Distance computation
#-----------------------------------------------------------------------------------------

A_pairs = pairUp(A_2)
D_pairs = pairUp(D_2)
#print(A_pairs)
#print(D_pairs)

#print(diff(A_pairs, D_pairs))
#print(diff(D_pairs, A_pairs))
#print(sym_diff(A_pairs, D_pairs))

Total_distance = sym_diff(A_pairs, D_pairs) 	
print(Total_distance)


											
										


