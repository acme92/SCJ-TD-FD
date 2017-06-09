from functools import reduce
from sys import argv

def negate(a):	
	assert(type(a)==str or type(a)==list)

	if type(a) == str:
		if len(a)==1:	
			# for a single character (+ or -)				
			return '-' if a == '+' else '+'
		else:
			# for a gene (eg +a, -b)
			return negate(a[0]) + a[1:]
	elif type(a) == list:
		a_negative = ["" for i in range(len(a))]
		for idx in range(len(a)):
			a_negative[idx] = negate(a[idx])
		return a_negative[::-1]

def charListToString(charList):
	return reduce(lambda x, y: x + y, charList, '')

def genesCompare(chromosome1, chromosome2):
	assert(type(chromosome1)==list and type(chromosome2==list))

	if charListToString(chromosome1) == charListToString(chromosome2):
		return True
	elif charListToString(negate(chromosome1)) == charListToString(chromosome2):
		return True
	else:
		return False

#constructs two lists of three genes each, with the duplicated gene in the centre
#one list in forward direction, one in reverse
def getTriplets(chromosome, gene_idx):
	triplet1 = []
	triplet2 = []
	if gene_idx > 0 and gene_idx < len(chromosome)-1:
		triplet1 = [chromosome[gene_idx-1], chromosome[gene_idx], chromosome[gene_idx+1]]
		triplet2 = negate(triplet1)
	return triplet1, triplet2

def updateSprLabel(S_pr, chromosome_idx, gene_idx, new_gene):
	'''
	If the gene has no label, we need to label it 
	'''
	if len(S_pr[chromosome_idx][gene_idx]):
		# check if unlabeled 
		if S_pr[chromosome_idx][gene_idx][-1][-1].isdigit():
			S_pr[chromosome_idx][gene_idx].append( new_gene )
		else:
			S_pr[chromosome_idx][gene_idx] = [new_gene]
	else:
		S_pr[chromosome_idx][gene_idx] = [new_gene]

	return S_pr

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

def updateIdxCount(S_idx_count, S_gene):
	if S_gene[1:] not in S_idx_count:
		S_idx_count[S_gene[1:]] = 1
	return S_idx_count

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
        A_double.append(
            negate(i[2]) + i[3] +
            negate(i[0]) + i[1]
        )
    return A_double

#A-B
def diff(A, B):
    B_dict = {i:None for i in duplicateOrder(B)}
    return [i for i in A if i not in B_dict]

def sym_diff(A, B):
    return len(diff(A, B)) + len(diff(B, A))

if len(argv) < 2:
	print('Usage: python d_SCjDup.py <genome_file>')
	exit(1)

genome_file = argv[1]

string = open(genome_file, "r").read()
genomes = string.split("\n")									#Genes S and D are on different lines in the input file
genomes = [i.strip().split("|") for i in genomes if len(i)]		#Each chromosome in a gene is split by |

#obtained in respective variables as a list of lists
S = [j.strip().split(" ") for j in genomes[0]]
D = [j.strip().split(" ") for j in genomes[1]]

print("S: ", S)
print("D: ", D)

FD = 0	#counter for number of free duplications

#Intermediate genomes S' and D'. These will be the input for the symmetric difference function above.
#Currently the genomes S and D have been copied as they are to form S' and D'.
# each gene in S_pr is a list, because they may have tandem duplicates later and we want to maintain the indices
S_pr = [ [[i] for i in chromosome] for chromosome in S]
D_pr = [ [i for i in chromosome] for chromosome in D]


search = [None]*2
for chromosome_idx in range(len(S)):
	for gene_idx in range(len(S_pr[chromosome_idx])):
		chromosome = S_pr[chromosome_idx]
		chromosome = [i[-1] for i in chromosome]
		print("Chromosome searched: ", chromosome)
	
		gene = chromosome[gene_idx]
		search[0] = gene 							#searches D for the occurence of the gene
		search[1] = negate(gene[0]) + gene[1:]		#or the the reverse

		found_indices = []							#initializes the list of co-ordinates of occurrences
		for sublist_chromosome_idx in range(len(D_pr)):	
			sublist = D_pr[sublist_chromosome_idx]
			for sublist_gene_idx in range(len(sublist)):
				sublist_gene = sublist[sublist_gene_idx]
				if search[0] == sublist_gene or search[1] == sublist_gene:
					found_indices.append( (sublist_chromosome_idx, sublist_gene_idx) )		#appends the co-ordinates (chromosome number, gene number) to found indices
		count = len(found_indices)

		if count > 1:								#count>1 implies duplicates exist in D_pr
			triplet1, triplet2 = getTriplets(chromosome, gene_idx)   						
			
			S_idx_count = {}
			D_idx_count = {}

			no_triplet_indices = []

			is_triplet_found = False
			if len(triplet1) and len(triplet2):
				#print ('triplet: ', triplet1, triplet2)

				#checks if the triplet from S matches any triplet in D.	
				for sublist_chromosome_idx, sublist_gene_idx in found_indices:
					sublist_chromosome = D_pr[sublist_chromosome_idx]
					D_triplet, _ = getTriplets(sublist_chromosome, sublist_gene_idx)
					#print ('D_triplet', D_triplet)

					if len(D_triplet):
						#converts lists to string in order to compare
						triplet1 = charListToString(triplet1)
						triplet2 = charListToString(triplet2)
						D_triplet = charListToString(D_triplet)
						S_gene = charListToString(S[chromosome_idx][gene_idx])
				
						if triplet1 == D_triplet:
							print('Found triplet: ', triplet1)
							is_triplet_found = True
							FD += 1
							
							S_pr = updateSprLabel( S_pr, chromosome_idx, gene_idx, S_gene + str(S_idx_count[S_gene[1:]]) )
							D_pr[adj_found_at[0][0]][adj_found_at[0][1]] += str(S_idx_count[S_gene[1:]])

							S_idx_count[S_gene[1:]] += 1

						elif triplet2 == D_triplet:
							print('Found triplet: ', triplet2)
							is_triplet_found = True
							FD += 1

							S_pr = updateSprLabel( S_pr, chromosome_idx, gene_idx, S_gene + str(S_idx_count[S_gene[1:]]) )
							D_pr[adj_found_at[0][0]][adj_found_at[0][1]] += str(S_idx_count[S_gene[1:]])

							S_idx_count[S_gene[1:]] += 1

						else: # elif not is_triplet_found:
							no_triplet_indices.append( (sublist_chromosome_idx, sublist_gene_idx) )
					else:
						no_triplet_indices.append( (sublist_chromosome_idx, sublist_gene_idx) )
			else:
				no_triplet_indices = [i for i in found_indices]	
			
			#If triplet is NOT found, we check for adjacencies (weak ...)
			if not is_triplet_found:
				# where the left and right adjacencies were found
				# index 0: left adjacency, index 1: right adjacency
				adj_found_at = [None, None]

				duplicate_adjacency = []

				for no_triplet_D_chromosome_idx, no_triplet_D_gene_idx in no_triplet_indices:
					no_triplet_D_chromosome = D_pr[no_triplet_D_chromosome_idx]
					
					D_adjacencies = getAdj(no_triplet_D_chromosome, no_triplet_D_gene_idx)
					S_adjacencies = getAdj(chromosome, gene_idx)

					for S_adj_idx in range(len(S_adjacencies)):
						S_adj = S_adjacencies[S_adj_idx]
						
						if len(S_adj):
							for D_adj_idx in range(len(D_adjacencies)):
								D_adj = D_adjacencies[D_adj_idx]
								
								if len(D_adj):
									if genesCompare(S_adj, D_adj):
										if adj_found_at[S_adj_idx] == None:
											adj_found_at[S_adj_idx] = (
												no_triplet_D_chromosome_idx, 
												no_triplet_D_gene_idx,
												D_adj_idx
											)
										else:
											# duplicate adjacency
											duplicate_adjacency.append((
												no_triplet_D_chromosome_idx, 
												no_triplet_D_gene_idx,
												D_adj_idx
											))

				print('Adjacencies of ', chromosome[gene_idx])
				print('at ', adj_found_at)

				# if left adjacency at D is found, it starts with one index less
				# eg: if S_adj = ab, is matched with D_adjacencies=['ab', 'bx'], the start index is one less for S
				if adj_found_at[0]:
					found_gene_idx = adj_found_at[0][1] - int(not adj_found_at[0][2])
					print('Left: ', D_pr[adj_found_at[0][0]][found_gene_idx:found_gene_idx+2])
				if adj_found_at[1]:
					found_gene_idx = adj_found_at[1][1] - int(not adj_found_at[1][2])
					print('Right: ', D_pr[adj_found_at[1][0]][found_gene_idx:found_gene_idx+2])

				S_gene = charListToString(S[chromosome_idx][gene_idx])
				if adj_found_at[0] != None:
					S_idx_count = updateIdxCount(S_idx_count, S_gene)

					# left tandem (according to S)
					S_pr = updateSprLabel( S_pr, chromosome_idx, gene_idx, S_gene + str(S_idx_count[S_gene[1:]]) )
					D_pr[adj_found_at[0][0]][adj_found_at[0][1]] += str(S_idx_count[S_gene[1:]])

					S_idx_count[S_gene[1:]] += 1

				if adj_found_at[1] != None:
					S_idx_count = updateIdxCount(S_idx_count, S_gene)

					# right tandem (according to S)
					S_pr = updateSprLabel( S_pr, chromosome_idx, gene_idx, S_gene + str(S_idx_count[S_gene[1:]]) )
					D_pr[adj_found_at[1][0]][adj_found_at[1][1]] += str(S_idx_count[S_gene[1:]])

					S_idx_count[S_gene[1:]] += 1

				if adj_found_at[0] == None and adj_found_at[1] == None:
					
					no_triplet_D_chromosome_idx, no_triplet_D_gene_idx = no_triplet_indices[0]
					S_idx_count = updateIdxCount(S_idx_count, S_gene)

					# left tandem (according to S)
					S_pr = updateSprLabel( S_pr, chromosome_idx, gene_idx, S_gene + str(S_idx_count[S_gene[1:]]) )
					D_pr[no_triplet_D_chromosome_idx][no_triplet_D_gene_idx] += str(S_idx_count[S_gene[1:]])

					S_idx_count[S_gene[1:]] += 1

					for no_triplet_D_chromosome_idx, no_triplet_D_gene_idx in no_triplet_indices[1:]:
						S_idx_count = updateIdxCount(S_idx_count, S_gene)

						S_pr.append( [[S_gene + str(S_idx_count[S_gene[1:]]) + S_gene + str(S_idx_count[S_gene[1:]])]] )
						D_pr[no_triplet_D_chromosome_idx][no_triplet_D_gene_idx] += str(S_idx_count[S_gene[1:]])
						
						S_idx_count[S_gene[1:]] += 1

				for duplicate_adjaceny_D_chromosome_idx, duplicate_adjaceny_D_gene_idx, _ in duplicate_adjacency:
					S_idx_count = updateIdxCount(S_idx_count, S_gene)

					S_pr.append( [[S_gene + str(S_idx_count[S_gene[1:]]) + S_gene + str(S_idx_count[S_gene[1:]])]] )
					D_pr[duplicate_adjaceny_D_chromosome_idx][duplicate_adjaceny_D_gene_idx] += str(S_idx_count[S_gene[1:]])
					
					S_idx_count[S_gene[1:]] += 1
				# no_triplet_indices = [i for i in found_indices]

			else:
				# floating duplicates
				S_gene = charListToString(S[chromosome_idx][gene_idx])

				for no_triplet_D_chromosome_idx, no_triplet_D_gene_idx in no_triplet_indices:
					S_idx_count = updateIdxCount(S_idx_count, S_gene)

					S_pr.append( [[S_gene + str(S_idx_count[S_gene[1:]]) + S_gene + str(S_idx_count[S_gene[1:]])]] )
					D_pr[no_triplet_D_chromosome_idx][no_triplet_D_gene_idx] += str(S_idx_count[S_gene[1:]])
					
					S_idx_count[S_gene[1:]] += 1

					print("added free dup: ", S_pr)

			print("S_pr: ", S_pr)
			print("D_pr: ", D_pr)
			#Needs work.	
			#Assigns appropriate labelling to tandem genes in S_pr and corresponding genes in D_pr		
			# for sublist_chromosome_idx, sublist_gene_idx in no_triplet_indices:
			# 	D_gene = D_pr[sublist_chromosome_idx][sublist_gene_idx]
			# 	print('No triplet: ', D_gene)
			# 	if D_gene[1:] in D_idx_count:
			# 		D_idx_count[D_gene[1:]] += 1
			# 	else:
			# 		D_idx_count[D_gene[1:]] = 1
					
			# 	D_pr[sublist_chromosome_idx][sublist_gene_idx] += str(D_idx_count[D_gene[1:]])	

# if no tandem duplicate, add the gene without label
for chromosome_idx in range(len(S_pr)):
	for gene_idx in range(len(S_pr[chromosome_idx])):
		if len(S_pr[chromosome_idx][gene_idx]) == 0:
			S_pr[chromosome_idx][gene_idx].append( charListToString(S[chromosome_idx][gene_idx]) )

print("\n\n\nEnd")
print('S_pr', S_pr)
print('D_pr', D_pr)
print('#Free Dups: ', FD)

# flatten S_pr so that you don't it's a single dimension
S_pr = [ [ gene for tandem in chromosome for gene in tandem ] for chromosome in S_pr  ]

print('S_pr unwrapped: ', S_pr)
A = pairUp(S_pr)
B = pairUp(D_pr)

#print( sorted([i for i in B]) ) 
#print( sorted(D) )
print(diff(A, B))
print(diff(B, A))
print(sym_diff(A, B))	#Note that S_pr and D_pr, now have all genes distinct.

Total_distance = sym_diff(A, B) + FD 	#Total dist is the sum of the symmetric difference and the number of cuts required for the free duplications.
print(Total_distance)