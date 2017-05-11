from functools import reduce

def negate(a):						
	return '-' if a == '+' else '+'

def charListToString(charList):
	return reduce(lambda x, y: x + y, charList, '')

#constructs two lists of three genes each, with the duplicated gene in the centre
#one list in forward direction, one in reverse
def getTriplets(chromosome, gene_idx):
	triplet1 = []
	triplet2 = []
	if gene_idx > 0 and gene_idx < len(chromosome)-1:
		triplet1.append(chromosome[gene_idx-1])
		triplet2.append(negate(chromosome[gene_idx-1][0]) +  chromosome[gene_idx-1][1:])
	
		triplet1.append(chromosome[gene_idx])
		triplet2.append(negate(chromosome[gene_idx][0]) +  chromosome[gene_idx][1:])
	
		triplet1.append(chromosome[gene_idx+1])
		triplet2.append(negate(chromosome[gene_idx+1][0]) +  chromosome[gene_idx+1][1:])
	
		# reverse triplet 2
		triplet2 = triplet2[::-1]
	return triplet1, triplet2

#def getAdj(chromosome, gene_idx): to be completed

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


string = open("genome.txt", "r").read()
genomes = string.split("\n")									#Genes S and D are on different lines in the input file
genomes = [i.strip().split("|") for i in genomes if len(i)]		#Each chromosome in a gene is split by |

#obtained in respective variables as a list of lists
S = [j.strip().split(" ") for j in genomes[0]]
D = [j.strip().split(" ") for j in genomes[1]]

FD = 0	#counter for number of free duplications

#Intermediate genomes S' and D'. These will be the input for the symmetric difference function above.
#Currently the genomes S and D have been copied as they are to form S' and D'.
S_pr = [ [i for i in chromosome] for chromosome in S]
D_pr = [ [i for i in chromosome] for chromosome in D]


search = [None]*2
for chromosome_idx in range(len(S)):
	chromosome = S[chromosome_idx]
	gene_idx_offset = 0
	for gene_idx in range(len(chromosome)):
		gene = chromosome[gene_idx]
		search[0] = gene 							#searches D for the occurence of the gene
		search[1] = negate(gene[0]) + gene[1:]		#or the the reverse

		found_indices = []							#initializes the list of co-ordinates of occurrences
		for sublist_chromosome_idx in range(len(D)):	
			sublist = D[sublist_chromosome_idx]
			for sublist_gene_idx in range(len(sublist)):
				sublist_gene = sublist[sublist_gene_idx]
				if search[0] == sublist_gene or search[1] == sublist_gene:
					found_indices.append( (sublist_chromosome_idx, sublist_gene_idx) )		#appends the co-ordinates (chromosome number, gene number) to found indices
		count = len(found_indices)

		if count > 1:								#count>1 implies duplicates exist in D
			triplet1, triplet2 = getTriplets(chromosome, gene_idx)   						
			
			S_idx_count = {}
			D_idx_count = {}

			no_triplet_indices = []

			is_triplet_found = False
			if len(triplet1) and len(triplet2):
				#print ('triplet: ', triplet1, triplet2)

				#checks if the triplet from S matches any triplet in D.	
				for sublist_chromosome_idx, sublist_gene_idx in found_indices:
					sublist_chromosome = D[sublist_chromosome_idx]
					D_triplet, _ = getTriplets(sublist_chromosome, sublist_gene_idx)
					#print ('D_triplet', D_triplet)

					if len(D_triplet):
						#converts lists to string in order to compare
						triplet1 = charListToString(triplet1)
						triplet2 = charListToString(triplet2)
						D_triplet = charListToString(D_triplet)
						
						if triplet1 == D_triplet:
							print('Found: ', triplet1)
							is_triplet_found = True
							FD += 1
						elif triplet2 == D_triplet:
							print('Found: ', triplet2)
							is_triplet_found = True
							FD += 1
						elif not is_triplet_found:
							no_triplet_indices.append( (sublist_chromosome_idx, sublist_gene_idx) )
			else:
				no_triplet_indices = [i for i in found_indices]	
			
			#If triplet is NOT found, we assume the TD model. 
			if not is_triplet_found:
				S_gene = S_pr[chromosome_idx][gene_idx + gene_idx_offset]
				print ('S_gene', S_gene)
				if S_gene[1:] in S_idx_count:
					S_idx_count[S_gene[1:]] += 1
				else:
					S_idx_count[S_gene[1:]] = 1
				S_pr[chromosome_idx][gene_idx+gene_idx_offset] += str(S_idx_count[S_gene[1:]])	

				S_idx_count[S_gene[1:]] += 1
				S_pr[chromosome_idx] = S_pr[chromosome_idx][0:gene_idx+gene_idx_offset+1] +  [S_gene + str(S_idx_count[S_gene[1:]])] + S_pr[chromosome_idx][gene_idx+gene_idx_offset+1:]
				gene_idx_offset += 1

				no_triplet_indices = [i for i in found_indices]	

			#Needs work.	
			#Assigns appropriate labelling to tandem genes in S_pr and corresponding genes in D_pr		
			for sublist_chromosome_idx, sublist_gene_idx in no_triplet_indices:
				D_gene = D_pr[sublist_chromosome_idx][sublist_gene_idx]
				print('No triplet: ', D_gene)
				if D_gene[1:] in D_idx_count:
					D_idx_count[D_gene[1:]] += 1
				else:
					D_idx_count[D_gene[1:]] = 1
					
				D_pr[sublist_chromosome_idx][sublist_gene_idx] += str(D_idx_count[D_gene[1:]])	

			
print('S_pr', S_pr)
print('D_pr', D_pr)
print('#Free Dups: ', FD)

A = pairUp(S_pr)
B = pairUp(D_pr)

#print( sorted([i for i in B]) ) 
#print( sorted(D) )
print(diff(A, B))
print(diff(B, A))
print(sym_diff(A, B))	#Note that S_pr and D_pr, now have all genes distinct.

Total_distance = sym_diff(A, B) + FD 	#Total dist is the sum of the symmetric difference and the number of cuts required for the free duplications.
print(Total_distance)
