from functools import reduce

def negate(a):
	return '-' if a == '+' else '+'

def charListToString(charList):
	return reduce(lambda x, y: x + y, charList, '')

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

#def getAdj(chromosome, gene_idx):


string = open("genome.txt", "r").read()
genomes = string.split("\n")
genomes = [i.strip().split("|") for i in genomes if len(i)]

S = [j.strip().split(" ") for j in genomes[0]]
D = [j.strip().split(" ") for j in genomes[1]]

FD = 0

search = [None]*2
S_pr = [ [i for i in chromosome] for chromosome in S]
D_pr = [ [i for i in chromosome] for chromosome in D]

for chromosome_idx in range(len(S)):
	chromosome = S[chromosome_idx]
	gene_idx_offset = 0
	for gene_idx in range(len(chromosome)):
		gene = chromosome[gene_idx]
		search[0] = gene
		search[1] = negate(gene[0]) + gene[1:]

		found_indices = []
		for sublist_chromosome_idx in range(len(D)):
			sublist = D[sublist_chromosome_idx]
			for sublist_gene_idx in range(len(sublist)):
				sublist_gene = sublist[sublist_gene_idx]
				if search[0] == sublist_gene or search[1] == sublist_gene:
					found_indices.append( (sublist_chromosome_idx, sublist_gene_idx) )
		count = len(found_indices)

		if count > 1:
			#append the gene in Sprime
			triplet1, triplet2 = getTriplets(chromosome, gene_idx)
			
			S_idx_count = {}
			D_idx_count = {}

			no_triplet_indices = []

			is_triplet_found = False
			if len(triplet1) and len(triplet2):
				print ('triplet: ', triplet1, triplet2)

				for sublist_chromosome_idx, sublist_gene_idx in found_indices:
					sublist_chromosome = D[sublist_chromosome_idx]
					D_triplet, _ = getTriplets(sublist_chromosome, sublist_gene_idx)
					print ('D_triplet', D_triplet)
					if len(D_triplet):
						# TODO: you can avoid doing four comparison
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

