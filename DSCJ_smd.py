from sys import argv
import random
import networkx as nx
import matplotlib.pyplot as plt


#Using python 3 interpreter

#Function definitions
#---------------------------------------------------------------------------
#Gene in opposite direction
def negate(gene):
	if gene:
		return gene[1:] if gene[0] == '-' else str('-' + gene)
	else:
		return None

#Forms a list of all genes in input genome
def listGene(genome):
	gene_list = []
	for chromosome in genome:
		for gene in chromosome:
			if gene[0] == '-':
				if negate(gene) not in gene_list:
					gene_list.append(negate(gene))
			else:
				if gene not in gene_list:
					gene_list.append(gene)
	return gene_list

#Forms adjacency list of input genome
def listAdj(genome):
	adj_list = []
	for chromosome in genome:
		for gene_idx in range(len(chromosome) - 1):
			if chromosome[gene_idx][0] == '-':
				left = (chromosome[gene_idx][1:], 't')
			else:
				left = (chromosome[gene_idx], 'h')
			if chromosome[gene_idx + 1][0] == '-':
				right = (chromosome[gene_idx + 1][1:], 'h')
			else:
				right = (chromosome[gene_idx + 1], 't')
			adj_list.append([left, right])
	return adj_list

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



#Main functions
#---------------------------------------------------------------------------
#Finds distance between the two given genomes in the input file
def distance(filename):
	string = open(filename, "r").read()
	substrings = string.split("\n")
	substrings = [line for line in substrings if line and line[0] != '#']	#Read line only if it is nonempty and not a comment.

	genomes = [[] for i in range(2)]										#A = genomes[0] and D = genomes[1]
	seen = set()															#Set of genes in A to check for trivialness
	trivial = True

	n_genes = [0,0]															#Maintain a count of number of genes in A and D, respectively
	i = -1
	for line in substrings:
		if line[-1] not in {')','|'}:										#Skip lines that do not end with | or ), as they aren't chromosomes.
			i += 1
		elif line[-1] == '|':												
			line = line.split(' ')
			line = [x for x in line if x != '|']
			if i == 0:														#Check trivialness only for A.
				for gene in line:
					if gene not in seen and negate(gene) not in seen:
						seen.add(gene)
					else:
						trivial = False
			genomes[i].append(line)
			n_genes[i] += len(line)
		else:
			line = line.split(' ')											#If circular chromosome, then add first gene after the last gene.
			line[-1] = line[0]												#For e.g., (a,b,c) = [a,b,c,a]. Required for tandem array check in circular chromosomes
			if i == 0:
				for gene in line[:-1]:
					if gene not in seen and negate(gene) not in seen:
						seen.add(gene)
					else:
						trivial = False
			genomes[i].append(line)
			n_genes[i] += len(line) - 1

	A = genomes[0]
	D = genomes[1]

	if trivial == False:													#If gene repeats, A is nontrivial. Terminate program.
		print("Error message: Ancestor genome must be trivial.")
		quit()

	if set(listGene(A)) != set(listGene(D)):								#If A and D on different set of gene families, terminate program.
		print("Error message: Ancestor and descendant genomes have different sets of gene families.")
		quit()

	TD_from_arrays = 0														#Maintain a count of TD from arrays.

	for chromosome in D:										
		for gene_idx in range(len(chromosome)):
			while gene_idx < len(chromosome) - 1:
				if chromosome[gene_idx] == chromosome[gene_idx + 1]:		#Remove tandem arrays, if any.
					del chromosome[gene_idx]	
					TD_from_arrays += 1
				else:
					gene_idx += 1

	n_genes[1] -= TD_from_arrays											#Number of genes in D after removing tandem arrays.
	n_duplicates = n_genes[1] - n_genes[0]									#Number of genes in D - number of genes in A

	A_adj = listAdj(A)
	D_adj = listAdj(D)

	preserved_adj = [adj for adj in A_adj if adj in D_adj or list(reversed(adj)) in D_adj]		#Intersection of adjacency sets, A and D
	n_cuts = len(A_adj) - len(preserved_adj)
	n_joins = len(D_adj) - len(preserved_adj)

	d_DSCJ = n_cuts + n_joins + 2*n_duplicates + TD_from_arrays				#d_DSCJ(A,D) = |A-D| + |D-A| + 2*n_d + TDA.

	print(d_DSCJ)
	print(n_cuts)
	print(n_joins)
	print(n_duplicates)
	print(TD_from_arrays)



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
	substrings = string.split("\n")
	substrings = [line for line in substrings if line and line[0] != '#']	#Read line only if it is nonempty and not a comment.

	genomes = [[] for i in range(2)]										#A = genomes[0] and D = genomes[1]
	seen = set()															#Set of genes in A to check for trivialness
	trivial = True

	chr_type = [[] for i in range(2)]										#A log of each chromosome being linear or circular

	i = -1
	for line in substrings:
		if line[-1] not in {')','|'}:										#Skip lines that do not end with | or ), as they aren't chromosomes.
			i += 1
		elif line[-1] == '|':
			line = line.split(' ')
			line = [x for x in line if x != '|']
			if i == 0:														#Check trivialness only for A.
				for gene in line:			
					if gene not in seen and negate(gene) not in seen:
						seen.add(gene)	
					else:
						trivial = False
			chr_type[i].append('L')
			genomes[i].append(line)	
		else:
			line = line.split(' ')											#If circular chromosome, then add first gene after the last gene.
			line[-1] = line[0]												#For e.g., (a,b,c) = [a,b,c,a]. Required for tandem array check in circular chromosomes
			if i == 0:
				for gene in line[:-1]:
					if gene not in seen and negate(gene) not in seen:
						seen.add(gene)
					else:
						trivial = False
			chr_type[i].append('C')											
			genomes[i].append(line)

	A = genomes[0]
	D = genomes[1]	

	if trivial == False:
		print("Error message: Ancestor genome must be trivial.")			#If gene repeats, A is nontrivial. Terminate program.
		quit()

	if set(listGene(A)) != set(listGene(D)):								#If A and D on different set of gene families, terminate program.
		print("Error message: Ancestor and descendant genomes have different sets of gene families.")
		quit()

	TD_from_arrays = 0														#Maintain a count of TD from arrays.

	for chromosome in D:				
		for gene_idx in range(len(chromosome)):
			while gene_idx < len(chromosome) - 1:
				if chromosome[gene_idx] == chromosome[gene_idx + 1]:
					del chromosome[gene_idx]								#Remove tandem arrays, if any.
					TD_from_arrays += 1
				else:
					gene_idx += 1
			
	A_dict = {}	#Dictionary for A. Key = gene. Value = (Idx, Sign, Left neighbor, Left neighbor index, Right neighbor, Right neighbor index, TD)		
	for i in range(len(A)):
		if chr_type[0][i] == 'L':									#If chromosome is circular, add appropriate neighbors for genes at both ends.
			for j in range(len(A[i])):								#For e.g.: In (a,b,c) LN of a is c and RN of c is a. 
				if j == 0:
					if A[i][j][0] == '-':							#If gene orientation is negative.
						A_dict[negate(A[i][j])] = {'Idx': (i,j), 'Sign': 'Neg', 'LN': None, 'LNIdx': None, 'RN': A[i][j+1], 'RNIdx': (i,j+1), 'TD': None}
					else:
						A_dict[(A[i][j])] = {'Idx': (i,j), 'Sign': 'Pos', 'LN': None, 'LNIdx': None, 'RN': A[i][j+1], 'RNIdx': (i,j+1), 'TD': None}
				elif j == len(A[i])-1:	
					if A[i][j][0] == '-':
						A_dict[negate(A[i][j])] = {'Idx': (i,j), 'Sign': 'Neg', 'LN': A[i][j-1], 'LNIdx': (i,j-1), 'RN': None, 'RNIdx': None, 'TD': None}
					else:
						A_dict[(A[i][j])] = {'Idx': (i,j), 'Sign': 'Pos', 'LN': A[i][j-1], 'LNIdx': (i,j-1), 'RN': None, 'RNIdx': None, 'TD': None}
				else:
					if A[i][j][0] == '-':
						A_dict[negate(A[i][j])] = {'Idx': (i,j), 'Sign': 'Neg', 'LN': A[i][j-1], 'LNIdx': (i,j-1), 'RN': A[i][j+1], 'RNIdx': (i,j+1), 'TD': None}
					else:
						A_dict[(A[i][j])] = {'Idx': (i,j), 'Sign': 'Pos', 'LN': A[i][j-1], 'LNIdx': (i,j-1), 'RN': A[i][j+1], 'RNIdx': (i,j+1), 'TD': None}					
		if chr_type[0][i] == 'C':
			A[i] = A[i][:-1]
			for j in range(len(A[i])):
				if j == 0:
					if A[i][j][0] == '-':
						A_dict[negate(A[i][j])] = {'Idx': (i,j), 'Sign': 'Neg', 'LN': A[i][-1], 'LNIdx': (i,len(A[i])-1), 'RN': A[i][j + 1], 'RNIdx': (i,j+1), 'TD': None}
					else:
						A_dict[(A[i][j])] = {'Idx': (i,j), 'Sign': 'Pos', 'LN': A[i][-1], 'LNIdx': (i,len(A[i])-1), 'RN': A[i][j + 1], 'RNIdx': (i,j+1), 'TD': None}
				elif j == len(A[i])-1:	
					if A[i][j][0] == '-':
						A_dict[negate(A[i][j])] = {'Idx': (i,j), 'Sign': 'Neg', 'LN': A[i][j-1], 'LNIdx': (i,j-1), 'RN': A[i][0], 'RNIdx': (i,0), 'TD': None}
					else:
						A_dict[(A[i][j])] = {'Idx': (i,j), 'Sign': 'Pos', 'LN': A[i][j-1], 'LNIdx': (i,j-1), 'RN': A[i][0], 'RNIdx': (i,0), 'TD': None}
				else:
					if A[i][j][0] == '-':
						A_dict[negate(A[i][j])] = {'Idx': (i,j), 'Sign': 'Neg', 'LN': A[i][j-1], 'LNIdx': (i,j-1), 'RN': A[i][j+1], 'RNIdx': (i,j+1), 'TD': None}
					else:
						A_dict[(A[i][j])] = {'Idx': (i,j), 'Sign': 'Pos', 'LN': A[i][j-1], 'LNIdx': (i,j-1), 'RN': A[i][j+1], 'RNIdx': (i,j+1), 'TD': None}										
	
	Idx_dict = {}	#Dictionary for indices. Key = Gene family, Value = List of positions of g in D. 
	D_dict = {}		#Dictionary for D. Key = Index. Value = (Sign, Left neighbor, Left neighbor index, Right neighbor, Right neighbor index)
	for i in range(len(D)):
		if chr_type[1][i] == 'L':									#If chromosome is circular, add appropriate neighbors for genes at both ends.
			for j in range(len(D[i])):								#For e.g.: In (a,b,c) LN of a is c and RN of c is a. 
				if D[i][j][0] == '-':								#If gene orientation is negative.
					try:
						Idx_dict[negate(D[i][j])].append((i,j))
					except KeyError:
						Idx_dict[negate(D[i][j])] = [(i,j)]
					if j == 0:
						D_dict[(i,j)] = {'Sign': 'Neg', 'LN': None, 'LNIdx': None, 'RN': D[i][j+1], 'RNIdx': (i,j+1)}
					elif j == len(D[i])-1:
						D_dict[(i,j)] = {'Sign': 'Neg', 'LN': D[i][j-1], 'LNIdx': (i,j-1), 'RN': None, 'RNIdx': None}
					else:
						D_dict[(i,j)] = {'Sign': 'Neg', 'LN': D[i][j-1], 'LNIdx': (i,j-1), 'RN': D[i][j+1], 'RNIdx': (i,j+1)}
				else:
					try:
						Idx_dict[D[i][j]].append((i,j))
					except KeyError:
						Idx_dict[D[i][j]] = [(i,j)]
					if j == 0:
						D_dict[(i,j)] = {'Sign': 'Pos', 'LN': None, 'LNIdx': None, 'RN': D[i][j+1], 'RNIdx': (i,j+1)}
					elif j == len(D[i])-1:
						D_dict[(i,j)] = {'Sign': 'Pos', 'LN': D[i][j-1], 'LNIdx': (i,j-1), 'RN': None, 'RNIdx': None}
					else:
						D_dict[(i,j)] = {'Sign': 'Pos', 'LN': D[i][j-1], 'LNIdx': (i,j-1), 'RN': D[i][j+1], 'RNIdx': (i,j+1)}
		if chr_type[1][i] == 'C':
			D[i] = D[i][:-1]
			for j in range(len(D[i])):
				if D[i][j][0] == '-':
					try:
						Idx_dict[negate(D[i][j])].append((i,j))
					except KeyError:
						Idx_dict[negate(D[i][j])] = [(i,j)]
					if j == 0:
						D_dict[(i,j)] = {'Sign': 'Neg', 'LN': D[i][-1], 'LNIdx': (i,len(D[i])-1), 'RN': D[i][j+1], 'RNIdx': (i,j+1)}
					elif j == len(D[i])-1:
						D_dict[(i,j)] = {'Sign': 'Neg', 'LN': D[i][j-1], 'LNIdx': (i,j-1), 'RN': D[i][0], 'RNIdx': (i,0)}
					else:
						D_dict[(i,j)] = {'Sign': 'Neg', 'LN': D[i][j-1], 'LNIdx': (i,j-1), 'RN': D[i][j+1], 'RNIdx': (i,j+1)}
				else:
					try:
						Idx_dict[D[i][j]].append((i,j))
					except KeyError:
						Idx_dict[D[i][j]] = [(i,j)]
					if j == 0:
						D_dict[(i,j)] = {'Sign': 'Pos', 'LN': D[i][-1], 'LNIdx': (i,len(D[i])-1), 'RN': D[i][j+1], 'RNIdx': (i,j+1)}
					elif j == len(D[i])-1:
						D_dict[(i,j)] = {'Sign': 'Pos', 'LN': D[i][j-1], 'LNIdx': (i,j-1), 'RN': D[i][0], 'RNIdx': (i,0)}
					else:
						D_dict[(i,j)] = {'Sign': 'Pos', 'LN': D[i][j-1], 'LNIdx': (i,j-1), 'RN': D[i][j+1], 'RNIdx': (i,j+1)}

	#outputfile = open("idx1.txt", "w")					
	#for key in Idx_dict:
	#	outputfile.write(str(key)+":"+str(len(Idx_dict[key]))+"\n")
					
	coords = []										#List of co-ordinates of A, shuffled for randomness							
	for i in range(len(A)):												
		for j in range(len(A[i])):
			coords.append([i,j])		
	random.shuffle(coords)

	FD = []											#List of FD created
	TD = []											#List of TD created
	n_TD = 0
	n_FD = 0

	for i, j in coords:
		gene = A[i][j]
		if gene[0] == '-':
			gene = negate(gene)

		if len(Idx_dict[gene]) > 1:
			strong_context = False 															#Checking if context strongly conserved
			for index in Idx_dict[gene]:
				if strong_context == False:			
					if A_dict[gene]['Sign'] == D_dict[index]['Sign']:
						if A_dict[gene]['LN'] and A_dict[gene]['LN'] == D_dict[index]['LN'] and A_dict[gene]['RN'] and A_dict[gene]['RN'] == D_dict[index]['RN']:
							strong_context = True
							LN, RN = A_dict[gene]['LN'], A_dict[gene]['RN']					#Updating A_dict
							if LN[0] == '-':
								A_dict[LN[1:]]['RN'] = A_dict[LN[1:]]['RN']+str('copy')+str(1)
							else:
								A_dict[LN]['RN'] = A_dict[LN]['RN']+str('copy')+str(1)
							if RN[0] == '-':
								A_dict[RN[1:]]['LN'] = A_dict[RN[1:]]['LN']+str('copy')+str(1)
							else:
								A_dict[RN]['LN'] = A_dict[RN]['LN']+str('copy')+str(1)
							A_dict[gene+str('copy')+str(1)] = A_dict[gene]
							del A_dict[gene]

							Idx_dict[gene+str('copy')+str(1)] = [index] 					#Updating Idx_dict and D_dict
							Idx_dict[gene].remove(index)
							LNIdx, RNIdx = D_dict[index]['LNIdx'], D_dict[index]['RNIdx']
							D_dict[LNIdx]['RN'] = D_dict[LNIdx]['RN']+str('copy')+str(1)
							D_dict[RNIdx]['LN'] = D_dict[RNIdx]['LN']+str('copy')+str(1)

					else:
						if A_dict[gene]['LN'] and A_dict[gene]['LN'] == negate(D_dict[index]['RN']) and A_dict[gene]['RN'] and A_dict[gene]['RN'] == negate(D_dict[index]['LN']):
							strong_context = True
							LN, RN = A_dict[gene]['LN'], A_dict[gene]['RN']					#Updating A_dict
							if LN[0] == '-':
								A_dict[LN[1:]]['RN'] = A_dict[LN[1:]]['RN']+str('copy')+str(1)
							else:
								A_dict[LN]['RN'] = A_dict[LN]['RN']+str('copy')+str(1)
							if RN[0] == '-':
								A_dict[RN[1:]]['LN'] = A_dict[RN[1:]]['LN']+str('copy')+str(1)
							else:
								A_dict[RN]['LN'] = A_dict[RN]['LN']+str('copy')+str(1)
							A_dict[gene+str('copy')+str(1)] = A_dict[gene]					#Introduce entry for relabeled gene in A_dict	
							del A_dict[gene]												#Delete the gene from A_dict since it is now relabeled

							Idx_dict[gene+str('copy')+str(1)] = [index] 					#Updating Idx_dict and D_dict
							Idx_dict[gene].remove(index)
							LNIdx, RNIdx = D_dict[index]['LNIdx'], D_dict[index]['RNIdx']
							D_dict[LNIdx]['RN'] = D_dict[LNIdx]['RN']+str('copy')+str(1)
							D_dict[RNIdx]['LN'] = D_dict[RNIdx]['LN']+str('copy')+str(1)

			if strong_context == True:														#If strong conserved context, remaining genes matched with FDs
				i = 2
				for index in Idx_dict[gene]:												#Updating list of FDs
					FD.append(gene+str('copy')+str(i))
					Idx_dict[gene+str('copy')+str(i)] = [index] 							#Updating D_dict and Idx_dict for each FD
					LNIdx, RNIdx = D_dict[index]['LNIdx'], D_dict[index]['RNIdx']
					if LNIdx:
						D_dict[LNIdx]['RN'] = D_dict[LNIdx]['RN']+str('copy')+str(i)
					if RNIdx:	
						D_dict[RNIdx]['LN'] = D_dict[RNIdx]['LN']+str('copy')+str(i)
					i += 1
					
			#Strong context (Case 1) checked.

			if strong_context == False:														#If not found, check for weak context.		
				left_adj_at, right_adj_at = None, None
				for index in Idx_dict[gene]:												#For every index in D having the gene get right and left neighbors
					if not left_adj_at:
						if A_dict[gene]['Sign'] == D_dict[index]['Sign']:					#Checking for left adjacency match
							if A_dict[gene]['LN'] and A_dict[gene]['LN'] == D_dict[index]['LN']:
								left_adj_at = index
						else:
							if A_dict[gene]['LN'] and A_dict[gene]['LN'] == negate(D_dict[index]['RN']):
								left_adj_at = index
					if not right_adj_at:
						if A_dict[gene]['Sign'] == D_dict[index]['Sign']:					#Checking for right adjacency match		
							if A_dict[gene]['RN'] and A_dict[gene]['RN'] == D_dict[index]['RN']:
								right_adj_at = index
						else:
							if A_dict[gene]['RN'] and A_dict[gene]['RN'] == negate(D_dict[index]['LN']):
								right_adj_at = index

				if left_adj_at and right_adj_at:
					TD.append(gene)
					Idx_dict[gene+str('copy')+str(1)] = [left_adj_at] 						#Updating Idx_dict
					Idx_dict[gene+str('copy')+str(2)] = [right_adj_at]
					Idx_dict[gene].remove(left_adj_at)
					Idx_dict[gene].remove(right_adj_at)				
					
					if A_dict[gene]['Sign'] == D_dict[left_adj_at]['Sign']:					#Updating D_dict
						LNIdx = D_dict[left_adj_at]['LNIdx']
						RNIdx = D_dict[left_adj_at]['RNIdx']
						if LNIdx:
							D_dict[LNIdx]['RN'] = D_dict[LNIdx]['RN']+str('copy')+str(1)
						if RNIdx:	
							D_dict[RNIdx]['LN'] = D_dict[RNIdx]['LN']+str('copy')+str(1)
					else:
						LNIdx = D_dict[left_adj_at]['RNIdx']
						RNIdx = D_dict[left_adj_at]['LNIdx']
						if LNIdx:
							D_dict[LNIdx]['LN'] = D_dict[LNIdx]['LN']+str('copy')+str(1)
						if RNIdx:
							D_dict[RNIdx]['RN'] = D_dict[RNIdx]['RN']+str('copy')+str(1)
					if A_dict[gene]['Sign'] == D_dict[right_adj_at]['Sign']:
						RNIdx = D_dict[right_adj_at]['RNIdx']
						LNIdx = D_dict[right_adj_at]['LNIdx']
						if RNIdx:
							D_dict[RNIdx]['LN'] = D_dict[RNIdx]['LN']+str('copy')+str(2)
						if LNIdx:	
							D_dict[LNIdx]['RN'] = D_dict[LNIdx]['RN']+str('copy')+str(2)
					else:
						RNIdx = D_dict[right_adj_at]['LNIdx']
						LNIdx = D_dict[right_adj_at]['RNIdx']
						if RNIdx:
							D_dict[RNIdx]['RN'] = D_dict[RNIdx]['RN']+str('copy')+str(2)
						if LNIdx:
							D_dict[LNIdx]['LN'] = D_dict[LNIdx]['LN']+str('copy')+str(2)
					
					LN = A_dict[gene]['LN']													#Updating A_dict
					RN = A_dict[gene]['RN']
					Idx = A_dict[gene]['Idx']
					Sign = A_dict[gene]['Sign']
					LNIdx = A_dict[gene]['LNIdx']
					RNIdx = A_dict[gene]['RNIdx']
					A_dict[gene+str('copy')+str(1)] = {'Idx': Idx, 'Sign': Sign, 'LN': LN, 'RN': RN, 'LNIdx': LNIdx, 'RNIdx': RNIdx}
					A_dict[gene+str('copy')+str(2)] = {'Idx': Idx, 'Sign': Sign, 'LN': LN, 'RN': RN, 'LNIdx': LNIdx, 'RNIdx': RNIdx}

					if LN[0] == '-':														#Introduce entry for TD in A_dict and relabel
						A_dict[LN[1:]]['RN'] = A_dict[LN[1:]]['RN']+str('copy')+str(1)		#Delete the gene from A_dict since it is now relabeled
						A_dict[gene+str('copy')+str(2)]['LN'] = A_dict[LN[1:]]['RN']
					else:
						A_dict[LN]['RN'] = A_dict[LN]['RN']+str('copy')+str(1)
						A_dict[gene+str('copy')+str(2)]['LN'] = A_dict[LN]['RN']
					if RN[0] == '-':
						A_dict[RN[1:]]['LN'] = A_dict[RN[1:]]['LN']+str('copy')+str(2)						
						A_dict[gene+str('copy')+str(1)]['RN'] = A_dict[RN[1:]]['LN']
					else:
						A_dict[RN]['LN'] = A_dict[RN]['LN']+str('copy')+str(2)
						A_dict[gene+str('copy')+str(1)]['RN'] = A_dict[RN]['LN']
					del A_dict[gene]

				#Context not conserved. Check for left match or right match or none	
				elif left_adj_at and not right_adj_at:										
					Idx_dict[gene+str('copy')+str(1)] = [left_adj_at]						#Updating Idx_dict
					Idx_dict[gene].remove(left_adj_at) 

					if A_dict[gene]['Sign'] == D_dict[left_adj_at]['Sign']:					#Updating D_dict
						LNIdx = D_dict[left_adj_at]['LNIdx']
						D_dict[LNIdx]['RN'] = D_dict[LNIdx]['RN']+str('copy')+str(1)
						RNIdx = D_dict[left_adj_at]['RNIdx']
						if RNIdx:
							D_dict[RNIdx]['LN'] = D_dict[RNIdx]['LN']+str('copy')+str(1)
					else:
						LNIdx = D_dict[left_adj_at]['RNIdx']
						D_dict[LNIdx]['LN'] = D_dict[LNIdx]['LN']+str('copy')+str(1) 
						RNIdx = D_dict[left_adj_at]['LNIdx']
						if RNIdx:
							D_dict[RNIdx]['RN'] = D_dict[RNIdx]['RN']+str('copy')+str(1)						

					LN = A_dict[gene]['LN']													#Updating A_dict
					RN = A_dict[gene]['RN']
						
					A_dict[gene+str('copy')+str(1)] = A_dict[gene]							#Introduce entry for relabeled gene in A_dict
					del A_dict[gene]														#Delete relabeled gene from A_dict

					if LN:
						if LN[0] == '-':
							A_dict[LN[1:]]['RN'] = A_dict[LN[1:]]['RN']+str('copy')+str(1)
						else:
							A_dict[LN]['RN'] = A_dict[LN]['RN']+str('copy')+str(1)
					if RN:
						if RN[0] == '-':
							A_dict[RN[1:]]['LN'] = A_dict[RN[1:]]['LN']+str('copy')+str(1)
						else:
							A_dict[RN]['LN'] = A_dict[RN]['LN']+str('copy')+str(1)

				elif right_adj_at and not left_adj_at:
					Idx_dict[gene+str('copy')+str(1)] = [right_adj_at]						#Updating Idx_dict and D_dict
					Idx_dict[gene].remove(right_adj_at)

					if A_dict[gene]['Sign'] == D_dict[right_adj_at]['Sign']:
						RNIdx = D_dict[right_adj_at]['RNIdx']
						D_dict[RNIdx]['LN'] = D_dict[RNIdx]['LN']+str('copy')+str(1)
						LNIdx = D_dict[right_adj_at]['LNIdx']
						if LNIdx:
							D_dict[LNIdx]['RN'] = D_dict[LNIdx]['RN']+str('copy')+str(1)					
					else:
						RNIdx = D_dict[right_adj_at]['LNIdx']
						D_dict[RNIdx]['RN'] = D_dict[RNIdx]['RN']+str('copy')+str(1)
						LNIdx = D_dict[right_adj_at]['RNIdx']
						if LNIdx:
							D_dict[LNIdx]['LN'] = D_dict[LNIdx]['LN']+str('copy')+str(1)
													
					LN = A_dict[gene]['LN']													#Updating A_dict
					RN = A_dict[gene]['RN']
						
					A_dict[gene+str('copy')+str(1)] = A_dict[gene]							#Introduce entry for relabeled gene in A_dict	
					del A_dict[gene]														#Delete relabeled gene from A_dict

					if LN:
						if LN[0] == '-':
							A_dict[LN[1:]]['RN'] = A_dict[LN[1:]]['RN']+str('copy')+str(1)
						else:
							A_dict[LN]['RN'] = A_dict[LN]['RN']+str('copy')+str(1)
					if RN:
						if RN[0] == '-':
							A_dict[RN[1:]]['LN'] = A_dict[RN[1:]]['LN']+str('copy')+str(1)
						else:
							A_dict[RN]['LN'] = A_dict[RN]['LN']+str('copy')+str(1)

				else:
					Idx_dict[gene+str('copy')+str(1)] = [index] 							#Updating Idx_dict and D_dict
					Idx_dict[gene].remove(index)

					if A_dict[gene]['Sign'] == D_dict[index]['Sign']:
						LNIdx = D_dict[index]['LNIdx']
						if LNIdx:
							D_dict[LNIdx]['RN'] = D_dict[LNIdx]['RN']+str('copy')+str(1)
						RNIdx = D_dict[index]['RNIdx']
						if RNIdx:
							D_dict[RNIdx]['LN'] = D_dict[RNIdx]['LN']+str('copy')+str(1)
					else:
						LNIdx = D_dict[index]['RNIdx']
						if LNIdx:
							D_dict[LNIdx]['LN'] = D_dict[LNIdx]['LN']+str('copy')+str(1)
						RNIdx = D_dict[index]['LNIdx']
						if RNIdx:
							D_dict[RNIdx]['RN'] = D_dict[RNIdx]['RN']+str('copy')+str(1)

					LN = A_dict[gene]['LN']													#Updating A_dict
					RN = A_dict[gene]['RN']
						
					A_dict[gene+str('copy')+str(1)] = A_dict[gene]							#Introduce entry for relabeled gene in A_dict
					del A_dict[gene]														#Delete relabeled gene from A_dict

					if LN:
						if LN[0] == '-':
							A_dict[LN[1:]]['RN'] = A_dict[LN[1:]]['RN']+str('copy')+str(1)
						else:
							A_dict[LN]['RN'] = A_dict[LN]['RN']+str('copy')+str(1)
					if RN:
						if RN[0] == '-':
							A_dict[RN[1:]]['LN'] = A_dict[RN[1:]]['LN']+str('copy')+str(1)
						else:
							A_dict[RN]['LN'] = A_dict[RN]['LN']+str('copy')+str(1)

				if left_adj_at and right_adj_at:
					i = 3
					for index in Idx_dict[gene]:											#Updating list of FDs in case of weakly conserved context
						FD.append(gene+str('copy')+str(i))
						Idx_dict[gene+str('copy')+str(i)] = [index] 						#Updating D_dict and Idx_dict for each FD
						LNIdx, RNIdx = D_dict[index]['LNIdx'], D_dict[index]['RNIdx']
						if LNIdx:
							D_dict[LNIdx]['RN'] = D_dict[LNIdx]['RN']+str('copy')+str(i)
						if RNIdx:	
							D_dict[RNIdx]['LN'] = D_dict[RNIdx]['LN']+str('copy')+str(i)
						i += 1
				else:
					i = 2
					for index in Idx_dict[gene]:											#Updating list of FDs when context not conserved
						FD.append(gene+str('copy')+str(i))
						Idx_dict[gene+str('copy')+str(i)] = [index] 						#Updating D_dict and Idx_dict for each FD
						LNIdx, RNIdx = D_dict[index]['LNIdx'], D_dict[index]['RNIdx']
						if LNIdx:
							D_dict[LNIdx]['RN'] = D_dict[LNIdx]['RN']+str('copy')+str(i)
						if RNIdx:	
							D_dict[RNIdx]['LN'] = D_dict[RNIdx]['LN']+str('copy')+str(i)
						i += 1

	D_adj = []													#Form adjacency set for relabeled genome D'							
	for x in sorted((k,v) for (k,v) in D_dict.items()):
		left, right = None, None
		if x[1]['RNIdx']:
			RNIdx = x[1]['RNIdx']
			if D_dict[RNIdx]['LN'][0] == '-':
				left = (D_dict[RNIdx]['LN'][1:], 't')
			else:
				left = (D_dict[RNIdx]['LN'], 'h')
			if D_dict[x[0]]['RN'][0] == '-': 
				right = (D_dict[x[0]]['RN'][1:], 'h')
			else:
				right = (D_dict[x[0]]['RN'], 't')
			D_adj.append([left, right])

	A_adj = []													#Form adjacency set for relabeled genome A'
	for x in sorted((v['Idx'],k) for (k,v) in A_dict.items()):
		left, right = None, None
		if A_dict[x[1]]['RN']:
			RN = A_dict[x[1]]['RN']
			if A_dict[x[1]]['Sign'] == 'Neg':
				left = (x[1], 't')
			else:
				left = (x[1], 'h')
			if RN[0] == '-':
				right = (RN[1:], 'h')
			else:
				right = (RN, 't')
			A_adj.append([left, right])			
	for x in FD:												#Append (g_h g_t) adjacencies for each FD
		A_adj.append([(x, 'h'),(x, 't')])

	preserved_adj = [adj for adj in A_adj if adj in D_adj or list(reversed(adj)) in D_adj]		#Intersection of adjacency sets, A' and D'
	n_cuts = len(A_adj) - len(preserved_adj)
	n_joins = len(D_adj) - len(preserved_adj)
	n_duplicates = len(FD) + len(TD)

	distance = n_cuts + n_joins + n_duplicates + TD_from_arrays	#d_DSCJ(A,D) = |A'-D'| + |D'-A'| + n_d + TDA.


	print(distance)
	print(n_cuts, n_joins)			
	print(len(FD), len(TD))

	#outputfile = open("output1.txt", "w")
	#outputfile.write(str(A_adj))
	#outputfile.write("\n")
	#outputfile.write(str(D_adj))
					


#---------------------------------------------------------------------------
#Finds the median of all given genomes in the input file
def median(filename):
	string = open(filename, "r").read()
	substrings = string.split("\n")
	substrings = [line for line in substrings if line and line[0] != '#']
	genomes = []
	i = -1
	for line in substrings:
		if line[-1] not in {')','|'}:
			genomes.append([])
			i += 1
		elif line[-1] == '|':
			line = line.split(' ')
			line = [x for x in line if x != '|']
			genomes[i].append(line)	
		else:
			line = line.split(' ')
			line[-1] = line[0]
			genomes[i].append(line)

	TD_from_arrays = [0] * len(genomes)
	i = 0
	for genome in genomes:
		for chromosome in genome:				
			for gene_idx in range(len(chromosome)):
				while gene_idx < len(chromosome) - 1:
					if chromosome[gene_idx] == chromosome[gene_idx + 1]:
						del chromosome[gene_idx]
						TD_from_arrays[i] += 1
						print(TD_from_arrays[i])
					else:
						gene_idx += 1
		i += 1

	adj_list = []
	total_gene_list = []	
	for genome in genomes:
		adj_list.append(listAdj(genome))
		total_gene_list = list(set(total_gene_list + listGene(genome)))

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