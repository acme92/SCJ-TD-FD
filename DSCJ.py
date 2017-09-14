__author__ = "Aniket Mane"
__email__ = "amane@sfu.ca"

"""
DSCJ.py implements the algorithm mentioned in the paper "A tractable variant of the Single Cut or Join distance with duplicated genes".
The following code consists of three main functions:
1.  Implementation: python DSCJ.py -d <inputfile>
    Given a genome A without duplicate genes and genome D having duplicate genes, 
    compute the SCJTDFD distance between the two.
2.  Implementation: python DSCJ.py -s <inputfile>
    Given a genome A without duplicate genes and genome D having duplicate genes, 
    compute one scenario that can transform A to D under the SCJTDFD model.
3.  Implementation: python DSCJ.py -m <inputfile>
    Given a set of genomes having duplicate genes,
    compute their median genome M having no duplicate genes.
"""

from sys import argv
import random
import networkx as nx
import matplotlib.pyplot as plt



#Using python 3 interpreter

#Auxiliary Functions
#---------------------------------------------------------------------------
#A genome is identified by its name and a list of chromosomes.
#Genome = [Name, [List of names of chromosomes]]. Data types: [String, [List of Strings]]

#Update genome list and initiate list for chromosome names. For the distance and scenario code, i = 0 => A, i = 1 => D
def update_genome_list(genome_list, genome_name, i):
    genome_list[i].append(genome_name)      
    genome_list[i].append([])
    return genome_list
#Update genome list by adding list of chromosome names for each genome.
def update_chr_list(genome_list, chr_name, i):
    genome_list[i][1].append(chr_name)
    return genome_list

#Return genome name.
def get_genome_name(genome):
    return genome[0]
#Return list of chromosomes in the genome.
def get_genome_chr_list(genome):
    return genome[1]


#Each chromosome has a type (linear or circular) and a set of ordered genes.
#Chromosome = [Name, Type, [List of ordered genes]]. Data types: [String, String, [List of Strings]]

#Update list of chromosomes by storing chromosome data.
def update_chr_data(chr_list, chromosome, chr_type, i):
    chr_list[i].append([chromosome[0], chr_type, chromosome[1:]])
    return chr_list

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


#A gene is characterised by its family, an orientation (BACKWARD or FORWARD) and its neighbors.
#The information about all the genes is captured in dictionaries (defined in the scenario function)
#Every entry in the gene dictionary has the following list of values for every key:
#[Gene posn, Orientation, Left neighbor, Right neighbor, Posn of left neighbor, Posn of right neighbor]
#Data type: [(Int,Int), String, String, String, (Int,Int), (Int,Int)]

#Each gene has an orientation. If it has orientation ('-g'/'g'), this function returns the opposite ('g'/'-g')  
def reverse(gene):
    if gene:
        return gene[1:] if gene[0] == '-' else str('-' + gene)
    else:
        return None

#Note: Since the distance function does not use specific data of the gene, 
#the functions to set/get gene information have been defined before the scenario function


#Check if genome A is trivial. Maintains a set of unique gene families. If same gene family repeats, the genome in non-trivial.
#The input is the set of 'seen' unique genes and a Boolean is_trivial (0 if trivial, 1 if not) and a chromosome.
#The output is the set 'seen' and the Boolean is_trivial after checking the chromosome for repeating genes.
def check_if_trivial(chromosome, seen, is_trivial):
    for gene in chromosome:
        if gene not in seen and reverse(gene) not in seen:
            seen.add(gene)
        else:
            is_trivial = False
        return (seen, is_trivial)   

#Forms a list of all gene families in input genome.
#The input is the genome (list of chromosomes).
#The output is a list of the names of gene families in the genome. 
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

#Forms a list of all gene families in input chromosome.
#The input is a chromosome (as defined above)
#The output is a set of the names of gene families in the chromosome.
def get_chr_gene_set(chromosome):
    gene_set = set()
    for gene in chromosome[2]:
        if gene[0] == '-':
            if reverse(gene) not in gene_set:
                gene_set.add(reverse(gene))
        else:
            if gene not in gene_set:
                gene_set.add(gene)
    return gene_set

#Retrieve genome data from text file
#Input is the list of strings read from the file (after removing comments and blank lines)
#Return a list of genomes, chromosomes (as defined above) and gene count for each genome.
def get_genome_data(string_list):
    genome_list = [[]*2 for i in range(2)]
    chr_list = [[]*3 for i in range(2)]
    seen = set()                                                            #Set of genes in A to check for trivialness
    is_trivial = True
    gene_count = [0,0]                                                      #Maintain a count of number of genes in A and D, respectively
     
    i = -1
    for line in string_list:
        if line[-1] not in {')','|'}:                                       #If not '|' or ')', then line provides genome name
            i += 1
            genome_list = update_genome_list(genome_list, line, i)
        elif line[-1] == '|':                                               #Linear chromosome
            line = line.split(' ')
            line = [x for x in line if x != '|']
            genome_list = update_chr_list(genome_list, line[0], i)
            chr_list = update_chr_data(chr_list, line, 'L', i)
            gene_count[i] += len(line[1:])
            if i == 0:                                                      #Check trivialness only for A.
                seen, trivial = check_if_trivial(line[1:], seen, is_trivial)            
        elif line[-1] == ')':                                               #Circular chromosome    
            line = line.split(' ')
            line[-1] = line[1]                                              #Replace the ')' by the first gene in the chromosome
            genome_list = update_chr_list(genome_list, line[0], i)
            chr_list = update_chr_data(chr_list, line, 'C', i)
            gene_count[i] += len(line[1:-1])
            if i == 0:                                                      #Check trivialness only for A.
                seen, trivial = check_if_trivial(line[1:-1], seen, is_trivial)
                
    if trivial == False:                                                    #If gene repeats, A is nontrivial. Terminate program.
        print("Error message: Ancestor genome must be trivial.")
        quit()

    if set(get_gene_list(chr_list[0])) != set(get_gene_list(chr_list[1])):  #If different set of gene families, terminate program.
        print("Error message: Ancestor and descendant genomes have different sets of gene families.")
        quit()      
    return (genome_list, chr_list, gene_count)  

#Counts TD from Arrays and removes them.
#Input is a genome (list of chromosomes)
#Output is the reduced genome and the count of TD from arrays
def remove_TDA(chr_list):
    TD_from_arrays = 0
    for chromosome in chr_list:
        gene_idx = 0
        if chromosome[1] == 'C':                                    #For circular genomes, len(list of genes) >/= 2,
            while gene_idx < len(chromosome[2]) - 1:                #as the first gene is appended after the last.
                if len(chromosome[2]) > 2:                           
                    if chromosome[2][gene_idx] == chromosome[2][gene_idx + 1]:
                        del chromosome[2][gene_idx]
                        TD_from_arrays += 1
                    else:
                        gene_idx += 1
                else:
                    gene_idx += 1
        elif chromosome[1] == 'L':
            for gene_idx in range(len(chromosome[2])):
                while gene_idx < len(chromosome[2]) - 1:
                    if chromosome[2][gene_idx] == chromosome[2][gene_idx + 1]:  #Remove tandem arrays, if any.
                        del chromosome[2][gene_idx] 
                        TD_from_arrays += 1
                    else:
                        gene_idx += 1
    return (chr_list, TD_from_arrays)

#Counts single-gene circular chromosomes (SGCC) that are duplicates and removes them
#Input is a genome (list of chromosomes)
#Output is the reduced genome and the count of SGCC that can be deleted
def remove_SGCC(chr_list):
    SGCC = []
    seen_set = set()
    i = 0
    count = 0
    for chromosome in chr_list:
        if chromosome[1] == 'C' and len(chromosome[2]) == 2:    #Check if chromosome is circular and of length 2 after removal of tandem arrays
            SGCC.append(i)                              
        else:
            seen_set = seen_set.union(get_chr_gene_set(chromosome)) #If not, add genes to a set of observed unique gene families not in SGCC
        i += 1
    for i in SGCC:
        gene = chr_list[i - count][2][0]                            #Gene found in SGCC 
        if gene in seen_set or reverse(gene) in seen_set:       #Check if gene found in set of observed unique gene families => duplicate
            del chr_list[i - count]
            count += 1
        else:                                                   #If not, => first instance of gene family, add to set of observed gene families 
            if gene[0] == '-':
                seen_set.add(gene[1:])
            else:
                seen_set.add(gene)  
    return chr_list, count          

#Forms adjacency list of input genome
#Obtains the genome as a list of chromosomes
#Iterates through each chromosome to form a list of adjacencies of neighboring genes
def get_adj_list(chr_list):
    adj_list = []                                                           #List of adjacencies where 
    for chromosome in chr_list:                                             #every adjacency is of the format:
        for gene_idx in range(len(chromosome[2]) - 1):                      #[(g1,'h'/'t'),(g2,'h'/'t')]
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



#Distance function
#---------------------------------------------------------------------------
#Finds distance between the two given genomes in the input file
def distance(filename):
    string = open(filename, "r").read()
    string_list = string.split("\n")
    string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.

    genome_list, chr_list, gene_count = get_genome_data(string_list)

    A = chr_list[0]         
    D = chr_list[1]
    
    D, TD_from_arrays = remove_TDA(D)                                       #Remove tandem arrays and duplicate single gene circular chromosomes
    D, SGCC = remove_SGCC(D)
    
    gene_count[1] -= TD_from_arrays + SGCC                                  #Number of genes in D after removing tandem arrays and SGCC
    n_duplicates = gene_count[1] - gene_count[0]                            #Number of genes in D - number of genes in A

    A_adj = get_adj_list(A)
    D_adj = get_adj_list(D)

    preserved_adj = [adj for adj in A_adj if adj in D_adj or list(reversed(adj)) in D_adj]      #Intersection of adjacency sets, A and D
    n_cuts = len(A_adj) - len(preserved_adj)                                #Adjacencies seen in A but NOT preserved in D
    n_joins = len(D_adj) - len(preserved_adj)                               #Adjacencies seen in D but NOT preserved from A

    d_DSCJ = n_cuts + n_joins + 2*n_duplicates + TD_from_arrays + SGCC      #d_DSCJ(A,D) = |A-D| + |D-A| + 2*n_d + TDA.

    print(d_DSCJ)
    print(n_cuts)
    print(n_joins)
    print(n_duplicates)
    print(TD_from_arrays)
    print(SGCC)



#Functions required for Scenario code
#---------------------------------------------------------------------------
#As discussed earlier, the following information is required for every gene:
#Gene family name, Orientation and its left and right neighbors.

#Assign gene family name
def set_gene_family(gene):
    if gene[0] == '-':
        return gene[1:]
    else:
        return gene
#Assign orientation of gene
def set_gene_orientation(gene):
    if gene[0] == '-':
        return 'BACKWARD'
    else:
        return 'FORWARD'
#Using the position of gene in the chromosome, find its left neighbor
def set_left_neighbor(chromosome, j):
    if j == 0:
        if get_chr_type(chromosome) == 'L':
            return None
        elif get_chr_type(chromosome) == 'C':
            return get_gene_by_posn(chromosome, -1)
    else:
        return get_gene_by_posn(chromosome, j-1)
#Using the position of gene in the chromosome, find its right neighbor
def set_right_neighbor(chromosome, j):
    if j == len(chromosome[2]) - 1:
        if get_chr_type(chromosome) == 'L':
            return None
        elif get_chr_type(chromosome) == 'C':
            return get_gene_by_posn(chromosome, 0)
    else:
        return get_gene_by_posn(chromosome, j+1)
#Using the position of gene in the genome, find the position of its left neighbor
def left_neighbor_posn(chromosome, posn):
    i, j = posn[0], posn[1]
    if j == 0:
        if get_chr_type(chromosome) == 'L':
            return None
        elif get_chr_type(chromosome) == 'C':
            return (i, len(chromosome[2]) - 1)
    else:
        return (i, j-1)
#Using the position of gene in the genome, find the position its right neighbor
def right_neighbor_posn(chromosome, posn):
    i, j = posn[0], posn[1]
    if j == len(chromosome[2]) - 1:
        if get_chr_type(chromosome) == 'L':
            return None
        elif get_chr_type(chromosome) == 'C':
            return (i, 0)
    else:
        return (i, j+1)

#The following functions return the features of the gene and information about its neighborhood from the dictionary entry
def get_gene_posn(g):
    return(g[0])
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

#Updates A_dict by relabeling gene matched with gene in A
def updateA(gene, A_dict):
    LN, RN = get_left_neighbor(A_dict[gene]), get_right_neighbor(A_dict[gene])                  #Left and right neighbors in A
    if LN:                                                                  #If Left neighbor exists, update corresponding entry.
        if LN[0] == '-':                                                    #(If LN exists, its right neighbor will be the gene itself, 
            A_dict[LN[1:]][3] = A_dict[LN[1:]][3]+str('copy')+str(1)        #so needs to be relabeled) 
        else: 
            A_dict[LN][3] = A_dict[LN][3]+str('copy')+str(1) 
    if RN:                                                                  #If Right neighbor exists, update corresponding entry.
        if RN[0] == '-':                                                    #(If RN exists, its left neighbor will be the gene itself,  
            A_dict[RN[1:]][2] = A_dict[RN[1:]][2]+str('copy')+str(1)        #so needs to be relabeled)
        else: 
            A_dict[RN][2] = A_dict[RN][2]+str('copy')+str(1)
    A_dict[gene+str('copy')+str(1)] = A_dict[gene]
    del A_dict[gene]
    return A_dict

#Updates Idx_dict and D_dict by relabeling gene matched with gene in A
def updateD(gene, Idx_dict, A_dict, D_dict, posn, i):
    Idx_dict[gene+str('copy')+str(i)] = [posn]                              #Update Idx_dict by introducing entry gcopy'i'
    Idx_dict[gene].remove(posn)                                             #and removing corresponding index from Idx_dict[g]

    LNIdx, RNIdx = get_left_neighbor_posn(D_dict[posn]), get_right_neighbor_posn(D_dict[posn])          #Left and right neighbors in D
    if LNIdx: D_dict[LNIdx][3] = D_dict[LNIdx][3]+str('copy')+str(i)        #If LN exists, relabel its right neighbor
    if RNIdx: D_dict[RNIdx][2] = D_dict[RNIdx][2]+str('copy')+str(i)        #If RN exists, relabel its left neighbor
        
    return Idx_dict, D_dict

#Updates list of Floating Duplicates and the labels of genes matched with floating duplicates in Idx_dict and D_dict
def updateFD(FD, i, gene, Idx_dict, D_dict):
    for posn in Idx_dict[gene]:
        FD.append(gene+str('copy')+str(i))                                  
        Idx_dict[gene+str('copy')+str(i)] = [posn]
        LNIdx, RNIdx = get_left_neighbor_posn(D_dict[posn]), get_right_neighbor_posn(D_dict[posn])          #Left and right neighbors in D
        if LNIdx: D_dict[LNIdx][3] = D_dict[LNIdx][3]+str('copy')+str(i)        #If LN exists, relabel its right neighbor
        if RNIdx: D_dict[RNIdx][2] = D_dict[RNIdx][2]+str('copy')+str(i)        #If RN exists, relabel its left neighbor
        i += 1
    return(FD, Idx_dict, D_dict)

#Scenario function
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
#   If g is from nontrivial family:
#       Check for context strongly conserved. Complexity O(copy number of g) 
#           If instance found: 
#               Update all the dictionaries.
#       If instance found: 
#           Update list of FDs. Complexity O(copy number of g)
#       
#       If not strongly conserved:
#           Check for conservation of adjacencies on either side. Complexity O(copy number of g)
#           Check if context weakly conserved: Complexity O(copy number of g)
#               If instance found:
#                   Update all the dictionaries.
#               Else check only if left adj conserved:  
#                   If so, update all the dictionaries accordingly.
#               Else check only if right adj conserved:  
#                   If so, update all the dictionaries accordingly.
#               Else: 
#                   Match with the first copy of g in D. Update all the dictionaries.
#
#           If weakly conserved:
#               Update list of FDs for remaining copies. Complexity O(copy number of g) 
#           Else not conserved:
#               Update list of FDs for remaining copies. Complexity O(copy number of g)
#
#Create an adjacency set of D using the dictionary. Complexity O(n_D)
#Create an adjacency set of A using the dictionary. Also add FDs. Complexity O(n_D)
#
#Find the distance: |A-D| + |D-A| + n_d + TDA

def scenario(filename):
    string = open(filename, "r").read()
    string_list = string.split("\n")
    string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.

    genome_list, chr_list, gene_count = get_genome_data(string_list)

    A = chr_list[0]         
    D = chr_list[1]
    
    D, TD_from_arrays = remove_TDA(D)                       #Remove tandem arrays and duplicate single gene circular chromosomes
    D, SGCC = remove_SGCC(D)
    
    A_dict = {}     #Dictionary for A. KEY = Gene Family Name. VALUE = (Sign, Left neighbor, Left neighbor index, Right neighbor, Right neighbor index)
    for i in range(len(A)):     
        if get_chr_type(A[i]) == 'C':
            A[i][2] = A[i][2][:-1]
        chromosome = get_chr_gene_list(A[i])
        for j in range(len(chromosome)):
            A_dict[set_gene_family(chromosome[j])]= [   (i,j),
                                                        set_gene_orientation(chromosome[j]),
                                                        set_left_neighbor(A[i],j),
                                                        set_right_neighbor(A[i],j)
                                                    ]

    Idx_dict = {}   #Dictionary for indices. KEY = Gene family, VALUE = List of positions of g in D. 
    D_dict = {}     #Dictionary for D. KEY = Index. VALUE = (Sign, Left neighbor, Left neighbor index, Right neighbor, Right neighbor index)    
    for i in range(len(D)):
        if get_chr_type(D[i]) == 'C':
            D[i][2] = D[i][2][:-1]
        chromosome = get_chr_gene_list(D[i])
        for j in range(len(chromosome)):
            try:
                Idx_dict[set_gene_family(chromosome[j])].append((i,j))
            except KeyError:
                Idx_dict[set_gene_family(chromosome[j])] = [(i,j)]
            D_dict[(i,j)] = [   (i,j), 
                                set_gene_orientation(chromosome[j]),
                                set_left_neighbor(D[i],j),
                                set_right_neighbor(D[i],j),
                                left_neighbor_posn(D[i],(i,j)),
                                right_neighbor_posn(D[i],(i,j))     
                            ]

    coords = []                                     #List of co-ordinates of genes in A, shuffled for randomness                            
    for i in range(len(A)):                                             
        for j in range(len(A[i][2])):
            coords.append([i,j])        
    random.shuffle(coords)

    FD = []                                         #List of FD created
    TD = []                                         #List of TD created

    for i, j in coords:                             #Iterate through all gene coordinates 
        gene = get_gene_by_posn(A[i],j)
        gene = set_gene_family(gene)
        if len(Idx_dict[gene]) > 1:                 #Enter block only if copy number of gene > 1
            #CASE 1: Context strongly conserved
            strong_context = 0                      
            weak_context = 0 
            for posn in Idx_dict[gene]:             #Iterate though all positions of the gene in D
                if strong_context == 0:
                    #If orientation of the original gene in A and the gene in current position in D matches,
                    #compare left (right) neighbor of the gene in A with left (right) neighbor of current gene in D 
                    if get_gene_orientation(A_dict[gene]) == get_gene_orientation(D_dict[posn]):    
                        if (get_left_neighbor(A_dict[gene]) and get_left_neighbor(A_dict[gene]) == get_left_neighbor(D_dict[posn]) and
                                get_right_neighbor(A_dict[gene]) and get_right_neighbor(A_dict[gene]) == get_right_neighbor(D_dict[posn])):
                            strong_context = 1
                            A_dict = updateA(gene, A_dict)                                  
                            Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, posn, 1)
                    #If orientation of the original gene in A and the gene in current position in D is opposite,
                    #compare left (right) neighbor of the gene in A with right (left) neighbor of current gene in D         
                    else:
                        if (get_left_neighbor(A_dict[gene]) and get_left_neighbor(A_dict[gene]) == reverse(get_right_neighbor(D_dict[posn])) and
                                get_right_neighbor(A_dict[gene]) and get_right_neighbor(A_dict[gene])==reverse(get_left_neighbor(D_dict[posn]))):
                            strong_context = 1
                            A_dict = updateA(gene, A_dict)
                            Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, posn, 1)

            if strong_context == 1:                 #Strong context has been found, remaining duplicates matched with FDs
                FD, Idx_dict, D_dict = updateFD(FD, 2, gene, Idx_dict, D_dict)

            if strong_context == 0:                 #Strong context not found
                left_adj_at, right_adj_at = None, None      #If adjacencies of the gene have been preserved in D, 
                for posn in Idx_dict[gene]:                 #find the position of the copy where they have been preserved
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
                if left_adj_at and right_adj_at:    #If both adjacencies preserved (and strong context not found) => weak context found 
                    weak_context = 1
                    TD.append(gene)
                    #Update Idx_dict and D by creating an entry for copy1 and copy2 and relabeling accordingly
                    Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, left_adj_at, 1)      
                    Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, right_adj_at, 2)

                    LN, RN = get_left_neighbor(A_dict[gene]), get_right_neighbor(A_dict[gene])      #Get information about gene in A
                    Posn, Orientation = get_gene_posn(A_dict[gene]), get_gene_orientation(A_dict[gene])

                    A_dict[gene+str('copy')+str(1)] = [Posn, Orientation, LN, RN]       #Create an entry for copy1 in A_dict
                    A_dict[gene+str('copy')+str(2)] = [Posn, Orientation, LN, RN]       #Create an entry for copy2 (tandem duplicate) in A_dict

                    if LN[0] == '-':                                                    #Relabel gene to genecopy1 and genecopy2 for relevant fields.
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
                    del A_dict[gene]                                                    #Delete entry corresponding to original gene in A

                #Case 3: Context not conserved, but at least one adjacency conserved.
                elif left_adj_at and not right_adj_at:
                    Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, left_adj_at, 1)  #Updating D by relabeling left match                        
                    A_dict = updateA(gene, A_dict)                                              #Updating A_dict

                elif right_adj_at and not left_adj_at:
                    Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, right_adj_at, 1) #Updating D by relabeling right match                       
                    A_dict = updateA(gene, A_dict)                                              #Updating A_dict

                #CASE 4: Context not conserved. No adjacency conserved.         
                else:
                    Idx_dict, D_dict = updateD(gene, Idx_dict, A_dict, D_dict, posn, 1) #Updating D by relabeling last copy (to original gene in A)                     
                    A_dict = updateA(gene, A_dict)                                      #Updating A_dict

                if left_adj_at and right_adj_at:                                    #If weakly conserved context, remaining genes matched with FDs
                    FD, Idx_dict, D_dict = updateFD(FD, 3, gene, Idx_dict, D_dict)                                  
                else:                                                               #If context not conserved, remaining genes matched with FDs
                    FD, Idx_dict, D_dict = updateFD(FD, 2, gene, Idx_dict, D_dict)                              

    #All Cases covered. Create adjacency lists from dictionaries

    #Logic:
    #Iterate through sorted dictionary. 
    #If right neighbor (RN) exists, adjacency is [left neighbor of RN, RN] (head or tail chosen according to orientation)
    D_adj = []                                                  #Form adjacency set for relabeled genome D'                         
    for x in sorted((k,v) for (k,v) in D_dict.items()):
        left, right = None, None                                #All adjacencies of the format: [(g1,'h'/'t'),(g2,'h'/'t')]
        if x[1][5]:
            RNIdx = x[1][5]
            if D_dict[RNIdx][2][0] == '-':                      #Left extremity of adjacency
                left = (D_dict[RNIdx][2][1:], 't')
            else:
                left = (D_dict[RNIdx][2], 'h')
            if D_dict[x[0]][3][0] == '-':                       #Right extremity of adjacency
                right = (D_dict[x[0]][3][1:], 'h')
            else:
                right = (D_dict[x[0]][3], 't')
            D_adj.append([left, right])

    #print(D_adj)
    A_adj = []                                                  #Form adjacency set for relabeled genome A'
    for x in sorted((v[0],k) for (k,v) in A_dict.items()):
        left, right = None, None                                #All adjacencies of the format: [(g1,'h'/'t'),(g2,'h'/'t')] 
        if A_dict[x[1]][3]:
            RN = A_dict[x[1]][3]
            if A_dict[x[1]][1] == 'BACKWARD':                   #Left extremity of adjacency
                left = (x[1], 't')
            else:
                left = (x[1], 'h')
            if RN[0] == '-':                                    #Right extremity of adjacency
                right = (RN[1:], 'h')
            else:
                right = (RN, 't')
            A_adj.append([left, right])         
    for x in FD:                                                #Append (g_h g_t) adjacencies for each FD
        A_adj.append([(x, 'h'),(x, 't')])   

    preserved_adj = [adj for adj in A_adj if adj in D_adj or list(reversed(adj)) in D_adj]      #Intersection of adjacency sets, A' and D'
    n_cuts = len(A_adj) - len(preserved_adj)                    #Adjacencies seen in A' but NOT preserved in D'
    n_joins = len(D_adj) - len(preserved_adj)                   #Adjacencies seen in D' but NOT preserved from A'
    n_duplicates = len(FD) + len(TD)

    distance = n_cuts + n_joins + n_duplicates + TD_from_arrays + SGCC  fget#d_DSCJ(A,D) = |A'-D'| + |D'-A'| + n_d + TDA + SGCC.

    print(distance)
    print(n_cuts, n_joins)          
    print(len(FD), len(TD), TD_from_arrays, SGCC)



#Functions required for Median code
#---------------------------------------------------------------------------

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

#Median function
#---------------------------------------------------------------------------
#Finds the median of all given genomes in the input file
def median(filename):
    string = open(filename, "r").read()
    string_list = string.split("\n")
    string_list = [line for line in string_list if line and line[0] != '#']
    genome_list = []
    chr_list = []

    i = -1
    for line in string_list:
        if line and line[-1] not in {')','|'}:
            i += 1
            genome_list.append([])
            chr_list.append([])
            genome_list = update_genome_list(genome_list, line, i)
        elif line[-1] == '|':
            line = line.split(' ')
            line = [x for x in line if x != '|']
            genome_list = update_chr_list(genome_list, line[0], i)
            chr_list = update_chr_data(chr_list, line, 'L', i)
        else:
            line = line.split(' ')
            line[-1] = line[1]
            genome_list = update_chr_list(genome_list, line[0], i)
            chr_list = update_chr_data(chr_list, line, 'C', i)      

    TD_from_arrays = [0] * len(genome_list)
    SGCC = [0] * len(genome_list)
    i = 0
    for genome in chr_list:
        genome, TD_from_arrays[i] = remove_TDA(genome)
        genome, SGCC[i] = remove_SGCC(genome)
        i += 1

    adj_list = [] #list of adjs per genome
    total_gene_list = []    
    for genome in chr_list:
        adj_list.append(get_adj_list(genome))
        total_gene_list = list(set(total_gene_list + get_gene_list(genome)))

    total_adj_list = [] #set of adj of all genomes
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
if len(argv) < 3:                                                           #Takes file with genomes as argument in command line        
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