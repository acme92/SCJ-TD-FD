from sys import argv
from functools import reduce
import os
import functools
import re

#Using python 3 interpreter

#Function definitions
#---------------------------------------------------------------------------
#Lists the extremities of each fragment of input genome
def listExtremities(genome):
	fragment_extremities = []
	for fragment in genome:
		if fragment[0][0] == '+':
			fragment_extremities.append(fragment[0][1:] + '_t')
		else:
			fragment_extremities.append(fragment[0][1:] + '_h')
		if fragment[-1][0] == '+':
			fragment_extremities.append(fragment[-1][1:] + '_h')
		else:
			fragment_extremities.append(fragment[-1][1:] + '_t')			
	return fragment_extremities	

#Forms adjacency list of input genome
def listAdj(genome):
	adj_list = []
	for fragment in genome:
		for idx in range(len(fragment) - 1):
			if fragment[idx][0] == '+':
				left = fragment[idx][1:] + '_h'
			else:
				left = fragment[idx][1:] + '_t'
			if fragment[idx + 1][0] == '+':
				right = fragment[idx + 1][1:] + '_t'
			else:
				right = fragment[idx + 1][1:] + '_h'
			adj_list.append([left, right])
	return adj_list


#Main function
#---------------------------------------------------------------------------
def DSCJ(foldername):
	outputfile = open("results_1.txt", "w")
	outputfile.write("filename \tnd \tcuts \tjoins \tTDA \td_DSCJ \tweak_cuts \tweak_joins \n")
	file_list = os.listdir(foldername)

	file_number = 0

	for filename in file_list:
		string = open(os.path.join(foldername, filename), "r").read()

		file_number += 1
		print(file_number,": ", filename)

		#Preprocessing
		#-------------------------------------------------------------------
		genomes = string.split("\n")										#Splits file into genomes and then splits genomes into individual fragments
		#genomes_test = [re.findall("#\d+_[AUX]", i) for i in genomes if len(i)]

		frag_A = []
		frag_U = []
		frag_X = []
		for seq in genomes:
			if len(seq) > 0:
				frag_A_curr = re.findall("#\d+_A[0-9 +-|]*", seq)
				frag_A.append(frag_A_curr)

				frag_U_curr = re.findall("\d+_U[0-9 +-|]*", seq)
				frag_U.append(frag_U_curr)

				frag_X_curr = re.findall("\d+_X[0-9 +-|]*", seq)
				frag_X.append(frag_X_curr)

		genomeA_A = [j.strip().split(" ") for j in frag_A[0]]								
		genomeA_A = [fragment[1:] for fragment in genomeA_A]
		genomeA_A = [fragment for fragment in genomeA_A if fragment]

		genomeD_A = [j.strip().split(" ") for j in frag_A[1]]
		genomeD_A = [fragment[1:] for fragment in genomeD_A]
		genomeD_A = [fragment for fragment in genomeD_A if fragment]

		genomeA_U = [j.strip().split(" ") for j in frag_U[0]]								
		genomeA_U = [fragment[1:] for fragment in genomeA_U]
		genomeA_U = [fragment for fragment in genomeA_U if fragment]

		genomeD_U = [j.strip().split(" ") for j in frag_U[1]]
		genomeD_U = [fragment[1:] for fragment in genomeD_U]
		genomeD_U = [fragment for fragment in genomeD_U if fragment]

		genomeA_X = [j.strip().split(" ") for j in frag_X[0]]								
		genomeA_X = [fragment[1:] for fragment in genomeA_X]
		genomeA_X = [fragment for fragment in genomeA_X if fragment]

		genomeD_X = [j.strip().split(" ") for j in frag_X[1]]
		genomeD_X = [fragment[1:] for fragment in genomeD_X]
		genomeD_X = [fragment for fragment in genomeD_X if fragment]

		TD_from_arrays = 0
		A = [genomeA_A, genomeA_U, genomeA_X]
		D = [genomeD_A, genomeD_U, genomeD_X]

		for chr_type in D:				
			for chr_idx in range(len(chr_type)):
				gene_idx = 0
				while gene_idx < len(chr_type[chr_idx]) - 1:
					if chr_type[chr_idx][gene_idx] == chr_type[chr_idx][gene_idx + 1]:
						
						print((chr_type[chr_idx]))
						del chr_type[chr_idx][gene_idx]
						print((chr_type[chr_idx]))
						TD_from_arrays += 1
					else:
						gene_idx = gene_idx + 1
						
		#Main Program
		#-------------------------------------------------------------------

		nA_A_genes = sum(len(fragment) for fragment in A[0])				#counts #genes in A and D' and #duplicates accordingly
		nD_A_genes = sum(len(fragment) for fragment in D[0])
		nA_U_genes = sum(len(fragment) for fragment in A[1])
		nD_U_genes = sum(len(fragment) for fragment in D[1])
		nA_X_genes = sum(len(fragment) for fragment in A[2])
		nD_X_genes = sum(len(fragment) for fragment in D[2])
		nA_genes = nA_A_genes + nA_U_genes + nA_X_genes
		nD_genes = nD_A_genes + nD_U_genes + nD_X_genes
		n_duplicates = nD_genes - nA_genes
		print(n_duplicates)

		

		A_A_extremities = listExtremities(A[0])								#lists extremities of A and D'
		A_U_extremities = listExtremities(A[1])
		A_X_extremities = listExtremities(A[2])
		D_A_extremities = listExtremities(D[0])
		D_U_extremities = listExtremities(D[1])
		D_X_extremities = listExtremities(D[2])

		A_A_adjset = listAdj(A[0])											#lists adjacency sets of A and D'
		A_U_adjset = listAdj(A[1])
		A_X_adjset = listAdj(A[2])
		D_A_adjset = listAdj(D[0])
		D_U_adjset = listAdj(D[1])
		D_X_adjset = listAdj(D[2])

		weak_cuts = []
		weak_joins = []

		for adjacency in A_adjset:											#lists weak cuts
			if adjacency[0] in D_extremities and adjacency[1] in D_extremities:
				weak_cuts.append(adjacency)

		for adjacency in D_adjset:											#lists weak joins
			if adjacency[0] in A_extremities and adjacency[1] in A_extremities:
				if not adjacency in weak_joins:
					weak_joins.append(adjacency)

		preserved_adj = [adj for adj in A_adjset if adj in D_adjset or list(reversed(adj)) in D_adjset]		#intersection of adjacency sets, A and D'
		n_cuts = len(A_adjset) - len(preserved_adj)
		n_joins = len(D_adjset) - len(preserved_adj)

		d_DSCJ = len(A_adjset) + len(D_adjset) - 2*len(preserved_adj) + 2*n_duplicates + TD_from_arrays

		outputfile.write(filename + "\t" + str(n_duplicates) + "\t" + str(n_cuts) + "\t" + str(n_joins) + "\t" + str(TD_from_arrays) + "\t" + str(d_DSCJ) + "\t" + str(len(weak_cuts)) + "\t\t" + str(len(weak_joins)) + "\n")
		
	outputfile.close()

#Input
#---------------------------------------------------------------------------
if len(argv) < 2:															#Takes file with genomes as argument in command line
	print('Usage: python d_SCjDup.py <genome_file>')
	exit(1)

foldername = argv[1]
DSCJ(foldername)