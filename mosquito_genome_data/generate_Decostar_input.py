__author__  = "Cedric Chauve"
__version__ = "1.0"
__email__   = "cedric.chauve@sfu.ca"
__status__  = "Development"

# This script takes as input the results of the DeCoStar pipeline and formats them into 
# pairs of genomes observed along the branches of a phylogeny
#
# Usage:
# python generate_Decostar_input.py <orthogroups_file> <scaffolds_file> <assignment_2_chromosomes_file> <threshold> <prefix_output_files>
#
# The format of the files is described in the relevant section of the code

import sys

# Orthogroups: An orthogroup is a gene family
# Each line represents an orthology relation between an ancestral gene and a descendant gene
# ancestral_species descendant_species ancestral_gene descendant_gene gene_tree_id
# If a gene has no descendant or ancestor, the missing gene is labelled as NIL
# Extant genes are named <species_name>@<gene_name>
# Ancestral genes are named by the corresponding node in the gene tree: <gene_tree_id>|<node_id>
orthogroups=open(sys.argv[1],"r").readlines()

# Scaffolds: Gene orders / genome fragments
# Each line represents a gene: species_id fragment_id gene_name orientation
# Where gene name follows the same convention than the orthogroups file
scaffolds=open(sys.argv[2],"r").readlines()

# Scaffold assignment: 
# Each line represents a scaffold / genome fragment: scaffold_name scaffold_length prob_X prob_2L prob_2R prob_3L prob_3R
# Scaffold name: ANC_scaffold_id (ancestral species) or EXT_scaffold_id (extant species)
# Scaffold length = number of genes
# prob_AA = probability ot be on chromosomal arm AA (can be X, 2L, 2R, 3L, 3R)
assignment=open(sys.argv[3],"r").readlines()
assign_threshold=float(sys.argv[4]) # A scaffold is assigned to chromosome X if prob_X >= 0.8

# Auxiliary function to extract the gene name of gene, filtering out the species name in case of an extant gene
def filter_extant_gene(g):
    if '@' in g:
        gene=g.split('@')[1]
    else:
        gene=g
    return(gene)
        
# Initializing genomes and families
GENOMES={} # Gene orders, indexed by species id
# Each entry is a quadruplet (index,scaffold_name,gene_name,sign)
# where index is a global index used for all genomes and giving ordering relations for genes on the same fragment
PARENT={}  # Parent of a given gene
NB_DESC={} # Number of descendants of a given gene
for l in scaffolds:
    if l[0]!="#":
        l1=l.rstrip().split()
        species=l1[0]
        GENOMES[species]=[]
        gene=filter_extant_gene(l1[2])
        PARENT[gene]="NIL"
        NB_DESC[gene]=0
        
# Reading branches and gene families
BRANCHES=[] # List of branches, each branch is a pair (ancestrl_species_id, descendant_spcies_id)
for l in orthogroups:
    l1=l.rstrip().split()
    ancestor=l1[0]
    descendant=l1[1]
    if not (ancestor,descendant) in BRANCHES:
        BRANCHES.append((ancestor,descendant))
        
# Recording genes in genomes
i=1
for l in scaffolds:
    if l[0]!="#":
        l1=l.rstrip().split()
        species=l1[0]
        scf=l1[1]
        gene=filter_extant_gene(l1[2])
        sign=l1[3]
        GENOMES[species].append((i,scf,gene,sign))
        i+=1

# Reading scaffold assignment
ASSGN={} # Assignment of each scaffold
for l in assignment:
    if l[0]!="#":
        l1=l.rstrip().split()
        species=l1[0].split("_")[0].replace("ANC","").replace("EXT","")
        scf=l1[0].split("_")[1]
        if float(l1[3])>=assign_threshold:
            ASSGN[(species,scf)]="X" # X chromosome
        elif float(l1[3])<=(1.0-assign_threshold):
            ASSGN[(species,scf)]="A" # Autosomes
        else:
            ASSGN[(species,scf)]="U" # Unassigned
            

ORTHO_LIST={}
for (anc,desc) in BRANCHES:
   ORTHO_LIST[(anc,desc)]=[] 
for l in orthogroups:
    l1=l.rstrip().split()
    ancestor=l1[0]
    descendant=l1[1]
    ORTHO_LIST[(ancestor,descendant)].append(l)

# Main body: Generating one file per branch
# Format:
# First line = ancestral genome
# Second line = descendant genome
# Fragments are marked by #<fragment_id>_<assignment=A_or_X_or_U>
# X is for X chromosome, A for autosomes, U is for unassigned
for (anc,desc) in BRANCHES:
    output=open(sys.argv[5]+"_"+anc+"_"+desc,"w")
    A = GENOMES[anc]
    D = GENOMES[desc]

    for l in ORTHO_LIST[(anc,desc)]:
        l1=l.rstrip().split()
        for i in range(3,len(l1)-1):
            gene=filter_extant_gene(l1[i])
            if gene in PARENT.keys() and l1[2] in NB_DESC.keys():
                PARENT[gene]=l1[2]
                NB_DESC[l1[2]]+=1
                   
    pscf=""
    nbgA=0
    sizeA=0
    for (i,scf,gene,sign) in A:
        if scf!=pscf:
            sizeA+=1
            output.write("#"+scf+"_"+ASSGN[(anc,scf)]+" ")
        if gene in NB_DESC.keys() and NB_DESC[gene]>0:
            sizeA+=1
            output.write(sign+gene+" ")
            NB_DESC[gene]=0
            nbgA+=1        
        pscf=scf
    output.write("\n")
    
    pscf=""
    nbgD=0
    sizeD=0
    for (i,scf,gene,sign) in D:
        if scf!=pscf:
            sizeD+=1
            output.write("#"+scf+"_"+ASSGN[(desc,scf)]+" ")
        if gene in PARENT.keys() and PARENT[gene]!="NIL":
            sizeD+=1
            output.write(sign+PARENT[gene]+" ")
            nbgD+=1
        pscf=scf
    output.write("\n")
    output.flush()
