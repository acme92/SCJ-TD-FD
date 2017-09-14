A = [['Chr1','L',['a','-b','c','d']]]
A.append(['Chr2','C',['-h','-h']])
A.append(['Chr3','C',['b','-d','f','-h']])
A.append(['Chr4','C',['b','b']])
A.append(['Chr5','C',['-g','-g']])
A.append(['Chr6','L',['e','f']])
A.append(['Chr7','C',['g','g']])




#print('Before:', A)

def reverse(gene):
    if gene:
        return gene[1:] if gene[0] == '-' else str('-' + gene)
    else:
        return None

def get_chr_gene_set(chromosome):
    gene_set = set()
    for gene in chromosome[2]:
        if gene[0] == '-':
            if reverse(gene) not in gene_set:
                gene_set.add(reverse(gene))
        else:
            if gene not in gene_set:
                gene_set.add(gene)
    #print(gene_set)
    return gene_set

def remove_FD(chr_list):
    SGCC = []
    seen_set = set()
    i = 0
    SFD = 0
    for chromosome in chr_list:
        if chromosome[1] == 'C' and len(chromosome[2]) == 2:
            SGCC.append(i)
        else:
            print(chromosome)
            print(get_chr_gene_set(chromosome))
            seen_set = seen_set.union(get_chr_gene_set(chromosome))
            print(seen_set)
        i += 1
    
    print(SGCC)
    for i in SGCC:
        #print(i)
        gene = chr_list[i - SFD][2][0]
        if gene in seen_set or reverse(gene) in seen_set:
            #print(gene)
            #print(seen_set)
            del chr_list[i - SFD]
            SFD += 1
        else:
            #print(gene)
            if gene[0] == '-':
                seen_set.add(gene[1:])
            else:
                seen_set.add(gene)  
    print(seen_set)
    return chr_list, SFD
    
print(remove_FD(A))