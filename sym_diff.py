def negate(a):
    return '-' if a == '+' else '+'

def pairUp(a):
    set_a = set()
    for chrom in a:
        for i in range(len(chrom) - 1):
            set_a.add(chrom[i] + chrom[i+1])
    return set_a

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

def diff(A, B):
    B_dict = {i:None for i in duplicateOrder(B)}
    return [i for i in A if i not in B_dict]

def sym_diff(A, B):
    return len(diff(A, B)) + len(diff(B, A))
    
string = open("genome.txt", "r").read()
genomes = string.split("\n")
# remove empty strings
genomes = [i.strip().split("|") for i in genomes if len(i)]
S = [j.strip().split(" ") for j in genomes[0]]
D = [j.strip().split(" ") for j in genomes[1]]
A = pairUp(S)
B = pairUp(D)

#print( sorted([i for i in B]) ) 
#print( sorted(D) )
print(diff(A, B))
print(diff(B, A))
print(sym_diff(A, B))
