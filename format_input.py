from sys import argv
import re

if len(argv) < 2:
	print('Usage: python format_input.py <input_file> <output_file>')
	exit(1)

filename = argv[1]
string = open(filename, "r").read()
substrings = string.split("\n")

outputfile = open(argv[2], "w")

i = 0
empty_chr = []	#Segments without genes

for line in substrings:
	if i == 2:		#Only two genomes allowed (A and D)	
		break
	else:
		outputfile.write("\n")
		outputfile.write("Genome " + str(i))
		line = line.replace('+', '')
		line += ' '
		chr_list = re.findall("#\d+_[A/U/X][0-9 +-|]*", line)	#Pattern followed by all segments (including the name and type A/U/X)
		for chromosome in chr_list:
			elem = chromosome.split(' ')
			if len(elem[1:-1]) > 0:				#Check if at least one gene exists in the segment
				outputfile.write("\n")
				outputfile.write(chromosome[1:] + "|")	#Assuming that all are linear (since no indication to the contrary)
			else:
				empty_chr.append(elem[0])
		i += 1		

		
		
		
		
		
