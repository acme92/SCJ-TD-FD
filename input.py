from sys import argv
import re

if len(argv) < 2:															#Takes file with genomes as argument in command line		
	print('Usage: python input.py <genome_file>')
	exit(1)

outputfile = open("inputfile.txt", "w")

filename = argv[1]
string = open(filename, "r").read()
substrings = string.split("\n")

i = 0

for line in substrings:
	if i == 2:
		break
	else:
		outputfile.write("\n")
		outputfile.write("Genome " + str(i) + ":\n")
		i += 1
		newline = line.replace('+', '')
		outputfile.write(newline)
		outputfile.write("|")
