#!/usr/bin/env python
import argparse
import sys

parser=argparse.ArgumentParser(description='Add compartment size file to *subcompartments_40kb.bed')
parser.add_argument('-i', help= 'input file (*subcompartments_40kb.bed)', type=str, required=True)
args=parser.parse_args()

	
def get_bin_size(d):
		line = d[0].strip().split('\t')
		start = int(line[1])
		end = int(line[2])
		b = (end - start)
		return b


def main():

	A1 = open('A1_size.txt', 'w')
	A2 = open('A2_size.txt', 'w')
	B1 = open('B1_size.txt', 'w')
	B2 = open('B2_size.txt', 'w')
	B3 = open('B3_size.txt', 'w')

	with open(args.i, 'r') as f, open(args.i[:-4] + '.size.bed', 'w') as OUT:
		f.readline()
		data = f.readlines()
		bin_size = get_bin_size(data)
		i = 0
		while i < len(data):
			splitline= data[i].strip().split('\t')
			chrom = splitline[0]
			sub = splitline[3]
			switch = False
			compartment_bins = 1
			while not switch:
				# Check if NA
				if data[i + (compartment_bins - 1)].strip().split('\t')[3] == 'NA':				
					break
				# Break loop if on last compartment
				if i+compartment_bins > len(data)-1:
					break
				nextline = data[i+compartment_bins].strip().split('\t')
				# Break loop if next sub is NA
				if nextline[3] == 'NA':
					break
				nextsub = nextline[3]
				# Break loop if on new chromosome
				nextchrom = nextline[0]
				if chrom != nextchrom:
					break
				# Check if subcompartment switch
				if sub != nextsub:
					switch = True
				else:
					compartment_bins += 1
				
			# Calculate size in bp
			compartment_size = bin_size * compartment_bins
			if sub == 'A1':
				A1.write(str(compartment_size) + '\n')
			elif sub == 'A2':
				A2.write(str(compartment_size) + '\n')
			elif sub == 'B1':
				B1.write(str(compartment_size) + '\n')
			elif sub == 'B2':
				B2.write(str(compartment_size) + '\n')
			elif sub == 'B3':
				B3.write(str(compartment_size) + '\n')

			for j in range(compartment_bins):
				if data[i+j].strip().split('\t')[3] == 'NA':
					printline = data[i+j].strip()
					OUT.write(printline + '\tNA\n')
				else:
					printline = data[i+j].strip()
					OUT.write(printline + '\t' + str(compartment_size) + '\n')

			i = i+compartment_bins


if __name__ == '__main__':
	main()
