#!/usr/bin/env python
import argparse
import sys

parser=argparse.ArgumentParser(description='Add compartment size file to *eigen1.sorted.bedGraph')
parser.add_argument('-i', help= 'input file (*eigen1.sorted.bedGraph)', type=str, required=True)
args=parser.parse_args()

	
def get_bin_size(d):
		line = d[0].strip().split('\t')
		start = int(line[1])
		end = int(line[2])
		b = (end - start) + 1	
		return b


def main():

	with open(args.i, 'r') as f, open(args.i[:-9] + '.size.bedGraph', 'w') as OUT:
		data = f.readlines()
		bin_size = get_bin_size(data)
		i = 0
		while i < len(data):
			splitline= data[i].strip().split('\t')
			chrom = splitline[0]
			eigen = float(splitline[3])
			switch = False
			compartment_bins = 1
			while not switch:
				# Break loop if on last compartment
				if i+compartment_bins > len(data)-1:
					break
				nextline = data[i+compartment_bins].strip().split('\t')
				nexteigen = float(nextline[3])
				# Break loop if on new chromosome
				nextchrom = nextline[0]
				if chrom != nextchrom:
					break
				# Check if compartment switch
				if eigen*nexteigen < 0 :
					switch = True
				else:
					compartment_bins += 1
			# Calc size in bp
			compartment_size = bin_size * compartment_bins
			for j in range(compartment_bins):
				printline = data[i+j].strip()
				OUT.write(printline + '\t' + str(compartment_size) + '\n')

			i = i+compartment_bins


if __name__ == '__main__':
	main()
