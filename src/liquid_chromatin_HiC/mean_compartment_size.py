#!/usr/bin/env python
import argparse
import sys
import numpy as np

	
def get_bin_size(d):
		line = d[0].strip().split('\t')
		start = int(line[1])
		end = int(line[2])
		b = (end - start)
		return b


def main():
	parser=argparse.ArgumentParser(description='Add compartment size file to *eigen1.sorted.bedGraph')
	parser.add_argument('-i', help= 'input file (*eigen1.sorted.bedGraph)', type=str, required=True)
	parser.add_argument('-c', help= 'calculate for specific chromosome (ex. chr19)', type=str, default='genome')
	args=parser.parse_args()

	with open(args.i, 'r') as f:
		data = f.readlines()
		bin_size = get_bin_size(data)
		i = 0
		A_sizes = []
		B_sizes = []
		if args.c == 'genome':
			while i < len(data):
				splitline= data[i].strip().split('\t')
				chrom = splitline[0]
				eigen = float(splitline[3])
				switch = False
				compartment_bins = 1
				while not switch:
					# Check if NA
					if data[i + (compartment_bins - 1)].strip().split('\t')[3] == 'nan':				
						break
					# Break loop if on last compartment
					if i+compartment_bins > len(data)-1:
						break
					nextline = data[i+compartment_bins].strip().split('\t')
					# Break loop in next eigen is NA
					if nextline[3] == 'nan':
						break
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
					
				# Calculate size in bp
				compartment_size = bin_size * compartment_bins
				if eigen > 0:
					A_sizes.append(compartment_size)
				elif eigen < 0: 
					B_sizes.append(compartment_size)

				i = i+compartment_bins

		else:
			while i < len(data):
				splitline= data[i].strip().split('\t')
				chrom = splitline[0]
				if chrom == args.c:
					eigen = float(splitline[3])
					switch = False
					compartment_bins = 1
					while not switch:
						# Check if NA
						if data[i + (compartment_bins - 1)].strip().split('\t')[3] == 'nan':				
							break
						# Break loop if on last compartment
						if i+compartment_bins > len(data)-1:
							break
						nextline = data[i+compartment_bins].strip().split('\t')
						# Break loop in next eigen is NA
						if nextline[3] == 'nan':
							break
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
						
					# Calculate size in bp
					compartment_size = bin_size * compartment_bins
					if eigen > 0:
						A_sizes.append(compartment_size)
					elif eigen < 0: 
						B_sizes.append(compartment_size)

					i = i+compartment_bins
				else:
					i = i + 1

	print(len(A_sizes))
	print(len(B_sizes))
	print(A_sizes[:29])
	print(B_sizes[:29])

	print 'Mean A size: '
	print np.mean(A_sizes)
	print 'Mean B size: '
	print np.mean(B_sizes)

	print 'Median A size: '
	print np.median(A_sizes)
	print 'Median B size: '
	print np.median(B_sizes)



if __name__ == '__main__':
	main()
