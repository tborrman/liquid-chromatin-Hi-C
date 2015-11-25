#!/usr/bin/env python
import argparse
import sys
import os
# Flush STOUT continuously
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

parser = argparse.ArgumentParser(description='Write list of distances to closest DpnII site for reads')
parser.add_argument('-s', help='SAM file with alignment reads (bowtie_output_mapped.sam)', type=str, dest='s', required=True)
args = parser.parse_args()

def main():

	# Store dpnII coordinates in dictionary
	dpnII_sites = {}
	chroms = map(str, range(1,23)) + ['X','Y']
	for chrom in chroms:
		DPN = open('../dpnII/dpnII_sites_chr'+chrom+'.bed', 'r')
		coords = []
		for line in DPN:
			coord1 = int(line.split()[1]) + 1
			coord2 = int(line.split()[2]) + 1
			coords.append(coord1)
			coords.append(coord2)
		dpnII_sites['chr'+chrom] = coords
		DPN.close()
	print 'Processed coordinate files'

	# Parse SAM file
	SAM = open(args.s, 'r')
	index = args.s[-3:]
	OUT = open('dpnII_distances_'+index,'w')
	counter = 0
	for line in SAM:
		counter += 1
		if counter%100 == 0:
			print 'On read: '+ str(counter)
		# Get chromosome
		chrom = line.split()[2]
		start = int(line.split()[3])
		# Find closest dpnII site
		closest = min(dpnII_sites[chrom], key=lambda x: abs(start-x))
		OUT.write(str(closest-start) + '\n')
	SAM.close()
	OUT.close()

if __name__ == '__main__':
	main()

