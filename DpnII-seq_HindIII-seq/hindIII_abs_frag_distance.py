#!/usr/bin/env python
import argparse
import sys
import os
# Flush STOUT continuously
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

parser = argparse.ArgumentParser(description='Write list of distances to closest HindIII cut site for fragment '+
	'(distance reported is the minimum distance from either fragment end)')
parser.add_argument('-s', help='SAM file with alignment reads (bowtie_output_mapped.sam)', type=str, dest='s', required=True)
args = parser.parse_args()

def main():

	# Store hindIII coordinates in dictionary
	hindIII_sites = {}
	chroms = map(str, range(1,23)) + ['X','Y']
	for chrom in chroms:
		DPN = open('../hindIII/hindIII_sites_chr'+chrom+'.bed', 'r')
		coords = []
		for line in DPN:
			# Change coordinates from BED file format
			coord1 = int(line.split()[1]) + 1
			coord2 = int(line.split()[2])
			# Get coordinates for nucleotides after cut site on 5' end and before cut site on 3' end
			coord1 = coord1 + 1
			coord2 = coord2 - 1 
			coords.append(coord1)
			coords.append(coord2)
		hindIII_sites['chr'+chrom] = coords
		DPN.close()
	print 'Processed coordinate files'

	# Parse SAM file
	SAM = open(args.s, 'r')
	index = args.s[-3:]
	OUT = open('hindIII_abs_frag_distances_'+index,'w')
	counter = 0
	while True:
		line1=SAM.readline()
		line2=SAM.readline()
		if not line2:
			break
		else:
			counter += 1
			if counter%100 == 0:
				print 'On fragment: '+ str(counter)
			# Get chromosome
			chrom1 = line1.split()[2]
			chrom2 = line2.split()[2]
			if chrom1 != chrom2:
				print 'ERROR: wrong chromosomes'
				quit()
			# Ignore mitochondrial DNA
			if chrom1 == 'chrM':
				continue
			# Check alignment has pair
			if (int(line1.split()[8]) < 0) or (int(line2.split()[8]) > 0):
				print 'ERROR: no pair'
				quit()
			# Get start and end of fragment
			start = int(line1.split()[3])
			end = int(line2.split()[3]) + 99
			# Find closest hindIII sites
			closest_start = min(hindIII_sites[chrom1], key=lambda x: abs(start-x))
			closest_end = min(hindIII_sites[chrom1], key=lambda x: abs(end-x))
			# Minimum distance to hindIII sites from either end
			min_distance = min(abs(start-closest_start), abs(end-closest_end))
			OUT.write(str(min_distance) + '\n')
	SAM.close()
	OUT.close()

if __name__ == '__main__':
	main()

