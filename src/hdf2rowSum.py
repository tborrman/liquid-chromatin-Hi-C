#!/usr/bin/env python2

import argparse
import matrix_functions as mf
import numpy as np
import h5py
import sys

parser = argparse.ArgumentParser(description='Sum rows of Hi-C hdf5 and output sums as a bedGraph')
parser.add_argument('-i', help='input hdf5 file', type=str, required=True)
args = parser.parse_args()


def main():
	# Get matrix of interactions
	m = mf.hdf5_2_numpy_matrix(args.i)
	rowSums = np.nansum(m, 1)
	f = h5py.File(args.i, 'r') 
	# Write output
	OUT = open(args.i[:-5] + '_rowSum.bedGraph', 'w')
	# Only using 22 autosomes and X
	y_chrom_bin_start =  f['chr_bin_range'][:][23][0]
	for i, b in enumerate(f['bin_positions'][:][:y_chrom_bin_start]):
		if b[0] == 22:
			chrom = 'chrX'
		else:
			chrom = 'chr' + str(b[0]+1)
		OUT.write(chrom+'\t'+str(b[1])+'\t'+str(b[2])+
		'\t'+ str(rowSums[i])+'\n')
	OUT.close()

if __name__ == '__main__':
	main()
