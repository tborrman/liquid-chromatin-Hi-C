#!/usr/bin/env python
import argparse
import numpy as np
import matrix_functions as mf 

parser=argparse.ArgumentParser(description='scalar multiply matrix by input factor')
parser.add_argument('-i', help= 'input matrix file', type=str, required=True)
parser.add_argument('-f', help= 'scalar factor', type=float, required=True)
args = parser.parse_args()

def main():

	# Filenaming
	if args.i[-10:] == '.matrix.gz':
		filename_prefix = args.i[:-10]
	elif args.i[-7:] == '.matrix':
		filename_prefix = args.i[:-7]
	else:
		print 'ERROR: not .matrix file'
		quit()

	X, xheader, yheader = mf.dekker_2_numpy_matrix(args.i)
	scaled_X = X*args.f
	mf.numpy_matrix_2_dekker(scaled_X, xheader, yheader, filename_prefix + '_scaleMatrixbyFactor_' + str(args.f) + '.matrix.gz')

if __name__ == '__main__':
	main()