#!/usr/bin/env python
import argparse
import matrix_functions as mf
import numpy as np

parser = argparse.ArgumentParser(description='Perform Rao et al. (2014) subcompartment analysis')
parser.add_argument('-i', help='input interaction matrix file', type=str, required=True)

args = parser.parse_args()

def main():
	X, xheader, yheader = mf.dekker_2_numpy_matrix(args.i)
	print X, xheader, yheader

	# Step 1
	#  Rows and columns for which more than 30% of the entries were either undefined or zeros
	#  were removed from the matrix

	X_step1, xheader, yheader = mf.remove_NA_zeros(X, xheader, yheader,.65)
	print X_step1, xheader, yheader







if __name__ == '__main__':
	main()
