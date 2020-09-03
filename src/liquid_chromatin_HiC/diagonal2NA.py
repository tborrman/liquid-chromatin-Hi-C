#!/usr/bin/env python
import argparse
import matrix_functions as mf
import shutil
import h5py
import numpy as np



def main():

	parser=argparse.ArgumentParser(description='Change diagonals of dekker format hdf5 interaction matrix to NA')
	parser.add_argument('-i', help='Input hdf5 file', type=str, required=True)
	parser.add_argument('-k', help='Number of diagonals to NA offset from main diagonal (use k = 0 to only NA main diagonal)',
		type=int, required=True)
	args=parser.parse_args()
	
	filename_prefix = args.i[:-5]
	nodiag_file = filename_prefix + '_diags_NA_' + str(args.k + 1) + '.hdf5'
	shutil.copy(args.i, nodiag_file)
	f = h5py.File(nodiag_file, 'r+')
	m = f['interactions'][:]
	for i in range(args.k + 1):
		nanarray = np.repeat(np.nan, np.diag(m,i).size)
		diag_matrix = np.diag(nanarray, i)
		m = m + diag_matrix
		if i > 0:
			diag_matrix = np.diag(nanarray, -i)
			m = m + diag_matrix
	blocksize = f['interactions'].chunks
	del f['interactions']
	f.create_dataset("interactions",  data=m, dtype='float64', compression='gzip', chunks=blocksize)
	f.close()

if __name__ == '__main__':
	main()