#!/usr/bin/env python
import shutil
import argparse
import numpy as np
import h5py
import matrix_functions as mf 


def main():

	parser=argparse.ArgumentParser(description='scalar multiply matrix by input factor')
	parser.add_argument('-i', help= 'input .matrix file or .hdf5 file', type=str, required=True)
	parser.add_argument('-f', help= 'scalar factor', type=float, required=True)
	args = parser.parse_args()
	
	# HDF5 format
	if args.i[-5:] == '.hdf5':
		filename_prefix = args.i[:-5]
		shutil.copy(args.i, filename_prefix + '_scaleMatrixbyFactor_' + str(args.f) + '.hdf5')
		f = h5py.File(filename_prefix + '_scaleMatrixbyFactor_' + str(args.f) + '.hdf5', 'r+')
		scaled_interactions = f['interactions'][:] * args.f
		blocksize = f['interactions'].chunks
		del f['interactions']
		f.create_dataset("interactions",  data=scaled_interactions, dtype='float64', compression='gzip', chunks=blocksize)
		f.close()
		
	else:
		if args.i[-10:] == '.matrix.gz':
			filename_prefix = args.i[:-10]
		elif args.i[-7:] == '.matrix':
			filename_prefix = args.i[:-7]
		else:
			print 'ERROR: not .matrix or .hdf5 file'
			quit()

		X, xheader, yheader = mf.dekker_2_numpy_matrix(args.i)
		scaled_X = X*args.f
		mf.numpy_matrix_2_dekker(scaled_X, xheader, yheader, filename_prefix + '_scaleMatrixbyFactor_' + str(args.f) + '.matrix.gz')

if __name__ == '__main__':
	main()