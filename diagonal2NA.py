#!/usr/bin/env python
import argparse
import matrix_functions as mf
import shutil
import h5py
import numpy as np

parser=argparse.ArgumentParser(description='Change diagonal of dekker format hdf5 interaction matrix to NA')
parser.add_argument('-i', help='Input hdf5 file', type=str, required=True)
args=parser.parse_args()

def main():
	filename_prefix = args.i[:-5]
	nodiag_file = filename_prefix + '_nodiag.hdf5'
	shutil.copy(args.i, nodiag_file)
	f = h5py.File(nodiag_file, 'r+')
	m = f['interactions'][:]
	np.fill_diagonal(m, np.nan)
	blocksize = f['interactions'].chunks
	del f['interactions']
	f.create_dataset("interactions",  data=m, dtype='float64', compression='gzip', chunks=blocksize)
	f.close()

if __name__ == '__main__':
	main()