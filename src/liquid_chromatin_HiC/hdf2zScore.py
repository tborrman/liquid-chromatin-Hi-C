#!/usr/bin/env python
import argparse
import matrix_functions as mf 
import numpy as np 
import shutil
import sys
import h5py

def main():
	
	parser = argparse.ArgumentParser(description='Return hdf5 with z-score cis matrices for input hdf5')
	parser.add_argument('-i', help='input hdf5 Hi-C file', type=str, required=True)
	args= parser.parse_args()
	
	# Generate z-score hdf5
	if args.i[-5:] == '.hdf5':
		filename_prefix = args.i[:-5]
		zscore_file = filename_prefix + '_zScore.hdf5'
		shutil.copy(args.i, zscore_file)
		f = h5py.File(zscore_file, 'r+')
		obs = f['interactions'][:]
		blocksize = f['interactions'].chunks
		# Loop through chromosomes:
		chroms = f['chrs'][:]
		for chrom in chroms:
			obs_cis = mf.get_cis_matrix(f, chrom)
			z = mf.z_score(obs_cis)
			# Replace observed cis matrix with zscore matrix
			chr_idx, = np.where(chroms == chrom)
			chr_idx = chr_idx[0]
			bins = f['chr_bin_range'][chr_idx]
			obs[bins[0]:bins[1] + 1, bins[0]: bins[1] + 1] = z

		del f['interactions']
		f.create_dataset("interactions",  data=obs, dtype='float64', 
			compression='gzip', chunks=blocksize)
		f.close()
	
	else:
		print 'ERROR: wrong file extension'
		sys.exit()

if __name__ == '__main__':
	main()
