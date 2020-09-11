#!/usr/bin/env python
import argparse
import matrix_functions as mf
import numpy as np 
import h5py


def main():

	parser=argparse.ArgumentParser(description='Calculate trans percentage for input hdf5')
	parser.add_argument('-i', help= 'input hdf5 file', type=str, required=True)
	args = parser.parse_args()

	h = h5py.File(args.i, 'r')
	t = mf.get_all_trans(h, False)
	trans_sum = np.sum(t) / 2.0
	
	x = mf.hdf5_2_numpy_matrix(args.i)
	total_reads = mf.total_reads(x)

	trans_pct = trans_sum/total_reads
	print 'Trans percent: ' + str(trans_pct)



if __name__ == '__main__':
	main()
