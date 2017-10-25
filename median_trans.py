#!/usr/bin/env python
import argparse
import matrix_functions as mf 
import numpy as np 
import h5py

parser = argparse.ArgumentParser(description='Calculate median trans signal from hdf5')
parser.add_argument('-i', help='input hdf5 matrix file', type=str, required=True)
args = parser.parse_args()


def main():

	h = h5py.File(args.i, 'r')
	t = mf.get_all_trans(h)
	mt = np.median(t)
	print str(args.i) + '\t' + str(round(mt, 4))


if __name__ == '__main__':
	main()