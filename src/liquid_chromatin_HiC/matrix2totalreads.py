#!/usr/bin/env python
import argparse
import sys
import matrix_functions as mf



def main():

	parser=argparse.ArgumentParser(description='Input a matrix and output count of total interactions in upper triangle including diagonal')
	parser.add_argument('-i', help='Input matrix', type=str, dest='i', required=True)
	args=parser.parse_args()

	if args.i[-5:] == '.hdf5':
		x = mf.hdf5_2_numpy_matrix(args.i)
		total_reads = mf.total_reads(x)
	elif args.i[-10:] == '.matrix.gz' or args.i[-7:] == '.matrix':
		x, colnames, rownames = mf.dekker_2_numpy_matrix(args.i)
		total_reads = mf.total_reads(x)
	else:
		print 'ERROR: not matrix or hdf format'
		sys.exit()
	print 'Total reads: ' + str(total_reads)


if __name__ == '__main__':
	main()




