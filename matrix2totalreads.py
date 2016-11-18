#!/usr/bin/env python
import argparse
import matrix_functions as mf

parser=argparse.ArgumentParser(description='Input a matrix and output count of total interactions in upper triangle including diagonal')
parser.add_argument('-i', help='Input matrix', type=str, dest='i', required=True)
args=parser.parse_args()

def main():

	x, colnames, rownames = mf.dekker_2_numpy_matrix(args.i)
	total_reads = mf.total_reads(x)
	print 'Total reads: ' + str(total_reads)


if __name__ == '__main__':
	main()




