#!/usr/bin/env python
import argparse
import sys
import os
import shutil
import h5py
import matrix_functions as mf


parser=argparse.ArgumentParser(description='Houda Tyler scaling: scale matrices such that ' +
	' allxall = billion total reads; ice matrices; scale again. Must be run in a resolution ' + 
	'directory from resulting Dekker mapping pipeline (ex. HiC_library/C-10000000/).')
parser.add_argument('-p', help='sample prefix before __hdf or __txt (ex. HBCRACKHiC-K562-MN-R1__hg19)',
	type=str, required=True)
parser.add_argument('-t', help="type of files: 'matrix' or 'hdf'", type=str, required=True)
args = parser.parse_args()


def all_x_all_total_reads(file):
	'''
	Calculate the total number of reads in the all_x_all matrix
	'''
	if file[-5:] == '.hdf5':
		x = mf.hdf5_2_numpy_matrix(file)
		total = mf.total_reads(x)
	elif file[-10:] == '.matrix.gz':
		x, colnames, rownames = mf.dekker_2_numpy_matrix(file)
		total = mf.total_reads(x)
	else:
		print 'ERROR: wrong file extension'
		sys.exit()
	return total

def scaleMatrixbyFactor(file, factor):
	# HDF5 format
	if file[-5:] == '.hdf5':
		filename_prefix = file[:-5]
		shutil.copy('raw/' + file, filename_prefix + '_scaleMatrixbyFactor_' + str(factor) + '.hdf5')
		f = h5py.File(filename_prefix + '_scaleMatrixbyFactor_' + str(factor) + '.hdf5', 'r+')
		scaled_interactions = f['interactions'][:] * factor
		blocksize = f['interactions'].chunks
		del f['interactions']
		f.create_dataset("interactions",  data=scaled_interactions, dtype='float64', compression='gzip', chunks=blocksize)
		f.close()
		
	else:
		if file[-10:] == '.matrix.gz':
			filename_prefix = file[:-10]
		elif file[-7:] == '.matrix':
			filename_prefix = file[:-7]
		else:
			print 'ERROR: not .matrix or .hdf5 file'
			quit()

		X, xheader, yheader = mf.dekker_2_numpy_matrix('raw/' + file)
		scaled_X = X*factor
		mf.numpy_matrix_2_dekker(scaled_X, xheader, yheader, filename_prefix + '_scaleMatrixbyFactor_' + str(factor) + '.matrix.gz')


def main():

	if not os.path.exists('hot'):
		os.makedirs('hot')
	###########################################################
	# Step 1. Scale 
	# Check if allxall matrix exits
	if args.t == 'matrix':
		ext = 'matrix.gz'
	elif args.t == 'hdf':
		ext = 'hdf5'
	else:
		print 'ERROR: incorrect file type'
		sys.exit()
	allxall = '../C-10000000/raw/' + args.p + '__genome__C-10000000-raw.' + ext
	try:
		open(allxall, 'r')
	except IOError as e:
		print 'File does not exist'
		print e
		sys.exit()
	# Get total reads of allxall matrix
	total_reads = all_x_all_total_reads(allxall)
	scale_factor = 1000000000.0 / total_reads

	# Scale all files
	for file in os.listdir('raw'):
		if file[-5:] == '.hdf5' or file[-10:] == '.matrix.gz' or file[-7:] =='.matrix':
			scaleMatrixbyFactor(file, scale_factor)


	







	




if __name__ == '__main__':
	main()
