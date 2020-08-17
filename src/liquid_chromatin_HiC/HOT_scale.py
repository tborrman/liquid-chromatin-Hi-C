#!/usr/bin/env python
import argparse
import sys
import os
import subprocess
import shutil
import h5py
import matrix_functions as mf


def all_x_all_total_reads(file):
	'''
	Calculate the total number of reads in the all_x_all matrix
	'''
	if file[-5:] == '.hdf5':
		x = mf.hdf5_2_numpy_matrix(file)
		total = mf.total_reads(x)
	else:
		print 'ERROR: wrong file extension'
		sys.exit()
	return total

def scaleMatrixbyFactor(file, factor):
	# HDF5 format
	if file[-5:] == '.hdf5':
		filename_prefix = file[:-5]
		scaled_file = filename_prefix + '_scaleBy_' + str(round(factor,2)) + '.hdf5'
		shutil.copy(file, scaled_file)
		f = h5py.File(scaled_file, 'r+')
		scaled_interactions = f['interactions'][:] * factor
		blocksize = f['interactions'].chunks
		del f['interactions']
		f.create_dataset("interactions",  data=scaled_interactions, dtype='float64', compression='gzip', chunks=blocksize)
		f.close()
		return scaled_file
		
	else:
		print 'ERROR: wrong file extension'
		sys.exit()

def main():

	parser=argparse.ArgumentParser(description='Houda Tyler scaling: scale matrices such that ' +
		' allxall = billion total reads; ice matrices; scale again. Must be run on hdf5 matrices. ' +
		'balance.py for icing step must be in PATH' )
	parser.add_argument('-i', help='input matrix hdf5 format (HBCRACKHiC-K562-DN__hg19__genome__C-10000000-raw.hdf5)',
		type=str, required=True)
	args = parser.parse_args()

	# Get total reads of allxall matrix
	total_reads = all_x_all_total_reads(args.i)
	scale_factor = 1000000000.0 / total_reads

	# Scale file
	scaled_file = scaleMatrixbyFactor(args.i, scale_factor)
	print 'scaling...'

	# Ice file
	p = subprocess.Popen(['balance.py', '-i', scaled_file])
	p.wait()
	iced_file = scaled_file[:-5] + '.balanced.hdf5'
	print 'icing...'

	# Rescale iced matrix
	total_reads = all_x_all_total_reads(iced_file)
	scale_factor = 1000000000.0 / total_reads
	hot_file = scaleMatrixbyFactor(iced_file, scale_factor)
	print 'rescaling...'
	
	# Cleanup
	os.remove(scaled_file)
	os.remove(scaled_file[:-5] + '.factors')
	os.remove(iced_file)
	


if __name__ == '__main__':
	main()
