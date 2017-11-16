#!/usr/bin/env python
import numpy as np
import gzip
import h5py
import sys

def dekker_2_numpy_matrix(filename):
	''' 
	Load Dekker format interaction matrix txt file into a numpy matrix
	'''
	x = []
	# Check zip
	if filename[-3:] == '.gz':
		IN = gzip.open(filename, 'rb')
	else:
		IN = open(filename, 'r')  
	
	is_xheader = True
	yheader = []
	for line in IN:
	    if '#' != line[0]:
	    	if is_xheader:
	    		xheader = line.strip().split('\t')[1:]
	    		is_xheader = False
	    	else:
	    		split_line = line.split('\t')
	    		yheader.append(split_line[0])
	    		x.append(map(lambda x: np.nan if 'nan' in x or 'NA' in x else float(x), split_line[1:]))

	IN.close()
	return np.array(x, dtype=float), xheader, yheader

def remove_NA_zeros(X, xheader, yheader, percent):
	'''
	Remove rows and columns with more than 'percent' NAs or zeros
	'''
	# Remove rows	
	#rows_to_keep = np.sum(np.logical_or(np.isnan(X), X==0), 1).astype(float)/len(X[0,:]) <= percent
	rows_to_keep = np.sum(np.isnan(X), 1).astype(float)/len(X[0,:]) <= percent
	X_removed_rows = X[rows_to_keep, :]
	# Update yheader
	yheader = np.array(yheader, dtype=str)
	yheader_update = yheader[rows_to_keep]

	# Remove columns
	#cols_to_keep = np.sum(np.logical_or(np.isnan(X), X==0), 0).astype(float)/len(X[:,0]) <= percent
	cols_to_keep = np.sum(np.isnan(X), 0).astype(float)/len(X[:,0]) <= percent
	X_removed_NA_zeros  =  X_removed_rows[:, cols_to_keep]
	# Update xheader
	xheader = np.array(xheader, dtype=str)
	xheader_update = xheader[cols_to_keep]
	return(X_removed_NA_zeros, xheader_update, yheader_update)

def numpy_matrix_2_dekker(X, colnames, rownames, filename):
	''' 
	Write numpy matrix into dekker format interaction matrix txt file
	'''
	shape = 'x'.join(map(str,X.shape))
	OUT = gzip.open(filename, 'wb')
	colnames.insert(0,shape)
	OUT.write('\t'.join(colnames) + '\n')
	for row_idx in range(len(X)):
		line = X[row_idx,:]
		fmt_line = ['{:.8f}'.format(x) for x in line]
		fmt_line.insert(0, rownames[row_idx])
		OUT.write('\t'.join(fmt_line) + '\n')
	OUT.close()

def total_reads(X):
	'''
	Input numpy matrix and output count of total interactions 
	in upper triangle including diagonal
	'''
	up_tri = np.triu(X)
	y = np.nansum(up_tri, axis=1)
	return np.nansum(y)

def hdf5_2_numpy_matrix(filename):
	'''
	Load hdf5 format interaction matrix into numpy matrix
	'''
	f = h5py.File(filename, 'r')
	x = f['interactions'][:]
	return x

def get_cis_matrix(f, chrom):
	'''
	Get cis interaction matrix
	
	Args:
		f: h5py Hi-C file object
	Returns:
		c: numpy cis interaciton matrix for 
			chrom
	'''
	chr_idx, = np.where(f['chrs'][:] == chrom)
	chr_idx = chr_idx[0]
	bins = f['chr_bin_range'][chr_idx]
	c = f['interactions'][bins[0]:bins[1] + 1, bins[0]: bins[1] + 1]
	return c

def get_all_trans(f):
	'''
	Get all trans interactions from hdf5

	NOTE: this is grabbing all trans from full matrix so double counting

	Args:
		f: h5py Hi-C file object
	Returns:
		trans : numpy array of all trans interactions
	'''

	# Test example
	# bins = np.array([[0,0],[1,2],[3,3],[4,5]])
	# x = np.array([[1,2,3,4,5,6],[5,6,7,8,9,1],[9,1,2,3,4,5],[4,5,6,7,8,9],[5,1,-2,-5,-1,-9], [-1,-3,4,-7,8,-9]])

	bins = f['chr_bin_range'][:]
	# Remove Y chrom since K562 and also mitochondrial chrom
	bins = bins[:-2,:]
	x = f['interactions'][:]
	

	trans = np.array([])
	row_length = bins[-1,:][1] + 1

	for chrom_bin in bins:
		chrom_length = chrom_bin[1] - chrom_bin[0] + 1
		trans_length = row_length - chrom_length
		rows = []
		cols = []
		for i in range(chrom_bin[0], chrom_bin[1] + 1):
			rows.append([i]*trans_length)
			cols.append(range(0,chrom_bin[0]) + range(chrom_bin[1] + 1, row_length))
		trans = np.append(trans, x[rows, cols])
	return trans


		



