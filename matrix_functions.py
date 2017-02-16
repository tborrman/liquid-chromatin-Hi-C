#!/usr/bin/env python
import numpy as np
import gzip
import h5py

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


