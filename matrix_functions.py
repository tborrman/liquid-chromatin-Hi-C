#!/usr/bin/env python
import numpy as np
import gzip

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