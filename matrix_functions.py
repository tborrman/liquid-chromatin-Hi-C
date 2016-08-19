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

def remove_NA_zeros(X, percent):
	'''
	Remove rows and columns with more than 'percent' NAs or zeros
	'''
	# Remove rows	
	X_removed_rows = X[np.sum(np.logical_or(np.isnan(X), X==0), 1).astype(float)/len(X[0,:]) <= percent, :]
	# Remove columns
	X_removed_NA_zeros =  X_removed_rows[:, np.sum(np.logical_or(np.isnan(X), X==0), 0).astype(float)/len(X[:,0]) <= percent]
	return(X_removed_NA_zeros)