#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import argparse
import matrix_functions as mf 
import h5py
import sys

parser = argparse.ArgumentParser(description='create boxplots of interactions for timecourse' +
	' at different distances')
parser.add_argument('-i', help= 'timecourse hdf5 files in order Mock -> Overnight digestion',
		nargs=7)
args = parser.parse_args()

def check_shape(x):
	'''
	Check if dimensions of numpy arrays are equal

	Args: 
		x : list of numpy ndarrays
	Returns:
		bool: True if dimensions are equal,
		false otherwise		
	'''
	dim_1 = x[0].shape
	flag = True
	for m in x[1:]:
		if m.shape != dim_1:
			flag = False 
	return flag


def main():

	# List of timecourse file objects in order
	f_obj_list = []
	for a in args.i:
		f_obj_list.append(h5py.File(a, 'r'))
	# List of chromosome 14 matrices for timecourse data
	# in order
	chr14_list = []
	for o in f_obj_list:
		chr14_list.append(mf.get_cis_matrix(o, 'chr14'))

	if not check_shape(chr14_list):
		print 'ERROR unequal dimensions for cis matrices'
		sys.exit()


	data = []
	for c in chr14_list:
		iu1 = np.triu_indices(c.shape[0])
		print sum(c[iu1])
		data.append(c[iu1])
	plt.figure()
	plt.boxplot(data, showmeans=True)
	plt.ylim(-10, 5000)
	plt.savefig('testboxplot.png', dpi=300)
	plt.close()
	




if __name__ == '__main__':
	main()
