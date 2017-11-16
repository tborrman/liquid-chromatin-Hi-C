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

def my_boxplot(d):
	'''
	Boxplots for timecourse interactions

	Args:
		d: list of numpy interaction arrays for
			timecourse timepoints
	'''
	# Full boxplot
	fig, ax = plt.subplots(figsize=(5,8))
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	flierpointprops = dict(marker='o', markeredgecolor='black',
		markerfacecolor='white')
	timelabels = ['MN', '5m', '1h', '2h', '3h', '4h', 'OVN']
	bplot = ax.boxplot(d, showmeans=True, flierprops=flierpointprops,
		meanline=True, labels=timelabels, patch_artist=True)
	plt.tick_params(axis= 'x',labelsize=15)
	ax.set_ylabel('Interactions', fontsize=15)
	plt.setp(bplot['whiskers'], color='k')
	plt.setp(bplot['medians'], color='k')
	plt.setp(bplot['means'], color = 'r')
	plt.setp(bplot['boxes'], facecolor='0.75')
	fig.tight_layout()
	plt.savefig('chr14_interactions_boxplot_timecourse.png', dpi=300)
	plt.close()
	
	# Zoom boxplot
	fig, ax = plt.subplots(figsize=(5,8))
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	bplot = ax.boxplot(d, showmeans=True, flierprops=flierpointprops,
		meanline=True, labels=timelabels, patch_artist=True)
	plt.tick_params(axis= 'x',labelsize=15)
	ax.set_ylabel('Interactions', fontsize=15)
	ax.set_ylim(-10,1000)
	# for box in bplot['boxes']:
	# 	box.set_facecolor('0.75')
	plt.setp(bplot['whiskers'], color='k')
	plt.setp(bplot['medians'], color='k')
	plt.setp(bplot['means'], color = 'r')
	plt.setp(bplot['boxes'], facecolor='0.75')
	fig.tight_layout()
	plt.savefig('chr14_interactions_boxplot_timecourse_zoom.png', dpi=300)
	plt.close()



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
		data.append(c[iu1])
	my_boxplot(data)
	




if __name__ == '__main__':
	main()
