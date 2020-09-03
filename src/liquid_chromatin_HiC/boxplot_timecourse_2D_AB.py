#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse
import matrix_functions as mf 
import h5py
import sys

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

def is_compartment_type(p, t, e):
	'''
	Args: 
		p: indices of pixel in heatmap
		t: type of compartment interaction AA, BB or AB
		e: eigenvector numpy array for chromosome
	Returns:
		b: boolean true if pixel is correct type of 
			compartment interaction
	'''
	
	c1 =  e[p[0], 2]
	c2 = e[p[1], 2]
	if c1 >= 0 and c2 >= 0:
		ptype = 'AA'
	elif c1 < 0 and c2 < 0:
		ptype = 'BB'
	elif (c1 * c2) <= 0:
		ptype = 'AB'
	else:
		print 'ERROR: cannot assign pixel compartment type'
		sys.exit()
	if ptype == t:
		b = True
	else:
		b = False
	return b
	

def get_compartment_interaction_mean(c, d, e, t):
	'''
	Get mean interaction at distance for type of compartment 
	interaction
	
	Args:
		c: numpy matrix for chromosome
		d: int index of diagonal
		e: eigenvector numpy array for chromosome
		t: string type of comparment interaction AA, BB, or AB
	Return:
		m: mean interaction
	'''
	# Compartment type interaction list
	c_inter = []

	# Check e and c dimensions
	n = c.shape[0]
	if n != e.shape[0]:
		print 'ERROR: unequal eigen and matrix dim'
		sys.exit()

	x_idx = np.array(range(0,n-d))
	y_idx = np.array(range(0,n-d)) + d

	if t == 'AA':
		for i in range(len(x_idx)):
			pixel = (x_idx[i], y_idx[i])
			if is_compartment_type(pixel, t, e):
				c_inter.append(c[pixel]) 
	elif t == 'BB':
		for i in range(len(x_idx)):
			pixel = (x_idx[i], y_idx[i])
			if is_compartment_type(pixel, t, e):
				c_inter.append(c[pixel]) 
	elif t == 'AB':
		for i in range(len(x_idx)):
			pixel = (x_idx[i], y_idx[i])
			if is_compartment_type(pixel, t, e):
				c_inter.append(c[pixel])
	else:
		print 'ERROR: unknown compartment interaction type'
		sys.exit()
	# If there are no compartment type interactions occuring at this distance
	# assign mean interactions to zero
	if len(c_inter) == 0:
		m = 0
	else:
		m = np.mean(c_inter)
	return m


def plot_mean_interaction_distance_AB(l, e, t):
	'''
	Plot mean interactions vs timepoints across different 
	distances for specific compartment interactions
	
	Args: 
		l: list of numpy matrices for each timepoint
		e: eigenvector numpy array for chromosome
		t: string type of compartment interaction AA, BB, or AB
	'''
	x = [1,2,3,4,5,6,7]
	timelabels = ['MN', '5m', '1h', '2h', '3h', '4h', 'OVN']
	ncol = l[0].shape[0]
	binsizeMb = .5
	distances = np.arange(ncol) * binsizeMb

	# For drawing colorbar
	# https://stackoverflow.com/questions/26545897/drawing-a-colorbar-aside-a-line-plot-using-matplotlib
	# this is identical to my_norm() function above
	norm = matplotlib.colors.Normalize(vmin=np.min(distances), vmax=np.max(distances))

	# choose a colormap
	c_m = plt.cm.cool

	# create a ScalarMappable and initialize a data structure
	s_m = plt.cm.ScalarMappable(cmap=c_m, norm=norm)
	s_m.set_array([])

	# Set figure parmaters 
	fig, ax = plt.subplots(figsize=(8,6))
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)	
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.set_ylabel('Mean Interactions', fontsize=12)

	# Get mean interactions at distances
	for i, d in enumerate(distances):
		interaction_means = []
		for c in l:
			cim = get_compartment_interaction_mean(c, i, e, t)
			interaction_means.append(cim)
		# Plot lines
		ax.plot(x, interaction_means, color=s_m.to_rgba(d))

	ax.set_xticklabels(timelabels)
	cbar = fig.colorbar(s_m)
	cbar.set_label('Distance (Mb)')
	fig.tight_layout()
	plt.savefig('chr14_mean_interactions_distance_timecourse_' +t+'.png', dpi=300)
	plt.close()
	

	# Zoom figure 
	fig, ax = plt.subplots(figsize=(8,6))
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)	
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.set_ylabel('Mean Interactions', fontsize=12)

	# Get mean interactions at distances
	for i, d in enumerate(distances):
		interaction_means = []
		for c in l:
			cim = get_compartment_interaction_mean(c, i, e, t)
			interaction_means.append(cim)
		ax.plot(x, interaction_means, color=s_m.to_rgba(d))
	ax.set_xticklabels(timelabels)
	cbar = fig.colorbar(s_m)
	cbar.set_label('Distance (Mb)')
	ax.set_ylim(0,500)
	fig.tight_layout()
	plt.savefig('chr14_mean_interactions_distance_timecourse_zoom_500_'+t+'.png', dpi=300)
	plt.close()


def get_chrom_eigen(IN, c):
	'''
	Get compartment eigenvector values
	for chromosome
	
	Args:
		IN: filehandle for compartment file
		c: string] chromosome
	Return:
		e: 3D numpy array of start, stop, eigen
	'''
	e_list = []
	for line in IN:
		spltline = line.split()
		if spltline[0] == c:
			e_list.append(map(float,spltline[1:]))
	e = np.array(e_list)
	return e


def main():
	parser = argparse.ArgumentParser(description='create mean interactions vs timecourse and color' +
	' by distance plots' )
	parser.add_argument('-i', help= 'timecourse hdf5 files in order Mock -> Overnight digestion',
		nargs=7)
	parser.add_argument('-c', help= 'compartment file (ex. HBHiC-K562-MN-Dp-1__hg19__genome__C-500000-raw_scaleBy_2.72.balanced_scaleBy_51.45__all.zScore.eigen1.sorted.bedGraph)', 
	type=str, required=True)
	args = parser.parse_args()

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

	# Make mean interaction plots by distane for 
	# A-A, B-B, A-B
	CF = open(args.c, 'r')
	eigen = get_chrom_eigen(CF, 'chr14')
	plot_mean_interaction_distance_AB(chr14_list, eigen, 'AA')
	plot_mean_interaction_distance_AB(chr14_list, eigen, 'BB')
	plot_mean_interaction_distance_AB(chr14_list, eigen, 'AB')



if __name__ == '__main__':
	main()
