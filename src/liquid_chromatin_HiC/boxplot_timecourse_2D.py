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

def my_norm(x, xmin, xmax):
    return (x - xmin) / float(xmax - xmin)


def plot_mean_interaction_distance(l):
	'''
	Plot mean interactions vs timepoints across different 
	distances
	
	Args: 
		c: list of numpy matrices for each timepoint
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
			interaction_means.append(np.mean(np.diag(c, k= i)))
		# Plot lines
		ax.plot(x, interaction_means, color=s_m.to_rgba(d))
	ax.set_xticklabels(timelabels)
	cbar = fig.colorbar(s_m)
	cbar.set_label('Distance (Mb)')
	fig.tight_layout()
	plt.savefig('chr14_mean_interactions_distance_timecourse.png', dpi=300)
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
			interaction_means.append(np.mean(np.diag(c, k= i)))
		# Plot lines
		# if i == 19:
		# 	ax.plot(x, interaction_means, color='k')
		# else:
		# 	ax.plot(x, interaction_means, color=s_m.to_rgba(d))
		ax.plot(x, interaction_means, color=s_m.to_rgba(d))
	ax.set_xticklabels(timelabels)
	cbar = fig.colorbar(s_m)
	cbar.set_label('Distance (Mb)')
	ax.set_ylim(0,2500)
	fig.tight_layout()
	plt.savefig('chr14_mean_interactions_distance_timecourse_zoom.png', dpi=300)
	plt.close()

def main():
	parser = argparse.ArgumentParser(description='create boxplots of interactions for timecourse' +
	' at different distances')
	parser.add_argument('-i', help= 'timecourse hdf5 files in order Mock -> Overnight digestion',
		nargs=7)
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

	# Max
	z = chr14_list[2]
	max_idx = np.unravel_index(np.argmax(z), z.shape)
	print max_idx
	print z[max_idx]
	quit()




	# Boxplots for all interactions
	data = []
	for c in chr14_list:
		iu1 = np.triu_indices(c.shape[0])
		data.append(c[iu1])
	my_boxplot(data)

	# Mean interaction plots by distance
	plot_mean_interaction_distance(chr14_list)
	




if __name__ == '__main__':
	main()
