#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import argparse
import matrix_functions as mf
import sys
import h5py
import shutil
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(description='Create half-life matrix from timecourse')
parser.add_argument('-i', help= 'timecourse hdf5 files in order Mock -> Overnight digestion', nargs=7)
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


def scatter_plot(i, t, x, y, out, hl, m):
	'''
	Plot a scatter_plot of interactions vs time

	Args:
		i: interactions from timecourse
		t: time in minutes
	'''
	fig, ax = plt.subplots(figsize=(8,4))
	ax.set_xlabel(r'Minutes', fontsize=10)
	ax.set_ylabel(r'Interactions', fontsize=10)
	#ax.set_xlim([-10,1200])
	#ax.set_ylim([-10,6500])
	ax.scatter(t, i)
	ax.plot(x, y)
	ax.vlines(hl, -1, 1, colors='r')
	ax.hlines(m, 0, 500, colors='r')
	#ax.text(800, 4000, '$t_{1/2}$' +' = ' + str(int(round(hl, 0))) + ' min')
	plt.savefig(out, dpi=300)
	plt.close()
	return

def half_life_equation(m, a, b, c):
	'''
	Get half life using exponential decay function
	solved for minutes at half loss of interactions

	Args:
		m: midpoint for delta interactions
		a,b,c: parameters for decay function

	Returns:
		hl: half-life

	'''
	if ((a-m)/b) <= 0 or c == 0:
		hl = np.nan
	else:
		hl = -(np.log((a-m)/b) / c)
	return hl

def exp_decay(minutes, a, b, c):
	'''
	Sigmoid function with adjustable parameters
	
	Args: 
	 	minutes: minutes for timecourse
	 	a,b,c: parameters for decay function
	
	Returns:
	r = range of function
	'''
	r = a - (b*np.exp(-c*minutes))
	return r




def get_half_life(d, t, m):
	'''
	Get half life of interactions for 2D position in
	Hi-C matrix

	Args: 
		d: delta interactions from timecourse
		t: time in minutes
		m: midpoint for delta interactions
	Returns:
		hl: half-life for pixel
		popt: list of fitted parameters
	'''
	# Curve fit
	try:
		popt, pcov = curve_fit(exp_decay, t, d, p0= (0.75, 0.5, 0.01))
	except RuntimeError as error:
		print error
		hl = np.nan
		popt = np.nan
	else:
		hl = half_life_equation(m, *popt)
	return hl, popt
	

def get_delta_intxns(n):
	'''
	Get delta interactions for given
	set of interactions for pixels from
	timecourse
	delta = (mock - digest) / mock
	'''
	n = np.array(n, dtype=float)
	mock = n[0]
	digest = n[1:]
	delta = (mock - digest) / mock
	return delta

def get_intxns_half(d):
	'''
	Get midpoint for delta interactions

	Args:
		d: delta interactions
	Returns:
		m: midpoint for delta interactions
	'''
	start = d[0]
	end= d[-1]
	m = start + ((end - start)/2)
	return m

def main():

	# Create copy of mock and use that to write over
	# with half life data
	shutil.copy(args.i[0], 'half_life_chr22_6Mb.hdf5')
	f = h5py.File('half_life_chr22_6Mb.hdf5', 'r+')

	# List of timecourse file objects in order
	f_obj_list = []
	for a in args.i:
		f_obj_list.append(h5py.File(a, 'r'))

	# List of chromosome 14 matrices for timecourse data
	# in order
	chr14_list = []
	for o in f_obj_list:
		chr14_list.append(mf.get_cis_matrix(o, 'chr22'))

	if not check_shape(chr14_list):
		print 'ERROR unequal dimensions for cis matrices'
		sys.exit()
	
	time = np.array([5, 60, 120, 180, 240, 960], dtype=float)
	
	# # 1 MB distance
	# interactions = []
	# for h in chr14_list:
	# 	interactions.append(h[1250, 1275])
	# # delta interactions
	# delta_intxns = get_delta_intxns(interactions)
	# mid = get_intxns_half(delta_intxns)
	# hl, params = get_half_life(delta_intxns, time, mid)
	# x = np.linspace(-10, 1000, 100, dtype=float)
	# y = exp_decay(x, *params)
	# scatter_plot(delta_intxns, time, x, y, 'test_1MB_distance_delta.png', hl, mid)
	

	# # 400kb distance
	# interactions = []
	# for h in chr14_list:
	# 	interactions.append(h[1000, 1010])
	# # delta interactions
	# delta_intxns = get_delta_intxns(interactions)
	# mid = get_intxns_half(delta_intxns)
	# hl, params = get_half_life(delta_intxns, time, mid)
	# x = np.linspace(-10, 1000, 100, dtype=float)
	# y = exp_decay(x, *params)
	# scatter_plot(delta_intxns, time, x, y, 'test_400kb_distance_delta.png', hl, mid)

	chrom = 'chr22'
	dist = (6000000/40000)/2
	#obs = f['interactions'][:]
	bin_positions = f['bin_positions'][:]
	num_bins = len(bin_positions)
	chr_idx, = np.where(f['chrs'][:] == chrom)
	chr_idx = chr_idx[0]
	bins = f['chr_bin_range'][chr_idx]
	for i in range(bins[0], bins[1] + 1):
	#for i in range(bins[0] + 1250, bins[1] + 1):
		print 'on row: ' + str(i)
		# 6Mb window
		# Check if row is all nan
		if np.all(np.isnan(f['interactions'][i])) or np.nansum(f['interactions'][i]) == 0:
			f['interactions'][i] = np.nan
		# Check if at start of chromosome (upstream range reaches trans)
		elif bin_positions[i,0] != bin_positions[i-dist,0]:
			f['interactions'][i] = np.nan
		# Check if at end of chromosome (downstream range reaches trans) 
		elif bin_positions[i,0] != bin_positions[i+(dist-1),0]:
			f['interactions'][i] = np.nan
		else:
			for j in range(bins[0], bins[1] + 1):
				# Check if in 6Mb window
				if j >= i-dist and j < i+dist:
					interactions = []
					for h in chr14_list:
						interactions.append(h[i - bins[0], j - bins[0]])
					if np.any(np.isnan(interactions)) or (interactions[0] == 0):
						f['interactions'][i, j] = np.nan
					else:
						delta_intxns = get_delta_intxns(interactions)
						mid = get_intxns_half(delta_intxns)
						hl, params = get_half_life(delta_intxns, time, mid)
						f['interactions'][i, j] = hl
						f['interactions'][j, i] = hl
				else:
					f['interactions'][i, j] = np.nan
	f.close()


if __name__ == '__main__':
	main()