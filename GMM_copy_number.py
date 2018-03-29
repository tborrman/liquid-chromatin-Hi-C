#!/usr/bin/env python
import argparse 
import numpy as np 
from sklearn import mixture


parser=argparse.ArgumentParser(description='Fit a Gaussian Mixture Model to the DpnII signal data to extract copy number states')
parser.add_argument('-d', help='DpnII signal bed file (ex. dpnII_coverage_sorted.bed)', type=str, required=True, dest='d')
args=parser.parse_args()

def get_outlier_cutoffs(r, s):
	'''
	Return min and max outlier read cutoffs
	using s standard deviations
	Args:
		r: numpy array of reads per bin
		s: number of standard deviations to 
			determine outlier cutoff
	Returns:
		minr: min read cutoff
		maxr: max read cutoff
	'''
	# Remove zeros
	rnz = r[r>0]
	# Get cutoffs
	minr = np.mean(rnz) - (np.std(rnz) * s)
	maxr = np.mean(rnz) + (np.std(rnz) * s)
	return minr, maxr

def remove_outliers(d, min_reads, max_reads):
	''' 
	Remove bins with < min_reads or > max_reads
	Args:
		d: numpy array with cols= [start, stop, dpnII_reads]
		min_reads: minimum read cutoff
		max_reads: maximum read cutoff
	Returns:
		f: numpy array of filtered data
			with cols = [start, stop, dpnII_reads, index]
	'''

	idx_keep = np.where((d[:,2] <= max_reads) & (d[:,2] >= min_reads))[0]
	d = d[idx_keep,]
	idx_keep = np.reshape(idx_keep, (len(idx_keep), 1))
	f = np.append(d, idx_keep, axis=1)
	return f

def get_copynumber_states(mu, p):
	'''
	Assign copy number states based on posterior 
	probabilites of each mixture component for 
	each observation
	Args: 
		mu: np.array of means for mixture
			compononets
		p: np.array of posterior probabilites
	Returns:
		c: list of copy number states 
	'''
	
	sort_mu = np.sort(mu)
	copy_number = np.array([np.where(sort_mu == x)[0][0] for x in mu]) + 2
	# Get copy number state
	c = []
	for row in p:
		if np.all(np.isnan(row)):
			c.append(np.nan)
		else:
			max_idx = np.argmax(row)
			c.append(copy_number[max_idx])
	return c


def main():

	dpnII_chr = np.genfromtxt(args.d, delimiter='\t', usecols=(0), dtype=str)
	dpnII = np.genfromtxt(args.d, delimiter='\t', usecols=(1,2,3))

	mino, maxo = get_outlier_cutoffs(dpnII[:,2], 2)

	df = remove_outliers(dpnII, mino, maxo)

	dpnII_reads = df[:,2]
	dpnII_reads = dpnII_reads.reshape((len(dpnII_reads),1))

	g=mixture.GMM(n_components=3, tol=0.0001)
	g.fit(dpnII_reads)

	print "converged"
	print g.converged_
	print "weights"
	print np.round(g.weights_, 2)
	print "means"
	print np.round(g.means_, 2)
	print "covariance"
	print np.round(g.covars_, 2)

	
	np.set_printoptions(suppress=True)
	means = np.round(g.means_, 2).flatten()
	probs = g.score_samples(dpnII_reads)[1]
	probs_full = np.zeros((len(dpnII), 3))
	probs_full[:] = np.nan
	kept_idx = df[:, 3].astype(int)
	probs_full[kept_idx, :] = probs

	copynum = get_copynumber_states(means, probs_full)

	# Write copy number results to file
	OUT=open('GMM_post_probs.txt', 'w')
	OUT.write('chrom\tstart\tend\treads\t' + str(means[0]) +'\t'+ str(means[1]) +
		'\t'+ str(means[2]) + '\tcopy\n')
	for i in range(len(dpnII)):
		OUT.write('\t'.join([dpnII_chr[i], str(int(dpnII[i,0])), str(int(dpnII[i,1])), 
			str(int(dpnII[i,2])), '{:.5f}'.format(probs_full[i][0]),
			'{:.5f}'.format(probs_full[i][1]),
			'{:.5f}'.format(probs_full[i][2]),
			str(copynum[i])])
			 + '\n')
	OUT.close()
	
if __name__ == '__main__':
	main()




