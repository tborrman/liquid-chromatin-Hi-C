#!/usr/bin/env python
import argparse 
import numpy as np 
from sklearn import mixture

parser=argparse.ArgumentParser(description='Fit a Gaussian Mixture Model to the DpnII signal data to extract copy number states')
parser.add_argument('-d', help='DpnII signal bed file (ex. dpnII_coverage_sorted.bed)', type=str, required=True, dest='d')
args=parser.parse_args()

def main():

	dpnII_chr = np.genfromtxt(args.d, delimiter='\t', usecols=(0), dtype=str)
	dpnII = np.genfromtxt(args.d, delimiter='\t', usecols=(1,2,3))

	# Remove outlier peaks
	dpnII_chr = dpnII_chr[(dpnII[:,2] <= 3500) & (dpnII[:,2] >= 500)]
	dpnII = dpnII[(dpnII[:,2] <= 3500) & (dpnII[:,2] >= 500)]
	dpnII_reads = dpnII[:,2]
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

	# Write posterior probabilities output file
	np.set_printoptions(suppress=True)
	means = np.round(g.means_, 2)
	probs = g.score_samples(dpnII_reads)[1]
	OUT=open('GMM_post_probs.txt', 'w')
	OUT.write('\t\t\t\t' + str(means[0,0]) +'\t'+ str(means[1,0]) +'\t'+ str(means[2,0]) + '\n')
	for i in range(len(dpnII)):
		OUT.write('\t'.join([dpnII_chr[i], str(int(dpnII[i,0])), str(int(dpnII[i,1])), 
			str(int(dpnII[i,2])), '{:.5f}'.format(probs[i][0]),
			'{:.5f}'.format(probs[i][1]),
			'{:.5f}'.format(probs[i][2])])
			 + '\n')
	OUT.close()
	
if __name__ == '__main__':
	main()




