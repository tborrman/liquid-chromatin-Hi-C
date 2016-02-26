#!/usr/bin/env python
import argparse 
import numpy as np 
from sklearn import mixture

parser=argparse.ArgumentParser(description='Fit a Gaussian Mixture Model to the DpnII signal data to extract copy number states')
parser.add_argument('-d', help='DpnII signal bed file (ex. dpnII_coverage_sorted.bed)', type=str, required=True, dest='d')
args=parser.parse_args()

dpnII = np.genfromtxt(args.d, delimiter='\t', usecols=3)
# Remove outlier peaks

dpnII = dpnII[(dpnII <= 3500) & (dpnII >= 500)]

dpnII = dpnII.reshape((len(dpnII),1))
dpnII = dpnII

g=mixture.GMM(n_components=3)
g.fit(dpnII)

print "converged"
print np.round(g.converged_,2)
print "weights"
print np.round(g.weights_, 2)
print "means"
print np.round(g.means_, 2)
print "covaraiance"
print np.round(g.covars_, 2)
print "score_samples"
print len(g.score_samples(dpnII))


