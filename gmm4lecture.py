#!/usr/bin/env python
import numpy as np 
from sklearn import mixture
data = np.genfromtxt(args.d, delimiter='\t', usecols=3)
mino, maxo = get_outlier_cutoffs(data)
dpnII_reads = remove_outliers(data, mino, maxo)
dpnII_reads = dpnII_reads.reshape((len(dpnII_reads),1))
g=mixture.GMM(n_components=3)
g.fit(dpnII_reads)
