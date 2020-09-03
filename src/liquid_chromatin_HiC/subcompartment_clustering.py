#!/usr/bin/env python
import argparse
import matrix_functions as mf
import numpy as np
from scipy import stats
from hmmlearn import hmm

def main():

	parser = argparse.ArgumentParser(description='Perform Rao et al. (2014) HMM subcompartment analysis')
	parser.add_argument('-i', help='input interaction matrix file', type=str, required=True)
	parser.add_argument('-c', help='run HMM on columns/even chromosomes)', type=bool, default=False)

	args = parser.parse_args()
	
	X, xheader, yheader = mf.dekker_2_numpy_matrix(args.i)
	#print X, xheader, yheader
	print 'Matrix shape: ' + str(np.shape(X))

	# Step 1
	#  Rows and columns for which more than 30% of the entries were either undefined or zeros
	#  were removed from the matrix
	print 'Removing NAs...'
	X_step1, xheader, yheader = mf.remove_NA_zeros(X, xheader, yheader,.65)
	#print X_step1, xheader, yheader
	print 'Matrix shape: ' + str(np.shape(X_step1))

	# Step 2
	# Determine whether to transpose matrix
	if args.c:
		# Transpose matrix
		print 'Transposing matrix...'
		X_step2 = np.transpose(X_step1)
	else:
		X_step2 = X_step1


	# Step 3
	# Transform matrix to zscores by row
	print 'Transforming rows to zscores...'
	X_step3 = stats.zscore(X_step2, axis = 1)
	

	# Step 4
	# Train hmm model
	print 'Training HMM'
	model = hmm.GaussianHMM(n_components=5, covariance_type='diag', n_iter=1000)
	model.fit(X_step3)
	print 'HMM output'
	classes= model.predict(X_step3)
	if args.c:
		print '# of columns: ' + str(len(xheader))
		print '# of class entries: ' + str(len(classes))
		OUT = open('even_compartments.tab', 'w')
		for i in range(len(xheader)):
			OUT.write(xheader[i] + '\t' + str(classes[i]) + '\n')
		OUT.close()
	else:
		print '# of rows: ' + str(len(yheader))
		print '# of class entries: ' + str(len(classes))
		OUT = open('odd_compartments.tab', 'w')
		for i in range(len(yheader)):
			OUT.write(yheader[i] + '\t' + str(classes[i]) + '\n')
		OUT.close()

if __name__ == '__main__':
	main()
