#!/usr/bin/env python
import argparse
import matrix_functions as mf
import numpy as np


parser = argparse.ArgumentParser(description='Calculate contact enrichment to merge odd/even subcompartment calls')
parser.add_argument('-i', help='input interaction matrix file', type=str, required=True)
parser.add_argument('-o', help='odd chrom subcompartment calls (ex. odd_compartments.tab)', type=str, required=True)
parser.add_argument('-e', help='even chrom subcompartment calls (ex. even_compartments.tab)', type=str, required=True)
args = parser.parse_args()

def main():

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

	#Calculate contact enrichment
	contact_enrich = np.zeros((5,5), dtype=object)
	for i in range(len(contact_enrich)):
		for j in range(len(contact_enrich[i])):
			contact_enrich[i,j] = []
	ODD_FH = open(args.o, 'r')
	EVEN_FH = open(args.e, 'r')
	odd = []
	even = []
	for line in ODD_FH:
		odd.append(int(line.split()[1]))
	ODD_FH.close()
	for line in EVEN_FH:
		even.append(int(line.split()[1]))
	EVEN_FH.close()
	# Size check
	if len(odd) != len(X_step1) or len(even) != len(X_step1[0,:]):
		print 'ERROR: size mismatches exist'
		quit()
 
	for i in range(len(X_step1)):
		if i % 1000 == 0:
			print 'On row: ' + str(i)
		for j in range(len(X_step1[0,:])):
			odd_class = odd[i]
			even_class = even[j] 
			contact_enrich[odd_class, even_class].append(X_step1[i,j])

	for i in range(len(contact_enrich)):
		for j in range(len(contact_enrich[0,:])):
			contact_enrich[i,j] = np.mean(contact_enrich[i,j])

	np.savetxt('contact_enrichment.tab', contact_enrich, delimiter='\t', fmt='%1.2f')


if __name__ == '__main__':
	main()