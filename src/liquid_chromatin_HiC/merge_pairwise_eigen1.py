#!/usr/bin/env python

import argparse
import gzip
import pandas as pd
import re

def load_eigen1(file):
	'''
	Load eigen1 file as pandas dataframe
	'''
	eigen_df = pd.read_csv(file, sep='\t', header=None, names= ['chrom', 'start', 'stop', 'eigen'], skiprows=1)
	return eigen_df


def main():

	parser=argparse.ArgumentParser(description="merge '.pairwise-score' files and '.eigen1' files")
	parser.add_argument('-p', help='pairwise-score file (ex. HBCRACKHiC-K562-MN-R1__hg19__genome__C-100000-iced__chr1.scaled-1000000--ic--it.pairwise-score.txt.gz)',
		type=str, required=True)
	parser.add_argument('-e', help='eigen1 file (ex. HBCRACKHiC-K562-MN-R1__hg19__genome__C-100000-iced__chr1.zScore.eigen1.bedGraph)', 
		type=str, required=True)
	parser.add_argument('-o', help='output file', type=str, required=True)
	args=parser.parse_args()

	# Open pairwise-file
	if args.p[-3:] == '.gz':
		PF = gzip.open(args.p, 'rb')
	else:
		PF = open(args.p, 'r')
	# Load eigen file
	eigen_df = load_eigen1(args.e)

	OUT = open(args.o, 'w')

	PF.readline()
	chrom = eigen_df.loc[0,'chrom']
	counter = 0
	for line in PF:
		if counter%1000 == 0:
			print 'On row: ' + str(counter)
		splitline = line.strip().split('\t')
		# Grab start positions
		search_obj = re.search('chr\d+:(\d+)', splitline[0]) 
		I_start = search_obj.group(1)
		search_obj = re.search('chr\d+:(\d+)', splitline[1])
		J_start = search_obj.group(1)

		eigen_I = str(float(eigen_df[eigen_df['start'] == int(I_start)]['eigen']))
		eigen_J = str(float(eigen_df[eigen_df['start'] == int(J_start)]['eigen']))

		
		IF = splitline[2]
		OUT.write('\t'.join([chrom, I_start, J_start, eigen_I, eigen_J, IF]) + '\n')
		counter +=1
	PF.close()
	OUT.close()


if __name__ == '__main__':
	main()
