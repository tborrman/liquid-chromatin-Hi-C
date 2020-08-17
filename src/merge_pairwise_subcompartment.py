#!/usr/bin/env python

import argparse
import gzip
import pandas as pd
import re

parser=argparse.ArgumentParser(description="merge '.pairwise-score' files and  GSE63525_GM12878_subcompartments_sorted_100kb.bed)")
parser.add_argument('-p', help='pairwise-score file (ex. HBCRACKHiC-K562-MN-R1__hg19__genome__C-100000-iced__chr1.scaled-1000000--ic--it.pairwise-score.txt.gz)',
	type=str, required=True)
parser.add_argument('-s', help='subcompartment file file (ex.  GSE63525_GM12878_subcompartments_sorted_100kb.bed)', 
	type=str, required=True)
parser.add_argument('-o', help='output file', type=str, required=True)

args=parser.parse_args()


def load_subcomp(file):
	'''
	Load subcompartment file as pandas dataframe
	'''
	subcomp_df = pd.read_csv(file, sep='\t', header=None, names= ['chrom', 'start', 'stop', 'subcompartment'])
	return subcomp_df

def main():
	# Open pairwise-file
	if args.p[-3:] == '.gz':
		PF = gzip.open(args.p, 'rb')
	else:
		PF = open(args.p, 'r')
	# Load subcompartment file
	subcomp_df = load_subcomp(args.s)

	OUT = open(args.o, 'w')

	PF.readline()
	chrom = subcomp_df.loc[0,'chrom']
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

		subcomp_I = str(list(subcomp_df[subcomp_df['start'] == int(I_start)]['subcompartment'])[0])
		subcomp_J = str(list(subcomp_df[subcomp_df['start'] == int(J_start)]['subcompartment'])[0])

		
		IF = splitline[2]
		OUT.write('\t'.join([chrom, I_start, J_start, subcomp_I, subcomp_J, IF]) + '\n')
		counter +=1
	PF.close()
	OUT.close()


if __name__ == '__main__':
	main()