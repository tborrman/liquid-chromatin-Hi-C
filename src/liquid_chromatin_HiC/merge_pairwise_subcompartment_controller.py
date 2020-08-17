#!/usr/bin/env python
import subprocess
import argparse

def main():

	parser=argparse.ArgumentParser(description='merge pairwise and subcompartment files for all chromosomes given two directories')
	parser.add_argument('-p', help='path to pairwise-files directory (ex./home/tb37w/project/Research/digest/HiC/HBCRACKHiC-K562-MN-R1__hg19__txt/C-100000/iced/scaled_matrices/pairwise/)',
	 type=str, required=True)
	parser.add_argument('-s', help='path to subcompartment files directory (ex. /home/tb37w/project/Research/digest/subcompartment/barplots/)',
	 type=str, required=True)
	parser.add_argument('-o', help='prefix for input and output files (ex. HBCRACKHiC-K562-MN-R1__hg19__genome__C-100000-iced__chr)', type=str, required=True)
	args = parser.parse_args()

	for chrom in map(str,range(1, 23)):
		subprocess.call("bsub -q short -W 4:00 -R 'rusage[mem=10000]' -o merge.out -e merge.err /home/tb37w/project/Research/digest/digest-Hi-C/merge_pairwise_subcompartment.py " +
			"-p " + args.p + args.o + chrom + ".scaled-1000000--ic--it.pairwise-score.txt.gz " + 
			"-s " + args.s + "GSE63525_GM12878_subcompartments_sorted_100kb_chr" + chrom + ".bed " + 
			"-o " + args.o + chrom + "_merge_pairwise_subcompartment.txt", 
			shell=True)

if __name__ == '__main__':
	main()

