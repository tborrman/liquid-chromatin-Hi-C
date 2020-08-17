#!/usr/bin/env python
import subprocess
import argparse

parser=argparse.ArgumentParser(description='merge pairwise and eigen1 files for all chromosomes given two directories')
parser.add_argument('-p', help='path to pairwise-files directory (ex./home/tb37w/project/Research/digest/HiC/HBCRACKHiC-K562-MN-R1__hg19__txt/C-100000/iced/scaled_matrices/pairwise/)',
 type=str, required=True)
parser.add_argument('-e', help='path to eigen1 files directory (ex. /home/tb37w/project/Research/digest/HiC/HBCRACKHiC-K562-MN-R1__hg19__txt/C-100000/iced/eigen/)',
 type=str, required=True)
parser.add_argument('-o', help='prefix for input and output files (ex. HBCRACKHiC-K562-MN-R1__hg19__genome__C-100000-iced__chr)', type=str, required=True)
parser.add_argument('-s', help='scaled matrix total (ex. 1000000)', type=str, default=None)
args = parser.parse_args()

if args.s:	
	for chrom in map(str,range(1, 23)) + ['X']:
		subprocess.call("bsub -q short -W 4:00 -R 'rusage[mem=10000]' -o merge.out -e merge.err /home/tb37w/project/Research/digest/digest-Hi-C/merge_pairwise_eigen1.py " +
			"-p " + args.p + args.o + chrom + ".scaled-" + args.s + "--ic--it.pairwise-score.txt.gz " + 
			"-e " + args.e + args.o + chrom + ".zScore.eigen1.bedGraph " + 
			"-o " + args.o + chrom + ".scaled-" + args.s + "_merge_pairwise_eigen1.txt", 
			shell=True)

else:
	for chrom in map(str,range(1, 23)) + ['X']:
		subprocess.call("bsub -q short -W 4:00 -R 'rusage[mem=10000]' -o merge.out -e merge.err /home/tb37w/project/Research/digest/digest-Hi-C/merge_pairwise_eigen1.py " +
			"-p " + args.p + args.o + chrom + "--ic--it.pairwise-score.txt.gz " + 
			"-e " + args.e + args.o + chrom + ".zScore.eigen1.bedGraph " + 
			"-o " + args.o + chrom + "_merge_pairwise_eigen1.txt", 
			shell=True)