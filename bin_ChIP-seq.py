#!/usr/bin/env python

import argparse
import subprocess

parser = argparse.ArgumentParser(description='Clean and bin ChIP-seq bigWig into bedGraph')
parser.add_argument('-bw', help='bigWig file', type=str, required=True)
parser.add_argument('-b', help='binning bed file (ex. hg19_500kb.bed)', type=str, required=True)
args = parser.parse_args()




def main():
	f = args.bw[:-7]
	# Convert to bedGraph
	print 'Converting to bedGraph...'
	p = subprocess.Popen("./bigWigToBedGraph " + args.bw + " " + f + ".bedGraph", shell=True)
	p.wait()
	print 'Success'
	# Remove unwanted chromosomes and sort 
	print 'Removing chrY and chrM and sorting...'
	p = subprocess.Popen("grep -v 'chrM' " + f +".bedGraph  | " +
		"grep -v 'chrY' | sort -k1,1V -k2,2n > " + f + "_clean.bedGraph", shell=True)
	p.wait() 
	print 'Success'
	# bin
	print 'Binning...'
	p = subprocess.Popen("bedtools map -a " + args.b + " -b "+ f + "_clean.bedGraph " + 
	 "-c 4 -o mean -null NA > " + f + "_500kb.bedGraph", shell=True)
	print 'Success'
	p.wait()
	print 'Removing tmp file...'
	p = subprocess.Popen("rm " + f + "_clean.bedGraph", shell=True)
	print 'Success'


if __name__ == '__main__':
	main()