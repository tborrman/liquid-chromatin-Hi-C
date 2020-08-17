#!/usr/bin/env python

import argparse
import subprocess

parser = argparse.ArgumentParser(description='Clean and bin WGBS bed into bedGraph')
parser.add_argument('-bw', help='WGBS bed file', type=str, required=True)
parser.add_argument('-b', help='binning bed file (ex. hg19_500kb.bed)', type=str, required=True)
parser.add_argument('-r', help='resolution (ex. 500kb)', type=str, default='500kb')
args = parser.parse_args()

def main():
	f = args.bw[:-4]
	# Convert to bedGraph
	# Remove unwanted chromosomes and sort 
	print 'Removing chrY and chrM and sorting...'
	p = subprocess.Popen("grep -v 'chrM' " + f +".bed  | " +
		"grep -v 'chrY' | sort -k1,1V -k2,2n > " + f + "_clean.bed", shell=True)
	p.wait() 
	print 'Success'
	# bin
	print 'Binning...'
	p = subprocess.Popen("bedtools map -a " + args.b + " -b "+ f + "_clean.bed " + 
	 "-c 11 -o mean -null NA > " + f + "_" + args.r + ".bedGraph", shell=True)
	print 'Success'
	p.wait()
	print 'Removing tmp file...'
	p = subprocess.Popen("rm " + f + "_clean.bed", shell=True)
	print 'Success'


if __name__ == '__main__':
	main()
