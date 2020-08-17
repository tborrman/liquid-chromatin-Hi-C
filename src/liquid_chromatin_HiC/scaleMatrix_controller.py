#!/usr/bin/env python
import subprocess
import argparse

def main(): 

	parser=argparse.ArgumentParser(description='run scaleMatrix.pl on all chromosome maps in current directory')
	parser.add_argument('-f', help='file prefix (ex. HBCRACKHiC-K562-DN__hg19__genome__C-500000-iced__chr)', type=str, required=True)
	args = parser.parse_args()

	for chrom in map(str,range(1, 23)) + ['X']:
		subprocess.call("bsub -q short -W 1:00 -R 'rusage[mem=5000]' perl /home/tb37w/project/Research/ENCODE/cworld-dekker/scripts/perl/scaleMatrix.pl " +
			"-i " + "./" + args.f + chrom + ".matrix.gz", 
			shell=True)
		# subprocess.call("bsub -q short -W 1:00 -R 'rusage[mem=5000]' perl /home/tb37w/project/Research/ENCODE/cworld-dekker/scripts/perl/scaleMatrix.pl " +
		# 	"-i " + MN_path + "HBCRACKHiC-K562-MN-R1__hg19__genome__C-500000-iced__chr" + chrom + ".matrix", 
		# 	shell=True)
if __name__ == '__main__':
	main()
	