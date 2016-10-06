#!/usr/bin/env python
import subprocess
import argparse

parser=argparse.ArgumentParser(description='run matrix2compartment.pl on all chromosome maps in current directory')
parser.add_argument('-f', help='file prefix (ex. HBCRACKHiC-K562-DN__hg19__genome__C-500000-iced__chr)', type=str, required=True)
args = parser.parse_args()

for chrom in map(str,range(1, 23)) + ['X']:
	subprocess.call("bsub -q short -W 1:00 -R 'rusage[mem=5000]' perl /home/tb37w/project/Research/ENCODE/cworld-dekker/scripts/perl/matrix2compartment.pl " +
		"-i " + "./" + args.f + chrom + ".matrix.gz", 
		shell=True)