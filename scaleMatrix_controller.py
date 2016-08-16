#!/usr/bin/env python
import subprocess

DN_path = '/home/tb37w/project/Research/digest/HiC/HBCRACKHiC-K562-DN-R1__hg19__txt/C-500000/iced/'
MN_path = '/home/tb37w/project/Research/digest/HiC/HBCRACKHiC-K562-MN-R1__hg19__txt/C-500000/iced/'

for chrom in map(str,range(1, 23)) + ['X']:
	# subprocess.call("bsub -q short -W 1:00 -R 'rusage[mem=5000]' perl /home/tb37w/project/Research/ENCODE/cworld-dekker/scripts/perl/scaleMatrix.pl " +
	# 	"-i " + DN_path + "HBCRACKHiC-K562-DN-R1__hg19__genome__C-500000-iced__chr" + chrom + ".matrix", 
	# 	shell=True)
	subprocess.call("bsub -q short -W 1:00 -R 'rusage[mem=5000]' perl /home/tb37w/project/Research/ENCODE/cworld-dekker/scripts/perl/scaleMatrix.pl " +
		"-i " + MN_path + "HBCRACKHiC-K562-MN-R1__hg19__genome__C-500000-iced__chr" + chrom + ".matrix", 
		shell=True)