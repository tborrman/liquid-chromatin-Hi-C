#!/usr/bin/env python
import subprocess

def main():
	DN_path = '/home/tb37w/project/Research/digest/HiC/HBCRACKHiC-K562-DN-R1__hg19__txt/C-500000/iced/scaled_matrices/'
	MN_path = '/home/tb37w/project/Research/digest/HiC/HBCRACKHiC-K562-MN-R1__hg19__txt/C-500000/iced/scaled_matrices/'

	for chrom in map(str,range(1, 23)) + ['X']:
		subprocess.call("bsub -q short -W 1:00 -R 'rusage[mem=5000]' perl /home/tb37w/project/Research/ENCODE/cworld-dekker/scripts/perl/compareMatrices.pl " +
			"-1 " + DN_path + "HBCRACKHiC-K562-DN-R1__hg19__genome__C-500000-iced__chr" + chrom + ".scaled-1000000.matrix.gz " +
			"-2 " + MN_path + "HBCRACKHiC-K562-MN-R1__hg19__genome__C-500000-iced__chr" + chrom + ".scaled-1000000.matrix.gz " +
			"--cm subtract",
			shell=True)

if __name__ == '__main__':
	main()