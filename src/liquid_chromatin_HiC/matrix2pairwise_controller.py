#!/usr/bin/env python
import subprocess

def main():
	process = subprocess.Popen(['ls'], shell=True, stdout=subprocess.PIPE)
	out = process.communicate()[0]
	files = out.split('\n')
	for f in files:
		if 'HBCRACKHiC' in f:
			subprocess.call("bsub -q short -W 1:00 -R 'rusage[mem=5000]' " +
			"perl ~/project/Research/ENCODE/cworld-dekker/scripts/perl/matrix2pairwise.pl -i " + f, shell=True)

if __name__ == '__main__':
	main()
