#!/usr/bin/env python
import subprocess

def main():

	for i in range(293):
		command = 'bsub -q short -R "rusage[mem=5000]" -n 1 -W 2:00 -o ' \
			+ str(i) + '.out -e ' + str(i) + '.err ./get_gene_positions.py -i RNAseq_' \
			'{:03d}'.format(i)
		subprocess.call(command, shell=True)

if __name__ == '__main__':
	main()
	
