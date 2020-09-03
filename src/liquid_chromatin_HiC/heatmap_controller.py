#!/usr/bin/env python
import subprocess

def main():
	p = subprocess.Popen('ls', shell=True, stdout= subprocess.PIPE)
	files = p.communicate()[0].split('\n')

	for f in files:
		subprocess.call("bsub -q short -W 1:00 -R 'rusage[mem=5000]' perl /home/tb37w/project/Research/ENCODE/cworld-dekker/scripts/perl/heatmap.pl " +
			"-i " + f, shell=True)

if __name__ == '__main__':
	main()
	