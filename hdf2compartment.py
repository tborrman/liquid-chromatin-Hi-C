#!/usr/bin/env python

import argparse
import os
import sys
import subprocess
import re


parser = argparse.ArgumentParser(description='Extract all cis maps from hdf5 and '+
	'run matrix2compartment.pl PCA analysis on extracted maps')
parser.add_argument('-i', help='input hdf5 interaction file', type=str, required=True)
parser.add_argument('-cp', help='path to cworld-dekker (eg. /home/username/cworld-dekker',
	type=str, required=True)
parser.add_argument('-hp', help='path to hdf2tab (eg. /home/username/hdf2tab',
	type=str, required=True)
args = parser.parse_args()


def sleep(attempt):
	q = subprocess.Popen('sleep 2', shell=True)
	q.wait()
	attempt += 1
	return attempt


def job_success(jobID, attempt):
	'''
	Check if job successfully completed
	'''
	# bjobs STAT does not immediately switch to DONE even when job is complete
	# retry check if STAT is RUN or JOBID is not found yet
	# (figure out a better way to handle job checks later)
	p=subprocess.Popen('bjobs ' + jobID, shell=True, stdout=subprocess.PIPE)
	jobStat = p.communicate()[0]
	if attempt > 10:
		success = False
	elif len(jobStat.strip().split('\n')) != 2: 
		attempt = sleep(attempt)
		success = job_success(jobID, attempt)
	elif jobStat.split('\n')[1].split()[2] != 'DONE':
		attempt = sleep(attempt)
		success = job_success(jobID, attempt)
	else:
		success = True
	return success


def extract_cis_matrices(hdf, hpath):
	'''
	Extract cis interaction matrices from hdf5 file
	'''
	print 'Extracting cis matrices into directory matrices/....'
	if not os.path.exists('matrices'):
		os.makedirs('matrices')

	p=subprocess.Popen('bsub -K -q short -W 4:00 -R "rusage[mem=10000]" -o hdf2tab.out '+
		'-e hdf2tab.err python ' + hpath + '/scripts/hdf2tab.py -i ' + hdf + ' -v -wm cis '+
		'-o matrices/' + hdf[:-5], shell=True, stdout=subprocess.PIPE)
	p.wait()

	myID = re.search('<(\d+)>', p.communicate()[0]).group(1)
	if job_success(myID, 1):
		print 'Cis matrices succesfully extracted'
	else:
		print 'ERROR in extracting matrices'


def main():
	# Extract cis matrices from hdf5
	extract_cis_matrices(args.i, args.hp)




if __name__ == '__main__':
	main()
