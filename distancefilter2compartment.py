#!/usr/bin/env python2
import argparse
import gzip
import multiprocessing as mp
import time

parser = argparse.ArgumentParser(description='Count A-A and B-B interactions for distance filtered dedups')
parser.add_argument('-i', help='file prefix before *_cisdistance*', type=str, required=True)
parser.add_argument('-t', help='total interactions for timepoint', type=float, required=True)
args = parser.parse_args()


def get_compartment_counts(e):
	'''
	Count total A bins and total B bins

	Args: 
		e: string eigenvector file
	Returns: 
		A_count: float A bin counts
		B_count: float B bin counts
	'''
	A_count = 0
	B_count = 0
	zero_count = 0
	with open(e, 'r') as IN:
		for line in IN:
			splitline = line.strip().split('\t')
			if splitline[3] != 'nan':
				eigen = float(splitline[3])
				if eigen > 0:
					A_count += 1
				elif eigen < 0:
					B_count += 1
	return A_count, B_count


def test(my_i, my_q):
	time.sleep(5)
	print my_i
	print my_q
	print mp.current_process()
	my_q.put((my_i*my_i, my_i))

	



def main():
	distances = ['1000-2000', '2000-3000', '3000-4000', '4000-5000',
				 '5000-6000', '6000-7000', '7000-8000', '8000-9000',
				 '9000-10000']
	eigenfile = '/home/tb37w/project/Research/digest/feature_analysis/eigen/eigen1_40kb.bedGraph'
	total_intxns = args.t
	A_bins, B_bins = get_compartment_counts(eigenfile)

	x = [1,2,3,4]
	q = mp.Queue()
	processes = [mp.Process(target=test, args=(i, q)) for i in x]
	for p in processes:
		p.start()
	for p in processes:
		p.join()

	print 'len: ' + str(len(processes))
	result = [q.get() for i in range(len(processes))]
	print 'result: ' + str(result)
	
if __name__ == '__main__':
	main()



