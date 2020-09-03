#!/usr/bin/env python2
import argparse
import gzip
import multiprocessing as mp
import time
import sys




def get_compartment_counts(e):
	'''
	Count total A bins and total B bins

	Args: 
		e: string eigenvector file
	Returns: 
		A_count: float A bin counts
		B_count: float B bin counts
	'''
	A_count = 0.0
	B_count = 0.0

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


def make_eigen_dict(e):
	'''
	Make dictionary of eigenvector data using chromosome keys
	
	Args:
		e: string eigenvector bedGraph file
	Returns:
		ed: dictionary of eigenvector data
	'''

	# Initialize eigenvector dictionary 
	ed = {}
	chroms = ['chr' + x for x in map(str, range(1,23))] + ['chrX']
	for chrom in chroms:
		ed[chrom] = []
	with open(e, 'r') as EIGEN:
		for line in EIGEN:
			splitline = line.strip().split('\t')
			linechrom = splitline[0]
			ed[linechrom].append(line)
	return ed
	
def split_by_compartment(f, t, a, b, ed, q):
	'''
	Split dedup validpair distance filtered
	file by compartment interactions
	
	Args:
		f: string dedup validpair file
		t: float total interactions for sample
		a: float A bin counts
		b: float B bin counts
		ed: eigenvector dictionary
		q: mp Queue object for storing results
	'''
	
	# Specific for gzip files
	IN = gzip.open(f, 'rb')
	A_intxns = 0.0
	B_intxns = 0.0
	for i, line in enumerate(IN):
		if i%1000 == 0:
			print 'On row: ' + str(i)
		splitline = line.strip().split('\t')
		chrom = splitline[1]
		validchroms = ['chr' + x for x in map(str, range(1,23))] + ['chrX']
		if chrom not in validchroms:
			continue
		p1 = int(splitline[2])
		p2 = int(splitline[8])
		# Check which compartment reads belong to and 
		# throw out if reads belong to separate compartments
		foundPos1Eigen = False
		foundPos2Eigen = False
		for eline in ed[chrom]:
			esplit = eline.strip().split()
			e1 = int(esplit[1])
			e2 = int(esplit[2])
			if p1 < e2 and p1 >= e1:
				p1_eigen = esplit[3]
				foundPos1Eigen = True
			if p2 < e2 and p2 >= e1:
				p2_eigen = esplit[3]
				foundPos2Eigen = True
			if foundPos1Eigen and foundPos2Eigen:
				break
		if foundPos1Eigen and foundPos2Eigen:
			if p1_eigen == 'nan' or p2_eigen == 'nan':
				continue
			elif float(p1_eigen) > 0 and float(p2_eigen) > 0:
				A_intxns += 1
			elif float(p1_eigen) < 0 and float(p2_eigen) < 0:
				B_intxns += 1
			else:
				pass
	IN.close()
	# Normalize interactions
	norm_A = A_intxns/(t*a)
	norm_B = B_intxns/(t*b)
	# Save in queue
	result_dict = {f: [("A", norm_A), ("B", norm_B)]}
	q.put(result_dict)



def write_compartment_table(d, f, results, c, p):
	'''
	Write results to table

	Args:
		d: list of distances
		f: list of distance files
		results: list of dictionary results
		c: string compartment ('A' or 'B')
		p: prefix for sample
	'''
	with open(p + '_cisdistance_' + c  + '.txt', 'w') as OUT:
		OUT.write('\t' + '\t'.join(d) + '\n')
		OUT.write(p)
		for df in f:
			for r in results: 
				if df in r: 
					if c == 'A':
						OUT.write('\t' + str(r[df][0][1]))
					elif c == 'B':
						OUT.write('\t' + str(r[df][1][1]))
					else:
						print 'ERROR: unknown compartment type'
						sys.exit()
					break
		OUT.write('\n')
		


def main():

	parser = argparse.ArgumentParser(description='Count A-A and B-B interactions for distance filtered dedups')
	parser.add_argument('-i', help='file prefix before *_cisdistance*', type=str, required=True)
	parser.add_argument('-t', help='total interactions for timepoint', type=float, required=True)
	args = parser.parse_args()

	distances = ['1000-2000', '2000-3000', '3000-4000', '4000-5000',
				 '5000-6000', '6000-7000', '7000-8000', '8000-9000',
				 '9000-10000']
	eigenfile = '/home/tb37w/project/Research/digest/feature_analysis/eigen/eigen1_40kb.bedGraph'
	eigen_dict = make_eigen_dict(eigenfile)
	A_bins, B_bins = get_compartment_counts(eigenfile)
	total_intxns = args.t
	distance_files = [args.i + '_cisdistance' + d + '.gz' for d in distances]


	# x = [1,2,3,4]
	# q = mp.Queue()
	# processes = [mp.Process(target=test, args=(i, q)) for i in x]
	# for p in processes:
	# 	p.start()
	# for p in processes:
	# 	p.join()

	# print 'len: ' + str(len(processes))
	# result = [q.get() for i in range(len(processes))]
	# print 'result: ' + str(result)

	q = mp.Queue()
	processes = [mp.Process(target=split_by_compartment, 
		args=(d, total_intxns, A_bins, B_bins, eigen_dict, q)) for d in distance_files]
	for p in processes:
		p.start()
	for p in processes:
		p.join()
	results = [q.get() for i in range(len(processes))]
	# Write results to compartment tables
	write_compartment_table(distances, distance_files, results, 'A', args.i)
	write_compartment_table(distances, distance_files, results, 'B', args.i)

if __name__ == '__main__':
	main()



