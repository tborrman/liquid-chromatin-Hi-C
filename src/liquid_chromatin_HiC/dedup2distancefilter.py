#!/usr/bin/env python2

import argparse
import gzip
import sys


def check_zip(filename):
	'''
	Check if input file is gzipped or not
	and return filehandle

	Args:
		filename: string
	Returns:
		FH: output filehandle 
	'''


	if filename[-3:] == '.gz':
		FH = gzip.open(filename, 'rb')
	else:
		FH = open(filename, 'r')
	return FH

def filter_by_distance(f, d):
	'''
	Filter dedup file for cis interactions using d distances
	
	Args: 
		f: string filename
		d: list distancest
	Returns:
		count_dict: dictionary counts for interactions within d distances
		linecount: int number of lines of f
	'''

	IN = check_zip(f)
	OUT_list = []
	count_dict = {}
	linecount = 0

	for i in range(len(d)-1):
		# Open output files
		if f[-3:] == '.gz':
			OUT_list.append(gzip.open(f[:-3] + '_cisdistance' + str(d[i]) + '-' + str(d[i+1]) + '.gz', 'wb'))
		else:
			OUT_list.append(gzip.open(f + '_cisdistance' + str(d[i]) + '-' + str(d[i+1]) + '.gz', 'wb'))
		# Initialize counts
		count_dict[str(d[i]) + '-' + str(d[i+1])] = 0
	for line in IN:
		linecount += 1
		splitline = line.split('\t')
		if len(splitline) != 12:
			print 'Error with input file! Must have 12 columns!'
			print splitline
			sys.exit()
		chromosome_1=splitline[1]
		readPos_1=int(splitline[2])
		chromosome_2=splitline[7]
		readPos_2=int(splitline[8])

		if chromosome_1 == chromosome_2:
			cis_dist =abs(readPos_2 - readPos_1)
			for i in range(len(d)-1):
				if cis_dist > d[i] and cis_dist < d[i+1]:
					OUT_list[i].write(line)
					count_dict[str(d[i]) + '-' + str(d[i+1])] += 1

	IN.close()
	for OUT in OUT_list:
		OUT.close()

	return count_dict, linecount

def main():	

	parser=argparse.ArgumentParser(description='filter *_dedup files by distance')
	parser.add_argument('-i', help='dedup file', type=str, required=True)
	parser.add_argument('-d', help='distances', type=int, nargs='+', default=[1000, 3000, 5000, 10000, 25000, 50000])
	args=parser.parse_args()

	distances = args.d
	c, l = filter_by_distance(args.i, distances)
	print '***************************'
	print 'Total interactons: ' + str(l)
	print '***************************'
	print '*****DISTANCE COUNTS*******'
	for key in c.keys():
		print key + ': ' + str(c[key]) + ' (' + str(round((float(c[key])/float(l))*100, 2)) + '%)'  
	print '***************************'

if __name__ == '__main__':
	main()