#!/usr/bin/env python2
import argparse
import gzip

parser = argparse.ArgumentParser(description='filter validPairs.txt.gz file for valid pairs greater than d apart')
parser.add_argument('-i', help='input validPairs.txt.gz file', type=str, required=True)
parser.add_argument('-d', help='minimum bp distance between valid pairs (default 1kb)', type=int, default=1000)
args = parser.parse_args()

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

def filter_valid_pairs(NOFILTER, FILTERED, d):
	'''
	Filter valid pairs file for valid pairs that
	have a bp distance greater than d
	
	Args:
		NOFILTER: filehandle for unfiltered valid pairs
		FILTERED: filenandle for filtered valid pairs
		d: minimum bp distance between valid pairs
	'''
	counter = 0
	for line in NOFILTER:
		splitline = line.split('\t')
		# Check if cis
		if splitline[1] == splitline [7]:
			# Check valid distance
			if abs(int(splitline[2]) - int(splitline[8])) < d:
				counter += 1
			else:
				FILTERED.write(line)
		else:
			FILTERED.write(line)
	return counter
				
def main():

	IN = check_zip(args.i)
	suffix_idx = args.i[::-1].find('__') + 2
	new_name = args.i[:-suffix_idx] + '-filter' + str(args.d) + args.i[-suffix_idx:]
	OUT = gzip.open(new_name, 'wb')
	c = filter_valid_pairs(IN, OUT, args.d)
	print 'Filtered ' + str(c) + ' valid pairs from ' + args.i

if __name__ == '__main__':
	main()