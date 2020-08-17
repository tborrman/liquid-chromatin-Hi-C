#!/usr/bin/env python
import sys

def main():
	IN = open('test_sub.txt', 'r')

	lines = IN.readlines()
	splitline = lines[0].strip().split('\t')
	current_chrom = splitline[0]
	current_start = splitline[1]
	current_end = splitline[2]
	current_sub = splitline[3]
	for i in range(len(lines) - 1):
		nextline = lines[i+1].strip().split('\t')
		next_chrom = nextline[0]
		next_start = nextline[1]
		next_end = nextline[2]
		next_sub = nextline[3]
		if (current_chrom != next_chrom) or (current_sub != next_sub):
			print '\t'.join([current_chrom, current_start, current_end, current_sub])
			# Update
			current_chrom = next_chrom
			current_start = next_start
			current_end = next_end
			current_sub = next_sub
			if i == (len(lines) - 2):
				print '\t'.join([current_chrom, current_start, current_end, current_sub])
		elif current_sub == next_sub:
			current_end = next_end
			if i == (len(lines) - 2):
				print '\t'.join([current_chrom, current_start, current_end, current_sub])
		else:
			print 'WTF!?'
			sys.exit()

if __name__ == '__main__':
	main()

