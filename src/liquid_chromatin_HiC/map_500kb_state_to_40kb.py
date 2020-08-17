#!/usr/bin/env python
import argparse


def main():
	parser=argparse.ArgumentParser(description='Assign 500kb copy number states to 40kb bins')
	parser.add_argument('-i', help='bedtools intersect file (ex. intersect_40kb_500kb.txt)', type=str, required=True)
	args = parser.parse_args()
	print 'chrom\tstart\tend\tcopy'
	i = 0
	IN = open(args.i, 'r')
	lines = IN.readlines()
	while i < len(lines):
		if (i + 1) < len(lines):
			start = lines[i].split()[1]
			nextstart = lines[i+1].split()[1]
			if start == nextstart:
				copy = lines[i].split()[6]
				nextcopy = lines[i+1].split()[6]
				if copy == nextcopy:
					print '\t'.join(lines[i].split()[:3] + [copy])
				else:
					print '\t'.join(lines[i].split()[:3] + ['NA'])
				i += 2
			else:
				copy = lines[i].split()[6]
				print '\t'.join(lines[i].split()[:3] + [copy])
				i +=1
		# Last line
		else:
			copy = lines[i].split()[6]
			print '\t'.join(lines[i].split()[:3] + [copy])
			i +=1

if __name__ == '__main__':
	main()