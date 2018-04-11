#!/usr/bin/env python
import argparse
import sys

parser=argparse.ArgumentParser(description='Correct final dpnII/hindIII-seq bedfile by K562 copy number')
parser.add_argument('-i', help='input dpnII/hindIII-seq file (ex. dpnII_3_distance_coverage_sorted_500kb_R1.bed)',
	type=str, required=True)
parser.add_argument('-c', help='copy number file (ex. K562_copynumber_500kb.bed)', type=str, required=True)
args = parser.parse_args()


def check_total_bins(s, c):
	SEQ = open(s, 'r')
	COP = open(c, 'r')
	seqTotal = len(SEQ.readlines())
	copyTotal = len(COP.readlines())
	if seqTotal != copyTotal:
		print 'ERROR: unequal number of bins'
		sys.exit()
	SEQ.close()
	COP.close()
	return

def check_bins_match(s, c):
	
	if s.split()[:3] != c.split()[:3]:
		print 'ERROR: mismatching bins'
		sys.exit()
	return

def correct_reads(r, c):
	if c == 'NA':
		correct = 'NA'
	elif c == '2':
		correct = round(int(r)/1.0, 2)
	elif c == '3':
		correct = round(int(r)/1.5, 2)
	elif c == '4':
		correct = round(int(r)/2.0, 2)
	else:
		print 'ERROR: unknown copy number'
		sys.exit()
	return correct

	

def main():
	# Check files have the same number of bins
	seqfile = args.i
	copyfile = args.c
	check_total_bins(seqfile, copyfile)
	SEQFH = open(seqfile, 'r')
	seqlines = SEQFH.readlines()
	COPYFH = open(copyfile, 'r')
	clines = COPYFH.readlines()
	for i in range(len(seqlines)):
		check_bins_match(seqlines[i], clines[i])
		reads = seqlines[i].split()[3]
		copynum = clines[i].split()[3]
		corr_reads = correct_reads(reads, copynum)
		print '\t'.join(seqlines[i].split()[:3] + [str(corr_reads)])

	SEQFH.close()
	COPYFH.close()

if __name__ == '__main__':
	main()