#!/usr/bin/env python
import argparse
import sys

parser=argparse.ArgumentParser(description='Add NAs to bins based on NAs assigned from half-life analysis')
parser.add_argument('-i', help='bedGraph file to remove outliers from', type=str, required=True)
parser.add_argument('-hd', help='Header exists for input file', type=bool, default=False)
args = parser.parse_args()

def check_coordinates(i, j):
	# Check if coordinates match
	coord1 = i.split('\t')[:3]
	coord2 = j.split('\t')[:3]
	if coord1 == coord2:
		return True
	else:
		return False

def main():

	HL = open('/home/Tyler/Research/digest/cis_percent/timecourse/C-40000/half-life_LOS/half-life_exponential_40kb_removed_outliers.bed', 'r')
	halflifes = HL.readlines()
	IN = open(args.i, 'r')
	OUT = open(args.i[:-9] + '_removed_outliers.bedGraph', 'w')
	if args.hd:
		IN.readline()
	inlines = IN.readlines()
	for i in range(len(inlines)):
		if check_coordinates(inlines[i], halflifes[i]):
			if halflifes[i].strip().split('\t')[3] == 'NA':
				OUT.write(halflifes[i])
			else:
				OUT.write(inlines[i])
		else:
			print 'ERROR: mismatching coordinates in files'
			sys.exit()

if __name__ == '__main__':
	main()

