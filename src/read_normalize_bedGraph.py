#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='normalize input bedGraph by input number of total reads')
parser.add_argument('-i', help='input bedGraph', type=str, required=True)
parser.add_argument('-t', help='total number of reads', type=float, required=True)
args = parser.parse_args()

def main():
	with open(args.i, 'r') as f:
		for line in f:
			splitline = line.split('\t')
			signal = float(splitline[3])
			norm_signal = signal / args.t
			print '\t'.join(splitline[0:3] + [str(norm_signal)])

if __name__ == '__main__':
	main()
	