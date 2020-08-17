#!/usr/bin/env python
import argparse

parser=argparse.ArgumentParser(description='transform *std.bedGraph to *variance.bedGraph')
parser.add_argument('-i', help='input *std.bedGraph file', type=str, required=True)
args=parser.parse_args()



def main():
	with open(args.i, 'r') as f:
		for line in f:
			s = line.strip().split('\t')
			std = s[3]
			if std == 'nan':
				print line.strip()
			else:
				std = float(std)
				var = std*std
				l = s[:3] + [str(var)]
				print '\t'.join(l)

if __name__ == '__main__':
	main()