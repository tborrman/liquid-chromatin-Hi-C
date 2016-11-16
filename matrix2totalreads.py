#!/usr/bin/env python
import argparse
import numpy as np
import gzip

parser=argparse.ArgumentParser(description='Input a matrix and output count of total interactions in upper triangle including diagonal')
parser.add_argument('-i', help='Input matrix', type=str, dest='i', required=True)
args=parser.parse_args()

if args.i[-3:] == '.gz':
	IN = gzip.open(args.i, 'rb')
else:
	IN = open(args.i, 'r')
# Count lines to skip in header
counter=0
for line in IN:
    if '#' in line:
        counter += 1

x=np.genfromtxt(args.i, comments='#', delimiter='\t', skip_header=counter + 1)
up_tri = np.triu(x)
y = np.nansum(up_tri, axis=1)
print 'Total reads: ' + str(np.nansum(y))



