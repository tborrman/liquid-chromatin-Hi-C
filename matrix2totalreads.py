#!/usr/bin/env python
import argparse
import numpy as np

parser=argparse.ArgumentParser(description='Input a matrix and output count of total reads')
parser.add_argument('-i', help='Input matrix', type=str, dest='i', required=True)
args=parser.parse_args()

# Count lines to skip in header
IN = open(args.i, 'r')
counter=0
for line in IN:
    if '#' in line:
        counter += 1

x=np.genfromtxt(args.i, comments='#', delimiter='\t', skip_header=counter + 1, missing_values='nan')
y = np.nansum(x, axis=1)
print 'Total reads: ' + str(np.nansum(y))



