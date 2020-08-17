#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Convert from Lieberman Aiden matrix format (Rao et al 2014) to my5C matrix format')
parser.add_argument('-l', help='Lieberman Aiden format matrix', type=str, required=True, dest='l')
parser.add_argument('-n', help='name (ex. K562)', type=str, required=True, dest='n')
parser.add_argument('-a', help='assembly (ex. hg19)', type=str, required=True, dest='a')
parser.add_argument('-b', help='bin size in kb (ex. 500)', type= int, required=True, dest='b')
parser.add_argument('-c', help='chromosome (ex. 1 or Y)', type=str, required=True, dest='c')
parser.add_argument('-o', help='output file prefix', type=str, required=True, dest='o')



args = parser.parse_args()

# Open lieberman matrix
l_df = pd.read_csv(args.l, sep='\t', header=None)
l_df  = l_df.sort_values(by=[0,1])

last_bin =  int(l_df.loc[len(l_df[0]) - 1,0])
first_bin = int(l_df.loc[0,0])
bins = range(first_bin, last_bin + 1, args.b * 1000)

columns = []
for bin in bins:
	colname = '|'.join([args.n, args.a, 'chr'+args.c+':'+str(bin + 1)+'-'+str(bin+(args.b * 1000))])
	columns.append(colname)

my5c_df = pd.DataFrame(0, columns=columns, index=columns)

# Transform data
counter = 0
for i, row in l_df.iterrows():
	if counter % 10 == 0:
		print 'On row: ' + str(counter)
	x = '|'.join([args.n, args.a, 'chr'+args.c+':'+str(int(row[0]) + 1)+'-'+str(int(row[0])+(args.b * 1000))])
	y = '|'.join([args.n, args.a, 'chr'+args.c+':'+str(int(row[1]) + 1)+'-'+str(int(row[1])+(args.b * 1000))])
	my5c_df.ix[x, y] = row[2]
	# for symmetry
	my5c_df.ix[y, x] = row[2]
	counter += 1

# Write to file
my5c_df.to_csv(args.o + '.matrix', sep='\t', na_rep='nan')



