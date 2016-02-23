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

# Get all bin indices
bin_ind = list(l_df.loc[l_df[0] == list(l_df[0])[0], 1])

# Initialize my5C matrix
# dim = str(len(bin_ind)) + 'x' + str(len(bin_ind))

columns = []
for bin in bin_ind:
	colname = '|'.join([args.n, args.a, 'chr'+args.c+':'+str((bin + 1)*1000)+'-'+str((bin+args.b)*1000)])
	columns.append(colname)

my5c_df = pd.DataFrame(columns=columns, index=columns)

# Transform data
for i, row in l_df.iterrows():
	my5c_df.iloc[bin_ind.index(row[0]), bin_ind.index(row[1])] = row[2]
	# for symmetry
	my5c_df.iloc[bin_ind.index(row[1]), bin_ind.index(row[0])] = row[2]

# Write to file
	my5c_df.to_csv(args.o + '.matrix', sep='\t')



