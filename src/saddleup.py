#!/usr/bin/env python
import cooler
import my_saddles
import pandas as pd
import numpy as np


def format_loc(df, bp_col, res):
	'''
	args
	df: pandas df with cols chrom, bp
	bp_col: name of bp column
	res: int resolution

	return
	list of loci formatted chrom:start-end
	'''
	chrom = list(df['chrom'])
	end = df[bp_col] + (res -1)
	start = map(str, list(df[bp_col]))
	end = map(str, list(end))
	loc = [c + ':' +  start[i] + '-' + end[i] for i, c in enumerate(chrom)]
	return loc




def main():
	res = 500000
	my_file = 'eigenI_eigenJ_MN_interactome_chr14.txt'

	# bins
	df  = pd.read_csv(my_file, sep='\t')
	print df.tail()
	my_bins = df.loc[:,['chrom', 'I']]
	my_bins['end'] = my_bins['I'] + (res - 1)
	my_bins.columns = ['chrom', 'start', 'end']
	my_bins = my_bins.sort_values('start')
	my_bins.index = range(my_bins.shape[0])

	# pixels
	#df.columns = ['chrom', 'start', 'start2', 'count', 'Ei', 'Ej']
	df = df.sort_values(['chrom', 'I', 'J'])
	df.index = range(df.shape[0])
	loc1 = format_loc(df, 'I', res)
	loc2 = format_loc(df, 'J', res)
	x = np.array([np.array(loc1), np.array(loc2), df['cScore']]).T
	y = pd.DataFrame(x, columns = ['bin1_id', 'bin2_id', 'count'])

	print my_bins.iloc[1234]
	print y.iloc[1234]
	
	
	z = my_bins.astype(str)
	my_bins['chrom'] = my_bins['chrom'].astype('|S')
	print my_bins['chrom'].dtype

	cooler.io.create('test.cool', my_bins, y)
	



if __name__ == '__main__':
	main()


	