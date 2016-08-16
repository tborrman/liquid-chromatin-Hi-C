#!/usr/bin/env python
import argparse

parser=argparse.ArgumentParser(description='bin ChIP-seq bed file by summing signal')
parser.add_argument('-i', help='ChIP-seq bed file (ex. H3K9me3_sorted_clean.bed)', type=str, dest='i', required=True)
parser.add_argument('-g', help='genome bin file (ex. hg19_500kb.bed)', type=str, dest='g', required=True)

args=parser.parse_args()

def main():
	
	GENOME = open(args.g, 'r')
	for gen_line in GENOME:
		chrom = gen_line.split()[0]
		start = int(gen_line.split()[1])
		end = int(gen_line.split()[2])
		bin_signals = []
		ChIP = open(args.i, 'r')
		for line in ChIP:
			splitline = line.split()
			if (splitline[0] == chrom) and (int(splitline[1]) >= start) and (int(splitline[1]) < end):
				bin_signals.append(int(splitline[3]))
		ChIP.close()
		signal_sum = sum(bin_signals)
		print '\t'.join([chrom, str(start), str(end), str(signal_sum)])

	GENOME.close()






if __name__ == '__main__':
	main()