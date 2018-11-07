#!/usr/bin/env python
import sys

IN = open('../ensGene_hg19.txt', 'r')
OUT = open('exons_per_kb.bed', 'w')

IN.readline()

for line in IN:
	splitline=line.split('\t')
	chrom = splitline[2]
	exonCount = float(splitline[8])
	txStart = int(splitline[4])
	txEnd = int(splitline[5])
	gene_length = int(txEnd) - int(txStart)
	if gene_length < 0:
		print 'ERROR: gene length'
		sys.exit()
	EPK = (exonCount/gene_length) * 1000
	OUT.write('\t'.join(map(str, [chrom, txStart, txEnd, EPK])) + '\n')
IN.close()
OUT.close()