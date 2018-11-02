#!/usr/bin/env python
import argparse
parser=argparse.ArgumentParser(description='Get positions of genes for RNAseq')
parser.add_argument('-i', help='RNAseqfile', type=str, required=True)
args=parser.parse_args()

index = args.i[-3:]
OUT=open('bedFPKM_' + index, 'w')
RNA_FH = open(args.i, 'r')
for i, RNAline in enumerate(RNA_FH):
	if i % 100 == 0:
		print 'On row: ' + str(i)
	if RNAline[:4] == 'ENSG':
		ensgid = RNAline[:15]
		GENES = open('ensGene_hg19.txt', 'r')
		for GENEline in GENES:
			splitline = GENEline.split('\t')
			ensgid2 = splitline[12]
			if ensgid == ensgid2:
				FPKM = RNAline.split('\t')[6]
				OUT.write(splitline[2]+'\t'+splitline[4]+'\t'+splitline[5]+'\t'+FPKM+'\n')
				break
		GENES.close()
RNA_FH.close()
OUT.close()
		
	
