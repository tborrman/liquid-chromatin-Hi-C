#!/usr/bin/env python

import argparse
from Bio import SeqIO, Restriction
from Bio.Alphabet import IUPAC

parser=argparse.ArgumentParser(description='record hindIII sites for given chromosome FASTA in bed file format')
parser.add_argument('-f', help='chromsome FASTA file (ex. chr21.fa)', type=str, dest='f', required=True)
args = parser.parse_args()

seq_record = SeqIO.read(args.f, "fasta", IUPAC.ambiguous_dna)
coords = Restriction.HindIII.search(seq_record.seq)

chrom = seq_record.id
OUT=open('hindIII_sites_' + chrom + '.bed', 'w')
for start in coords:
	# Note: compensate for search function finding first base after the 
	# position the enzyme will cut.
	OUT.write('\t'.join([chrom, str(start-2), str((start-2) + 6)]) + '\n')
OUT.close()

