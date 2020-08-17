#!/usr/bin/env python

import argparse
from Bio import SeqIO, Restriction
from Bio.Alphabet import IUPAC

parser=argparse.ArgumentParser(description='record FatI sites for given chromosome FASTA in bed file format')
parser.add_argument('-f', help='chromsome FASTA file (ex. chr21.fa)', type=str, dest='f', required=True)
args = parser.parse_args()

seq_record = SeqIO.read(args.f, "fasta", IUPAC.ambiguous_dna)
coords = Restriction.FatI.search(seq_record.seq)

chrom = seq_record.id
OUT=open('FatI_sites_' + chrom + '.bed', 'w')
for start in coords:
	OUT.write('\t'.join([chrom, str(start-1), str((start-1) + 4)]) + '\n')
OUT.close()

