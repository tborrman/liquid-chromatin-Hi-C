#!/usr/bin/env python2
# Identical to deduplication step of cMapping pipeline found here:
# https://github.com/dekkerlab/cMapping/blob/4126432e1130b9a78f895dfa1e0f5c90c08cd406/scripts/perl/collapseHiCValidPairs.pl
import argparse
import gzip
import re
import sys


def check_zip(filename):
	'''
	Check if input file is gzipped or not
	and return filehandle 
	
	Args:
		filename: string
	Returns:
		FH: output filehandle 
	'''
	
	if filename[-3:] == '.gz':
		FH = gzip.open(filename, 'rb')
	else:
		FH = open(filename, 'r')
	return FH


def dedup(f):
	''' 
	Run deduplication of validPairs file
	
	Args:
		f: validPair file
	Returns:
		validPairCount: original number of valid pairs
		interactionTotal: number of interactions after duplicates removed
		moleculeTotal: total molecules (should equal validPairCount)
		totalRedundantMolecules: total duplicates removed
	'''
	previousInteractionKey=''
	previousMoleculeKey=''
	previousLine = ''
	interactionLines = []

	interactionCount=0
	moleculeCount=0
	validPairCount=0
	interactionTotal=0
	moleculeTotal=0
	totalRedundantMolecules=0

	IN = check_zip(f)
	if f[-3:] == '.gz':
		OUT = gzip.open(f[:-3] + '_dedup.gz', 'wb')
	else:
		OUT = gzip.open(f + '_dedup.gz', 'wb')

	lineNum = 1
	for line in IN:

		if line[0] == '#' or line == '':
			continue
		splitline = line.strip().split('\t')
		if len(splitline) != 12:
			print 'Error with input file! Must have 12 columns!'
			print splitline
			sys.exit()

		chromosome_1=splitline[1]
		readPos_1=splitline[2]
		strand_1=splitline[3]
		readID_1=splitline[4]
		fragmentIndex_1=int(splitline[5])
		
		chromosome_2=splitline[7]
		readPos_2=splitline[8]
		strand_2=splitline[9]
		readID_2=splitline[10];
		fragmentIndex_2=int(splitline[11])
		
		validPairCount += 1

		if fragmentIndex_1 > fragmentIndex_2:
			print 'ERROR: line ' + str(lineNum) + ': frag1 > frag2'
			print str(fragmentIndex_1) + ' > ' + str(fragmentIndex_2)
			print line
			sys.exit()

		interactionKey=str(fragmentIndex_1) + '\t' + str(fragmentIndex_2)
		moleculeKey=str(fragmentIndex_1) + '@' + readPos_1 + '\t' + str(fragmentIndex_2) + '@' + readPos_2


		if previousInteractionKey == '':
			previousInteractionKey = interactionKey

		if interactionKey == previousInteractionKey:
			if moleculeKey != previousMoleculeKey:
				interactionCount += 1
				interactionLines.append(line)
			moleculeCount += 1
		else:
			OUT.write(''.join(interactionLines))
			interactionLines = []
			interactionLines.append(line)
			interactionTotal += interactionCount
			moleculeTotal += moleculeCount
			numRedundantMolecules = moleculeCount-interactionCount
			interactionCount = 1
			moleculeCount = 1
			totalRedundantMolecules += numRedundantMolecules


		previousInteractionKey=interactionKey
		previousMoleculeKey=moleculeKey
		lineNum += 1

	OUT.write(''.join(interactionLines))
	interactionTotal += interactionCount
	moleculeTotal += moleculeCount
	numRedundantMolecules = moleculeCount - interactionCount
	totalRedundantMolecules += numRedundantMolecules

	IN.close()
	OUT.close()
	
	return(validPairCount, interactionTotal, moleculeTotal, totalRedundantMolecules)

def main():

	parser=argparse.ArgumentParser('Remove duplicates from cMapping validPair file')
	parser.add_argument('-i', help='validPair file', type=str, required=True)
	args=parser.parse_args()

	v, i, m, t, = dedup(args.i)
	print '*****STATISTICS******'
	print 'Valid Pairs: ' + str(v)
	print 'Interaction Total: ' + str(i)
	print 'Molecule Total: ' + str(m)
	print 'Total Redundant Molecules: ' + str(t)
	print '*********************'

if __name__ == '__main__':
	main()

