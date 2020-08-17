#!/usr/bin/env python
import linecache
import sys

def check_coordinates(l, f):
	coord1 = l.split()[:3]
	coord2 = f.split()[:3]
	if coord1 == coord2:
		return True
	else:
		return False

def main():

	features = [
		('half-life/half-life_500kb_NA.bed', 'half-life'),
		('ChIP-seq/ENCFF223BKS_H3K36me3_R1_500kb.bedGraph', 'H3K36me3_R1'),
		('ChIP-seq/ENCFF013QMA_H3K36me3_R2_500kb.bedGraph', 'H3K36me3_R2'),
		('ChIP-seq/ENCFF006RIO_H3K27ac_R1_500kb.bedGraph', 'H3K27ac_R1'),
		('ChIP-seq/ENCFF685PCO_H3K27ac_R2_500kb.bedGraph', 'H3K27ac_R2'),
		('ChIP-seq/ENCFF463AQS_H3K4me1_R1_500kb.bedGraph', 'H3K4me1_R1'),
		('ChIP-seq/ENCFF624RNZ_H3K4me1_R2_500kb.bedGraph', 'H3K4me1_R2'),
		('ChIP-seq/ENCFF242XLB_H4K20me1_R1_500kb.bedGraph', 'H4K20me1_R1'),
		('ChIP-seq/ENCFF474KOM_H4K20me1_R2_500kb.bedGraph', 'H4K20me1_R2'),
		('ChIP-seq/ENCFF430FSQ_H2AFZ_R1_500kb.bedGraph', 'H2AFZ_R1'),
		('ChIP-seq/ENCFF551PZO_H2AFZ_R2_500kb.bedGraph', 'H2AFZ_R2'),
		('ChIP-seq/ENCFF700FQH_H3K9me3_R1_500kb.bedGraph', 'H3K9me3_R1'),
		('ChIP-seq/ENCFF742JJQ_H3K9me3_R2_500kb.bedGraph', 'H3K9me3_R2'),
		('ChIP-seq/ENCFF936BVT_H3K27me3_R1_500kb.bedGraph', 'H3K27me3_R1'),
		('ChIP-seq/ENCFF274JTE_H3K27me3_R2_500kb.bedGraph', 'H3K27me3_R2'),
		('ChIP-seq/ENCFF545UAE_CTCF_R1_500kb.bedGraph', 'CTCF_R1'),
		('ChIP-seq/ENCFF061IVP_CTCF_R2_500kb.bedGraph', 'CTCF_R2'),
		('DNase-seq/ENCFF111KJD_DNase_R1_500kb.bedGraph', 'DNase-seq_R1'),
		('DNase-seq/ENCFF456JVK_DNase_R2_500kb.bedGraph', 'DNase-seq_R2'),
		('DpnII-seq/DpnII-seq_copy_correct_500kb_NA.bed', 'DpnII-seq'),
		('lamina/Lamina_mean_OE_Clone.14.1N.OE_500kb.bedGraph', 'LAD_clone14'),
		('lamina/Lamina_mean_OE_Clone.5-5.1N.OE_500kb.bedGraph', 'LAD_clone5-5'),
		('loops/houda_control_DpnII_K562_anchor_list_500kb.bedGraph', 'loops'),
		('loops/GSE63525_K562_HiCCUPS_anchor_list_500kb.bedGraph', 'loops_Rao'),
		('eigen/HBHiC-K562-MN-Dp-1__hg19__genome__C-500000-raw_scaleBy_2.72.balanced_scaleBy_51.45__all.zScore.eigen1.sorted_bpfix.bedGraph', 'PCA_eigen1')
		]
	IN = open('hg19_500kb.bed', 'r')
	OUT = open('feature_matrix.txt', 'w')
	header = ['chrom', 'start', 'end']
	for f in features:
		header.append(f[1])
	OUT.write('\t'.join(header) + '\n')
	for i, line in enumerate(IN):
		OUT.write(line.strip() + '\t')
		for f in features:
			fline = linecache.getline(f[0], i + 1)
			if check_coordinates(line, fline):
				if features[-1] == f:
					OUT.write(fline.split()[3] + '\n')
				else:
					OUT.write(fline.split()[3] + '\t')		
			else:
				print 'ERROR in coordinates'
				sys.exit()

	IN.close()
	OUT.close()
		



if __name__ == '__main__':
	main()

