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
		('half-life/half-life_exponential_40kb.bed', 'half-life_LOS'),
		('half-life_std/half-life_std_exponential_40kb.bed', 'half-life_std'),
		('LOS/HBHiCK562DN10-5m-DPnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph_noheader', 'LOS_5m'),
		('LOS/HBHiCK562DN10-1h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph_noheader', 'LOS_1h'),
		('LOS/HBHiCK562DN10-2h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph_noheader', 'LOS_2h'),
		('LOS/HBHiCK562DN10-3h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph_noheader', 'LOS_3h'),
		('LOS/HBHiCK562DN10-4h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph_noheader', 'LOS_4h'),
		('LOS/HBHiCK562DN10-ON-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph_noheader', 'LOS_ON'),
		('ChIP-seq/ENCFF223BKS_H3K36me3_R1_40kb.bedGraph', 'H3K36me3_R1'),
		('ChIP-seq/ENCFF013QMA_H3K36me3_R2_40kb.bedGraph', 'H3K36me3_R2'),
		('ChIP-seq/ENCFF006RIO_H3K27ac_R1_40kb.bedGraph', 'H3K27ac_R1'),
		('ChIP-seq/ENCFF685PCO_H3K27ac_R2_40kb.bedGraph', 'H3K27ac_R2'),
		('ChIP-seq/ENCFF463AQS_H3K4me1_R1_40kb.bedGraph', 'H3K4me1_R1'),
		('ChIP-seq/ENCFF624RNZ_H3K4me1_R2_40kb.bedGraph', 'H3K4me1_R2'),
		('ChIP-seq/ENCFF242XLB_H4K20me1_R1_40kb.bedGraph', 'H4K20me1_R1'),
		('ChIP-seq/ENCFF474KOM_H4K20me1_R2_40kb.bedGraph', 'H4K20me1_R2'),
		('ChIP-seq/ENCFF430FSQ_H2AFZ_R1_40kb.bedGraph', 'H2AFZ_R1'),
		('ChIP-seq/ENCFF551PZO_H2AFZ_R2_40kb.bedGraph', 'H2AFZ_R2'),
		('ChIP-seq/ENCFF700FQH_H3K9me3_R1_40kb.bedGraph', 'H3K9me3_R1'),
		('ChIP-seq/ENCFF742JJQ_H3K9me3_R2_40kb.bedGraph', 'H3K9me3_R2'),
		('ChIP-seq/ENCFF936BVT_H3K27me3_R1_40kb.bedGraph', 'H3K27me3_R1'),
		('ChIP-seq/ENCFF274JTE_H3K27me3_R2_40kb.bedGraph', 'H3K27me3_R2'),
		('ChIP-seq/ENCFF545UAE_CTCF_R1_40kb.bedGraph', 'CTCF_R1'),
		('ChIP-seq/ENCFF061IVP_CTCF_R2_40kb.bedGraph', 'CTCF_R2'),
		('DNase-seq/ENCFF111KJD_DNase_R1_40kb.bedGraph', 'DNase-seq_R1'),
		('DNase-seq/ENCFF456JVK_DNase_R2_40kb.bedGraph', 'DNase-seq_R2'),
		('DpnII-seq/HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_copy_correct_coverage_40kb.bed', 'DpnII-seq'),
		('lamina/Lamina_mean_OE_Clone.14.1N.OE_40kb.bedGraph', 'LAD_clone14'),
		('lamina/Lamina_mean_OE_Clone.5-5.1N.OE_40kb.bedGraph', 'LAD_clone5-5'),
		('loops/houda_control_DpnII_K562_anchor_list_40kb.bedGraph', 'loops'),
		('loops/GSE63525_K562_HiCCUPS_anchor_list_40kb.bedGraph', 'loops_Rao'),
		('eigen/eigen1_40kb.bedGraph', 'PCA_eigen1'),
		('ChIP-seq/remodellers/ENCFF358INE_CBX3_R1_40kb.bedGraph','CBX3_Bernstein'),
		('ChIP-seq/remodellers/ENCFF630YDI_CBX3_R1_R2_40kb.bedGraph','CBX3_Myers'),		
		('ChIP-seq/remodellers/ENCFF014XHT_CBX5_R1_R2_40kb.bedGraph', 'CBX5'),
		('ChIP-seq/remodellers/ENCFF982RMW_EHMT2_R1_R2_40kb.bedGraph','EHMT2'),
		('ChIP-seq/remodellers/ENCFF206WVX_CBX8_R1_R2_40kb.bedGraph','CBX8'),
		('ChIP-seq/remodellers/ENCFF071CIY_RNF2_R1_R2_40kb.bedGraph','RNF2'),
		('ChIP-seq/remodellers/ENCFF514QBY_BMI1_R1_R2_40kb.bedGraph','BMI1'),
		('ChIP-seq/remodellers/ENCFF363DWX_SUZ12_R1_R2_40kb.bedGraph','SUZ12'),
		('ChIP-seq/remodellers/ENCFF058KTS_RBBP5_R1_R2_40kb.bedGraph', 'RBBP5'),
		('ChIP-seq/remodellers/ENCFF112LWM_CTBP1_R1_R2_40kb.bedGraph','CTBP1'),
		('ChIP-seq/remodellers/ENCFF164NLF_KAT2B_R1_40kb.bedGraph','KAT2B'),
		('ChIP-seq/remodellers/ENCFF260JHC_BRD4_R1_R2_40kb.bedGraph','BRD4'),
		('ChIP-seq/remodellers/ENCFF292YLV_NCOR1_R1_R2_40kb.bedGraph','NCOR1'),
		('ChIP-seq/remodellers/ENCFF293III_KDM5B_R1_R2_40kb.bedGraph','KDM5B'),
		('ChIP-seq/remodellers/ENCFF532SIS_HDAC2_R1_R2_40kb.bedGraph','HDAC2_Bernstein'),
		('ChIP-seq/remodellers/ENCFF915GWT_HDAC2_R1_R2_40kb.bedGraph','HDAC2_Snyder'),
		('ChIP-seq/remodellers/ENCFF816KCQ_SAP30_R1_R2_40kb.bedGraph','SAP30'),
		('ChIP-seq/remodellers/ENCFF864WOR_WHSC1_R1_40kb.bedGraph','WHSC1'),
		('ChIP-seq/remodellers/ENCFF427QTV_PHF8_R1_R2_40kb.bedGraph','PHF8'),
		('ChIP-seq/remodellers/ENCFF518QUW_REST_R1_R2_40kb.bedGraph','REST'),
		('tad/Houda_Ctrl_DpnII_K562.40000_allchr.txt--is480000--nt0.1--ids320000--ss40000--immean.insulation.bpfix.bed', 'TAD_insulation')
		]
	IN = open('hg19_40kb.bed', 'r')
	OUT = open('feature_matrix_40kb.txt', 'w')
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
				print f
				sys.exit()

	IN.close()
	OUT.close()
		



if __name__ == '__main__':
	main()

