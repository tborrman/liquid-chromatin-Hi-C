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
		('half-life/half-life_exponential_40kb_removed_outliers_range6Mb.bedGraph', 'half-life_LOS'),
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
		('ChIP-seq/ENCFF778DNU_H3K4me2_R1_40kb.bedGraph', 'H3K4me2_R1'),
		('ChIP-seq/ENCFF622YGE_H3K4me2_R2_40kb.bedGraph', 'H3K4me2_R2'),
		('ChIP-seq/ENCFF330TJF_H3K4me3_R1_40kb.bedGraph', 'H3K4me3_R1'),
		('ChIP-seq/ENCFF803FOA_H3K4me3_R2_40kb.bedGraph', 'H3K4me3_R2'),
		('ChIP-seq/ENCFF242XLB_H4K20me1_R1_40kb.bedGraph', 'H4K20me1_R1'),
		('ChIP-seq/ENCFF474KOM_H4K20me1_R2_40kb.bedGraph', 'H4K20me1_R2'),
		('ChIP-seq/ENCFF430FSQ_H2AFZ_R1_40kb.bedGraph', 'H2AFZ_R1'),
		('ChIP-seq/ENCFF551PZO_H2AFZ_R2_40kb.bedGraph', 'H2AFZ_R2'),
		('ChIP-seq/ENCFF526UWC_H3K9me1_R1_40kb.bedGraph', 'H3K9me1_R1'),
		('ChIP-seq/ENCFF700FQH_H3K9me3_R1_40kb.bedGraph', 'H3K9me3_R1'),
		('ChIP-seq/ENCFF742JJQ_H3K9me3_R2_40kb.bedGraph', 'H3K9me3_R2'),
		('ChIP-seq/ENCFF527JUP_H3K9ac_R1_40kb.bedGraph', 'H3K9ac_R1'),
		('ChIP-seq/ENCFF098LBU_H3K9ac_R2_40kb.bedGraph', 'H3K9ac_R2'),
		('ChIP-seq/ENCFF936BVT_H3K27me3_R1_40kb.bedGraph', 'H3K27me3_R1'),
		('ChIP-seq/ENCFF274JTE_H3K27me3_R2_40kb.bedGraph', 'H3K27me3_R2'),
		('ChIP-seq/ENCFF619JRY_H3K79me2_R1_40kb.bedGraph', 'H3K79me2_R1'),
		('ChIP-seq/ENCFF453XIK_H3K79me2_R2_40kb.bedGraph', 'H3K79me2_R2'),
		('ChIP-seq/ENCFF545UAE_CTCF_R1_40kb.bedGraph', 'CTCF_R1'),
		('ChIP-seq/ENCFF061IVP_CTCF_R2_40kb.bedGraph', 'CTCF_R2'),
		('ChIP-seq/ENCFF253VJL_SMC3_R1_40kb.bedGraph', 'SMC3_R1'),
		('ChIP-seq/ENCFF003UFU_SMC3_R2_40kb.bedGraph', 'SMC3_R2'),
		('ChIP-seq/ENCFF000YXZ_RAD21_40kb.bedGraph', 'RAD21'),
		('DNase-seq/ENCFF111KJD_DNase_R1_40kb.bedGraph', 'DNase-seq_R1'),
		('DNase-seq/ENCFF456JVK_DNase_R2_40kb.bedGraph', 'DNase-seq_R2'),
		('DpnII-seq/HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_copy_correct_coverage_40kb.bed', 'DpnII-seq'),
		('lamina/KBM7/Lamina_mean_OE_Clone.14.1N.OE_40kb.bedGraph', 'LAD_KBM7_clone14'),
		('lamina/KBM7/Lamina_mean_OE_Clone.5-5.1N.OE_40kb.bedGraph', 'LAD_KBM7_clone5-5'),
		('lamina/K562/K562_LmnB1_DamID_log2ratios_hg19_40kb.bedGraph', 'LAD_K562'),
		('loops/houda_control_DpnII_K562_anchor_list_40kb.bedGraph', 'loops'),
		('loops/GSE63525_K562_HiCCUPS_anchor_list_40kb.bedGraph', 'loops_Rao'),
		('eigen/eigen1_40kb.bedGraph', 'PCA_eigen1'),
		('ChIP-seq/ENCFF985MBQ_CBX1_R1_40kb.bedGraph', 'CBX1_R1'),
		('ChIP-seq/ENCFF824OMP_CBX1_R2_40kb.bedGraph', 'CBX1_R2'),
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
		('ChIP-seq/remodellers/ENCFF734YKJ_KDM1A_R1_R2_40kb.bedGraph', 'KDM1A_Snyder1'),
		('ChIP-seq/remodellers/ENCFF758MEL_KDM1A_R1_R2_40kb.bedGraph', 'KDM1A_Snyder2'),
		('ChIP-seq/PolII/ENCFF749YKR_POLR2A_R1_40kb.bedGraph', 'POLR2A'),
		('ChIP-seq/PolII/ENCFF452WGO_POLR2B_R1_40kb.bedGraph', 'POLR2B'),
		('ChIP-seq/PolII/ENCFF015NSS_POLR2G_R1_40kb.bedGraph', 'POLR2G'),
		('tad/Houda_Ctrl_DpnII_K562.40000_allchr.txt--is480000--nt0.1--ids320000--ss40000--immean.insulation.bpfix.bed', 'TAD_insulation'),
		('Repli-Seq/ENCFF001GRX_G1_40kb.bedGraph', 'G1_Repli-seq'),
		('Repli-Seq/ENCFF001GSF_S1_40kb.bedGraph', 'S1_Repli-seq'),
		('Repli-Seq/ENCFF001GSJ_S2_40kb.bedGraph', 'S2_Repli-seq'),
		('Repli-Seq/ENCFF001GSN_S3_40kb.bedGraph', 'S3_Repli-seq'),
		('Repli-Seq/ENCFF001GSP_S4_40kb.bedGraph', 'S4_Repli-seq'),
		('Repli-Seq/ENCFF001GSB_G2_40kb.bedGraph', 'G2_Repli-seq'),
		('WGBS/ENCFF867JRG_WGBS_R1_liftOverhg19_40kb.bedGraph', 'WGBS_R1'),
		('WGBS/ENCFF721JMB_WGBS_R2_liftOverhg19_40kb.bedGraph', 'WGBS_R2'),
		('NADs/NADs_HeLa_2010_40kb.bedGraph', 'NADs_HeLA'),
		('NADs/NADs_IMR90_2017_40kb.bedGraph', 'NADs_IMR90'),
		('RNA-seq/ENCFF039CTL_RNA-seq_polyA_plus_R1_40kb.bedGraph', 'RNA-seq_polyA_+_R1'),
		('RNA-seq/ENCFF809GEF_RNA-seq_polyA_plus_R2_40kb.bedGraph', 'RNA-seq_polyA_+_R2'),
		('RNA-seq/ENCFF960DJT_RNA-seq_polyA_minus_R1_40kb.bedGraph', 'RNA-seq_polyA_-_R1'),
		('RNA-seq/ENCFF511CNS_RNA-seq_polyA_minus_R2_40kb.bedGraph', 'RNA-seq_polyA_-_R2'),
		('RNA-seq/ENCFF717RYC_RNA-seq_total_plus_R1_40kb.bedGraph', 'RNA-seq_total_+_R1'),
		('RNA-seq/ENCFF091RAW_RNA-seq_total_plus_R2_40kb.bedGraph', 'RNA-seq_total_+_R2'),
		('RNA-seq/ENCFF633AHL_RNA-seq_total_minus_R1_40kb.bedGraph', 'RNA-seq_total_-_R1'),
		('RNA-seq/ENCFF078XGZ_RNA-seq_total_minus_R2_40kb.bedGraph', 'RNA-seq_total_-_R2'),
		('gene_density/gene_density_40kb.bedGraph', 'gene_density'),
		('exons/exons_per_kb_40kb.bedGraph', 'exons_per_kb'),
		('TSA-seq/GSE81553_SON_TSA-Seq_Condition2_40kb.bedGraph', 'SON_TSA-seq'),
		('TSA-seq/GSE81553_LaminAC_TSA-Seq_Condition2_40kb.bedGraph', 'LaminAC_TSA-seq'),
		('TSA-seq/GSE81553_LaminB_TSA-Seq_Condition2_40kb.bedGraph', 'LaminB_TSA-seq'),
		('TSA-seq/GSE81553_Pol2_TSA-Seq_Condition1_40kb.bedGraph', 'Pol2_TSA-seq'),
		('PML/ENCFF157IUB_PML_R1_40kb.bedGraph', 'PML_R1'),
		('PML/ENCFF412XML_PML_R2_40kb.bedGraph', 'PML_R2')
		]
	IN = open('hg19_40kb.bed', 'r')
	OUT = open('feature_matrix_v5_40kb.txt', 'w')
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
				print line
				print fline
				print f
				sys.exit()

	IN.close()
	OUT.close()
		



if __name__ == '__main__':
	main()

