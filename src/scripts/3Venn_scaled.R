#!/usr/bin/env Rscript
library(argparse)
library(VennDiagram)

parser <- ArgumentParser(description= "Make scaled Venn Diagram plots for Rao_K562, DpnII_K562, HindIII_K562")
parser$add_argument("-r", help= "radius for merging in kb", type="integer", default=25)
args <- parser$parse_args()


count_overlap_pair <- function(F1_row, F2, radius) {
  # Get data
  chr_1 <- as.character(F1_row$chr1)
  chr_2 <- as.character(F1_row$chr2)
  x1 <- F1_row$x1
  y1 <- F1_row$y1
  
  # Count overlaps
  F1_F2_overlap <- any((chr_1 == F2$chr1) & (chr_2 == F2$chr2) & ((abs(x1 - F2$x1) <= radius)) & (abs(y1 - F2$y1) <= radius))

  if (F1_F2_overlap) {
    return(1)
  }
  else {
    return(0)
  }
}

count_overlap_all <- function(F1_row, F2, F3, radius) {
  # Get data
  chr_1 <- as.character(F1_row$chr1)
  chr_2 <- as.character(F1_row$chr2)
  x1 <- F1_row$x1
  y1 <- F1_row$y1
  # Count overlaps
  F1_F2_overlap <- any((chr_1 == F2$chr1) & (chr_2 == F2$chr2) & ((abs(x1 - F2$x1) <= radius)) & (abs(y1 - F2$y1) <= radius))
  F1_F3_overlap <- any((chr_1 == F3$chr1) & (chr_2 == F3$chr2) & ((abs(x1 - F3$x1) <= radius)) & (abs(y1 - F3$y1) <= radius))
  if (F1_F2_overlap & F1_F3_overlap) {
    return(1)
  }
  else {
    return(0)
  }
}


DpnII <- read.table('../houda_control_DpnII_K562/merged_loops', header=TRUE, sep="\t")
HindIII <- read.table('../houda_control_HindIII_K562/merged_loops', header=TRUE, sep="\t")
MboI <- read.table('../Rao_MboI_K562/merged_loops', header=TRUE, sep="\t")

total_DpnII <- length(DpnII$x1)
total_HindIII <- length(HindIII$x1)
total_MboI <- length(MboI$x1)


######################################################################################################
#Fig 1. DpnII, HindIII, MboI

DpnII_HindIII <- 0
DpnII_MboI <- 0
HindIII_MboI <- 0
all_3 <- 0

for (row_idx in 1:nrow(DpnII)) {
  DpnII_HindIII <- DpnII_HindIII + count_overlap_pair(DpnII[row_idx,], HindIII, args$r * 1000)
}



for (row_idx in 1:nrow(DpnII)) {
  DpnII_MboI <- DpnII_MboI + count_overlap_pair(DpnII[row_idx,], MboI, args$r * 1000)
}

for (row_idx in 1:nrow(HindIII)) {
  HindIII_MboI <- HindIII_MboI + count_overlap_pair(HindIII[row_idx,], MboI,  args$r * 1000)
}

for (row_idx in 1:nrow(DpnII)) {
  all_3 <- all_3 + count_overlap_all(DpnII[row_idx,], HindIII, MboI, args$r * 1000)
}
overrideTriple = 1
png("DpnII_HindIII_MboI.png", width=2500, height=2500, res=300)
draw.triple.venn(area1 = total_DpnII, area2 = total_HindIII, area3= total_MboI, 
                 n12=DpnII_HindIII, n23=HindIII_MboI, n13 = DpnII_MboI, n123=all_3,
                 category = c("Houda_DpnII_K562", "Houda_HindIII_K562", "Rao_MboI_K562"), 
                 fill = c("seagreen2", "orangered", "mediumorchid4"), scaled = TRUE)
dev.off()