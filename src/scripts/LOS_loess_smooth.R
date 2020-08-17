#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Loess smooth LOS data")
parser$add_argument("-i", help="LOS bedGraph", type="character", dest="i", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()


get_loess_prediction <- function(x, y) {
  # Calculate loess predictions to smooth signal
  ## Args:
  #  x: x-axis data (predictor variable)
  #  y: y-axis data (response variable)
  # Returns:
  #  d: data.frame of x and p
  #     x: x-axis for predictions (includes NAs)
  #     p: loess predictions
  x[is.na(y)] <- NA
  lo_obj <- loess(y ~ x, span = 0.01)
  p <- predict(lo_obj, x)
  d <- data.frame(x, p)
  return(d)
}

d <- read.table(paste(wkdir, args$i, sep="/"),
                       sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS"))

p_LOS <- c()

for (chrom in c(paste("chr", 1:22,sep=""), "chrX")) {
  chrom_start <- d$start[d$chrom == chrom]
  chrom_LOS <- d$LOS[d$chrom  == chrom]
  # Skip smoothing if less than 100 valid LOS in chrom
  if (sum(!is.na(chrom_LOS)) <= 100) {
    p_LOS <- c(p_LOS, chrom_LOS)
  }
  else {
    # Smooth LOS data
    chrom_pred_df <- get_loess_prediction(chrom_start, chrom_LOS)
    p_LOS <- c(p_LOS, chrom_pred_df$p)
  }
}

loess_d <- cbind(d[c("chrom", "start", "end")], p_LOS)

outfile <- paste(substr(args$i, 0, nchar(args$i)-9), "_loess_smooth.bedGraph", sep="")

write.table(loess_d, paste(wkdir, outfile, sep="/"), quote=FALSE, sep="\t", 
            row.names=FALSE, col.names=FALSE)
