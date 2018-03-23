#!/usr/bin/env Rscript
library(argparse)
parser <- ArgumentParser(description="Create half-life bed file using exponential decay to fit LOS data")
parser$add_argument("-i", help="table of LOS values for timecourse (ex. LOS_500kb.txt)", type="character", 
                    required=TRUE)
parser$add_argument("-f", help="plot figures", type="logical", default=FALSE)
args <- parser$parse_args()

# Estimate half lives from LOS using exponential decay 
options(scipen=100)
df <- read.table(args$i, sep="\t", header=TRUE)
minutes <- c(5, 60, 120, 180, 240, 960)

exp_decay <- function(minutes, a, b, c) {
  # Sigmoid function with adjustable parameters
  return(a - (b*exp(-c*minutes)))
}

half_life_equation <- function(LOS_half, a, b, c) {
  # Get half-life using exponential decay function 
  # solved for minutes at half-Loss of structure
  hl <- -(log((a-LOS_half)/b, base=exp(1)) / c)
  return(hl)
}

get_half_life <- function(LOS, minutes, LOS_half) {
  # Fit LOS data with exponential decay curve and return 
  # half-life
  fitM <- try(nls(LOS ~ exp_decay(minutes, a, b, c), start = list(a=0.75, b=0.5, c=0.01), 
              trace=TRUE, control = nls.control(warnOnly = TRUE)))
  if(class(fitM) == "try-error") {
    return(NA)
  }
  else{
    params <- coef(fitM)
    half_life <- half_life_equation(LOS_half, params["a"], params["b"], params["c"])
    return(half_life)
  }
}


get_LOS_half <- function(row){
  LOS_start <- row$min5_LOS
  LOS_end <- row$OVN_LOS
  LOS_half <- LOS_start + ((LOS_end - LOS_start)/2)
  return(LOS_half)
}
if (args$f) {
  # Ex.chr1 49.5Mb-50Mb
  LOS <- as.numeric(df[100,4:9])
  # LOS_start <- LOS[1]
  # LOS_end <- LOS[length(LOS)]
  # LOS_half <- LOS_start + ((LOS_end - LOS_start)/2)
  LOS_half <- get_LOS_half(df[100,])
  
  # Exponential decay fit
  fitModel <- nls(LOS ~ exp_decay(minutes, a, b, c), start = list(a=0.75, b=0.5, c=0.01), trace=TRUE)
  params <- coef(fitModel)
  x <- seq(-50,1000, 0.5)
  p <- predict(fitModel, list(minutes = x))
  #hl_idx <- which.min(abs(p-LOS_half))
  #half_life <- x[hl_idx]
  half_life <- half_life_equation(LOS_half, params["a"], params["b"], params["c"])
  
  
  
  png("chr1_49.5Mb_exponential.png", height=1800, width=2500, res=300)
  plot(minutes, LOS, pch=20, col = "darkgreen", xlab="Minutes of DpnII Digestion",
       ylab= "Loss of Structure",main= "Chr1:49.5Mb-50Mb", ylim=c(-0.3,0.7))
  lines(x,p, col="chartreuse4")
  h_linex <- seq(-50,half_life, 0.5)
  v_linex <- seq(-0.5, LOS_half, 0.001)
  lines(h_linex,rep(LOS_half, length(h_linex)), lty=2, col="blue")
  lines(rep(half_life, length(v_linex)), v_linex, lty=2, col="blue")
  text(700,0.2, bquote("t"[1/2] ~ "=" ~ .(round(half_life,2)) ~ "min"), cex=1.5)
  dev.off()
}

# Create half-life vector
half_lives <- c()
for (row in 1:nrow(df)) {
  if (sum(is.na(df[row,])) == 0) {
    LOS_half <- get_LOS_half(df[row,])
    half_life <- get_half_life(as.numeric(df[row,4:9]), minutes, LOS_half)
    half_lives <- c(half_lives, half_life)
  }
  else {
    half_lives <- c(half_lives, NA)
  }  
}

dfHL <- cbind(df[,1:3], half_lives)

# Write half-life to file
write.table(dfHL, "half-life_exponential_500kb.bed", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

if (args$f) {
# Chromosome tracks
  for (chrom in c(paste("chr", 1:22, sep=""), "chrX")) {
  
    png(paste("plot_half_life_exponential_", chrom, ".png",sep=""), width=8000, height=2000, res=400)
    par(mar=c(5,6,4,8) + 0.1, lwd=2)
  
    mid <- dfHL$start + (abs(dfHL$end - dfHL$start)/2.0)
    chrom_mid <- mid[dfHL$chrom == chrom]
    chrom_hl <- dfHL$half_lives[dfHL$chrom == chrom]
    plot(chrom_mid, chrom_hl, type='l', col="chartreuse4", xlab= chrom,
         ylab=bquote("t"[1/2] ~ "(minutes)"), cex.lab = 2, ylim=c(0,300),
         main= "Loss of Structure", cex.main=2)
    dev.off()
  }
}
