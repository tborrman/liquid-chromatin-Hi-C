# Estimate half lives using loess
df <- read.table("LOS_500kb.txt", sep="\t", header=TRUE)

minutes <- c(5, 60, 120, 180, 240, 960)

LOS_half <- round(median(df$OVN_LOS, na.rm=TRUE)/2.0, 2)

# Ex.chr1 32Mb-32.5Mb
LOS <- as.numeric(df[65,4:9])
x <- seq(-50,1000, 0.5)
lo_obj <- loess(LOS ~ minutes)
p <- predict(lo_obj, x)
hl_idx <- which.min(abs(p-LOS_half))
half_life <- x[hl_idx]

png("plot_half_life/chr1_32Mb_hl_loess.png", height=1800, width=2500, res=300)
plot(minutes, LOS, pch=20, col = "darkgreen", xlab="Minutes of DpnII Digestion",
     ylab= "Loss of Structure",main= "Chr1:32Mb-32.5Mb", ylim=c(-0.03,0.7))
lines(x,p, col="chartreuse4")
h_linex <- seq(-50,half_life, 0.5)
v_linex <- seq(-0.3, LOS_half, 0.001)
lines(h_linex,rep(LOS_half, length(h_linex)), lty=2, col="blue")
lines(rep(half_life, length(v_linex)), v_linex, lty=2, col="blue")
text(700,0.2, bquote("t"[1/2] ~ "=" ~ .(half_life) ~ "min"), cex=1.5)
dev.off()

# Ex.chr1 49.5Mb-50Mb
LOS <- as.numeric(df[100,4:9])
x <- seq(-50,1000, 0.5)
lo_obj <- loess(LOS ~ minutes, span= 1.0)
p <- predict(lo_obj, x)
hl_idx <- which.min(abs(p-LOS_half))
half_life <- x[hl_idx]

png("plot_half_life/chr1_49.5Mb_hl_loess_span_1.0.png", height=1800, width=2500, res=300)
plot(minutes, LOS, pch=20, col = "darkgreen", xlab="Minutes of DpnII Digestion",
     ylab= "Loss of Structure",main= "Chr1:49.5Mb-50Mb", ylim=c(-0.03,0.7))
lines(x,p, col="chartreuse4")
h_linex <- seq(-50,half_life, 0.5)
v_linex <- seq(-0.3, LOS_half, 0.001)
lines(h_linex,rep(LOS_half, length(h_linex)), lty=2, col="blue")
lines(rep(half_life, length(v_linex)), v_linex, lty=2, col="blue")
text(700,0.2, bquote("t"[1/2] ~ "=" ~ .(half_life) ~ "min"), cex=1.5)
dev.off()

getHalfLife <- function(r, minutes, LOS_half) {
  LOS <- as.numeric(r[4:9])
  x <- seq(-50,1000, 0.5)
  lo_obj <- loess(LOS ~ minutes, span=1.0)
  p <- predict(lo_obj, x)
  hl_idx <- which.min(abs(p-LOS_half))
  hl <- x[hl_idx]
  return(hl)
}

# Create half-life vector
half_lives <- c()
for (row in 1:nrow(df)) {

  if (sum(is.na(df[row,])) == 0) {
    half_life <- getHalfLife(df[row,], minutes, LOS_half)
    half_lives <- c(half_lives, half_life)
  }
  else {
    half_lives <- c(half_lives, NA)
  }  
}

df <- cbind(df, half_lives)


# Chromosome tracks
for (chrom in c(paste("chr", 1:22, sep=""), "chrX", "chrY")) {

  png(paste("plot_half_life/plot_half_life_span_1.0_", chrom, ".png",sep=""), width=8000, height=2000, res=400)
  par(mar=c(5,6,4,8) + 0.1, lwd=2)

  mid <- df$start + (abs(df$end - df$start)/2.0)
  chrom_mid <- mid[df$chrom == chrom]
  chrom_hl <- df$half_lives[df$chrom == chrom]
  plot(chrom_mid, chrom_hl, type='l', col="chartreuse4", xlab= chrom,
       ylab=bquote("t"[1/2] ~ "(minutes)"), cex.lab = 2, ylim=c(0,500))
  dev.off()
}

