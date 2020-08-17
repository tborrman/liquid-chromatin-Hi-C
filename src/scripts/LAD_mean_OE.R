options(scipen=100)
# Lamin analysis
OE <- read.table("C:/cygwin64/home/Tyler/Research/digest/lamina/mmc2/DataS2_Clone.5-5.1N.OE.txt", header=TRUE, sep="\t")
OE$start <- OE$start - 1
OE_data <- OE[4:ncol(OE)]
means <- apply(OE_data, 1, mean)
OE_means <- cbind(OE[1:3], means)
write.table(OE_means,"Lamina_mean_OE_Clone.5-5.1N.OE.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
