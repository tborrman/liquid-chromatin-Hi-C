t1 <- read.table("timecourse1/cis_distances_1-10kb_timecourse1.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
t2 <- read.table("timecourse2/cis_distances_1-10kb_timecourse2.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
t3 <- read.table("timecourse3/cis_distances_1-10kb_timecourse3.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)

x <- seq(1.5, 9.5, 1) 
newx <- seq(1,10, .01)
timepoints <- rownames(t1)
mycolors <- c("dodgerblue", "orange", "forestgreen", "pink", "purple", "cyan", "red")

png("cis_distance_1-10kb.png", width=1500, height=2000, res=300)
# MN
plot(x,t1["MN",], pch=20, col="dodgerblue", cex=0.25, ylim=c(0,0.31),
     ylab= "% of total interactions", xlab= "cis interaction distance (kb)")
axis(side=1, at=1:10)
points(x, t2["MN",], pch=20, col="dodgerblue", cex=0.25)
points(x, t3["MN",], pch=20, col="dodgerblue", cex=0.25)
MN <- rbind(t1["MN",], t2["MN",], t3["MN",])
meanMN <- apply(MN, 2, mean)
loess.model <- loess(meanMN ~ x)
p <- predict(loess.model, newdata=newx)
lines(newx, p, col="dodgerblue")

# Timecourse
for (i in 2:length(timepoints)) {
  points(x, t1[timepoints[i],], pch=20, col=mycolors[i], cex=0.25)
  points(x, t2[timepoints[i],], pch=20, col=mycolors[i], cex=0.25)
  points(x, t3[timepoints[i],], pch=20, col=mycolors[i], cex=0.25)
  DN <- rbind(t1[timepoints[i],], t2[timepoints[i],], t3[timepoints[i],])
  meanDN <- apply(DN, 2, mean)
  loess.model <- loess(meanDN ~ x)
  p <- predict(loess.model, newdata=newx)
  lines(newx, p, col=mycolors[i])
}
legend("topright", timepoints, fill=mycolors, cex=0.75)
dev.off()

# semi-Log plot
png("cis_distance_1-10kb_semi-log10.png", width=1500, height=2000, res=300)
# MN
plot(x,log10(t1["MN",]), pch=20, col="dodgerblue", cex=0.25, ylim=c(-2.6, -0.2),
     ylab= "log10(% of total interactions)", xlab= "cis interaction distance (kb)")
axis(side=1, at=1:10)
points(x, log10(t2["MN",]), pch=20, col="dodgerblue", cex=0.25)
points(x, log10(t3["MN",]), pch=20, col="dodgerblue", cex=0.25)
MN <- rbind(log10(t1["MN",]), log10(t2["MN",]), log10(t3["MN",]))
meanMN <- apply(MN, 2, mean)
loess.model <- loess(meanMN ~ x)
p <- predict(loess.model, newdata=newx)
lines(newx, p, col="dodgerblue")

# Timecourse
for (i in 2:length(timepoints)) {
  points(x, log10(t1[timepoints[i],]), pch=20, col=mycolors[i], cex=0.25)
  points(x, log10(t2[timepoints[i],]), pch=20, col=mycolors[i], cex=0.25)
  points(x, log10(t3[timepoints[i],]), pch=20, col=mycolors[i], cex=0.25)
  DN <- rbind(log10(t1[timepoints[i],]), log10(t2[timepoints[i],]), log10(t3[timepoints[i],]))
  meanDN <- apply(DN, 2, mean)
  loess.model <- loess(meanDN ~ x)
  p <- predict(loess.model, newdata=newx)
  lines(newx, p, col=mycolors[i])
}
legend("topright", timepoints, fill=mycolors, cex=0.75)
dev.off()

# log-log plot
x <- log10(x*1000)
newx <- log10(seq(1000,10000, 10))
png("cis_distance_1-10kb_log10-log10.png", width=1500, height=2000, res=300)
# MN
plot(x,log10(t1["MN",]), pch=20, col="dodgerblue", cex=0.25, ylim=c(-2.6, -0.2),
     ylab= "log10(% of total interactions)", xlab= "log10(cis interaction distance)")
points(x, log10(t2["MN",]), pch=20, col="dodgerblue", cex=0.25)
points(x, log10(t3["MN",]), pch=20, col="dodgerblue", cex=0.25)
MN <- rbind(log10(t1["MN",]), log10(t2["MN",]), log10(t3["MN",]))
meanMN <- apply(MN, 2, mean)
loess.model <- loess(meanMN ~ x)
p <- predict(loess.model, newdata=newx)
lines(newx, p, col="dodgerblue")

# Timecourse
for (i in 2:length(timepoints)) {
  points(x, log10(t1[timepoints[i],]), pch=20, col=mycolors[i], cex=0.25)
  points(x, log10(t2[timepoints[i],]), pch=20, col=mycolors[i], cex=0.25)
  points(x, log10(t3[timepoints[i],]), pch=20, col=mycolors[i], cex=0.25)
  DN <- rbind(log10(t1[timepoints[i],]), log10(t2[timepoints[i],]), log10(t3[timepoints[i],]))
  meanDN <- apply(DN, 2, mean)
  loess.model <- loess(meanDN ~ x)
  p <- predict(loess.model, newdata=newx)
  lines(newx, p, col=mycolors[i])
}
legend("topright", timepoints, fill=mycolors, cex=0.75)
dev.off()

# Derivative plot
png("cis_distance_1-10kb_log10-log10_dydx.png", width=1500, height=2000, res=300)
# MN
MN <- rbind(log10(t1["MN",]), log10(t2["MN",]), log10(t3["MN",]))
meanMN <- apply(MN, 2, mean)
loess.model <- loess(meanMN ~ x)
p <- predict(loess.model, newdata=newx)
dY <- diff(p)/diff(newx)
midx <- rowMeans(embed(newx, 2))
plot(midx, dY, type='l', col="dodgerblue", xlab= "log10(cis interaction distance)",
     ylab="dy/dx", ylim=c(-2.15,0.15), xlim=c(3.17, 4))
# Timecourse
for (i in 2:length(timepoints)) {
  DN <- rbind(log10(t1[timepoints[i],]), log10(t2[timepoints[i],]), log10(t3[timepoints[i],]))
  meanDN <- apply(DN, 2, mean)
  loess.model <- loess(meanDN ~ x)
  p <- predict(loess.model, newdata=newx)
  dY <- diff(p)/diff(newx)
  midx <- rowMeans(embed(newx, 2))
  lines(midx, dY, type='l', col=mycolors[i])
}
legend("topright", timepoints, fill=mycolors, cex=0.75)
dev.off()

