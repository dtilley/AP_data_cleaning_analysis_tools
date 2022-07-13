source("~/projects/code-share-112018/jgpplotter.R")
source("~/projects/code-share-112018/APSetReader.R")
library(plyr)

get.TriShapeCoef <- function(APDsAPVs) {
+     ## T_shape_coef: ~1 is triangulated, ~2.5 dome-like
+     tsc  <- (APDsAPVs$apd40-APDsAPVs$apd30)/(APDsAPVs$apd80-APDsAPVs$apd70)
+     return(tsc)
+ }

## Read Cell 1 Data
cell1.cntrl <- read.table("../Data/cell_1_122420_3_drew_drp/cntrl_ishihara/APDsAPVs_cntrl_ishi_122420_3_drew.txt", header=T)
cell1.cal.15 <- read.table("../Data/cell_1_122420_3_drew_drp/ical_-0.15/APDsAPVs_ical_decrease_122420_3_drew.txt", header=T)
cell1.cal.7 <- read.table("../Data/cell_1_122420_3_drew_drp/ical_0.7/APDsAPVs_ical_increase_122420_3_drew.txt", header=T)
cell1.kr.25 <- read.table("../Data/cell_1_122420_3_drew_drp/ikr_-0.25/APDsAPVs_ikr_decrease_122420_3_drew.txt", header=T)
cell1.kr.9 <- read.table("../Data/cell_1_122420_3_drew_drp/ikr_0.9/APDsAPVs_ikr_increase_122420_3_drew.txt", header=T)
cell1.ks.10 <- read.table("../Data/cell_1_122420_3_drew_drp/iks_10x/APDsAPVs_iks_10x_122420_3_drew.txt", header=T)
cell1.ks.4 <- read.table("../Data/cell_1_122420_3_drew_drp/iks_4x/APDsAPVs_iks_4x_122420_3_drew.txt", header=T)
cell1.to.9 <- read.table("../Data/cell_1_122420_3_drew_drp/ito_-0.9/APDsAPVs_ito_decrease_122420_3_drew.txt", header=T)
cell1.to.1.5 <- read.table("../Data/cell_1_122420_3_drew_drp/ito_1.5/APDsAPVs_ito_increase_122420_3_drew.txt", header=T)

cell1.APDsAPVs <- list(cell1.cntrl, cell1.cal.15, cell1.cal.7,
                       cell1.kr.25, cell1.kr.9, cell1.ks.10,
                       cell1.ks.4, cell1.to.9, cell1.to.1.5)

labels <- c("cntrl", "ICaL -0.15", "ICaL 0.7",
            "IKr -0.25", "IKr 0.9","IKs 10.0",
            "IKs 4.0", "Ito -0.9","Ito 1.5")
cols <- c("black", "blue", "blue",
          "red", "red", "magenta",
          "magenta", "darkorange", "darkorange")
ltys <- c(1, 4, 1,
          4, 1, 1,
          4, 4, 1)

cell1.TSC <- lapply(cell1.APDsAPVs, get.TriShapeCoef)
for (i in seq(length(cell1.TSC))) {
    filename <- paste(prefix,gsub(" ", "", labels[i]),suffix, sep="", collapse="")
    df <- as.data.frame(cell1.TSC[[i]])
    names(df) <- "TriSC"
    write.table(x=df, file=filename, quote=FALSE, sep="", row.names=FALSE)
}

## Find min and max for binning
apd90s <- list(cell1.cntrl$apd90, cell1.cal.15$apd90, cell1.cal.7$apd90,
               cell1.kr.25$apd90, cell1.kr.9$apd90, cell1.ks.10$apd90,
               cell1.ks.4$apd90, cell1.to.9$apd90, cell1.to.1.5$apd90)

apd90.extrema <- c(min(sapply(apd90s, min)), max(sapply(apd90s, max)))
apd90.bounds <- round_any(apd90.extrema, 20)

apd50s <- list(cell1.cntrl$apd50, cell1.cal.15$apd50, cell1.cal.7$apd50,
               cell1.kr.25$apd50, cell1.kr.9$apd50, cell1.ks.10$apd50,
               cell1.ks.4$apd50, cell1.to.9$apd50, cell1.to.1.5$apd50)

apd50.extrema <- c(min(sapply(apd50s, min)), max(sapply(apd50s, max)))
apd50.bounds <- round_any(apd50.extrema, 20)

## Plot histograms
breaks.apd90 <- seq(from=apd90.bounds[1], to=apd90.bounds[2], by=20)


relprob.apd90 <- list()
for (i in seq(length(apd90s))) {
    h <- hist(apd90s[[i]], breaks=breaks.apd90, plot=FALSE)
    relprob.apd90[[i]] <- c(h$counts/sum(h$counts),0)
}

labels <- c("cntrl", "ICaL -0.15", "ICaL 0.7",
            "IKr -0.25", "IKr 0.9","IKs 10.0",
            "IKs 4.0", "Ito -0.9","Ito 1.5")
cols <- c("black", "blue", "blue",
          "red", "red", "magenta",
          "magenta", "darkorange", "darkorange")
ltys <- c(1, 4, 1,
          4, 1, 1,
          4, 4, 1)
y.max <- round_any(max(sapply(relprob.apd90, max)), accuracy=0.05, f=ceiling)

jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
for (i in seq(length(relprob.apd90))) {
    lines(breaks.apd90,relprob.apd90[[i]], type="s",col=cols[i],lwd=1.5,lty=ltys[i])
}

jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[2]], type="s",col=cols[2],lwd=1.5,lty=ltys[2])
lines(breaks.apd90,relprob.apd90[[3]], type="s",col=cols[3],lwd=1.5,lty=ltys[3])
lines(breaks.apd90,relprob.apd90[[4]], type="s",col=cols[4],lwd=1.5,lty=ltys[4])
lines(breaks.apd90,relprob.apd90[[5]], type="s",col=cols[5],lwd=1.5,lty=ltys[5])
lines(breaks.apd90,relprob.apd90[[6]], type="s",col=cols[6],lwd=1.5,lty=ltys[6])
lines(breaks.apd90,relprob.apd90[[7]], type="s",col=cols[7],lwd=1.5,lty=ltys[7])
lines(breaks.apd90,relprob.apd90[[8]], type="s",col=cols[8],lwd=1.5,lty=ltys[8])
lines(breaks.apd90,relprob.apd90[[9]], type="s",col=cols[9],lwd=1.5,lty=ltys[9])

legend("topright", legend=labels, bty="n", text.col=cols,lty=ltys, col=cols)

## Plot ICaL shift
pdf(file="cell1_cntrl_icalblck.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[2]], type="s",col=cols[2],lwd=1.5,lty=ltys[2])
legend("topright", legend=labels[1:2], bty="n", text.col=cols[1:2],lty=ltys[1:2], col=cols[1:2])
dev.off()

# pdftoppm -png cell1_cntrl_icalblck.pdf cell1_cntrl_icalblck

pdf(file="cell1_cntrl_icalaug.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[3]], type="s",col=cols[3],lwd=1.5,lty=ltys[3])
legend("topright", legend=labels[c(1,3)], bty="n", text.col=cols[c(1,3)],lty=ltys[c(1,3)], col=cols[c(1,3)])
dev.off()

pdf(file="cell1_cntrl_ical.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[2]], type="s",col=cols[2],lwd=1.5,lty=ltys[2])
lines(breaks.apd90,relprob.apd90[[3]], type="s",col=cols[3],lwd=1.5,lty=ltys[3])
legend("topright", legend=labels[1:3], bty="n", text.col=cols[1:3],lty=ltys[1:3], col=cols[1:3])
dev.off()

## Plot IKr shift
pdf(file="cell1_cntrl_ikrblck.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[4]], type="s",col=cols[4],lwd=1.5,lty=ltys[4])
legend("topright", legend=labels[c(1,4)], bty="n", text.col=cols[c(1,4)],lty=ltys[c(1,4)], col=cols[c(1,4)])
dev.off()

pdf(file="cell1_cntrl_ikraug.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[5]], type="s",col=cols[5],lwd=1.5,lty=ltys[5])
legend("topright", legend=labels[c(1,5)], bty="n", text.col=cols[c(1,5)],lty=ltys[c(1,5)], col=cols[c(1,5)])
dev.off()

pdf(file="cell1_cntrl_ikr.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[4]], type="s",col=cols[4],lwd=1.5,lty=ltys[4])
lines(breaks.apd90,relprob.apd90[[5]], type="s",col=cols[5],lwd=1.5,lty=ltys[5])
legend("topright", legend=labels[c(1,4,5)], bty="n", text.col=cols[c(1,4,5)],lty=ltys[c(1,4,5)], col=cols[c(1,4,5)])
dev.off()

## Plot Iks shift
pdf(file="cell1_cntrl_iks10.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[6]], type="s",col=cols[6],lwd=1.5,lty=ltys[6])
legend("topright", legend=labels[c(1,6)], bty="n", text.col=cols[c(1,6)],lty=ltys[c(1,6)], col=cols[c(1,6)])
dev.off()

pdf(file="cell1_cntrl_iks4.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[7]], type="s",col=cols[7],lwd=1.5,lty=ltys[7])
legend("topright", legend=labels[c(1,7)], bty="n", text.col=cols[c(1,7)],lty=ltys[c(1,7)], col=cols[c(1,7)])
dev.off()

pdf(file="cell1_cntrl_iks.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[6]], type="s",col=cols[6],lwd=1.5,lty=ltys[6])
lines(breaks.apd90,relprob.apd90[[7]], type="s",col=cols[7],lwd=1.5,lty=ltys[7])
legend("topright", legend=labels[c(1,6,7)], bty="n", text.col=cols[c(1,6,7)],lty=ltys[c(1,6,7)], col=cols[c(1,6,7)])
dev.off()


## Plot Ito shift
pdf(file="cell1_cntrl_itoblck.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[8]], type="s",col=cols[8],lwd=1.5,lty=ltys[8])
legend("topright", legend=labels[c(1,8)], bty="n", text.col=cols[c(1,8)],lty=ltys[c(1,8)], col=cols[c(1,8)])
dev.off()

pdf(file="cell1_cntrl_itoaug.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[9]], type="s",col=cols[9],lwd=1.5,lty=ltys[9])
legend("topright", legend=labels[c(1,9)], bty="n", text.col=cols[c(1,9)],lty=ltys[c(1,9)], col=cols[c(1,9)])
dev.off()

pdf(file="cell1_cntrl_ito.pdf", width=7, height=7)
jgpplotter(breaks.apd90,relprob.apd90[[1]], type="n", xlab="APD90 (ms)",
           ylab="rel. probability", ylim=c(0,y.max))
lines(breaks.apd90,relprob.apd90[[1]], type="s",col=cols[1],lwd=1.5,lty=ltys[1])
lines(breaks.apd90,relprob.apd90[[8]], type="s",col=cols[8],lwd=1.5,lty=ltys[8])
lines(breaks.apd90,relprob.apd90[[9]], type="s",col=cols[9],lwd=1.5,lty=ltys[9])
legend("topright", legend=labels[c(1,8,9)], bty="n", text.col=cols[c(1,8,9)],lty=ltys[c(1,8,9)], col=cols[c(1,8,9)])
dev.off()

