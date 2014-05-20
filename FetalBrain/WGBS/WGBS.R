# read coordinates for all CpGs: 
CG <- read.delim("/projects/epigenomics/kraghavan/CpGData/Full_O2Files/CG.BED", as.is = T, head = F)
# integrate coverage and fractional methylation calls: ignore coverage < 5
setwd("~/FetalBrain/WGBS/")
cortex02m <- read.delim("cortex02.fractional.bedGraph", head = F)
rownames(cortex02m) <- paste0(cortex02m$V1, ":", cortex02m$V2)
ge02m <- read.delim("GE02.fractional.bedGraph", head = F)
rownames(ge02m) <- paste0(ge02m$V1, ":", ge02m$V2)
cortex04m <- read.delim("cortex04.fractional.bedGraph", head = F)
rownames(cortex04m) <- paste0(cortex04m$V1, ":", cortex04m$V2)
ge04m <- read.delim("GE04.fractional.bedGraph", head = F)
rownames(ge04m) <- paste0(ge04m$V1, ":", ge04m$V2)
cortex02c <- read.delim("cortex02.coverage.bedGraph", head = F)
rownames(cortex02c) <- paste0(cortex02c$V1, ":", cortex02c$V2)
ge02c <- read.delim("GE02.coverage.bedGraph", head = F)
rownames(ge02c) <- paste0(ge02c$V1, ":", ge02c$V2)
cortex04c <- read.delim("cortex04.coverage.bedGraph", head = F)
rownames(cortex04c) <- paste0(cortex04c$V1, ":", cortex04c$V2)
ge04c <- read.delim("GE04.coverage.bedGraph", head = F)
rownames(ge04c) <- paste0(ge04c$V1, ":", ge04c$V2)
WGBS <- data.frame(chr = rep(CG$V1, times = 2), start = c(CG$V2, CG$V2+1), end = c(CG$V3-1, CG$V3))
WGBS$chr <- factor(WGBS$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))
rownames(WGBS) <- paste0(rep(CG$V1, times = 2), ":", c(CG$V2, CG$V2+1))
WGBS[,c("cortex02m", "ge02m", "cortex04m", "ge04m")] <- -10000
WGBS[,c("cortex02c", "ge02c", "cortex04c", "ge04c")] <- 0
WGBS[rownames(cortex02m),]$cortex02m <- cortex02m$V4
WGBS[rownames(cortex02c),]$cortex02c <- abs(cortex02c$V4)
WGBS[rownames(ge02m),]$ge02m <- ge02m$V4
WGBS[rownames(ge02c),]$ge02c <- abs(ge02c$V4)
WGBS[rownames(cortex04m),]$cortex04m <- cortex04m$V4
WGBS[rownames(cortex04c),]$cortex04c <- abs(cortex04c$V4)
WGBS[rownames(ge04m),]$ge04m <- ge04m$V4
WGBS[rownames(ge04c),]$ge04c <- abs(ge04c$V4)
CpGs <- paste0(rep(CG$V1, times = 2), ":", c(CG$V2, CG$V2+1))
WGBS <- WGBS[CpGs,]
WGBS <- WGBS[order(WGBS$chr, WGBS$start),]
WGBS[WGBS$cortex02c < 5,]$cortex02m <- -10000
WGBS[WGBS$ge02c < 5,]$ge02m <- -10000
WGBS[WGBS$cortex04c < 5,]$cortex04m <- -10000
WGBS[WGBS$ge04c < 5,]$ge04m <- -10000
save(WGBS, file = "WGBS_56M.Rdata")

# symmetry or asymmetry on different strands
setwd("~/FetalBrain/WGBS/")
load("WGBS_56M.Rdata")
WGBS_pos <- WGBS[seq(1, nrow(WGBS)-1, by = 2), ]
WGBS_neg <- WGBS[seq(2, nrow(WGBS), by = 2), ]
rm(WGBS)
pdf("strand_symmetry.pdf")
smoothScatter(WGBS_pos$cortex02m, WGBS_neg$cortex02m, main = "Methylation on opposite strands - cortex02", xlab = "pos", ylab = "neg")
smoothScatter(WGBS_pos$ge02m, WGBS_neg$ge02m, main = "Methylation on opposite strands - ge02", xlab = "pos", ylab = "neg")
smoothScatter(WGBS_pos$cortex04m, WGBS_neg$cortex04m, main = "Methylation on opposite strands - cortex04", xlab = "pos", ylab = "neg")
smoothScatter(WGBS_pos$ge04m, WGBS_neg$ge04m, main = "Methylation on opposite strands - ge04", xlab = "pos", ylab = "neg")
plot(c(-10, 10), c(0, 0.5), type = "n", main = "Methylation asymmetry: pos-neg", xlab = "pos-neg", ylab = "density")
lines(density(na.omit(WGBS_pos$cortex02m-WGBS_neg$cortex02m), adjust = 15), col = 1, lwd = 2)
lines(density(na.omit(WGBS_pos$ge02m-WGBS_neg$ge02m), adjust = 15), col = 2, lwd = 2)
lines(density(na.omit(WGBS_pos$cortex04m-WGBS_neg$cortex04m), adjust = 15), col = 3, lwd = 2)
lines(density(na.omit(WGBS_pos$ge04m-WGBS_neg$ge04m), adjust = 15), col = 4, lwd = 2)
legend("topright", c("cortex02", "GE02", "cortex04", "GE04"), col = c(1:4), lwd = 5)
dev.off()

###############################################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R
# take weighted average on coverage to convert 56M to 28M
setwd("~/FetalBrain/WGBS/")
CG <- read.delim("/projects/epigenomics/kraghavan/CpGData/Full_O2Files/CG.BED", as.is = T, head = F)
load("WGBS_56M.Rdata")
WGBS_pos <- WGBS[seq(1, nrow(WGBS)-1, by = 2), ]
WGBS_neg <- WGBS[seq(2, nrow(WGBS), by = 2), ]
WGBS_28M <- data.frame(chr = CG$V1, start = CG$V2, end = CG$V3)
rm(CG, WGBS)
rownames(WGBS_28M) <- paste0(CG$V1, ":", CG$V2, "-", CG$V3)
rownames(WGBS_pos) <- rownames(WGBS_28M)
rownames(WGBS_neg) <- rownames(WGBS_28M)

rows <- intersect(rownames(WGBS_pos[WGBS_pos$cortex02m > -1,]), rownames(WGBS_neg[WGBS_neg$cortex02m > -1,]))
non <- rownames(WGBS_28M)[!rownames(WGBS_28M) %in% rows]
WGBS_28M[rows,]$cortex02 <- (WGBS_pos[rows,]$cortex02m * WGBS_pos[rows,]$cortex02c + WGBS_neg[rows,]$cortex02m * WGBS_neg[rows,]$cortex02c) / (WGBS_pos[rows,]$cortex02c + WGBS_neg[rows,]$cortex02c) # both strand have enough coverage
WGBS_28M[non,]$cortex02 <- WGBS_pos[non,]$cortex02m + WGBS_neg[non,]$cortex02m + 10000 # at least one strand have not enough coverage
WGBS_28M$cortex02 <- WGBS_28M$cortex02/10       # scale to [0,1]
WGBS_28M[WGBS_28M$cortex02 < 0,]$cortex02 <- NA # not enough coverage on both strand

rows <- intersect(rownames(WGBS_pos[WGBS_pos$ge02m > -1,]), rownames(WGBS_neg[WGBS_neg$ge02m > -1,]))
non <- rownames(WGBS_28M)[!rownames(WGBS_28M) %in% rows]
WGBS_28M[rows,]$ge02 <- (WGBS_pos[rows,]$ge02m * WGBS_pos[rows,]$ge02c + WGBS_neg[rows,]$ge02m * WGBS_neg[rows,]$ge02c) / (WGBS_pos[rows,]$ge02c + WGBS_neg[rows,]$ge02c)
WGBS_28M[non,]$ge02 <- WGBS_pos[non,]$ge02m + WGBS_neg[non,]$ge02m + 10000 
WGBS_28M$ge02 <- WGBS_28M$ge02/10
WGBS_28M[WGBS_28M$ge02 < 0,]$ge02 <- NA

rows <- intersect(rownames(WGBS_pos[WGBS_pos$cortex04m > -1,]), rownames(WGBS_neg[WGBS_neg$cortex04m > -1,]))
non <- rownames(WGBS_28M)[!rownames(WGBS_28M) %in% rows]
WGBS_28M[rows,]$cortex04 <- (WGBS_pos[rows,]$cortex04m * WGBS_pos[rows,]$cortex04c + WGBS_neg[rows,]$cortex04m * WGBS_neg[rows,]$cortex04c) / (WGBS_pos[rows,]$cortex04c + WGBS_neg[rows,]$cortex04c)
WGBS_28M[non,]$cortex04 <- WGBS_pos[non,]$cortex04m + WGBS_neg[non,]$cortex04m + 10000 
WGBS_28M$cortex04 <- WGBS_28M$cortex04/10
WGBS_28M[WGBS_28M$cortex04 < 0,]$cortex04 <- NA

rows <- intersect(rownames(WGBS_pos[WGBS_pos$ge04m > -1,]), rownames(WGBS_neg[WGBS_neg$ge04m > -1,]))
non <- rownames(WGBS_28M)[!rownames(WGBS_28M) %in% rows]
WGBS_28M[rows,]$ge04 <- (WGBS_pos[rows,]$ge04m * WGBS_pos[rows,]$ge04c + WGBS_neg[rows,]$ge04m * WGBS_neg[rows,]$ge04c) / (WGBS_pos[rows,]$ge04c + WGBS_neg[rows,]$ge04c)
WGBS_28M[non,]$ge04 <- WGBS_pos[non,]$ge04m + WGBS_neg[non,]$ge04m + 10000 
WGBS_28M$ge04<- WGBS_28M$ge04/10
WGBS_28M[WGBS_28M$ge04 < 0,]$ge04 <- NA

WGBS_28M <- WGBS_28M[WGBS_28M$chr %in% c(paste0("chr", 1:22), "chrX"),]
WGBS_28M$chr <- factor(WGBS_28M$chr, levels = c(paste0("chr", 1:22), "chrX"))
summary(WGBS_28M)
save(WGBS_28M, file = "WGBS_28M.Rdata")

###############################################################################################################
# global patterns
setwd("~/FetalBrain/WGBS/")
# load("WGBS_28M.Rdata")
load("~/FetalBrain/MeDIPMRE/MethylCRF.Rdata")

pdf("ECDF_WGBS_MethylCRF.pdf")
plot(c(0, 1), c(0, 1), type = "n", xlab = "Methylation level", ylab = "Ecdf", main = "ECDF of WGBS and MethylCRF")
lines(ecdf(MethylCRF$cortex02), lwd = 2, col = 1)
lines(ecdf(MethylCRF$ge02), lwd = 2, col = 2)
lines(ecdf(na.omit(WGBS_28M$cortex02)), lwd = 2, col = 3, lty = 2)
lines(ecdf(na.omit(WGBS_28M$ge02)), lwd = 2, col = 4, lty = 2)
lines(ecdf(na.omit(WGBS_28M$cortex04)), lwd = 2, col = 5, lty = 2)
lines(ecdf(na.omit(WGBS_28M$ge04)), lwd = 2, col = 6, lty = 2)
legend("topleft", c("CRFcortex02", "CRFge02", "WGBScortex02", "WGBSge02", "WGBScortex04", "WGBSge04"), col = 1:6, lwd = 5)
dev.off()
pdf("ECDF_WGBS_MethylCRF_diff.pdf")
plot(c(-1, 1), c(0, 1), type = "n", xlab = "Difference in methylation level", ylab = "Ecdf", main = "ECDF of cortex-GE WGBS vs MethylCRF")
lines(ecdf(MethylCRF$cortex02 - MethylCRF$ge02), lwd = 2, col = 1)
lines(ecdf(na.omit(WGBS_28M$cortex02 - WGBS_28M$ge02)), lwd = 2, col = 2, lty = 2)
lines(ecdf(na.omit(WGBS_28M$cortex04 - WGBS_28M$ge04)), lwd = 2, col = 3, lty = 2)
legend("bottomright", c("CRF", "WGBS-HuFNSC02", "WGBS-HuFNSC04"), col = 1:3, lwd = 5)
dev.off()

## boxplot of methylation level per chr
# install.packages("ggplot2", lib = "/home/lli/R")
require(ggplot2, lib.loc = "/home/lli/R")
require(plyr, lib.loc = "/home/lli/R")
WGBS_28M$method <- "WGBS" 
MethylCRF$method  <- "MethylCRF"

cortex02_WGBS_CRF <- rbind(WGBS_28M[,c("chr", "start", "end", "cortex02", "method")], MethylCRF[,c("chr", "start", "end", "cortex02", "method")])
cortex02_WGBS_CRF$method <- factor(cortex02_WGBS_CRF$method)
cortex02_WGBS_CRF <- droplevels(cortex02_WGBS_CRF[!(cortex02_WGBS_CRF$chr %in% c("chrM", "chrY")), ])
cortex02_quantile <- ddply(cortex02_WGBS_CRF, .(chr, method), summarize, "lower" = quantile(cortex02, 0.25, na.rm = T), "middle" = median(cortex02, na.rm = T), "upper" = quantile(cortex02, 0.75, na.rm = T), 
                           "ymin" = median(cortex02, na.rm = T) - 1.5*(quantile(cortex02, 0.75, na.rm = T) - quantile(cortex02, 0.25, na.rm = T)), 
                           "ymax" = median(cortex02, na.rm = T) + 1.5*(quantile(cortex02, 0.75, na.rm = T) - quantile(cortex02, 0.25, na.rm = T)))
cortex02_quantile$ymax <- 1
cortex02_quantile$chr <- factor(gsub("chr", "", as.character(cortex02_quantile$chr)), levels = c(as.character(1:22), "X"))
cortex02_boxplot <- ggplot(cortex02_quantile, aes(x = chr, lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax, fill = chr)) + 
  geom_boxplot(stat = "identity") +   
  facet_wrap(~ method) + 
  xlab("chr") + 
  ylab("methylation level") + 
  ggtitle("WGBS vs MethylCRF - cortex02") + 
  guides(fill = F) + 
  theme_bw()
ggsave(cortex02_boxplot, file = "cortex02_WGBS-CRF_boxplot.pdf", width = 10, height = 5)

ge02_WGBS_CRF <- rbind(WGBS_28M[,c("chr", "start", "end", "ge02", "method")], MethylCRF[,c("chr", "start", "end", "ge02", "method")])
ge02_WGBS_CRF$method <- factor(ge02_WGBS_CRF$method)
ge02_WGBS_CRF <- droplevels(ge02_WGBS_CRF[!(ge02_WGBS_CRF$chr %in% c("chrM", "chrY")), ])
ge02_quantile <- ddply(ge02_WGBS_CRF, .(chr, method), summarize, "lower" = quantile(ge02, 0.25, na.rm = T), "middle" = median(ge02, na.rm = T), "upper" = quantile(ge02, 0.75, na.rm = T), 
                       "ymin" = median(ge02, na.rm = T) - 1.5*(quantile(ge02, 0.75, na.rm = T) - quantile(ge02, 0.25, na.rm = T)), 
                       "ymax" = median(ge02, na.rm = T) + 1.5*(quantile(ge02, 0.75, na.rm = T) - quantile(ge02, 0.25, na.rm = T)))
ge02_quantile$ymax <- 1
ge02_quantile$chr <- factor(gsub("chr", "", as.character(ge02_quantile$chr)), levels = c(as.character(1:22), "X"))
ge02_boxplot <- ggplot(ge02_quantile, aes(x = chr, lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax, fill = chr)) + 
  geom_boxplot(stat = "identity") +   
  facet_wrap(~ method) + 
  xlab("chr") + 
  ylab("methylation level") + 
  ggtitle("WGBS vs MethylCRF - GE02") + 
  guides(fill = F) + 
  theme_bw()
ggsave(ge02_boxplot, file = "ge02_WGBS-CRF_boxplot.pdf", width = 10, height = 5)

## unmathylated: 0-0.3, intermediate: 0.3-0.7, methylated: 0.7-1. compare No. of CpGs in each category - need to normalize against total No. of CpGs


# bin into 200bp bins
require(plyr, lib.loc = "/home/lli/R")
chrlen <- read.csv("~/hg19/chrlen_hg19.csv",row.names=1)
bin200 <- read.csv("~/FetalBrain/MeDIPMRE/bin200.csv")
chrlen$startbin <- chrlen$X.1%/%200 + 2       # 1st bin of the chr. 
chrlen$startbin[1] <- 1

MethylCRF <- MethylCRF[!(MethylCRF$chr %in% c("chrM", "chrY")), ]
MethylCRF$startbin <- chrlen[MethylCRF$chr,]$startbin  #first bin of the chr
MethylCRF$bin <- MethylCRF$startbin + MethylCRF$start%/%200    
MethylCRF$count <- 1

crfave200 <- ddply(MethylCRF, ~bin, "count" = sum(count), "brain01" = mean(brain01), "brain02" = mean(brain02), "cortex01" = mean(cortex01), "cortex02" = mean(cortex02), "ge01" = mean(ge01), "ge02" = mean(ge02))

crfave200 <- aggregate(cbind(count, brain01, brain02, cortex01, cortex02, ge01, ge02) ~ bin, data = MethylCRF, sum)
crfave200$brain01 <- crfave200$brain01/crfave200$count
crfave200$brain02 <- crfave200$brain02/crfave200$count
crfave200$cortex01 <- crfave200$cortex01/crfave200$count
crfave200$cortex02 <- crfave200$cortex02/crfave200$count
crfave200$ge01 <- crfave200$ge01/crfave200$count
crfave200$ge02 <- crfave200$ge02/crfave200$count

crfave200$chr <- bin200[crfave200$bin,]$chr
crfave200$start <- bin200[crfave200$bin,]$start
crfave200$end <- bin200[crfave200$bin,]$end
save(crfave200, file = "crfave200.Rdata")


cortex02 <- cortex02[!(cortex02$chr %in% c("chrM", "chrY")), ]
cortex02$startbin <- chrlen[cortex02$chr,]$startbin  #first bin of the chr
cortex02$bin <- cortex02$startbin + cortex02$start%/%200    
cortex02$count <- 1

cortex02ave200 <- ddply(cortex02, ~bin, "count" = sum(count), "cortex02" = mean(cortex02))

cortex02ave200 <- aggregate(cbind(count, cortex02) ~ bin, data = cortex02, sum)
cortex02ave200$cortex02 <- cortex02ave200$cortex02/cortex02ave200$count

cortex02ave200$chr <- bin200[cortex02ave200$bin,]$chr
cortex02ave200$start <- bin200[cortex02ave200$bin,]$start
cortex02ave200$end <- bin200[cortex02ave200$bin,]$end
save(cortex02ave200, file = "cortex02ave200.Rdata")


ge02 <- ge02[!(ge02$chr %in% c("chrM", "chrY")), ]
ge02$startbin <- chrlen[ge02$chr,]$startbin  #first bin of the chr
ge02$bin <- ge02$startbin + ge02$start%/%200    
ge02$count <- 1

ge02ave200 <- ddply(ge02, ~bin, "count" = sum(count), "ge02" = mean(ge02))

ge02ave200 <- aggregate(cbind(count, ge02) ~ bin, data = ge02, sum)
ge02ave200$ge02 <- ge02ave200$ge02/ge02ave200$count

ge02ave200$chr <- bin200[ge02ave200$bin,]$chr
ge02ave200$start <- bin200[ge02ave200$bin,]$start
ge02ave200$end <- bin200[ge02ave200$bin,]$end
save(ge02ave200, file = "ge02ave200.Rdata")

rownames(cortex02ave200) <- cortex02ave200$bin
rownames(ge02ave200) <- ge02ave200$bin
commonBin <- intersect(cortex02ave200$bin, ge02ave200$bin)
WGBS.cortex02_ge02.ave200 <- cbind(cortex02ave200[commonBin,], ge02ave200[commonBin,]$count, ge02ave200[commonBin,]$ge02)
colnames(WGBS.cortex02_ge02.ave200) <- c("bin", "cortex02.count", "cortex02.m", "chr", "start", "end", "ge02.count", "ge02.m")
WGBS.cortex02_ge02.ave200$cortex02_ge02 <- WGBS.cortex02_ge02.ave200$cortex02.m - WGBS.cortex02_ge02.ave200$ge02.m
WGBS.cortex02_ge02.ave200 <- na.omit(WGBS.cortex02_ge02.ave200)
save(WGBS.cortex02_ge02.ave200, file = "WGBS.cortex02_ge02.ave200.Rdata") 

# global correlation
rownames(crfave200) <- crfave200$bin
rownames(WGBS.cortex02_ge02.ave200) <- WGBS.cortex02_ge02.ave200$bin
commonBin <- intersect(crfave200$bin, WGBS.cortex02_ge02.ave200$bin)
smoothScatter(crfave200[commonBin,]$cortex02, WGBS.cortex02_ge02.ave200[commonBin,]$cortex02.m)
smoothScatter(crfave200[commonBin,]$ge02, WGBS.cortex02_ge02.ave200[commonBin,]$ge02.m)
smoothScatter(crfave200[commonBin,]$cortex02_ge02, WGBS.cortex02_ge02.ave200[commonBin,]$cortex02_ge02)

# identify DMRs from 200bp bins 
setwd("~/FetalBrain/MeDIPMRE/")
load("crfave200.Rdata")
setwd("~/FetalBrain/WGBS/")
load("WGBS.cortex02_ge02.ave200.Rdata")

crfave200$brain01_brain02 <- crfave200$brain01 - crfave200$brain02
crfave200$cortex01_cortex02 <- crfave200$cortex01 - crfave200$cortex02
crfave200$ge01_ge02 <- crfave200$ge01 - crfave200$ge02
crfave200$cortex01_ge01 <- crfave200$cortex01 - crfave200$ge01
crfave200$cortex02_ge02 <- crfave200$cortex02 - crfave200$ge02

x <- seq(-1, 1, length = 10000)
plot(c(-0.5, 0.5), c(0, 25), type = "n", xlab = "Difference in methylation", ylab = "Density")
lines(x, dnorm(x, mean = mean(crfave200$brain01_brain02), sd = sd(crfave200$brain01_brain02)), lwd = 1, col = 1, lty = 2)
lines(density(crfave200$brain01_brain02, adjust = 20), lwd = 1, col = 2)
lines(density(crfave200$cortex01_cortex02, adjust = 20), lwd = 2, col = 3)
lines(density(crfave200$ge01_ge02, adjust = 20), lwd = 2, col = 4)
lines(density(crfave200$cortex01_ge01, adjust = 20), lwd = 2, col = 5)
lines(density(crfave200$cortex02_ge02, adjust = 20), lwd = 2, col = 6)
lines(density(WGBS.cortex02_ge02.ave200$cortex02_ge02, adjust = 20), lwd = 2, col = 7)
legend("topright", c("Normal", "brain01_brain02", "cortex01_cortex02", "ge01_ge02", "cortex01_ge01", "cortex02_ge02", "WGBS_cortex02_ge02"), col = 1:7, lty = c(2, rep(1, times=6)), lwd = 2)
plot(c(-1, 1), c(0, 1), type = "n", xlab = "Difference in methylation", ylab = "Ecdf")
lines(x, pnorm(x, mean = mean(crfave200$brain01_brain02), sd = sd(crfave200$brain01_brain02)), lwd = 1, col = 1, lty = 2)
lines(ecdf(crfave200$brain01_brain02), lwd = 2, col = 2)
lines(ecdf(crfave200$cortex01_cortex02), lwd = 2, col = 3)
lines(ecdf(crfave200$ge01_ge02), lwd = 2, col = 4)
lines(ecdf(crfave200$cortex01_ge01), lwd = 2, col = 5)
lines(ecdf(crfave200$cortex02_ge02), lwd = 2, col = 6)
lines(ecdf(WGBS.cortex02_ge02.ave200$cortex02_ge02), lwd = 2, col = 7)
legend("bottomright", c("Normal", "brain01_brain02", "cortex01_cortex02", "ge01_ge02", "cortex01_ge01", "cortex02_ge02", "WGBS_cortex02_ge02"), 
       cex = 0.5, col = 1:7, lty = c(2, rep(1, times=6)), lwd = 2)

# use 0.5% and 99.5% quantile as cutoff for DMRs 
(BScutneg <- quantile(WGBS.cortex02_ge02.ave200$cortex02_ge02, 0.005))
(BScutpos <- quantile(WGBS.cortex02_ge02.ave200$cortex02_ge02, 0.995))
(CRFcutneg <- quantile(crfave200$cortex02_ge02, 0.005))
(CRFcutpos <- quantile(crfave200$cortex02_ge02, 0.995))
BShypo <- WGBS.cortex02_ge02.ave200[WGBS.cortex02_ge02.ave200$cortex02_ge02 < BScutneg,]
BShyper <- WGBS.cortex02_ge02.ave200[WGBS.cortex02_ge02.ave200$cortex02_ge02 > BScutpos,]
CRFhypo <- crfave200[crfave200$cortex02_ge02 < CRFcutneg,]
CRFhyper <- crfave200[crfave200$cortex02_ge02 > CRFcutpos,]


