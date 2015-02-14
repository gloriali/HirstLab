# join MeDIP only fractional calls from different samples and chrs into one file
# /gsc/software/linux-x86_64-centos5/R-2.15.2/bin/R
setwd("/home/lli/FetalBrain/MeDIP/CG_25_around_chr/")
names <- c("HS2788.MeDIP.Brain01.q5.F1028.SET_174", 
           "HS2790.MeDIP.Brain02.q5.F1028.SET_174", 
           "HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174", 
           "HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174", 
           "HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157", 
           "HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166")
chrs <- c(paste0("chr", as.character(1:22)), "chrX")
brain01 <- read.delim(paste0("./", names[1], "/chr1/chr1.", names[1], ".dip"), head = F, as.is = T)[,1:2]
brain02 <- read.delim(paste0("./", names[2], "/chr1/chr1.", names[2], ".dip"), head = F, as.is = T)[,1:2]
cortex01 <- read.delim(paste0("./", names[3], "/chr1/chr1.", names[3], ".dip"), head = F, as.is = T)[,1:2]
cortex02 <- read.delim(paste0("./", names[4], "/chr1/chr1.", names[4], ".dip"), head = F, as.is = T)[,1:2]
ge01 <- read.delim(paste0("./", names[5], "/chr1/chr1.", names[5], ".dip"), head = F, as.is = T)[,1:2]
ge02 <- read.delim(paste0("./", names[6], "/chr1/chr1.", names[6], ".dip"), head = F, as.is = T)[,1:2]
for(chr in chrs[2:length(chrs)]){
  print(chr)
  brain01 <- rbind(brain01, read.delim(paste0("./", names[1], "/", chr,"/", chr,".", names[1], ".dip"), head = F, as.is = T)[,1:2])
  brain02 <- rbind(brain02, read.delim(paste0("./", names[2], "/", chr,"/", chr,".", names[2], ".dip"), head = F, as.is = T)[,1:2])
  cortex01 <- rbind(cortex01, read.delim(paste0("./", names[3], "/", chr,"/", chr,".", names[3], ".dip"), head = F, as.is = T)[,1:2])
  cortex02 <- rbind(cortex02, read.delim(paste0("./", names[4], "/", chr,"/", chr,".", names[4], ".dip"), head = F, as.is = T)[,1:2])
  ge01 <- rbind(ge01, read.delim(paste0("./", names[5], "/", chr,"/", chr,".", names[5], ".dip"), head = F, as.is = T)[,1:2])
  ge02 <- rbind(ge02, read.delim(paste0("./", names[6], "/", chr,"/", chr,".", names[6], ".dip"), head = F, as.is = T)[,1:2])
}
print(all.equal(brain01$V1, brain02$V1))
print(all.equal(brain01$V1, cortex01$V1))
print(all.equal(brain01$V1, cortex02$V1))
print(all.equal(brain01$V1, ge01$V1))
print(all.equal(brain01$V1, ge02$V1))
MeDIP <- data.frame(brain01 = brain01$V2, brain02 = brain02$V2, cortex01 = cortex01$V2, cortex02 = cortex02$V2, ge01 = ge01$V2, ge02 = ge02$V2)
rownames(MeDIP) <- brain01$V1
dim(MeDIP)
MeDIP$chr <- paste0("chr", gsub("_[0-9]+", "", rownames(MeDIP)))
MeDIP$start <- as.numeric(gsub("[0-9X]+_", "", rownames(MeDIP))) + 23
MeDIP$end <- MeDIP$start + 2
save(MeDIP, file = "~/FetalBrain/MeDIP/MeDIP_FetalBrain_fractional.Rdata")
fractional <- data.frame(chr = MeDIP$chr, start = MeDIP$start, end = MeDIP$end, ID = rownames(MeDIP), brain01 = brain01$V2, brain02 = brain02$V2, cortex01 = cortex01$V2, cortex02 = cortex02$V2, ge01 = ge01$V2, ge02 = ge02$V2)
write.table(fractional, file = "~/FetalBrain/MeDIP/MeDIP_FetalBrain_fractional.txt", sep = "\t", col.names = F, row.names = F, quote = F)

####################################################################################################################################################
# compare distribution with WGBS and MethylCRF
## ECDF plots
setwd("/home/lli/FetalBrain/MeDIP/")
load("MeDIP_FetalBrain_fractional.Rdata")
load("~/FetalBrain/MeDIPMRE/MethylCRF.Rdata")
load("~/FetalBrain/WGBS/WGBS_28M.Rdata")
pdf("ECDF_MeDIP.pdf")
plot(c(0, 1), c(0, 1), type = "n", xlab = "Methylation level", ylab = "Ecdf", main = "ECDF of MeDIP")
lines(ecdf(na.omit(MeDIP$brain01)), lwd = 2, col = 1, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$brain02)), lwd = 2, col = 1, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$cortex01)), lwd = 2, col = 2, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$cortex02)), lwd = 2, col = 2, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$ge01)), lwd = 2, col = 3, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$ge02)), lwd = 2, col = 3, lty = 1, cex = 0.1)
legend("topleft", c("Brain", "Cortex", "GE"), col = 1:3, lwd = 5)
plot(c(0, 1), c(0, 1), type = "n", xlab = "Methylation level", ylab = "Ecdf", main = "ECDF of MeDIP")
lines(ecdf(na.omit(MeDIP$brain01)), lwd = 2, col = 1, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$brain02)), lwd = 2, col = 2, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$cortex01)), lwd = 2, col = 1, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$cortex02)), lwd = 2, col = 2, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$ge01)), lwd = 2, col = 1, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$ge02)), lwd = 2, col = 2, lty = 1, cex = 0.1)
legend("topleft", c("HuFNSC01", "HuFNSC02"), col = 1:2, lwd = 5)
dev.off()
cor <- cor(MeDIP)
write.table(cor, file = "correlation between samples.txt", sep = "\t", quote = F)
MeDIPnX <- MeDIP[!grepl("X", rownames(MeDIP)),]
heatmap(cor(MeDIPnX))

pdf("ECDF_MeDIP_MethylCRF_WGBS.pdf")
plot(c(0, 1), c(0, 1), type = "n", xlab = "Methylation level", ylab = "Ecdf", main = "ECDF of MeDIP vs MethylCRF vs WGBS")
lines(ecdf(na.omit(MeDIP$brain01)), lwd = 2, col = 1, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MeDIP$brain02)), lwd = 2, col = 1, lty = 2, cex = 0.1)
lines(ecdf(na.omit(MeDIP$cortex01)), lwd = 2, col = 1, lty = 3, cex = 0.1)
lines(ecdf(na.omit(MeDIP$cortex02)), lwd = 2, col = 1, lty = 4, cex = 0.1)
lines(ecdf(na.omit(MeDIP$ge01)), lwd = 2, col = 1, lty = 5, cex = 0.1)
lines(ecdf(na.omit(MeDIP$ge02)), lwd = 2, col = 1, lty = 6, cex = 0.1)
lines(ecdf(na.omit(MethylCRF$brain01)), lwd = 2, col = 2, lty = 1, cex = 0.1)
lines(ecdf(na.omit(MethylCRF$brain02)), lwd = 2, col = 2, lty = 2, cex = 0.1)
lines(ecdf(na.omit(MethylCRF$cortex01)), lwd = 2, col = 2, lty = 3, cex = 0.1)
lines(ecdf(na.omit(MethylCRF$cortex02)), lwd = 2, col = 2, lty = 4, cex = 0.1)
lines(ecdf(na.omit(MethylCRF$ge01)), lwd = 2, col = 2, lty = 5, cex = 0.1)
lines(ecdf(na.omit(MethylCRF$ge02)), lwd = 2, col = 2, lty = 6, cex = 0.1)
lines(ecdf(na.omit(WGBS_28M$cortex02)), lwd = 2, col = 3, lty = 4, cex = 0.1)
lines(ecdf(na.omit(WGBS_28M$cortex04)), lwd = 2, col = 3, lty = 7, cex = 0.1)
lines(ecdf(na.omit(WGBS_28M$ge02)), lwd = 2, col = 3, lty = 6, cex = 0.1)
lines(ecdf(na.omit(WGBS_28M$ge04)), lwd = 2, col = 3, lty = 8, cex = 0.1)
legend("topleft", c("MeDIP only", "MethylCRF", "WGBS"), col = 1:3, lwd = 5)
plot(c(-1, 1), c(0, 1), type = "n", xlab = "Difference in methylation level", ylab = "Ecdf", main = "ECDF of cortex02-GE02 MeDIP vs MethylCRF vs WGBS")
lines(ecdf(na.omit(MeDIP$cortex02 - MeDIP$ge02)), lwd = 2, col = 1, lty = 1)
lines(ecdf(na.omit(MethylCRF$cortex02 - MethylCRF$ge02)), lwd = 2, col = 2, lty = 2)
lines(ecdf(na.omit(WGBS_28M$cortex02 - WGBS_28M$ge02)), lwd = 2, col = 3, lty = 3)
legend("bottomright", c("MeDIP only", "MethylCRF", "WGBS"), col = 1:3, lty = 1:3, lwd = 5)
dev.off()

## boxplot of methylation levels per chromosome
require(ggplot2, lib.loc = "/home/lli/R")
require(plyr, lib.loc = "/home/lli/R")
MeDIP$method <- "MeDIP_only"
WGBS_28M$method <- "WGBS" 
MethylCRF$method  <- "MethylCRF"

cortex02_3methods <- rbind(MeDIP[,c("chr", "start", "cortex02", "method")], WGBS_28M[,c("chr", "start", "cortex02", "method")], MethylCRF[,c("chr", "start", "cortex02", "method")])
cortex02_3methods$method <- factor(cortex02_3methods$method)
cortex02_3methods <- droplevels(cortex02_3methods[!(cortex02_3methods$chr %in% c("chrM", "chrY")), ])
cortex02_quantile <- ddply(cortex02_3methods, .(chr, method), summarize, "lower" = quantile(cortex02, 0.25, na.rm = T), "middle" = median(cortex02, na.rm = T), "upper" = quantile(cortex02, 0.75, na.rm = T), 
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
ggsave(cortex02_boxplot, file = "cortex02_3methods_boxplot.pdf", width = 15, height = 5)
ge02_3methods <- rbind(MeDIP[,c("chr", "start", "ge02", "method")], WGBS_28M[,c("chr", "start", "ge02", "method")], MethylCRF[,c("chr", "start", "ge02", "method")])
ge02_3methods$method <- factor(ge02_3methods$method)
ge02_3methods <- droplevels(ge02_3methods[!(ge02_3methods$chr %in% c("chrM", "chrY")), ])
ge02_quantile <- ddply(ge02_3methods, .(chr, method), summarize, "lower" = quantile(ge02, 0.25, na.rm = T), "middle" = median(ge02, na.rm = T), "upper" = quantile(ge02, 0.75, na.rm = T), 
                       "ymin" = median(ge02, na.rm = T) - 1.5*(quantile(ge02, 0.75, na.rm = T) - quantile(ge02, 0.25, na.rm = T)), 
                       "ymax" = median(ge02, na.rm = T) + 1.5*(quantile(ge02, 0.75, na.rm = T) - quantile(ge02, 0.25, na.rm = T)))
ge02_quantile$ymax <- 1
ge02_quantile$chr <- factor(gsub("chr", "", as.character(ge02_quantile$chr)), levels = c(as.character(1:22), "X"))
ge02_boxplot <- ggplot(ge02_quantile, aes(x = chr, lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax, fill = chr)) + 
  geom_boxplot(stat = "identity") +   
  facet_wrap(~ method) + 
  xlab("chr") + 
  ylab("methylation level") + 
  ggtitle("WGBS vs MethylCRF - ge02") + 
  guides(fill = F) + 
  theme_bw()
ggsave(ge02_boxplot, file = "ge02_3methods_boxplot.pdf", width = 15, height = 5)

####################################################################################################################################################

# bin into 200bp bins
options(scipen = 10)
bin200 <- read.csv("~/hg19/bin200.csv")
bin200$id <- paste0(bin200$chr, ":", bin200$start, "-", bin200$end)
write.table(bin200[, c("chr", "start", "end", "id")], file = "~/hg19/bin200.bed", sep = "\t", col.names = F, row.names = F, quote = F)
CG <- read.delim("~/hg19/CG.BED", head = F, as.is = T)
# CG$id <- paste0(CG$V1, ":", CG$V2, "-", CG$V3)
# write.table(CG, file = "~/hg19/CG.BED", sep = "\t", col.names = F, row.names = F, quote = F)
# less ~/hg19/CG.BED | awk '{print $0"\t"$1":"$2"-"$3}' > ~/hg19/CpG.bed
# /Applications/bedtools-2.17.0/bin/intersectBed -a ~/hg19/CpG.bed -b ~/hg19/bin200.bed -wa -wb > ~/hg19/CG_200bin.txt
# /Applications/bedtools-2.17.0/bin/intersectBed -a ~/hg19/CpG.bed -b ~/hg19/bin200.bed -v > ~/hg19/CG_nobin.txt
# less ~/hg19/CG_200bin.txt | awk '{print $4"\t"$8}' > ~/hg19/CpG_200bin.txt

# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R
require(plyr, lib.loc = "/home/lli/R")
setwd("/home/lli/FetalBrain/MeDIP/")
load("MeDIP_FetalBrain_fractional.Rdata")
load("~/FetalBrain/MeDIPMRE/MethylCRF.Rdata")
load("~/FetalBrain/WGBS/WGBS_28M.Rdata")
load("~/hg19/CpG_200bin.Rdata")
# CG_bin <- read.delim("~/hg19/CpG_200bin.txt", head = F, as.is = T)
# rownames(CG_bin) <- CG_bin$V1
# save(CG_bin, file = "~/hg19/CpG_200bin.Rdata")
# rownames(MeDIP) <- paste0(MeDIP$chr, ":", MeDIP$start, "-", MeDIP$end)
# save(MeDIP, file = "MeDIP_FetalBrain_fractional.Rdata")
# rownames(MethylCRF) <- paste0(MethylCRF$chr, ":", MethylCRF$start, "-", MethylCRF$end)
# save(MethylCRF, file = "~/FetalBrain/MeDIPMRE/MethylCRF.Rdata")

# MeDIP
MeDIP$bin <- CG_bin[rownames(MeDIP), ]$V2 
MeDIP_bin <- ddply(MeDIP, ~ bin, summarise, brain01 = mean(brain01), brain02 = mean(brain02), cortex01 = mean(cortex01), cortex02 = mean(cortex02), ge01 = mean(ge01), ge02 = mean(ge02))
rownames(MeDIP_bin) <- MeDIP_bin$bin
cortex02_ge02_MeDIP <- ddply(MeDIP, ~bin, summarize, p = try(t.test(cortex02, ge02, paired = T)$p.value, silent = T), cortex02_ge02 = abs(mean(cortex02)-mean(ge02)))
cortex02_ge02_MeDIP <- cortex02_ge02_MeDIP[!grepl("Error", cortex02_ge02_MeDIP$p), ]
cortex02_ge02_MeDIP$adjusted.p <- p.adjust(cortex02_ge02_MeDIP$p, method = "fdr")
rownames(cortex02_ge02_MeDIP) <- cortex02_ge02_MeDIP$bin

# MethylCRF
MethylCRF$bin <- CG_bin[rownames(MethylCRF), ]$V2 
MethylCRF_bin <- ddply(MethylCRF, ~ bin, summarise, brain01 = mean(brain01), brain02 = mean(brain02), cortex01 = mean(cortex01), cortex02 = mean(cortex02), ge01 = mean(ge01), ge02 = mean(ge02))
rownames(MethylCRF_bin) <- MethylCRF_bin$bin
cortex02_ge02_MethylCRF <- ddply(MethylCRF, ~bin, summarize, p = try(t.test(cortex02, ge02, paired = T)$p.value, silent = T), cortex02_ge02 = abs(mean(cortex02)-mean(ge02)))
cortex02_ge02_MethylCRF <- cortex02_ge02_MethylCRF[!grepl("Error", cortex02_ge02_MethylCRF$p), ]
cortex02_ge02_MethylCRF$adjusted.p <- p.adjust(cortex02_ge02_MethylCRF$p, method = "fdr")
rownames(cortex02_ge02_MethylCRF) <- cortex02_ge02_MethylCRF$bin

# WGBS_28M
WGBS_28M$bin <- CG_bin[rownames(WGBS_28M), ]$V2 
WGBS_28M_bin <- ddply(WGBS_28M, ~ bin, summarise, cortex02 = mean(cortex02, na.rm = T), ge02 = mean(ge02, na.rm = T), cortex04 = mean(cortex04, na.rm = T), ge04 = mean(ge04, na.rm = T))
rownames(WGBS_28M_bin) <- WGBS_28M_bin$bin
cortex02_ge02_WGBS_28M <- ddply(WGBS_28M, ~bin, summarize, p = try(t.test(cortex02, ge02, paired = T)$p.value, silent = T), cortex02_ge02 = abs(mean(cortex02, na.rm = T)-mean(ge02, na.rm = T)))
cortex02_ge02_WGBS_28M <- cortex02_ge02_WGBS_28M[!grepl("Error", cortex02_ge02_WGBS_28M$p), ]
cortex02_ge02_WGBS_28M$adjusted.p <- p.adjust(cortex02_ge02_WGBS_28M$p, method = "fdr")
rownames(cortex02_ge02_WGBS_28M) <- cortex02_ge02_WGBS_28M$bin

save(MeDIP_bin, cortex02_ge02_MeDIP, MethylCRF_bin, cortex02_ge02_MethylCRF, WGBS_28M_bin, cortex02_ge02_WGBS_28M, file = "compare3_bin200.Rdata")

# correlation between MeDIP, MethylCRF and WGBS
# install.packages("hexbin", lib = "~/R/")
require(hexbin, lib.loc = "/home/lli/R")
inx <- intersect(rownames(MeDIP_bin), intersect(rownames(MethylCRF_bin), rownames(WGBS_28M_bin)))
cortex02_bin200 <- data.frame(MeDIP = MeDIP_bin[inx, "cortex02"], MethylCRF = MethylCRF_bin[inx, "cortex02"], WGBS = WGBS_28M_bin[inx, "cortex02"])
ge02_bin200 <- data.frame(MeDIP = MeDIP_bin[inx, "ge02"], MethylCRF = MethylCRF_bin[inx, "ge02"], WGBS = WGBS_28M_bin[inx, "ge02"])
pdf("bin200_3methods_correlation.pdf")
splom(cortex02_bin200, panel = panel.hexbinplot, main = "correlation between methods - cortex02")
splom(ge02_bin200, panel = panel.hexbinplot, main = "correlation between methods - ge02")
dev.off()

# cortex02 vs ge02 DMRs 
cutoff <- 0.05
cortex02_ge02_MeDIP <- cortex02_ge02_MeDIP[cortex02_ge02_MeDIP$adjusted.p < cutoff,]
cortex02_ge02_MethylCRF <- cortex02_ge02_MethylCRF[cortex02_ge02_MethylCRF$adjusted.p < cutoff,]
cortex02_ge02_WGBS_28M <- cortex02_ge02_WGBS_28M[cortex02_ge02_WGBS_28M$adjusted.p < cutoff,]
save(cortex02_ge02_MeDIP, cortex02_ge02_MethylCRF, cortex02_ge02_WGBS_28M, file = "DMR_3methods_bin200.Rdata")

######################################################################################
# DNA methylation asymmetry between MZ twins
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
setwd("~/FetalBrain/MeDIP/")
load("MeDIP_FetalBrain_fractional.Rdata")

e <- 1e-6
MeDIP$logfold_brain <- log2((MeDIP$brain01 + e) / (MeDIP$brain02 + e))
MeDIP$logfold_cortex <- log2((MeDIP$cortex01 + e) / (MeDIP$cortex02 + e))
MeDIP$logfold_ge <- log2((MeDIP$ge01 + e) / (MeDIP$ge02 + e))
# MeDIP$logfold_brain <- log2((MeDIP$brain01) / (MeDIP$brain02))
# MeDIP$logfold_cortex <- log2((MeDIP$cortex01) / (MeDIP$cortex02))
# MeDIP$logfold_ge <- log2((MeDIP$ge01) / (MeDIP$ge02))
# fold_brain <- MeDIP[is.finite(MeDIP$logfold_brain), ]
# fold_cortex <- MeDIP[is.finite(MeDIP$logfold_cortex), ]
# fold_ge <- MeDIP[is.finite(MeDIP$logfold_ge), ]
fold_brain <- MeDIP[abs(MeDIP$logfold_brain) > 1, ]
fold_cortex <- MeDIP[abs(MeDIP$logfold_cortex) > 1, ]
fold_ge <- MeDIP[abs(MeDIP$logfold_ge) > 1, ]

# xrange <- range(fold_brain$logfold_brain, fold_cortex$logfold_cortex, fold_ge$logfold_ge)
# yrange <- range(density(fold_brain$logfold_brain)$y, density(fold_cortex$logfold_cortex)$y, density(fold_ge$logfold_ge)$y)
xrange <- c(-10, 10)
yrange <- c(0, 0.05)
pdf(file="Asymmetry.pdf")
plot(xrange, yrange, type="n", main="Global DNA Methylation", xlab="log2(MeDIP01/MeDIP02)", ylab="density")
lines(density(fold_brain$logfold_brain, adjust = 0.8), col=1, lty=1, lwd=3)
lines(density(fold_cortex$logfold_cortex, adjust = 0.8), col=2, lty=1, lwd=3)
lines(density(fold_ge$logfold_ge, adjust = 0.8), col=3, lty=1, lwd=3)
legend("topright", c("brain","cortex","GE"), col=c(1:3), lty=1, lwd=5)
dev.off()


