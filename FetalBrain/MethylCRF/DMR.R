# read in MethylCRF data and compute pairwise differences in methylation levels between samples
setwd("~/FetalBrain/MeDIPMRE/")
load("MethylCRF.Rdata")
MethylCRF$brain01_brain02 <- MethylCRF$brain01 - MethylCRF$brain02
MethylCRF$cortex01_cortex02 <- MethylCRF$cortex01 - MethylCRF$cortex02
MethylCRF$ge01_ge02 <- MethylCRF$ge01 - MethylCRF$ge02
MethylCRF$cortex01_ge01 <- MethylCRF$cortex01 - MethylCRF$ge01
MethylCRF$cortex02_ge02 <- MethylCRF$cortex02 - MethylCRF$ge02
save(MethylCRF, file= "MethylCRF.Rdata")
summary <- data.frame(mean = c(mean(MethylCRF$brain01_brain02), mean(MethylCRF$cortex01_cortex02), mean(MethylCRF$ge01_ge02), mean(MethylCRF$cortex01_ge01), mean(MethylCRF$cortex02_ge02)), 
                      sd = c(sd(MethylCRF$brain01_brain02), sd(MethylCRF$cortex01_cortex02), sd(MethylCRF$ge01_ge02), sd(MethylCRF$cortex01_ge01), sd(MethylCRF$cortex02_ge02)), 
                      min = c(min(MethylCRF$brain01_brain02), min(MethylCRF$cortex01_cortex02), min(MethylCRF$ge01_ge02), min(MethylCRF$cortex01_ge01), min(MethylCRF$cortex02_ge02)), 
                      max = c(max(MethylCRF$brain01_brain02), max(MethylCRF$cortex01_cortex02), max(MethylCRF$ge01_ge02), max(MethylCRF$cortex01_ge01), max(MethylCRF$cortex02_ge02)))
rownames(summary) <- c("brain01_brain02", "cortex01_cortex02", "ge01_ge02", "cortex01_ge01", "cortex02_ge02")
summary

# plot differences in methylation levels distribution 
x <- seq(-1, 1, length = 10000)
pdf("DMdensity.pdf")
plot(c(-0.5, 0.5), c(0, 10), type = "n", xlab = "Difference in methylation", ylab = "Density")
lines(x, dnorm(x, mean = mean(summary$mean), sd = mean(summary$sd)), lwd = 1, col = 1, lty = 2)
lines(density(MethylCRF$brain01_brain02, adjust = 20), lwd = 1, col = 2)
lines(density(MethylCRF$cortex01_cortex02, adjust = 20), lwd = 2, col = 3)
lines(density(MethylCRF$ge01_ge02, adjust = 20), lwd = 2, col = 4)
lines(density(MethylCRF$cortex01_ge01, adjust = 20), lwd = 2, col = 5)
lines(density(MethylCRF$cortex02_ge02, adjust = 20), lwd = 2, col = 6)
aline(v = 0)
legend("topright", c("Normal", "brain01_brain02", "cortex01_cortex02", "ge01_ge02", "cortex01_ge01", "cortex02_ge02"), col = 1:6, lty = c(2, rep(1, times=5)), lwd = 2)
dev.off()
pdf("DMecdf.pdf")
# cutoffs for differential methylation
plot(c(-1, 1), c(0, 1), type = "n", xlab = "Difference in methylation", ylab = "Ecdf")
lines(x, pnorm(x, mean = mean(summary$mean), sd = mean(summary$sd)), lwd = 1, col = 1, lty = 2)
lines(ecdf(MethylCRF$brain01_brain02), lwd = 1, col = 2)
lines(ecdf(MethylCRF$cortex01_cortex02), lwd = 2, col = 3)
lines(ecdf(MethylCRF$ge01_ge02), lwd = 2, col = 4)
lines(ecdf(MethylCRF$cortex01_ge01), lwd = 2, col = 5)
lines(ecdf(MethylCRF$cortex02_ge02), lwd = 2, col = 6)
legend("bottomright", c("Normal", "brain01_brain02", "cortex01_cortex02", "ge01_ge02", "cortex01_ge01", "cortex02_ge02"), col = 1:6, lty = c(2, rep(1, times=5)), lwd = 2)
dev.off()
pdf("DMcutoffneg.pdf")
plot(c(-0.8, -0.2), c(0, 0.005), type = "n", xlab = "Difference in methylation", ylab = "Ecdf")
lines(ecdf(MethylCRF$brain01_brain02), lwd = 1, col = 2)
lines(ecdf(MethylCRF$cortex01_cortex02), lwd = 2, col = 3)
lines(ecdf(MethylCRF$ge01_ge02), lwd = 2, col = 4)
lines(ecdf(MethylCRF$cortex01_ge01), lwd = 2, col = 5)
lines(ecdf(MethylCRF$cortex02_ge02), lwd = 2, col = 6)
abline(v = mean(summary$mean) - 6*mean(summary$sd))
legend("bottomright", c("brain01_brain02", "cortex01_cortex02", "ge01_ge02", "cortex01_ge01", "cortex02_ge02"), col = 2:6, lty = 1, lwd = 2)
dev.off()
pdf("DMcutoffpos.pdf")
plot(c(0.2, 0.8), c(0.993, 1), type = "n", xlab = "Difference in methylation", ylab = "Ecdf")
lines(ecdf(MethylCRF$brain01_brain02), lwd = 1, col = 2)
lines(ecdf(MethylCRF$cortex01_cortex02), lwd = 2, col = 3)
lines(ecdf(MethylCRF$ge01_ge02), lwd = 2, col = 4)
lines(ecdf(MethylCRF$cortex01_ge01), lwd = 2, col = 5)
lines(ecdf(MethylCRF$cortex02_ge02), lwd = 2, col = 6)
abline(v = mean(summary$mean) + 6*mean(summary$sd))
legend("bottomright", c("brain01_brain02", "cortex01_cortex02", "ge01_ge02", "cortex01_ge01", "cortex02_ge02"), col = 2:6, lty = 1, lwd = 2)
dev.off()

setwd("~/FetalBrain/MeDIPMRE/")
load("MethylCRF.Rdata")
cutoff = 6 # mean +/- cutoff*SD, from Ecdf plot, p < 0.01
MethylCRF$DM.brain01_brain02 <- factor("same", levels = c("hypo", "same", "hyper")) 
MethylCRF[MethylCRF$brain01_brain02 > mean(MethylCRF$brain01_brain02)+cutoff*sd(MethylCRF$brain01_brain02),]$DM.brain01_brain02 <- "hyper"
MethylCRF[MethylCRF$brain01_brain02 < mean(MethylCRF$brain01_brain02)-cutoff*sd(MethylCRF$brain01_brain02),]$DM.brain01_brain02 <- "hypo"
summary(MethylCRF$DM.brain01_brain02)
# MethylCRF[MethylCRF$DM.brain01_brain02!="same"
#           &(mapply(min, MethylCRF$brain01, MethylCRF$brain02, simplify=T)>0.4
#             |mapply(max, MethylCRF$brain01, MethylCRF$brain02, simplify=T)<0.6),]$DM.brain01_brain02 <- "same"
# summary(MethylCRF$DM.brain01_brain02)
# MethylCRF[MethylCRF$DM.brain01_brain02!="same"
#           &(mapply(min, MethylCRF$brain01, MethylCRF$brain02, simplify=T)>0.3
#             |mapply(max, MethylCRF$brain01, MethylCRF$brain02, simplify=T)<0.7),]$DM.brain01_brain02 <- "same"
# summary(MethylCRF$DM.brain01_brain02)
# range(mapply(min, MethylCRF[MethylCRF$DM.brain01_brain02!="same", ]$brain01, MethylCRF[MethylCRF$DM.brain01_brain02!="same", ]$brain02))
# range(mapply(max, MethylCRF[MethylCRF$DM.brain01_brain02!="same", ]$brain01, MethylCRF[MethylCRF$DM.brain01_brain02!="same", ]$brain02))

# collapse adjacent DM CpGs into DMRs
setwd("~/FetalBrain/MeDIPMRE/")
load("MethylCRF.Rdata")
cutoff = 6 # mean +/- cutoff*SD, from Ecdf plot, p < 0.01
MethylCRF$DM.brain01_brain02 <- factor("same", levels = c("hypo", "same", "hyper")) 
MethylCRF[MethylCRF$brain01_brain02 > mean(MethylCRF$brain01_brain02)+cutoff*sd(MethylCRF$brain01_brain02),]$DM.brain01_brain02 <- "hyper"
MethylCRF[MethylCRF$brain01_brain02 < mean(MethylCRF$brain01_brain02)-cutoff*sd(MethylCRF$brain01_brain02),]$DM.brain01_brain02 <- "hypo"
DM.brain01_brain02 <- data.frame(chr = MethylCRF$chr, start = MethylCRF$start, end = MethylCRF$end, id = MethylCRF$id, DM = MethylCRF$DM.brain01_brain02)
write.table(DM.brain01_brain02, "DM.brain01_brain02.1.txt", sep="\t", row.names=F, col.names=F, quote=F)
# ~/FetalBrain/scripts/collapseDM.sh
require(plyr, lib.loc = "/home/lli/R")
require(ggplot2, lib.loc = "/home/lli/R")
collapse.brain01_brain02 <- read.delim("collapse.brain01_brain02.txt", head = F, as.is = T)
colnames(collapse.brain01_brain02) <- c("id", "DM", "DMR")
all.equal(collapse.brain01_brain02$id, MethylCRF$id)
MethylCRF$DMR.brain01_brain02 <- collapse.brain01_brain02$DMR
DMR.brain01_brain02 <- MethylCRF[MethylCRF$DMR.brain01_brain02 != 0,]
DMR.brain01_brain02$DMR.brain01_brain02 <- factor(as.character(DMR.brain01_brain02$DMR.brain01_brain02))
DMR.brain01_brain02$count <- 1
DMregion.brain01_brain02 <- ddply(DMR.brain01_brain02, ~DMR.brain01_brain02, summarize, chr = chr[1], start = min(start), end = max(end), DM = DM.brain01_brain02[1], count = sum(count))
DMregion.brain01_brain02$id <- paste0(DMregion.brain01_brain02$chr, ":", DMregion.brain01_brain02$start, "-", DMregion.brain01_brain02$end)
DMregion.brain01_brain02$length <- DMregion.brain01_brain02$end - DMregion.brain01_brain02$start
DMregion.brain01_brain02$chr <- factor(DMregion.brain01_brain02$chr, levels = c(paste0("chr",as.character(1:22)), "chrX"))
DMregion.brain01_brain02 <- DMregion.brain01_brain02[order(DMregion.brain01_brain02$chr, DMregion.brain01_brain02$start),]
DMregion.brain01_brain02$pos <- (DMregion.brain01_brain02$start + DMregion.brain01_brain02$end)/2
DMregion.brain01_brain02$dis <- c(0, DMregion.brain01_brain02[2:nrow(DMregion.brain01_brain02),]$pos - DMregion.brain01_brain02[1:nrow(DMregion.brain01_brain02)-1,]$pos)
DMregion.brain01_brain02[DMregion.brain01_brain02$dis<0,]$dis <- 0
# should collapse nearby DMRs before proceed but no time
DMregion.brain01_brain02 <- DMregion.brain01_brain02[DMregion.brain01_brain02$count>1,]
DMregion.brain01_brain02 <- na.omit(DMregion.brain01_brain02)
write.table(DMregion.brain01_brain02, file = "DMR_brain12_dynamic.txt", sep = "\t", quote = F, row.names = F)
write.table(DMregion.brain01_brain02[,c("chr", "start", "end", "id")], file = "DMR_brain12_dynamic.bed", sep = "\t", quote = F, row.names = F, col.names = F)
DMR.length <- ggplot(DMregion.brain01_brain02, aes(chr, length, fill = chr)) + 
  geom_boxplot() + 
  geom_abline(intercept = mean(DMregion.brain01_brain02$length), show_guide = T) + 
  xlab("chr") + 
  ylab("length (bp)") + 
  ggtitle("Dynamic growth approach - DMR length") + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR.length, file = "DMRlength_dynamic.pdf", width = 10, height = 5)

DMR.dis <- ggplot(DMregion.brain01_brain02, aes(chr, dis, fill = chr)) + 
  geom_boxplot() + 
  geom_abline(intercept = mean(DMregion.brain01_brain02$dis), show_guide = T) + 
  xlab("chr") + 
  ylab("distance between adjacent DMRs (bp)") + 
  ggtitle("Dynamic growth approach - distance between adjacent DMRs") + 
  guides(fill = F) + 
  coord_cartesian(ylim = c(0,10000)) + 
  theme_bw()
ggsave(DMR.dis, file = "DMRdis_dynamic.pdf", width = 10, height = 5)

pdf("DMRdis_dynamic.pdf")
plot(x = c(0, 10000), y = c(0, 1), type= "n", main = "Dynamic growth approach - distance between adjacent DMRs", xlab = "distance between adjacent DMRs (bp)", ylab = "Ecdf")
lines(ecdf(DMregion.brain01_brain02$dis))
dev.off()

# Asymmetry?
# DMregion.brain01_brain02 <- read.delim("DMR_brain12_dynamic.txt", as.is = T)
DMregion.brain01_brain02$DM <- factor(DMregion.brain01_brain02$DM)
DMregion.brain01_brain02$chr <- factor(DMregion.brain01_brain02$chr, levels = c("chrX", paste0("chr",as.character(22:1))))
(DMregion.brain01_brain02_asymmetry <- ggplot(DMregion.brain01_brain02, aes(x = chr, fill = factor(DM))) + 
   geom_bar(width = .7) + 
   coord_flip() + 
   ggtitle("Global DNA methylation between MZ twins") + 
   ylab("No. of DMRs") +
   xlab("chr") +
   theme_bw())
ggsave(DMregion.brain01_brain02_asymmetry, file = "DMregion.brain01_brain02_asymmetry.pdf")

chrlength <- read.csv("~/快盘/hg19/chrlen_hg19.csv", as.is = T)[1:23,]
chrlength$chr <- factor(chrlength$chr, levels = c(paste0("chr",as.character(1:22)), "chrX"))
(DMregion.brain01_brain02_coord <- ggplot(DMregion.brain01_brain02) + 
   geom_linerange(aes(factor(chr, levels = c("chrX", paste0("chr",as.character(22:1)))), ymin = 0, ymax = length), data = chrlength, alpha = 0.5) + 
   geom_point(aes(x = as.numeric(chr) + 0.5*(as.numeric(as.factor(DM))-1.5), y = pos, color = factor(DM)), position = position_jitter(width = 0.05), size = 0.5, alpha = 0.5) +  
   ggtitle("DMR positions on chromosomes") + 
   ylab("Position of DMRs") +
   xlab("chr") +
   coord_flip() + 
   theme_bw())
ggsave(DMregion.brain01_brain02_coord, file = "DMregion.brain01_brain02_coord.pdf")

# intersecting with genomic regions: bedtools  
# /home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic_promoter.txt
# /home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic_genebody.txt
# /home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic_CGI.txt
setwd("~/快盘/FetalBrain/MeDIPMRE/")
brain12_promoter <- read.delim("DMR_brain12_dynamic_promoter.txt", head = F, as.is = T)
brain12_promoter$Ensembl <- gsub("_[a-zA-Z0-9_]+", "", brain12_promoter$V8)
brain12_promoter_pc <- brain12_promoter[grepl("protein_coding", brain12_promoter$V8),]
write.table(brain12_promoter, file = "DMR_brain12_dynamic_promoter.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(brain12_promoter_pc, file = "DMR_brain12_dynamic_promoter_pc.txt", sep = "\t", quote = F, row.names = F, col.names = F)

brain12_genebody <- read.delim("DMR_brain12_dynamic_genebody.txt", head = F, as.is = T)
brain12_genebody$Ensembl <- gsub("_[a-zA-Z0-9_]+", "", brain12_genebody$V8)
brain12_genebody_pc <- brain12_genebody[grepl("protein_coding", brain12_genebody$V8),]
write.table(brain12_genebody, file = "DMR_brain12_dynamic_genebody.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(brain12_genebody_pc, file = "DMR_brain12_dynamic_genebody_pc.txt", sep = "\t", quote = F, row.names = F, col.names = F)

brain12_CGI <- read.delim("DMR_brain12_dynamic_CGI.txt", head = F, as.is = T)


