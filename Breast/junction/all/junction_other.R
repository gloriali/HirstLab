# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# read junction files
setwd("/projects/mbilenky/REMC/breast/RNA-seq/junctions_new/")
lum084 <- read.delim("A17918.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
myo084 <- read.delim("A17919.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
stem084 <- read.delim("A17920.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
lum080 <- read.delim("A01029.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
myo080 <- read.delim("A01030.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
stem080 <- read.delim("A01031.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
lum035 <- read.delim("HS1187.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
myo035 <- read.delim("HS1188.q5.F516.jb.bedGraph.gz", head = F, as.is = T)

setwd("~/REMC/junction/all/")

# compare the new junction files to old ones 
load("~/REMC/junction/junction_RM084.Rdata")
# load("~/快盘/REMC/junction/junction_valid.Rdata")
# save(lum084_myo084_isoform_all, lum084_stem084_isoform_all, myo084_stem084_isoform_all, lum084_myo084_isoform_all_gene, lum084_stem084_isoform_all_gene, myo084_stem084_isoform_all_gene, file = "RM084_isoform_all.Rdata")
compare <- data.frame(coverage = c(junction[junction$lum084N > 0,]$lum084N, junction[junction$myo084N > 0,]$myo084N, junction[junction$stem084N > 0,]$stem084N, 
                               lum084$V3, myo084$V3, stem084$V3), 
                  sample = factor(c(rep("lum.old", nrow(junction[junction$lum084N > 0,])), rep("myo.old", nrow(junction[junction$myo084N > 0,])), rep("stem.old", nrow(junction[junction$stem084N > 0,])), 
                             rep("lum.new", nrow(lum084)), rep("myo.new", nrow(myo084)), rep("stem.new", nrow(stem084))), levels = c("lum.old", "lum.new", "myo.old", "myo.new", "stem.old", "stem.new")), 
                  version = c(rep("old", (nrow(junction[junction$lum084N > 0,]) + nrow(junction[junction$myo084N > 0,]) + nrow(junction[junction$stem084N > 0,]))), 
                              rep("new", (nrow(lum084) + nrow(myo084) + nrow(stem084)))))
require(ggplot2, lib.loc = "~/R-3.0.2/")
require(labeling, lib.loc = "~/R-3.0.2/")
(comparebox <- ggplot(compare, aes(x = sample, y = coverage, fill = version)) + 
   geom_boxplot() + 
   theme_bw() + 
   ggtitle("Compare new version junction coverage with the old") + 
   coord_cartesian(ylim = c(0, 500)))
ggsave(comparebox, file = "compare.pdf")

junction <- data.frame(id = unique(c(lum084$V1, myo084$V1, stem084$V1, lum080$V1, myo080$V1, stem080$V1, lum035$V1, myo035$V1)))
rownames(junction) <- junction$id
junction[, c("lum084N", "myo084N", "stem084N", "lum080N", "myo080N", "stem080N", "lum035N", "myo035N")] <- 0
junction[lum084$V1, "lum084N"] <- lum084$V3
junction[myo084$V1, "myo084N"] <- myo084$V3
junction[stem084$V1, "stem084N"] <- stem084$V3
junction[lum080$V1, "lum080N"] <- lum080$V3
junction[myo080$V1, "myo080N"] <- myo080$V3
junction[stem080$V1, "stem080N"] <- stem080$V3
junction[lum035$V1, "lum035N"] <- lum035$V3
junction[myo035$V1, "myo035N"] <- myo035$V3

# calculate RPKM: 10^3*10^6 / (N*read_length)  
Nlum084 <- 191290094   # from /projects/epigenomics/ep50/internal/jqc.1.7.6/A17918/A17918.report: Total number of exonic reads for RPKM (protein coding; no MT, no ribo proteins, top expressed 0.005 exons excluded) 
Nmyo084 <- 215406797
Nstem084 <- 220214822
Nlum080 <- 54482020
Nmyo080 <- 58851034
Nstem080 <- 55793105
Nlum035 <- 82129610
Nmyo035 <- 103893357
# check read length: /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools view /projects/epigenomics/ep50/internal/jqc.1.7.6/A17918/A17918.bam -- 75M
read_length <- 75
junction$lum084rpkm <- junction$lum084N*10^3*10^6 / (Nlum084*read_length) 
junction$myo084rpkm <- junction$myo084N*10^3*10^6 / (Nmyo084*read_length) 
junction$stem084rpkm <- junction$stem084N*10^3*10^6 / (Nstem084*read_length) 
junction$lum080rpkm <- junction$lum080N*10^3*10^6 / (Nlum080*read_length) 
junction$myo080rpkm <- junction$myo080N*10^3*10^6 / (Nmyo080*read_length) 
junction$stem080rpkm <- junction$stem080N*10^3*10^6 / (Nstem080*read_length) 
junction$lum035rpkm <- junction$lum035N*10^3*10^6 / (Nlum035*read_length) 
junction$myo035rpkm <- junction$myo035N*10^3*10^6 / (Nmyo035*read_length) 
save(junction, file = "junction_all.Rdata")

# ECDF of coverage & RPKM
pdf("coverage_distribution.pdf")
plot(c(0, 500), c(0, 1), type = "n", main = "Ecdf of junction read count", xlab = "coverage", ylab = "Ecdf")
lines(ecdf(junction$lum084N), col = 2, lty = 1, lwd = 2)
lines(ecdf(junction$myo084N), col = 3, lty = 1, lwd = 2)
lines(ecdf(junction$stem084N), col = 4, lty = 1, lwd = 2)
lines(ecdf(junction$lum080N), col = 2, lty = 2, lwd = 2)
lines(ecdf(junction$myo080N), col = 3, lty = 2, lwd = 2)
lines(ecdf(junction$stem080N), col = 4, lty = 2, lwd = 2)
lines(ecdf(junction$lum035N), col = 2, lty = 3, lwd = 2)
lines(ecdf(junction$myo035N), col = 3, lty = 3, lwd = 2)
legend("bottomright", c("lum", "myo", "stem-like"), col = 2:4, lty = 1, lwd = 5)
legend("topleft", c("RM084", "RM080", "RM035"), col = 1, lty = 1:3, lwd = 5)
plot(c(0, 50), c(0, 1), type = "n", main = "Ecdf of junction read count", xlab = "coverage", ylab = "Ecdf")
lines(ecdf(junction$lum084N), col = 2, lty = 1, lwd = 2)
lines(ecdf(junction$myo084N), col = 3, lty = 1, lwd = 2)
lines(ecdf(junction$stem084N), col = 4, lty = 1, lwd = 2)
lines(ecdf(junction$lum080N), col = 2, lty = 2, lwd = 2)
lines(ecdf(junction$myo080N), col = 3, lty = 2, lwd = 2)
lines(ecdf(junction$stem080N), col = 4, lty = 2, lwd = 2)
lines(ecdf(junction$lum035N), col = 2, lty = 3, lwd = 2)
lines(ecdf(junction$myo035N), col = 3, lty = 3, lwd = 2)
abline(v = 2) 
legend("bottomright", c("lum", "myo", "stem-like"), col = 2:4, lty = 1, lwd = 5)
legend("topleft", c("RM084", "RM080", "RM035"), col = 1, lty = 1:3, lwd = 5)
dev.off()
pdf("RPKM_distribution.pdf")
plot(c(0, 50), c(0, 1), type = "n", main = "Ecdf of junction RPKM", xlab = "RPKM", ylab = "Ecdf")
lines(ecdf(junction$lum084rpkm), col = 2, lty = 1, lwd = 2)
lines(ecdf(junction$myo084rpkm), col = 3, lty = 1, lwd = 2)
lines(ecdf(junction$stem084rpkm), col = 4, lty = 1, lwd = 2)
lines(ecdf(junction$lum080rpkm), col = 2, lty = 2, lwd = 2)
lines(ecdf(junction$myo080rpkm), col = 3, lty = 2, lwd = 2)
lines(ecdf(junction$stem080rpkm), col = 4, lty = 2, lwd = 2)
lines(ecdf(junction$lum035rpkm), col = 2, lty = 3, lwd = 2)
lines(ecdf(junction$myo035rpkm), col = 3, lty = 3, lwd = 2)
abline(v = 0.1) 
legend("bottomright", c("lum", "myo", "stem-like"), col = 2:4, lty = 1, lwd = 5)
legend("topleft", c("RM084", "RM080", "RM035"), col = 1, lty = 1:3, lwd = 5)
plot(c(0, 3), c(0, 1), type = "n", main = "Ecdf of junction RPKM", xlab = "RPKM", ylab = "Ecdf")
lines(ecdf(junction$lum084rpkm), col = 2, lty = 1, lwd = 2)
lines(ecdf(junction$myo084rpkm), col = 3, lty = 1, lwd = 2)
lines(ecdf(junction$stem084rpkm), col = 4, lty = 1, lwd = 2)
lines(ecdf(junction$lum080rpkm), col = 2, lty = 2, lwd = 2)
lines(ecdf(junction$myo080rpkm), col = 3, lty = 2, lwd = 2)
lines(ecdf(junction$stem080rpkm), col = 4, lty = 2, lwd = 2)
lines(ecdf(junction$lum035rpkm), col = 2, lty = 3, lwd = 2)
lines(ecdf(junction$myo035rpkm), col = 3, lty = 3, lwd = 2)
abline(v = 0.1) 
legend("bottomright", c("lum", "myo", "stem-like"), col = 2:4, lty = 1, lwd = 5)
legend("topleft", c("RM084", "RM080", "RM035"), col = 1, lty = 1:3, lwd = 5)
dev.off()

require(ggplot2, lib.loc = "~/R-3.0.2/")
require(labeling, lib.loc = "~/R-3.0.2/")
coverage <- data.frame(Coverage = c(junction$lum084N, junction$myo084N, junction$stem084N, junction$lum080N, junction$myo080N, junction$stem080N, junction$lum035N, junction$myo035N), 
                       Sample = factor(rep(c("lum084", "myo084", "stem084", "lum080", "myo080", "stem080", "lum035", "myo035"), each = nrow(junction)), levels = c("lum035", "myo035", "lum080", "myo080", "stem080", "lum084", "myo084", "stem084")), 
                       Donor = c(rep(c("RM084", "RM080"), each = 3*nrow(junction)), rep("RM035", 2*nrow(junction))))
(coverage_box <- ggplot(coverage, aes(x = Sample, y = Coverage, fill = Donor)) + 
   geom_boxplot() + 
   theme_bw() + 
   ggtitle("Boxplot of junction coverage") + 
   coord_cartesian(ylim = c(0, 200)))
ggsave(coverage_box, file = "coverage_boxplot.pdf")
rpkm <- data.frame(RPKM = c(junction$lum084rpkm, junction$myo084rpkm, junction$stem084rpkm, junction$lum080rpkm, junction$myo080rpkm, junction$stem080rpkm, junction$lum035rpkm, junction$myo035rpkm), 
                       Sample = factor(rep(c("lum084", "myo084", "stem084", "lum080", "myo080", "stem080", "lum035", "myo035"), each = nrow(junction)), levels = c("lum035", "myo035", "lum080", "myo080", "stem080", "lum084", "myo084", "stem084")), 
                       Donor = c(rep(c("RM084", "RM080"), each = 3*nrow(junction)), rep("RM035", 2*nrow(junction))))
(rpkm_box <- ggplot(rpkm, aes(x = Sample, y = RPKM, fill = Donor)) + 
   geom_boxplot() + 
   theme_bw() + 
   ggtitle("Boxplot of junction RPKM") + 
   coord_cartesian(ylim = c(0, 20)))
ggsave(rpkm_box, file = "RPKM_boxplot.pdf")

################################################################################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# validate previously identified isoform exons with junction RPKM 
junction_valid <- matrix(NA, nrow = 7, ncol = 6, dimnames = list(c("lum084_myo084", "lum084_stem084", "myo084_stem084", "lum080_myo080", "lum080_stem080", "myo080_stem080", "lum035_myo035"), c("No.isoform.exons", "No.isoform.genes", "No.exons.with.junction.cov", "No.genes.with.junction.cov", "No.exons.with.junction.support", "No.genes.with.junction.support")))

# DE genes FDR = 0.015
setwd("~/REMC/breast/gene/")
lum084_myo084_DEup <- read.delim("UP.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DEdn <- read.delim("DN.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DE <- c(lum084_myo084_DEup$V1, lum084_myo084_DEdn$V1)
lum084_stem084_DEup <- read.delim("UP.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DEdn <- read.delim("DN.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DE <- c(lum084_stem084_DEup$V1, lum084_stem084_DEdn$V1)
myo084_stem084_DEup <- read.delim("UP.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DEdn <- read.delim("DN.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DE <- c(myo084_stem084_DEup$V1, myo084_stem084_DEdn$V1)
rm(lum084_myo084_DEup, lum084_myo084_DEdn, lum084_stem084_DEup, lum084_stem084_DEdn, myo084_stem084_DEup, myo084_stem084_DEdn)
lum080_myo080_DEup <- read.delim("UP.lumRM080_myoRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum080_myo080_DEdn <- read.delim("DN.lumRM080_myoRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum080_myo080_DE <- c(lum080_myo080_DEup$V1, lum080_myo080_DEdn$V1)
lum080_stem080_DEup <- read.delim("UP.lumRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum080_stem080_DEdn <- read.delim("DN.lumRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum080_stem080_DE <- c(lum080_stem080_DEup$V1, lum080_stem080_DEdn$V1)
myo080_stem080_DEup <- read.delim("UP.myoRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo080_stem080_DEdn <- read.delim("DN.myoRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo080_stem080_DE <- c(myo080_stem080_DEup$V1, myo080_stem080_DEdn$V1)
rm(lum080_myo080_DEup, lum080_myo080_DEdn, lum080_stem080_DEup, lum080_stem080_DEdn, myo080_stem080_DEup, myo080_stem080_DEdn)
lum035_myo035_DEup <- read.delim("UP.lumRM035_myoRM035.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum035_myo035_DEdn <- read.delim("DN.lumRM035_myoRM035.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum035_myo035_DE <- c(lum035_myo035_DEup$V1, lum035_myo035_DEdn$V1)
rm(lum035_myo035_DEup, lum035_myo035_DEdn)

# previously identified isoform exons & junction RPKM 
setwd("~/REMC/junction/all/")
# load("junction_RM084_new.Rdata")
load("junction_all.Rdata")
load("~/hg19/hg19v65_genes.Rdata")
junction_exon <- read.delim("~/hg19/hg19v65_junctions_vs_exons_for_genes.relations.unique", head = T, as.is = T)
junction_exon$junctionID <- paste0("chr", junction_exon$junctionID)

# isoform exons with junction support: junction RPKM in one sample < 0.1(cov <= 1), > 0.1 in the other(cov >= 2) 
high_fold <- function(exon, up, sample1, sample2, junctions = junction, junction_exons = junction_exon){
  e <- 1e-6
  covcut <- 1 # cutoff for sum junction coverage of two samples
  junctions <- na.omit(junctions[junction_exons[junction_exons$exonID == exon, "junctionID"], ])
  junctions <- junctions[(junctions[, paste0(sample1, "N")] + junctions[, paste0(sample2, "N")]) > covcut,]
  if(nrow(junctions) == 0){
    return(c(NA, NA, NA)) # not enough coverage for junctionss of this gene
  }else{
    junctions$fold <- log2((junctions[, paste0(sample1, "rpkm")] + e) / (junctions[, paste0(sample2, "rpkm")] + e))
    if(up > 0){
      highest <- junctions[which.max(junctions$fold),]
    }else{
      highest <- junctions[which.min(junctions$fold),]
    }
    return(c(highest[, paste0(sample1, "rpkm")], highest[, paste0(sample2, "rpkm")], highest[, "fold"]))
  }
}
Vhigh_fold <- Vectorize(high_fold, vectorize.args = c("exon", "up"), SIMPLIFY = T)

logfoldcut <- 1 # cutoff for abs(log2(fold change junction rpkm))
jcut <- 0.1 # cutoff for junction RPKM, i.e.one sample > 0.1, the other < 0.1

load("RM084_isoform_all.Rdata")
# lum084 vs myo084
lum084_myo084_isoform_valid <- lum084_myo084_isoform_all
(junction_valid["lum084_myo084", "No.isoform.exons"] <- nrow(lum084_myo084_isoform_valid))
(junction_valid["lum084_myo084", "No.isoform.genes"] <- length(unique(lum084_myo084_isoform_valid$id)))
sum(is.element(lum084_myo084_isoform_valid$V1, junction_exon$exonID))
lum084_myo084_isoform_valid <- data.frame(exon = lum084_myo084_isoform_valid$V1, gene = lum084_myo084_isoform_valid$id, DEfine.p = lum084_myo084_isoform_valid$V5, 
                                          lum084gene = lum084_myo084_isoform_valid$ave3, myo084gene = lum084_myo084_isoform_valid$ave2, 
                                          lum084exon = lum084_myo084_isoform_valid$V3, myo084exon = lum084_myo084_isoform_valid$V2, 
                                          t(Vhigh_fold(exon = lum084_myo084_isoform_valid$V1, up = lum084_myo084_isoform_valid$V3 - lum084_myo084_isoform_valid$V2, sample1 = "lum084", sample2 = "myo084")))
colnames(lum084_myo084_isoform_valid)[(ncol(lum084_myo084_isoform_valid)-2):ncol(lum084_myo084_isoform_valid)] <- c("lum084junction", "myo084junction", "junction_logfold")
(junction_valid["lum084_myo084", "No.isoform.exons"] <- nrow(lum084_myo084_isoform_valid))
(junction_valid["lum084_myo084", "No.isoform.genes"] <- length(unique(lum084_myo084_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
lum084_myo084_isoform_valid <- na.omit(lum084_myo084_isoform_valid)
(junction_valid["lum084_myo084", "No.exons.with.junction.cov"] <- nrow(lum084_myo084_isoform_valid))
(junction_valid["lum084_myo084", "No.genes.with.junction.cov"] <- length(unique(lum084_myo084_isoform_valid$gene)))
# No. of exons/genes with junction support
#lum084_myo084_isoform_valid <- lum084_myo084_isoform_valid[abs(lum084_myo084_isoform_valid$junction_logfold) > logfoldcut,]
lum084_myo084_isoform_valid <- lum084_myo084_isoform_valid[(lum084_myo084_isoform_valid$lum084junction >= jcut & lum084_myo084_isoform_valid$myo084junction < jcut) | (lum084_myo084_isoform_valid$lum084junction < jcut & lum084_myo084_isoform_valid$myo084junction >= jcut),]
lum084_myo084_isoform_valid$junction_logfold <- NULL
lum084_myo084_isoform_valid <- lum084_myo084_isoform_valid[order((lum084_myo084_isoform_valid$lum084gene + lum084_myo084_isoform_valid$myo084gene), decreasing = T), ]
lum084_myo084_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum084_myo084_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["lum084_myo084", "No.exons.with.junction.support"] <- nrow(lum084_myo084_isoform_valid))
(junction_valid["lum084_myo084", "No.genes.with.junction.support"] <- length(unique(lum084_myo084_isoform_valid$gene)))
write.table(lum084_myo084_isoform_valid, file = "lum084_myo084_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# lum084 vs stem084
lum084_stem084_isoform_valid <- lum084_stem084_isoform_all
(junction_valid["lum084_stem084", "No.isoform.exons"] <- nrow(lum084_stem084_isoform_valid))
(junction_valid["lum084_stem084", "No.isoform.genes"] <- length(unique(lum084_stem084_isoform_valid$id)))
sum(is.element(lum084_stem084_isoform_valid$V1, junction_exon$exonID))
lum084_stem084_isoform_valid <- data.frame(exon = lum084_stem084_isoform_valid$V1, gene = lum084_stem084_isoform_valid$id, DEfine.p = lum084_stem084_isoform_valid$V5, 
                                           lum084gene = lum084_stem084_isoform_valid$ave2, stem084gene = lum084_stem084_isoform_valid$ave3, 
                                           lum084exon = lum084_stem084_isoform_valid$V2, stem084exon = lum084_stem084_isoform_valid$V3, 
                                           t(Vhigh_fold(exon = lum084_stem084_isoform_valid$V1, up = lum084_stem084_isoform_valid$V2 - lum084_stem084_isoform_valid$V3, sample1 = "lum084", sample2 = "stem084")))
colnames(lum084_stem084_isoform_valid)[(ncol(lum084_stem084_isoform_valid)-2):ncol(lum084_stem084_isoform_valid)] <- c("lum084junction", "stem084junction", "junction_logfold")
(junction_valid["lum084_stem084", "No.isoform.exons"] <- nrow(lum084_stem084_isoform_valid))
(junction_valid["lum084_stem084", "No.isoform.genes"] <- length(unique(lum084_stem084_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
lum084_stem084_isoform_valid <- na.omit(lum084_stem084_isoform_valid)
(junction_valid["lum084_stem084", "No.exons.with.junction.cov"] <- nrow(lum084_stem084_isoform_valid))
(junction_valid["lum084_stem084", "No.genes.with.junction.cov"] <- length(unique(lum084_stem084_isoform_valid$gene)))
# No. of exons/genes with junction support
#lum084_stem084_isoform_valid <- lum084_stem084_isoform_valid[abs(lum084_stem084_isoform_valid$junction_logfold) > logfoldcut,]
lum084_stem084_isoform_valid <- lum084_stem084_isoform_valid[(lum084_stem084_isoform_valid$lum084junction >= jcut & lum084_stem084_isoform_valid$stem084junction < jcut) | (lum084_stem084_isoform_valid$lum084junction < jcut & lum084_stem084_isoform_valid$stem084junction >= jcut),]
lum084_stem084_isoform_valid$junction_logfold <- NULL
lum084_stem084_isoform_valid <- lum084_stem084_isoform_valid[order((lum084_stem084_isoform_valid$lum084gene + lum084_stem084_isoform_valid$stem084gene), decreasing = T), ]
lum084_stem084_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum084_stem084_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["lum084_stem084", "No.exons.with.junction.support"] <- nrow(lum084_stem084_isoform_valid))
(junction_valid["lum084_stem084", "No.genes.with.junction.support"] <- length(unique(lum084_stem084_isoform_valid$gene)))
write.table(lum084_stem084_isoform_valid, file = "lum084_stem084_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# myo084 vs stem084
myo084_stem084_isoform_valid <- myo084_stem084_isoform_all
(junction_valid["myo084_stem084", "No.isoform.exons"] <- nrow(myo084_stem084_isoform_valid))
(junction_valid["myo084_stem084", "No.isoform.genes"] <- length(unique(myo084_stem084_isoform_valid$id)))
sum(is.element(myo084_stem084_isoform_valid$V1, junction_exon$exonID))
myo084_stem084_isoform_valid <- data.frame(exon = myo084_stem084_isoform_valid$V1, gene = myo084_stem084_isoform_valid$id, DEfine.p = myo084_stem084_isoform_valid$V5, 
                                           myo084gene = myo084_stem084_isoform_valid$ave2, stem084gene = myo084_stem084_isoform_valid$ave3, 
                                           myo084exon = myo084_stem084_isoform_valid$V2, stem084exon = myo084_stem084_isoform_valid$V3, 
                                           t(Vhigh_fold(exon = myo084_stem084_isoform_valid$V1, up = myo084_stem084_isoform_valid$V2 - myo084_stem084_isoform_valid$V3, sample1 = "myo084", sample2 = "stem084")))
colnames(myo084_stem084_isoform_valid)[(ncol(myo084_stem084_isoform_valid)-2):ncol(myo084_stem084_isoform_valid)] <- c("myo084junction", "stem084junction", "junction_logfold")
(junction_valid["myo084_stem084", "No.isoform.exons"] <- nrow(myo084_stem084_isoform_valid))
(junction_valid["myo084_stem084", "No.isoform.genes"] <- length(unique(myo084_stem084_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
myo084_stem084_isoform_valid <- na.omit(myo084_stem084_isoform_valid)
(junction_valid["myo084_stem084", "No.exons.with.junction.cov"] <- nrow(myo084_stem084_isoform_valid))
(junction_valid["myo084_stem084", "No.genes.with.junction.cov"] <- length(unique(myo084_stem084_isoform_valid$gene)))
# No. of exons/genes with junction support
#myo084_stem084_isoform_valid <- myo084_stem084_isoform_valid[abs(myo084_stem084_isoform_valid$junction_logfold) > logfoldcut,]
myo084_stem084_isoform_valid <- myo084_stem084_isoform_valid[(myo084_stem084_isoform_valid$myo084junction >= jcut & myo084_stem084_isoform_valid$stem084junction < jcut) | (myo084_stem084_isoform_valid$myo084junction < jcut & myo084_stem084_isoform_valid$stem084junction >= jcut),]
myo084_stem084_isoform_valid$junction_logfold <- NULL
myo084_stem084_isoform_valid <- myo084_stem084_isoform_valid[order((myo084_stem084_isoform_valid$myo084gene + myo084_stem084_isoform_valid$stem084gene), decreasing = T), ]
myo084_stem084_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(myo084_stem084_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["myo084_stem084", "No.exons.with.junction.support"] <- nrow(myo084_stem084_isoform_valid))
(junction_valid["myo084_stem084", "No.genes.with.junction.support"] <- length(unique(myo084_stem084_isoform_valid$gene)))
write.table(myo084_stem084_isoform_valid, file = "myo084_stem084_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
# save.image(file = "junction_valid_RM084_new.Rdata")

# lum080 vs myo080
lum080_myo080up_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84up_isoform.csv", head = T, as.is = T)
lum080_myo080dn_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84dn_isoform.csv", head = T, as.is = T)
lum080_myo080_isoform_all <- rbind(lum080_myo080up_isoform_valid, lum080_myo080dn_isoform_valid)
lum080_myo080_isoform_all <- lum080_myo080_isoform_all[!is.element(lum080_myo080_isoform_all$id, lum080_myo080_DE),]
rm(lum080_myo080up_isoform_valid, lum080_myo080dn_isoform_valid)
lum080_myo080_isoform_valid <- lum080_myo080_isoform_all
(junction_valid["lum080_myo080", "No.isoform.exons"] <- nrow(lum080_myo080_isoform_valid))
(junction_valid["lum080_myo080", "No.isoform.genes"] <- length(unique(lum080_myo080_isoform_valid$id)))
sum(is.element(lum080_myo080_isoform_valid$V1, junction_exon$exonID))
lum080_myo080_isoform_valid <- data.frame(exon = lum080_myo080_isoform_valid$V1, gene = lum080_myo080_isoform_valid$id, DEfine.p = lum080_myo080_isoform_valid$V5, 
                                          lum080gene = lum080_myo080_isoform_valid$ave3, myo080gene = lum080_myo080_isoform_valid$ave2, 
                                          lum080exon = lum080_myo080_isoform_valid$V3, myo080exon = lum080_myo080_isoform_valid$V2, 
                                          t(Vhigh_fold(exon = lum080_myo080_isoform_valid$V1, up = lum080_myo080_isoform_valid$V3 - lum080_myo080_isoform_valid$V2, sample1 = "lum080", sample2 = "myo080")))
colnames(lum080_myo080_isoform_valid)[(ncol(lum080_myo080_isoform_valid)-2):ncol(lum080_myo080_isoform_valid)] <- c("lum080junction", "myo080junction", "junction_logfold")
(junction_valid["lum080_myo080", "No.isoform.exons"] <- nrow(lum080_myo080_isoform_valid))
(junction_valid["lum080_myo080", "No.isoform.genes"] <- length(unique(lum080_myo080_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
lum080_myo080_isoform_valid <- na.omit(lum080_myo080_isoform_valid)
(junction_valid["lum080_myo080", "No.exons.with.junction.cov"] <- nrow(lum080_myo080_isoform_valid))
(junction_valid["lum080_myo080", "No.genes.with.junction.cov"] <- length(unique(lum080_myo080_isoform_valid$gene)))
# No. of exons/genes with junction support
#lum080_myo080_isoform_valid <- lum080_myo080_isoform_valid[abs(lum080_myo080_isoform_valid$junction_logfold) > logfoldcut,]
lum080_myo080_isoform_valid <- lum080_myo080_isoform_valid[(lum080_myo080_isoform_valid$lum080junction >= jcut & lum080_myo080_isoform_valid$myo080junction < jcut) | (lum080_myo080_isoform_valid$lum080junction < jcut & lum080_myo080_isoform_valid$myo080junction >= jcut),]
lum080_myo080_isoform_valid$junction_logfold <- NULL
lum080_myo080_isoform_valid <- lum080_myo080_isoform_valid[order((lum080_myo080_isoform_valid$lum080gene + lum080_myo080_isoform_valid$myo080gene), decreasing = T), ]
lum080_myo080_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum080_myo080_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["lum080_myo080", "No.exons.with.junction.support"] <- nrow(lum080_myo080_isoform_valid))
(junction_valid["lum080_myo080", "No.genes.with.junction.support"] <- length(unique(lum080_myo080_isoform_valid$gene)))
write.table(lum080_myo080_isoform_valid, file = "lum080_myo080_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# lum080 vs stem080
lum080_stem080up_isoform_valid <- read.csv("~/REMC/breast/tissue/lum_stem_84up_isoform.csv", head = T, as.is = T)
lum080_stem080dn_isoform_valid <- read.csv("~/REMC/breast/tissue/lum_stem_84dn_isoform.csv", head = T, as.is = T)
lum080_stem080_isoform_all <- rbind(lum080_stem080up_isoform_valid, lum080_stem080dn_isoform_valid)
lum080_stem080_isoform_all <- lum080_stem080_isoform_all[!is.element(lum080_stem080_isoform_all$id, lum080_stem080_DE),]
rm(lum080_stem080up_isoform_valid, lum080_stem080dn_isoform_valid)
lum080_stem080_isoform_valid <- lum080_stem080_isoform_all
(junction_valid["lum080_stem080", "No.isoform.exons"] <- nrow(lum080_stem080_isoform_valid))
(junction_valid["lum080_stem080", "No.isoform.genes"] <- length(unique(lum080_stem080_isoform_valid$id)))
sum(is.element(lum080_stem080_isoform_valid$V1, junction_exon$exonID))
lum080_stem080_isoform_valid <- data.frame(exon = lum080_stem080_isoform_valid$V1, gene = lum080_stem080_isoform_valid$id, DEfine.p = lum080_stem080_isoform_valid$V5, 
                                           lum080gene = lum080_stem080_isoform_valid$ave2, stem080gene = lum080_stem080_isoform_valid$ave3, 
                                           lum080exon = lum080_stem080_isoform_valid$V2, stem080exon = lum080_stem080_isoform_valid$V3, 
                                           t(Vhigh_fold(exon = lum080_stem080_isoform_valid$V1, up = lum080_stem080_isoform_valid$V2 - lum080_stem080_isoform_valid$V3, sample1 = "lum080", sample2 = "stem080")))
colnames(lum080_stem080_isoform_valid)[(ncol(lum080_stem080_isoform_valid)-2):ncol(lum080_stem080_isoform_valid)] <- c("lum080junction", "stem080junction", "junction_logfold")
(junction_valid["lum080_stem080", "No.isoform.exons"] <- nrow(lum080_stem080_isoform_valid))
(junction_valid["lum080_stem080", "No.isoform.genes"] <- length(unique(lum080_stem080_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
lum080_stem080_isoform_valid <- na.omit(lum080_stem080_isoform_valid)
(junction_valid["lum080_stem080", "No.exons.with.junction.cov"] <- nrow(lum080_stem080_isoform_valid))
(junction_valid["lum080_stem080", "No.genes.with.junction.cov"] <- length(unique(lum080_stem080_isoform_valid$gene)))
# No. of exons/genes with junction support
#lum080_stem080_isoform_valid <- lum080_stem080_isoform_valid[abs(lum080_stem080_isoform_valid$junction_logfold) > logfoldcut,]
lum080_stem080_isoform_valid <- lum080_stem080_isoform_valid[(lum080_stem080_isoform_valid$lum080junction >= jcut & lum080_stem080_isoform_valid$stem080junction < jcut) | (lum080_stem080_isoform_valid$lum080junction < jcut & lum080_stem080_isoform_valid$stem080junction >= jcut),]
lum080_stem080_isoform_valid$junction_logfold <- NULL
lum080_stem080_isoform_valid <- lum080_stem080_isoform_valid[order((lum080_stem080_isoform_valid$lum080gene + lum080_stem080_isoform_valid$stem080gene), decreasing = T), ]
lum080_stem080_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum080_stem080_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["lum080_stem080", "No.exons.with.junction.support"] <- nrow(lum080_stem080_isoform_valid))
(junction_valid["lum080_stem080", "No.genes.with.junction.support"] <- length(unique(lum080_stem080_isoform_valid$gene)))
write.table(lum080_stem080_isoform_valid, file = "lum080_stem080_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# myo080 vs stem080
myo080_stem080up_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_stem_84up_isoform.csv", head = T, as.is = T)
myo080_stem080dn_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_stem_84dn_isoform.csv", head = T, as.is = T)
myo080_stem080_isoform_all <- rbind(myo080_stem080up_isoform_valid, myo080_stem080dn_isoform_valid)
myo080_stem080_isoform_all <- myo080_stem080_isoform_all[!is.element(myo080_stem080_isoform_all$id, myo080_stem080_DE),]
rm(myo080_stem080up_isoform_valid, myo080_stem080dn_isoform_valid)
myo080_stem080_isoform_valid <- myo080_stem080_isoform_all
(junction_valid["myo080_stem080", "No.isoform.exons"] <- nrow(myo080_stem080_isoform_valid))
(junction_valid["myo080_stem080", "No.isoform.genes"] <- length(unique(myo080_stem080_isoform_valid$id)))
sum(is.element(myo080_stem080_isoform_valid$V1, junction_exon$exonID))
myo080_stem080_isoform_valid <- data.frame(exon = myo080_stem080_isoform_valid$V1, gene = myo080_stem080_isoform_valid$id, DEfine.p = myo080_stem080_isoform_valid$V5, 
                                           myo080gene = myo080_stem080_isoform_valid$ave2, stem080gene = myo080_stem080_isoform_valid$ave3, 
                                           myo080exon = myo080_stem080_isoform_valid$V2, stem080exon = myo080_stem080_isoform_valid$V3, 
                                           t(Vhigh_fold(exon = myo080_stem080_isoform_valid$V1, up = myo080_stem080_isoform_valid$V2 - myo080_stem080_isoform_valid$V3, sample1 = "myo080", sample2 = "stem080")))
colnames(myo080_stem080_isoform_valid)[(ncol(myo080_stem080_isoform_valid)-2):ncol(myo080_stem080_isoform_valid)] <- c("myo080junction", "stem080junction", "junction_logfold")
(junction_valid["myo080_stem080", "No.isoform.exons"] <- nrow(myo080_stem080_isoform_valid))
(junction_valid["myo080_stem080", "No.isoform.genes"] <- length(unique(myo080_stem080_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
myo080_stem080_isoform_valid <- na.omit(myo080_stem080_isoform_valid)
(junction_valid["myo080_stem080", "No.exons.with.junction.cov"] <- nrow(myo080_stem080_isoform_valid))
(junction_valid["myo080_stem080", "No.genes.with.junction.cov"] <- length(unique(myo080_stem080_isoform_valid$gene)))
# No. of exons/genes with junction support
#myo080_stem080_isoform_valid <- myo080_stem080_isoform_valid[abs(myo080_stem080_isoform_valid$junction_logfold) > logfoldcut,]
myo080_stem080_isoform_valid <- myo080_stem080_isoform_valid[(myo080_stem080_isoform_valid$myo080junction >= jcut & myo080_stem080_isoform_valid$stem080junction < jcut) | (myo080_stem080_isoform_valid$myo080junction < jcut & myo080_stem080_isoform_valid$stem080junction >= jcut),]
myo080_stem080_isoform_valid$junction_logfold <- NULL
myo080_stem080_isoform_valid <- myo080_stem080_isoform_valid[order((myo080_stem080_isoform_valid$myo080gene + myo080_stem080_isoform_valid$stem080gene), decreasing = T), ]
myo080_stem080_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(myo080_stem080_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["myo080_stem080", "No.exons.with.junction.support"] <- nrow(myo080_stem080_isoform_valid))
(junction_valid["myo080_stem080", "No.genes.with.junction.support"] <- length(unique(myo080_stem080_isoform_valid$gene)))
write.table(myo080_stem080_isoform_valid, file = "myo080_stem080_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# lum035 vs myo035
lum035_myo035up_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84up_isoform.csv", head = T, as.is = T)
lum035_myo035dn_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84dn_isoform.csv", head = T, as.is = T)
lum035_myo035_isoform_all <- rbind(lum035_myo035up_isoform_valid, lum035_myo035dn_isoform_valid)
lum035_myo035_isoform_all <- lum035_myo035_isoform_all[!is.element(lum035_myo035_isoform_all$id, lum035_myo035_DE),]
rm(lum035_myo035up_isoform_valid, lum035_myo035dn_isoform_valid)
lum035_myo035_isoform_valid <- lum035_myo035_isoform_all
(junction_valid["lum035_myo035", "No.isoform.exons"] <- nrow(lum035_myo035_isoform_valid))
(junction_valid["lum035_myo035", "No.isoform.genes"] <- length(unique(lum035_myo035_isoform_valid$id)))
sum(is.element(lum035_myo035_isoform_valid$V1, junction_exon$exonID))
lum035_myo035_isoform_valid <- data.frame(exon = lum035_myo035_isoform_valid$V1, gene = lum035_myo035_isoform_valid$id, DEfine.p = lum035_myo035_isoform_valid$V5, 
                                          lum035gene = lum035_myo035_isoform_valid$ave3, myo035gene = lum035_myo035_isoform_valid$ave2, 
                                          lum035exon = lum035_myo035_isoform_valid$V3, myo035exon = lum035_myo035_isoform_valid$V2, 
                                          t(Vhigh_fold(exon = lum035_myo035_isoform_valid$V1, up = lum035_myo035_isoform_valid$V3 - lum035_myo035_isoform_valid$V2, sample1 = "lum035", sample2 = "myo035")))
colnames(lum035_myo035_isoform_valid)[(ncol(lum035_myo035_isoform_valid)-2):ncol(lum035_myo035_isoform_valid)] <- c("lum035junction", "myo035junction", "junction_logfold")
(junction_valid["lum035_myo035", "No.isoform.exons"] <- nrow(lum035_myo035_isoform_valid))
(junction_valid["lum035_myo035", "No.isoform.genes"] <- length(unique(lum035_myo035_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
lum035_myo035_isoform_valid <- na.omit(lum035_myo035_isoform_valid)
(junction_valid["lum035_myo035", "No.exons.with.junction.cov"] <- nrow(lum035_myo035_isoform_valid))
(junction_valid["lum035_myo035", "No.genes.with.junction.cov"] <- length(unique(lum035_myo035_isoform_valid$gene)))
# No. of exons/genes with junction support
#lum035_myo035_isoform_valid <- lum035_myo035_isoform_valid[abs(lum035_myo035_isoform_valid$junction_logfold) > logfoldcut,]
lum035_myo035_isoform_valid <- lum035_myo035_isoform_valid[(lum035_myo035_isoform_valid$lum035junction >= jcut & lum035_myo035_isoform_valid$myo035junction < jcut) | (lum035_myo035_isoform_valid$lum035junction < jcut & lum035_myo035_isoform_valid$myo035junction >= jcut),]
lum035_myo035_isoform_valid$junction_logfold <- NULL
lum035_myo035_isoform_valid <- lum035_myo035_isoform_valid[order((lum035_myo035_isoform_valid$lum035gene + lum035_myo035_isoform_valid$myo035gene), decreasing = T), ]
lum035_myo035_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum035_myo035_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["lum035_myo035", "No.exons.with.junction.support"] <- nrow(lum035_myo035_isoform_valid))
(junction_valid["lum035_myo035", "No.genes.with.junction.support"] <- length(unique(lum035_myo035_isoform_valid$gene)))
write.table(lum035_myo035_isoform_valid, file = "lum035_myo035_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rm(junction_exon)

save.image(file = "junction_valid.Rdata")

###########################################################################################################
setwd("~/快盘/REMC/junction/all/")
load("junction_valid_RM084_new.Rdata")
# gene RPKM
lum084 <- read.delim("~/REMC/gene/A17918.G.A.rpkm.pc", head = F, as.is =T)
myo084 <- read.delim("~/REMC/gene/A17919.G.A.rpkm.pc", head = F, as.is =T)
stem084 <- read.delim("~/REMC/gene/A17920.G.A.rpkm.pc", head = F, as.is =T)
RM084gene <- data.frame(id = lum084$V1, lum084 = lum084$V3, myo084 = myo084$V3, stem084 = stem084$V3)
rm(lum084, myo084, stem084)

# No. of exons for DE genes / isoform genes: use collapsed exons instead of Ensembl exons
exon <- read.delim("~/hg19/hg19v65_exons_collapsed.txt", head = F, as.is = T)
Nexon <- data.frame(id = levels(as.factor(exon$V2)), Nexon = sapply(levels(as.factor(exon$V2)), function(x) sum(exon$V2 == x)))
lum084_myo084_DE <- data.frame(id = lum084_myo084_DE, Nexon = Nexon[lum084_myo084_DE, ]$Nexon)
lum084_myo084_isoform_all_gene <- lum084_myo084_isoform_all[!duplicated(lum084_myo084_isoform_all$id), ]
lum084_myo084_isoform_all_gene$Nexon <- Nexon[lum084_myo084_isoform_all_gene$id, "Nexon"]
lum084_myo084_isoform_valid_gene$Nexon <- Nexon[as.character(lum084_myo084_isoform_valid_gene$gene), "Nexon"] 
lum084_stem084_DE <- data.frame(id = lum084_stem084_DE, Nexon = Nexon[lum084_stem084_DE, ]$Nexon)
lum084_stem084_isoform_all_gene <- lum084_stem084_isoform_all[!duplicated(lum084_stem084_isoform_all$id), ]
lum084_stem084_isoform_all_gene$Nexon <- Nexon[lum084_stem084_isoform_all_gene$id, "Nexon"]
lum084_stem084_isoform_valid_gene$Nexon <- Nexon[as.character(lum084_stem084_isoform_valid_gene$gene), "Nexon"] 
myo084_stem084_DE <- data.frame(id = myo084_stem084_DE, Nexon = Nexon[myo084_stem084_DE, ]$Nexon)
myo084_stem084_isoform_all_gene <- myo084_stem084_isoform_all[!duplicated(myo084_stem084_isoform_all$id), ]
myo084_stem084_isoform_all_gene$Nexon <- Nexon[myo084_stem084_isoform_all_gene$id, "Nexon"]
myo084_stem084_isoform_valid_gene$Nexon <- Nexon[as.character(myo084_stem084_isoform_valid_gene$gene), "Nexon"] 
rm(exon)
pdf("Nexon_DEvsIsoforms.pdf")
plot(c(0, 50), c(0, 0.08), type = "n", main = "No. of exons in DE genes and isoforms", xlab = "No. of exons", ylab = "density")
lines(density(na.omit(Nexon[Nexon$id %in% RM084gene[(RM084gene$lum084 + RM084gene$myo084 + RM084gene$stem084) > 0.005,]$id,]$Nexon)), col = 1, lty = 1, lwd = 5)
lines(density(na.omit(lum084_myo084_DE$Nexon)), col = 2, lty = 1, lwd = 3)
lines(density(na.omit(lum084_stem084_DE$Nexon)), col = 3, lty = 1, lwd = 3)
lines(density(na.omit(myo084_stem084_DE$Nexon)), col = 4, lty = 1, lwd = 3)
lines(density(na.omit(lum084_myo084_isoform_all_gene$Nexon)), col = 2, lty = 3, lwd = 3)
lines(density(na.omit(lum084_stem084_isoform_all_gene$Nexon)), col = 3, lty = 3, lwd = 3)
lines(density(na.omit(myo084_stem084_isoform_all_gene$Nexon)), col = 4, lty = 3, lwd = 3)
# lines(density(na.omit(lum084_myo084_isoform_valid_gene$Nexon)), col = 2, lty = 4, lwd = 3)
# lines(density(na.omit(lum084_stem084_isoform_valid_gene$Nexon)), col = 3, lty = 4, lwd = 3)
# lines(density(na.omit(myo084_stem084_isoform_valid_gene$Nexon)), col = 4, lty = 4, lwd = 3)
# legend("bottomright", c("DE genes", "all isoforms", "validated isoforms", "all genes"), col = 1, lty = c(1, 3, 4, 1), lwd = 5, cex = 0.8)
legend("topleft", c("all expressed genes", "DE genes", "isoforms"), col = 1, lty = c(1, 1, 3), lwd = 5, cex = 0.8)
legend("topright", c("all expressed genes", "lum084 vs myo084", "lum084 vs stem084", "myo084 vs stem084"), col = 1:4, lty = 1, lwd = 5, cex = 0.8)
dev.off()

# Isoform genes in DE genes with small No. of exons
# all DE exons from DEfine
lum084_myo084_exon <- read.delim("~/REMC/breast/exon/myo-RM084_lum-RM084.RPKM.corrected", head = F, as.is = T)
lum084_stem084_exon <- read.delim("~/REMC/breast/exon/lum-RM084_stem-RM084.RPKM.corrected", head = F, as.is = T)
myo084_stem084_exon <- read.delim("~/REMC/breast/exon/myo-RM084_stem-RM084.RPKM.corrected", head = F, as.is = T)
# DE genes with <= 5 exons
lum084_myo084_DE5 <- lum084_myo084_DE[lum084_myo084_DE$Nexon <= 5, ]
lum084_stem084_DE5 <- lum084_stem084_DE[lum084_stem084_DE$Nexon <= 5, ]
myo084_stem084_DE5 <- myo084_stem084_DE[myo084_stem084_DE$Nexon <= 5, ]
# No. of DE exons for each gene 
DEexons <- function(id, exons){
  cutoff <- 1
  e <- 1e-6
  id_exon <- exons[grep(id, exons$V1), ]
  id_exon$logfold <- log2((id_exon$V3 + e) / (id_exon$V5 + e))
  return(sum(abs(id_exon$logfold) >= cutoff))
}
vDEexons <- Vectorize(DEexons, vectorize.args = "id")
lum084_myo084_DE5$DEexons <- vDEexons(as.character(lum084_myo084_DE5$id), lum084_myo084_exon)
lum084_stem084_DE5$DEexons <- vDEexons(as.character(lum084_stem084_DE5$id), lum084_stem084_exon)
myo084_stem084_DE5$DEexons <- vDEexons(as.character(myo084_stem084_DE5$id), myo084_stem084_exon)
rm(lum084_myo084_exon, lum084_stem084_exon, myo084_stem084_exon)
pdf("PercentDEexons.pdf")
plot(c(0, 1), c(0, 30), type = "n", main = "DE exons in DE genes", xlab = "Proportion of DE exons in DE genes", ylab = "density")
lines(density(lum084_myo084_DE5$DEexons / lum084_myo084_DE5$Nexon), col = 2, lty = 1, lwd = 3)
lines(density(lum084_stem084_DE5$DEexons / lum084_stem084_DE5$Nexon), col = 3, lty = 1, lwd = 3)
lines(density(myo084_stem084_DE5$DEexons / myo084_stem084_DE5$Nexon), col = 4, lty = 1, lwd = 3)
legend("topleft", c("lum084 vs myo084", "lum084 vs stem084", "myo084 vs stem084"), col = 2:4, lty = 1, lwd = 5, cex = 0.8)
dev.off()

# Position of isoform exons on the gene
# lum 084 vs myo 084
# start / end point of the gene 
lum084_myo084_isoform_all$start <- ensembl[lum084_myo084_isoform_all$id, "start"]
lum084_myo084_isoform_all$end <- ensembl[lum084_myo084_isoform_all$id, "end"]
lum084_myo084_isoform_all$strand <- ensembl[lum084_myo084_isoform_all$id, "strand"]
# exon position: mid point of the exon
lum084_myo084_isoform_all$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", lum084_myo084_isoform_all$V1))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", lum084_myo084_isoform_all$V1))))/2
# exon postition on the gene: calculate + strand & - strand separately
lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == 1, ]$exon_pos <- (lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == 1, ]$exon_pos - lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == 1, ]$start) / (lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == 1, ]$end - lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == 1, ]$start) * 100
lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == -1, ]$exon_pos <- 100 - (lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == -1, ]$exon_pos - lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == -1, ]$start) / (lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == -1, ]$end - lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == -1, ]$start) * 100
lum084_myo084_isoform_valid$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", lum084_myo084_isoform_valid$exon))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", lum084_myo084_isoform_valid$exon))))/2
lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == 1, ]$exon_pos <- (lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == 1, ]$exon_pos - lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == 1, ]$start) / (lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == 1, ]$end - lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == 1, ]$start) * 100
lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == -1, ]$exon_pos <- 100 - (lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == -1, ]$exon_pos - lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == -1, ]$start) / (lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == -1, ]$end - lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == -1, ]$start) * 100
# lum 084 vs stem 084
lum084_stem084_isoform_all$start <- ensembl[lum084_stem084_isoform_all$id, "start"]
lum084_stem084_isoform_all$end <- ensembl[lum084_stem084_isoform_all$id, "end"]
lum084_stem084_isoform_all$strand <- ensembl[lum084_stem084_isoform_all$id, "strand"]
lum084_stem084_isoform_all$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", lum084_stem084_isoform_all$V1))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", lum084_stem084_isoform_all$V1))))/2
lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == 1, ]$exon_pos <- (lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == 1, ]$exon_pos - lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == 1, ]$start) / (lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == 1, ]$end - lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == 1, ]$start) * 100
lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == -1, ]$exon_pos <- 100 - (lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == -1, ]$exon_pos - lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == -1, ]$start) / (lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == -1, ]$end - lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == -1, ]$start) * 100
lum084_stem084_isoform_valid$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", lum084_stem084_isoform_valid$exon))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", lum084_stem084_isoform_valid$exon))))/2
lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == 1, ]$exon_pos <- (lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == 1, ]$exon_pos - lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == 1, ]$start) / (lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == 1, ]$end - lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == 1, ]$start) * 100
lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == -1, ]$exon_pos <- 100 - (lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == -1, ]$exon_pos - lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == -1, ]$start) / (lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == -1, ]$end - lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == -1, ]$start) * 100
# myo 084 vs stem 084
myo084_stem084_isoform_all$start <- ensembl[myo084_stem084_isoform_all$id, "start"]
myo084_stem084_isoform_all$end <- ensembl[myo084_stem084_isoform_all$id, "end"]
myo084_stem084_isoform_all$strand <- ensembl[myo084_stem084_isoform_all$id, "strand"]
myo084_stem084_isoform_all$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", myo084_stem084_isoform_all$V1))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", myo084_stem084_isoform_all$V1))))/2
myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == 1, ]$exon_pos <- (myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == 1, ]$exon_pos - myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == 1, ]$start) / (myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == 1, ]$end - myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == 1, ]$start) * 100
myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == -1, ]$exon_pos <- 100 - (myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == -1, ]$exon_pos - myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == -1, ]$start) / (myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == -1, ]$end - myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == -1, ]$start) * 100
myo084_stem084_isoform_valid$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", myo084_stem084_isoform_valid$exon))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", myo084_stem084_isoform_valid$exon))))/2
myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == 1, ]$exon_pos <- (myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == 1, ]$exon_pos - myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == 1, ]$start) / (myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == 1, ]$end - myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == 1, ]$start) * 100
myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == -1, ]$exon_pos <- 100 - (myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == -1, ]$exon_pos - myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == -1, ]$start) / (myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == -1, ]$end - myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == -1, ]$start) * 100
pdf("exon_position.pdf")
plot(c(0, 100), c(0, 0.015), type = "n", main = "Distribution of isoform exons along genes", xlab = "gene position", ylab = "density")
lines(density(na.omit(lum084_myo084_isoform_all$exon_pos)), col = 1, lty = 1, lwd = 3)
lines(density(na.omit(lum084_stem084_isoform_all$exon_pos)), col = 2, lty = 1, lwd = 3)
lines(density(na.omit(myo084_stem084_isoform_all$exon_pos)), col = 3, lty = 1, lwd = 3)
# lines(density(na.omit(lum084_myo084_isoform_valid$exon_pos)), col = 1, lty = 3, lwd = 3)
# lines(density(na.omit(lum084_stem084_isoform_valid$exon_pos)), col = 2, lty = 3, lwd = 3)
# lines(density(na.omit(myo084_stem084_isoform_valid$exon_pos)), col = 3, lty = 3, lwd = 3)
# legend("topleft", c("all", "validated"), col = 1, lty = c(1, 3), lwd = 5, cex = 0.8)
legend("bottomright", c("lum084 vs myo084", "lum084 vs stem084", "myo084 vs stem084"), col = 1:3, lty = 1, lwd = 5)
dev.off()

# Venn Diagram with average expression level,  No. of exons, and average exon length              
# all isoforms
isoform_all <- list(lum084_myo084 = lum084_myo084_isoform_all_gene$id, lum084_stem084 = lum084_stem084_isoform_all_gene$id, myo084_stem084 = myo084_stem084_isoform_all_gene$id)
# gene list of different sections of Venn diagram (lm: lum084_myo084, ls: lum084_stem084, ms: myo084_stem084)
all_lm_only <- RM084gene[RM084gene$id %in% setdiff(lum084_myo084_isoform_all_gene$id, union(lum084_stem084_isoform_all_gene$id, myo084_stem084_isoform_all_gene$id)), ]
all_ls_only <- RM084gene[RM084gene$id %in% setdiff(lum084_stem084_isoform_all_gene$id, union(lum084_myo084_isoform_all_gene$id, myo084_stem084_isoform_all_gene$id)), ]
all_ms_only <- RM084gene[RM084gene$id %in% setdiff(myo084_stem084_isoform_all_gene$id, union(lum084_myo084_isoform_all_gene$id, lum084_stem084_isoform_all_gene$id)), ]
all_lm_ls_not_ms <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_myo084_isoform_all_gene$id, lum084_stem084_isoform_all_gene$id), myo084_stem084_isoform_all_gene$id), ]
all_lm_ms_not_ls <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_myo084_isoform_all_gene$id, myo084_stem084_isoform_all_gene$id), lum084_stem084_isoform_all_gene$id), ]
all_ls_ms_not_lm <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_stem084_isoform_all_gene$id, myo084_stem084_isoform_all_gene$id), lum084_myo084_isoform_all_gene$id), ]
all_lm_ls_ms <- RM084gene[RM084gene$id %in% intersect(intersect(lum084_myo084_isoform_all_gene$id, lum084_stem084_isoform_all_gene$id), myo084_stem084_isoform_all_gene$id), ]
all_lm_only$Nexon <- Nexon[all_lm_only$id, ]$Nexon
all_ls_only$Nexon <- Nexon[all_ls_only$id, ]$Nexon
all_ms_only$Nexon <- Nexon[all_ms_only$id, ]$Nexon
all_lm_ls_not_ms$Nexon <- Nexon[all_lm_ls_not_ms$id, ]$Nexon
all_lm_ms_not_ls$Nexon <- Nexon[all_lm_ms_not_ls$id, ]$Nexon
all_ls_ms_not_lm$Nexon <- Nexon[all_ls_ms_not_lm$id, ]$Nexon
all_lm_ls_ms$Nexon <- Nexon[all_lm_ls_ms$id, ]$Nexon
write.table(all_lm_only, file = "./Venn/all_lm_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_ls_only, file = "./Venn/all_ls_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_ms_only, file = "./Venn/all_ms_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_lm_ls_not_ms, file = "./Venn/all_lm_ls_not_ms.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_lm_ms_not_ls, file = "./Venn/all_lm_ms_not_ls.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_ls_ms_not_lm, file = "./Venn/all_ls_ms_not_lm.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_lm_ls_ms, file = "./Venn/all_lm_ls_ms.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# validated isoforms
isoform_valid <- list(lum084_myo084 = as.character(lum084_myo084_isoform_valid_gene$gene), lum084_stem084 = as.character(lum084_stem084_isoform_valid_gene$gene), myo084_stem084 = as.character(myo084_stem084_isoform_valid_gene$gene))
# gene list of different sections of Venn diagram (lm: lum084_myo084, ls: lum084_stem084, ms: myo084_stem084)
valid_lm_only <- RM084gene[RM084gene$id %in% setdiff(lum084_myo084_isoform_valid_gene$gene, union(lum084_stem084_isoform_valid_gene$gene, myo084_stem084_isoform_valid_gene$gene)), ]
valid_ls_only <- RM084gene[RM084gene$id %in% setdiff(lum084_stem084_isoform_valid_gene$gene, union(lum084_myo084_isoform_valid_gene$gene, myo084_stem084_isoform_valid_gene$gene)), ]
valid_ms_only <- RM084gene[RM084gene$id %in% setdiff(myo084_stem084_isoform_valid_gene$gene, union(lum084_myo084_isoform_valid_gene$gene, lum084_stem084_isoform_valid_gene$gene)), ]
valid_lm_ls_not_ms <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_myo084_isoform_valid_gene$gene, lum084_stem084_isoform_valid_gene$gene), myo084_stem084_isoform_valid_gene$gene), ]
valid_lm_ms_not_ls <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_myo084_isoform_valid_gene$gene, myo084_stem084_isoform_valid_gene$gene), lum084_stem084_isoform_valid_gene$gene), ]
valid_ls_ms_not_lm <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_stem084_isoform_valid_gene$gene, myo084_stem084_isoform_valid_gene$gene), lum084_myo084_isoform_valid_gene$gene), ]
valid_lm_ls_ms <- RM084gene[RM084gene$id %in% intersect(intersect(lum084_myo084_isoform_valid_gene$gene, lum084_stem084_isoform_valid_gene$gene), myo084_stem084_isoform_valid_gene$gene), ]
valid_lm_only$Nexon <- Nexon[valid_lm_only$id, ]$Nexon
valid_ls_only$Nexon <- Nexon[valid_ls_only$id, ]$Nexon
valid_ms_only$Nexon <- Nexon[valid_ms_only$id, ]$Nexon
valid_lm_ls_not_ms$Nexon <- Nexon[valid_lm_ls_not_ms$id, ]$Nexon
valid_lm_ms_not_ls$Nexon <- Nexon[valid_lm_ms_not_ls$id, ]$Nexon
valid_ls_ms_not_lm$Nexon <- Nexon[valid_ls_ms_not_lm$id, ]$Nexon
valid_lm_ls_ms$Nexon <- Nexon[valid_lm_ls_ms$id, ]$Nexon
write.table(valid_lm_only, file = "./Venn/valid_lm_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_ls_only, file = "./Venn/valid_ls_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_ms_only, file = "./Venn/valid_ms_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_lm_ls_not_ms, file = "./Venn/valid_lm_ls_not_ms.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_lm_ms_not_ls, file = "./Venn/valid_lm_ms_not_ls.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_ls_ms_not_lm, file = "./Venn/valid_ls_ms_not_lm.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_lm_ls_ms, file = "./Venn/valid_lm_ls_ms.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# average exon length for isoforms on Venn diagram
lum084_myo084_isoform_all$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", lum084_myo084_isoform_all$V1)))
lum084_myo084_isoform_all$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", lum084_myo084_isoform_all$V1)))
lum084_myo084_isoform_all$exon_length <- abs(lum084_myo084_isoform_all$exon_end - lum084_myo084_isoform_all$exon_start)
lum084_stem084_isoform_all$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", lum084_stem084_isoform_all$V1)))
lum084_stem084_isoform_all$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", lum084_stem084_isoform_all$V1)))
lum084_stem084_isoform_all$exon_length <- abs(lum084_stem084_isoform_all$exon_end - lum084_stem084_isoform_all$exon_start)
myo084_stem084_isoform_all$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", myo084_stem084_isoform_all$V1)))
myo084_stem084_isoform_all$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", myo084_stem084_isoform_all$V1)))
myo084_stem084_isoform_all$exon_length <- abs(myo084_stem084_isoform_all$exon_end - myo084_stem084_isoform_all$exon_start)
all_lm_only$exon_length <- sapply(all_lm_only$id, function(x) mean(lum084_myo084_isoform_all[lum084_myo084_isoform_all$id == x, "exon_length"]))
all_ls_only$exon_length <- sapply(all_ls_only$id, function(x) mean(lum084_stem084_isoform_all[lum084_stem084_isoform_all$id == x, "exon_length"]))
all_ms_only$exon_length <- sapply(all_ms_only$id, function(x) mean(myo084_stem084_isoform_all[myo084_stem084_isoform_all$id == x, "exon_length"]))
all_lm_ls_not_ms$exon_length <- sapply(all_lm_ls_not_ms$id, function(x) mean(c(lum084_myo084_isoform_all[lum084_myo084_isoform_all$id == x, "exon_length"], lum084_stem084_isoform_all[lum084_stem084_isoform_all$id == x, "exon_length"])))
all_lm_ms_not_ls$exon_length <- sapply(all_lm_ms_not_ls$id, function(x) mean(c(lum084_myo084_isoform_all[lum084_myo084_isoform_all$id == x, "exon_length"], myo084_stem084_isoform_all[myo084_stem084_isoform_all$id == x, "exon_length"])))
all_ls_ms_not_lm$exon_length <- sapply(all_ls_ms_not_lm$id, function(x) mean(c(lum084_stem084_isoform_all[lum084_stem084_isoform_all$id == x, "exon_length"], myo084_stem084_isoform_all[myo084_stem084_isoform_all$id == x, "exon_length"])))
all_lm_ls_ms$exon_length <- sapply(all_lm_ls_ms$id, function(x) mean(c(lum084_myo084_isoform_all[lum084_myo084_isoform_all$id == x, "exon_length"], lum084_stem084_isoform_all[lum084_stem084_isoform_all$id == x, "exon_length"], myo084_stem084_isoform_all[myo084_stem084_isoform_all$id == x, "exon_length"])))

lum084_myo084_isoform_valid$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", lum084_myo084_isoform_valid$exon)))
lum084_myo084_isoform_valid$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", lum084_myo084_isoform_valid$exon)))
lum084_myo084_isoform_valid$exon_length <- abs(lum084_myo084_isoform_valid$exon_end - lum084_myo084_isoform_valid$exon_start)
lum084_stem084_isoform_valid$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", lum084_stem084_isoform_valid$exon)))
lum084_stem084_isoform_valid$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", lum084_stem084_isoform_valid$exon)))
lum084_stem084_isoform_valid$exon_length <- abs(lum084_stem084_isoform_valid$exon_end - lum084_stem084_isoform_valid$exon_start)
myo084_stem084_isoform_valid$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", myo084_stem084_isoform_valid$exon)))
myo084_stem084_isoform_valid$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", myo084_stem084_isoform_valid$exon)))
myo084_stem084_isoform_valid$exon_length <- abs(myo084_stem084_isoform_valid$exon_end - myo084_stem084_isoform_valid$exon_start)
valid_lm_only$exon_length <- sapply(valid_lm_only$id, function(x) mean(lum084_myo084_isoform_valid[as.character(lum084_myo084_isoform_valid$gene) == x, "exon_length"]))
valid_ls_only$exon_length <- sapply(valid_ls_only$id, function(x) mean(lum084_stem084_isoform_valid[as.character(lum084_stem084_isoform_valid$gene) == x, "exon_length"]))
valid_ms_only$exon_length <- sapply(valid_ms_only$id, function(x) mean(myo084_stem084_isoform_valid[as.character(myo084_stem084_isoform_valid$gene) == x, "exon_length"]))
valid_lm_ls_not_ms$exon_length <- sapply(valid_lm_ls_not_ms$id, function(x) mean(c(lum084_myo084_isoform_valid[as.character(lum084_myo084_isoform_valid$gene) == x, "exon_length"], lum084_stem084_isoform_valid[as.character(lum084_stem084_isoform_valid$gene) == x, "exon_length"])))
valid_lm_ms_not_ls$exon_length <- sapply(valid_lm_ms_not_ls$id, function(x) mean(c(lum084_myo084_isoform_valid[as.character(lum084_myo084_isoform_valid$gene) == x, "exon_length"], myo084_stem084_isoform_valid[as.character(myo084_stem084_isoform_valid$gene) == x, "exon_length"])))
valid_ls_ms_not_lm$exon_length <- sapply(valid_ls_ms_not_lm$id, function(x) mean(c(lum084_stem084_isoform_valid[as.character(lum084_stem084_isoform_valid$gene) == x, "exon_length"], myo084_stem084_isoform_valid[as.character(myo084_stem084_isoform_valid$gene) == x, "exon_length"])))
valid_lm_ls_ms$exon_length <- sapply(valid_lm_ls_ms$id, function(x) mean(c(lum084_myo084_isoform_valid[as.character(lum084_myo084_isoform_valid$gene) == x, "exon_length"], lum084_stem084_isoform_valid[as.character(lum084_stem084_isoform_valid$gene) == x, "exon_length"], myo084_stem084_isoform_valid[as.character(myo084_stem084_isoform_valid$gene) == x, "exon_length"])))

write.table(lum084_myo084_isoform_all_gene, file = "./Venn/lum084_myo084_isoform_all_gene.txt", sep = "\t", col.names = F, row.names = F)
write.table(lum084_stem084_isoform_all_gene, file = "./Venn/lum084_stem084_isoform_all_gene.txt", sep = "\t", col.names = F, row.names = F)
write.table(myo084_stem084_isoform_all_gene, file = "./Venn/myo084_stem084_isoform_all_gene.txt", sep = "\t", col.names = F, row.names = F)
write.table(lum084_myo084_isoform_valid_gene, file = "./Venn/lum084_myo084_isoform_valid_gene.txt", sep = "\t", col.names = F, row.names = F)
write.table(lum084_stem084_isoform_valid_gene, file = "./Venn/lum084_stem084_isoform_valid_gene.txt", sep = "\t", col.names = F, row.names = F)
write.table(myo084_stem084_isoform_valid_gene, file = "./Venn/myo084_stem084_isoform_valid_gene.txt", sep = "\t", col.names = F, row.names = F)

save.image(file = "junction_valid_RM084.Rdata")

library(VennDiagram)
isoform_all_N <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5), N = c(nrow(all_lm_only), nrow(all_ls_only), nrow(all_ms_only), nrow(all_lm_ls_not_ms), nrow(all_ls_ms_not_lm), nrow(all_lm_ms_not_ls), nrow(all_lm_ls_ms))) 
isoform_valid_N <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5), N = c(nrow(valid_lm_only), nrow(valid_ls_only), nrow(valid_ms_only), nrow(valid_lm_ls_not_ms), nrow(valid_ls_ms_not_lm), nrow(valid_lm_ms_not_ls), nrow(valid_lm_ls_ms))) 
pdf("venn_isoform.pdf")
(bubble_all_N <- ggplot(isoform_all_N) + 
   geom_point(aes(x, y, size = N), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(200, 1500), range = c(5, 20)) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_all <- venn.diagram(isoform_all, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "No. of isoform genes", main.cex = 1.5)
grid.draw(venn_all)
(bubble_valid_N <- ggplot(isoform_valid_N) + 
   geom_point(aes(x, y, size = N), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(200, 1500), range = c(5, 20)) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_valid <- venn.diagram(isoform_valid, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "No. of validated isoform genes", main.cex = 1.5)
grid.draw(venn_valid)
dev.off()

# Venn diagram with average gene rpkm
isoform_all_rpkm <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), RPKM = c(mean(c(all_lm_only$lum084, all_lm_only$myo084)), mean(c(all_ls_only$lum084, all_ls_only$stem084)), mean(c(all_ms_only$myo084, all_ms_only$stem084)), mean(c(all_lm_ls_not_ms$lum084, all_lm_ls_not_ms$myo084, all_lm_ls_not_ms$stem084)), mean(c(all_ls_ms_not_lm$lum084, all_ls_ms_not_lm$myo084, all_ls_ms_not_lm$stem084)), mean(c(all_lm_ms_not_ls$lum084, all_lm_ms_not_ls$myo084, all_lm_ms_not_ls$stem084)), mean(c(all_lm_ls_ms$lum084, all_lm_ls_ms$myo084, all_lm_ls_ms$stem084)), with(RM084gene, mean(c(lum084[lum084 > 0.1], myo084[myo084 > 0.1], stem084[stem084 > 0.1]))))) 
isoform_valid_rpkm <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), RPKM = c(mean(c(valid_lm_only$lum084, valid_lm_only$myo084)), mean(c(valid_ls_only$lum084, valid_ls_only$stem084)), mean(c(valid_ms_only$myo084, valid_ms_only$stem084)), mean(c(valid_lm_ls_not_ms$lum084, valid_lm_ls_not_ms$myo084, valid_lm_ls_not_ms$stem084)), mean(c(valid_ls_ms_not_lm$lum084, valid_ls_ms_not_lm$myo084, valid_ls_ms_not_lm$stem084)), mean(c(valid_lm_ms_not_ls$lum084, valid_lm_ms_not_ls$myo084, valid_lm_ms_not_ls$stem084)), mean(c(valid_lm_ls_ms$lum084, valid_lm_ls_ms$myo084, valid_lm_ls_ms$stem084)), with(RM084gene, mean(c(lum084[lum084 > 0.1], myo084[myo084 > 0.1], stem084[stem084 > 0.1]))))) 
pdf("venn_rpkm.pdf")
(bubble_all_rpkm <- ggplot(isoform_all_rpkm) + 
   geom_point(aes(x, y, size = RPKM), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size_area(max_size = 30, breaks = c(0.2, 2, 20)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_all <- venn.diagram(isoform_all, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average gene RPKM of isoforms", main.cex = 1.5)
grid.draw(venn_all)
(bubble_valid_rpkm <- ggplot(isoform_valid_rpkm) + 
   geom_point(aes(x, y, size = RPKM), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size_area(max_size = 30, breaks = c(0.2, 2, 20)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_valid <- venn.diagram(isoform_valid, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average gene RPKM of validated isoforms", main.cex = 1.5)
grid.draw(venn_valid)
dev.off()

# Venn diagram with average No. of exons
isoform_all_Nexon <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), No.exon = c(mean(all_lm_only$Nexon), mean(all_ls_only$Nexon), mean(all_ms_only$Nexon), mean(all_lm_ls_not_ms$Nexon), mean(all_ls_ms_not_lm$Nexon), mean(all_lm_ms_not_ls$Nexon), mean(all_lm_ls_ms$Nexon), mean(Nexon[RM084gene[(RM084gene$lum084 + RM084gene$myo084 + RM084gene$stem084) > 0.3, "id"], "Nexon"]))) 
isoform_valid_Nexon <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), No.exon = c(mean(valid_lm_only$Nexon), mean(valid_ls_only$Nexon), mean(valid_ms_only$Nexon), mean(valid_lm_ls_not_ms$Nexon), mean(valid_ls_ms_not_lm$Nexon), mean(valid_lm_ms_not_ls$Nexon), mean(valid_lm_ls_ms$Nexon), mean(Nexon[RM084gene[(RM084gene$lum084 + RM084gene$myo084 + RM084gene$stem084) > 0.3, "id"], "Nexon"]))) 
pdf("venn_Nexons.pdf")
(bubble_all_Nexon <- ggplot(isoform_all_Nexon) + 
   geom_point(aes(x, y, size = No.exon), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(9, 15), range = c(5, 15), breaks = c(10, 12)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_all <- venn.diagram(isoform_all, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average No. of exons of isoforms", main.cex = 1.5)
grid.draw(venn_all)
(bubble_valid_Nexon <- ggplot(isoform_valid_Nexon) + 
   geom_point(aes(x, y, size = No.exon), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(9, 15), range = c(5, 15), breaks = c(10, 12)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_valid <- venn.diagram(isoform_valid, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average No. of exons of validated isoforms", main.cex = 1.5)
grid.draw(venn_valid)
dev.off()

# Venn diagram with average isoform exon length
exon_length <- read.delim("~/hg19/hg19v65_exons_for_genes.length", head = F, as.is = T)
exon_length$gene <- gsub("chr[0-9XY:-]+<[-]*1_", "", exon_length$V1)
isoform_all_exon_length <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), exon_length = c(mean(all_lm_only$exon_length), mean(all_ls_only$exon_length), mean(all_ms_only$exon_length), mean(all_lm_ls_not_ms$exon_length), mean(all_ls_ms_not_lm$exon_length), mean(all_lm_ms_not_ls$exon_length), mean(all_lm_ls_ms$exon_length), mean(exon_length[exon_length$gene %in% RM084gene[(RM084gene$lum084 + RM084gene$myo084 + RM084gene$stem084) > 0.3, "id"], "V2"]))) 
isoform_valid_exon_length <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), exon_length = c(mean(valid_lm_only$exon_length), mean(valid_ls_only$exon_length), mean(valid_ms_only$exon_length), mean(valid_lm_ls_not_ms$exon_length), mean(valid_ls_ms_not_lm$exon_length), mean(valid_lm_ms_not_ls$exon_length), mean(valid_lm_ls_ms$exon_length), mean(exon_length[exon_length$gene %in% RM084gene[(RM084gene$lum084 + RM084gene$myo084 + RM084gene$stem084) > 0.3, "id"], "V2"]))) 
pdf("venn_exon_length.pdf")
(bubble_all_exon_length <- ggplot(isoform_all_exon_length) + 
   geom_point(aes(x, y, size = exon_length), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(100, 400), range = c(5, 20), breaks = c(200, 300, 400)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_all <- venn.diagram(isoform_all, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average length of isoform exons", main.cex = 1.5)
grid.draw(venn_all)
(bubble_valid_exon_length <- ggplot(isoform_valid_exon_length) + 
   geom_point(aes(x, y, size = exon_length), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(100, 400), range = c(5, 20), breaks = c(200, 300, 400)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_valid <- venn.diagram(isoform_valid, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average length of validated isoform exons", main.cex = 1.5)
grid.draw(venn_valid)
dev.off()

###########################################################################################################
setwd("~/快盘/REMC/junction/all/")
load("junction_valid.Rdata")
lum084_myo084_isoform_all_gene <- lum084_myo084_isoform_all[!duplicated(lum084_myo084_isoform_all$id), ]
lum084_stem084_isoform_all_gene <- lum084_stem084_isoform_all[!duplicated(lum084_stem084_isoform_all$id), ]
myo084_stem084_isoform_all_gene <- myo084_stem084_isoform_all[!duplicated(myo084_stem084_isoform_all$id), ]
lum084_myo084_isoform_valid_gene <- lum084_myo084_isoform_valid[!duplicated(lum084_myo084_isoform_valid$gene), ]
lum084_stem084_isoform_valid_gene <- lum084_stem084_isoform_valid[!duplicated(lum084_stem084_isoform_valid$gene), ]
myo084_stem084_isoform_valid_gene <- myo084_stem084_isoform_valid[!duplicated(myo084_stem084_isoform_valid$gene), ]
lum080_myo080_isoform_all_gene <- lum080_myo080_isoform_all[!duplicated(lum080_myo080_isoform_all$id), ]
lum080_stem080_isoform_all_gene <- lum080_stem080_isoform_all[!duplicated(lum080_stem080_isoform_all$id), ]
myo080_stem080_isoform_all_gene <- myo080_stem080_isoform_all[!duplicated(myo080_stem080_isoform_all$id), ]
lum080_myo080_isoform_valid_gene <- lum080_myo080_isoform_valid[!duplicated(lum080_myo080_isoform_valid$gene), ]
lum080_stem080_isoform_valid_gene <- lum080_stem080_isoform_valid[!duplicated(lum080_stem080_isoform_valid$gene), ]
myo080_stem080_isoform_valid_gene <- myo080_stem080_isoform_valid[!duplicated(myo080_stem080_isoform_valid$gene), ]
lum035_myo035_isoform_all_gene <- lum035_myo035_isoform_all[!duplicated(lum035_myo035_isoform_all$id), ]
lum035_myo035_isoform_valid_gene <- lum035_myo035_isoform_valid[!duplicated(lum035_myo035_isoform_valid$gene), ]
save.image("junction_valid.Rdata")

setwd("~/快盘/REMC/junction/all/Venn/")
library(VennDiagram)
# RM084
isoform_RM084_all <- list(lum084_myo084 = lum084_myo084_isoform_all_gene$id, lum084_stem084 = lum084_stem084_isoform_all_gene$id, myo084_stem084 = myo084_stem084_isoform_all_gene$id)
pdf("venn_RM084_all.pdf")
plot.new()
venn_RM084_all <- venn.diagram(isoform_RM084_all, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of all RM084 isoforms")
grid.draw(venn_RM084_all)
dev.off()
isoform_RM084_valid <- list(lum084_myo084 = as.character(lum084_myo084_isoform_valid_gene$gene), lum084_stem084 = as.character(lum084_stem084_isoform_valid_gene$gene), myo084_stem084 = as.character(myo084_stem084_isoform_valid_gene$gene))
pdf("venn_RM084_valid.pdf")
plot.new()
venn_RM084_valid <- venn.diagram(isoform_RM084_valid, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of validated RM084 isoforms")
grid.draw(venn_RM084_valid)
dev.off()
# RM080
isoform_RM080_all <- list(lum080_myo080 = lum080_myo080_isoform_all_gene$id, lum080_stem080 = lum080_stem080_isoform_all_gene$id, myo080_stem080 = myo080_stem080_isoform_all_gene$id)
pdf("venn_RM080_all.pdf")
plot.new()
venn_RM080_all <- venn.diagram(isoform_RM080_all, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of all RM080 isoforms")
grid.draw(venn_RM080_all)
dev.off()
isoform_RM080_valid <- list(lum080_myo080 = as.character(lum080_myo080_isoform_valid_gene$gene), lum080_stem080 = as.character(lum080_stem080_isoform_valid_gene$gene), myo080_stem080 = as.character(myo080_stem080_isoform_valid_gene$gene))
pdf("venn_RM080_valid.pdf")
plot.new()
venn_RM080_valid <- venn.diagram(isoform_RM080_valid, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of validated RM080 isoforms")
grid.draw(venn_RM080_valid)
dev.off()
# lum vs myo 
isoform_lum_myo_all <- list(RM084 = lum084_myo084_isoform_all_gene$id, RM080 = lum080_myo080_isoform_all_gene$id, RM035 = lum035_myo035_isoform_all_gene$id)
pdf("venn_lum_myo_all.pdf")
plot.new()
venn_lum_myo_all <- venn.diagram(isoform_lum_myo_all, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of all lum vs myo isoforms")
grid.draw(venn_lum_myo_all)
dev.off()
isoform_lum_myo_valid <- list(RM084 = lum084_myo084_isoform_valid_gene$gene, RM080 = lum080_myo080_isoform_valid_gene$gene, RM035 = lum035_myo035_isoform_valid_gene$gene)
pdf("venn_lum_myo_valid.pdf")
plot.new()
venn_lum_myo_valid <- venn.diagram(isoform_lum_myo_valid, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of valid lum vs myo isoforms")
grid.draw(venn_lum_myo_valid)
dev.off()
# lum vs stem 
isoform_lum_stem_all <- list(RM084 = lum084_stem084_isoform_all_gene$id, RM080 = lum080_stem080_isoform_all_gene$id)
pdf("venn_lum_stem_all.pdf")
plot.new()
venn_lum_stem_all <- venn.diagram(isoform_lum_stem_all, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of all lum vs stem isoforms")
grid.draw(venn_lum_stem_all)
dev.off()
isoform_lum_stem_valid <- list(RM084 = lum084_stem084_isoform_valid_gene$gene, RM080 = lum080_stem080_isoform_valid_gene$gene)
pdf("venn_lum_stem_valid.pdf")
plot.new()
venn_lum_stem_valid <- venn.diagram(isoform_lum_stem_valid, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of valid lum vs stem isoforms")
grid.draw(venn_lum_stem_valid)
dev.off()
# myo vs stem 
isoform_myo_stem_all <- list(RM084 = myo084_stem084_isoform_all_gene$id, RM080 = myo080_stem080_isoform_all_gene$id)
pdf("venn_myo_stem_all.pdf")
plot.new()
venn_myo_stem_all <- venn.diagram(isoform_myo_stem_all, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of all myo vs stem isoforms")
grid.draw(venn_myo_stem_all)
dev.off()
isoform_myo_stem_valid <- list(RM084 = myo084_stem084_isoform_valid_gene$gene, RM080 = myo080_stem080_isoform_valid_gene$gene)
pdf("venn_myo_stem_valid.pdf")
plot.new()
venn_myo_stem_valid <- venn.diagram(isoform_myo_stem_valid, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of valid myo vs stem isoforms")
grid.draw(venn_myo_stem_valid)
dev.off()

##############################################################################################################################
# sample clustering on junction RPKM  
setwd("~/快盘/REMC/junction/all/")
load("junction_all.Rdata")
c <- cor(junction[, 10:17], method="spearman")
d <- as.dist(1 - c)
hc <- as.dendrogram(hclust(d, method = "complete"))
pdf("junction_cluster.pdf")
plot(hc, main = "Junction RPKM clustering", horiz = TRUE)
dev.off()
write.table(c, file = "junction_cor.txt", sep = "\t", quote = F)





