# Epigenetic profile at exon boundaries
setwd("~/快盘/REMC/epiProfile/")
library(ggplot2)
# get lum-specific, myo-specific, expressed in both cell types exons  
isoform <- read.delim("lum084_myo084_isoform_all.txt", head = T, as.is = T)
isoform$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform$V1)
isoform$coord  <- gsub("chr", "", isoform$coord)
isoform$coord  <- gsub(":", "_", isoform$coord)
isoform$coord  <- gsub("-", "_", isoform$coord)
isoform$coord  <- paste0(isoform$id, "_", isoform$coord)
lum_specific <- isoform[isoform$V2 < isoform$V3, "coord"]
myo_specific <- isoform[isoform$V2 > isoform$V3, "coord"]
lum084 <- read.delim("~/REMC/epiProfile/A17918.G.exn.A.rpkm", head = F, as.is = T)
lum084$coord <- gsub("<[-]*1", "", lum084$V1)
lum084$coord  <- gsub("chr", "", lum084$coord)
lum084$coord  <- gsub(":", "_", lum084$coord)
lum084$coord  <- gsub("-", "_", lum084$coord)
lum084$coord  <- paste0(lum084$V2, "_", lum084$coord)
myo084 <- read.delim("~/REMC/epiProfile/A17919.G.exn.A.rpkm", head = F, as.is = T)
myo084$coord <- gsub("<[-]*1", "", myo084$V1)
myo084$coord  <- gsub("chr", "", myo084$coord)
myo084$coord  <- gsub(":", "_", myo084$coord)
myo084$coord  <- gsub("-", "_", myo084$coord)
myo084$coord  <- paste0(myo084$V2, "_", myo084$coord)
both <- intersect(lum084[lum084$V4 > 0.1, "coord"], myo084[myo084$V4 > 0.1, "coord"])
neither <- setdiff(lum084$coord, c(lum_specific, myo_specific, both))
rm(lum084, myo084, isoform)

######################################################################################################
# WGBS profile @ exon boundaries
lum_bismark_3p <- read.table("~/REMC/epiProfile/exons3p_200/lumRM066_bismark.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_bismark_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM045_bismark.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
lum_bismark_5p <- read.table("~/REMC/epiProfile/exons5p_200/lumRM066_bismark.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_bismark_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM045_bismark.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
WGBS_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), WGBS = -1)
WGBS_3p$group <- interaction(WGBS_3p$Cell_type, WGBS_3p$Expression)
WGBS_3p[WGBS_3p$group == "lum.lum-specific", "WGBS"] <- colMeans(lum_bismark_3p[lum_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "lum.myo-specific", "WGBS"] <- colMeans(lum_bismark_3p[myo_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "lum.expressed_in_both", "WGBS"] <- colMeans(lum_bismark_3p[both,], na.rm = T)
WGBS_3p[WGBS_3p$group == "lum.not_expressed", "WGBS"] <- colMeans(lum_bismark_3p[neither,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.lum-specific", "WGBS"] <- colMeans(myo_bismark_3p[lum_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.myo-specific", "WGBS"] <- colMeans(myo_bismark_3p[myo_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.expressed_in_both", "WGBS"] <- colMeans(myo_bismark_3p[both,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.not_expressed", "WGBS"] <- colMeans(myo_bismark_3p[neither,], na.rm = T)
WGBS_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), WGBS = -1)
WGBS_5p$group <- interaction(WGBS_5p$Cell_type, WGBS_5p$Expression)
WGBS_5p[WGBS_5p$group == "lum.lum-specific", "WGBS"] <- colMeans(lum_bismark_5p[lum_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "lum.myo-specific", "WGBS"] <- colMeans(lum_bismark_5p[myo_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "lum.expressed_in_both", "WGBS"] <- colMeans(lum_bismark_5p[both,], na.rm = T)
WGBS_5p[WGBS_5p$group == "lum.not_expressed", "WGBS"] <- colMeans(lum_bismark_5p[neither,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.lum-specific", "WGBS"] <- colMeans(myo_bismark_5p[lum_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.myo-specific", "WGBS"] <- colMeans(myo_bismark_5p[myo_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.expressed_in_both", "WGBS"] <- colMeans(myo_bismark_5p[both,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.not_expressed", "WGBS"] <- colMeans(myo_bismark_5p[neither,], na.rm = T)
WGBS_boundaries <- data.frame(rbind(WGBS_3p, WGBS_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(WGBS_3p)), levels = c("5-prime", "3-prime")))
(WGBS_boundaries_profile <- ggplot(WGBS_boundaries, aes(x = Position, y = WGBS, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type), size = 1.5) + 
   geom_point(aes(color = Expression, shape = Cell_type), size = 4) + 
   facet_wrap(~ End) + 
   ggtitle("WGBS profile around exon boundaries") + 
   ylab("Average DNA methylation level") + 
   theme_bw())
ggsave(WGBS_boundaries_profile, file = "WGBS_boundaries_profile.pdf")

######################################################################################################
# MeDIP profile @ exon boundaries
lum_MeDIP_3p <- read.table("~/REMC/epiProfile/exons3p_200/lumRM035_MeDIP.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_MeDIP_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM035_MeDIP.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
lum_MeDIP_5p <- read.table("~/REMC/epiProfile/exons5p_200/lumRM035_MeDIP.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_MeDIP_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM035_MeDIP.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
MeDIP_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), MeDIP = -1)
MeDIP_3p$group <- interaction(MeDIP_3p$Cell_type, MeDIP_3p$Expression)
MeDIP_3p[MeDIP_3p$group == "lum.lum-specific", "MeDIP"] <- colMeans(lum_MeDIP_3p[lum_specific,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "lum.myo-specific", "MeDIP"] <- colMeans(lum_MeDIP_3p[myo_specific,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "lum.expressed_in_both", "MeDIP"] <- colMeans(lum_MeDIP_3p[both,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "lum.not_expressed", "MeDIP"] <- colMeans(lum_MeDIP_3p[neither,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "myo.lum-specific", "MeDIP"] <- colMeans(myo_MeDIP_3p[lum_specific,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "myo.myo-specific", "MeDIP"] <- colMeans(myo_MeDIP_3p[myo_specific,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "myo.expressed_in_both", "MeDIP"] <- colMeans(myo_MeDIP_3p[both,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "myo.not_expressed", "MeDIP"] <- colMeans(myo_MeDIP_3p[neither,], na.rm = T)
MeDIP_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), MeDIP = -1)
MeDIP_5p$group <- interaction(MeDIP_5p$Cell_type, MeDIP_5p$Expression)
MeDIP_5p[MeDIP_5p$group == "lum.lum-specific", "MeDIP"] <- colMeans(lum_MeDIP_5p[lum_specific,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "lum.myo-specific", "MeDIP"] <- colMeans(lum_MeDIP_5p[myo_specific,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "lum.expressed_in_both", "MeDIP"] <- colMeans(lum_MeDIP_5p[both,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "lum.not_expressed", "MeDIP"] <- colMeans(lum_MeDIP_5p[neither,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "myo.lum-specific", "MeDIP"] <- colMeans(myo_MeDIP_5p[lum_specific,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "myo.myo-specific", "MeDIP"] <- colMeans(myo_MeDIP_5p[myo_specific,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "myo.expressed_in_both", "MeDIP"] <- colMeans(myo_MeDIP_5p[both,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "myo.not_expressed", "MeDIP"] <- colMeans(myo_MeDIP_5p[neither,], na.rm = T)
MeDIP_boundaries <- data.frame(rbind(MeDIP_3p, MeDIP_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(MeDIP_3p)), levels = c("5-prime", "3-prime")))
(MeDIP_boundaries_profile <- ggplot(MeDIP_boundaries, aes(x = Position, y = MeDIP, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type), size = 1.5) + 
   geom_point(aes(color = Expression, shape = Cell_type), size = 4) + 
   facet_wrap(~ End) + 
   ggtitle("MeDIP profile around exon boundaries") + 
   ylab("Average DNA methylation level") + 
   theme_bw())
ggsave(MeDIP_boundaries_profile, file = "MeDIP_boundaries_profile.pdf")

######################################################################################################
# H3K4me3 profile @ exon boundaries (no lum libraries available)
myo_H3K4me3_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM080_H3K4me3.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_H3K4me3_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM080_H3K4me3.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
H3K4me3_3p <- data.frame(Expression = rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K4me3 = -1)
H3K4me3_3p[H3K4me3_3p$Expression == "lum-specific", "H3K4me3"] <- colMeans(myo_bismark_3p[lum_specific,], na.rm = T)
H3K4me3_3p[H3K4me3_3p$Expression == "myo-specific", "H3K4me3"] <- colMeans(myo_bismark_3p[myo_specific,], na.rm = T)
H3K4me3_3p[H3K4me3_3p$Expression == "expressed_in_both", "H3K4me3"] <- colMeans(myo_bismark_3p[both,], na.rm = T)
H3K4me3_3p[H3K4me3_3p$Expression == "not_expressed", "H3K4me3"] <- colMeans(myo_bismark_3p[neither,], na.rm = T)
H3K4me3_5p <- data.frame(Expression = rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K4me3 = -1)
H3K4me3_5p[H3K4me3_5p$Expression == "lum-specific", "H3K4me3"] <- colMeans(myo_bismark_5p[lum_specific,], na.rm = T)
H3K4me3_5p[H3K4me3_5p$Expression == "myo-specific", "H3K4me3"] <- colMeans(myo_bismark_5p[myo_specific,], na.rm = T)
H3K4me3_5p[H3K4me3_5p$Expression == "expressed_in_both", "H3K4me3"] <- colMeans(myo_bismark_5p[both,], na.rm = T)
H3K4me3_5p[H3K4me3_5p$Expression == "not_expressed", "H3K4me3"] <- colMeans(myo_bismark_5p[neither,], na.rm = T)
H3K4me3_boundaries <- data.frame(rbind(H3K4me3_3p, H3K4me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K4me3_3p)), levels = c("5-prime", "3-prime")))
(H3K4me3_boundaries_profile <- ggplot(H3K4me3_boundaries, aes(x = Position, y = H3K4me3, group = Expression)) + 
   geom_line(aes(color = Expression), size = 1.5) + 
   geom_point(aes(color = Expression), size = 4) + 
   facet_wrap(~ End) + 
   ggtitle("myo H3K4me3 profile around exon boundaries") + 
   ylab("Average H3K4me3 signal") + 
   theme_bw())
ggsave(H3K4me3_boundaries_profile, file = "H3K4me3_boundaries_profile.pdf")

######################################################################################################
# H3K4me1 profile @ exon boundaries
lum_H3K4me1_3p <- read.table("~/REMC/epiProfile/exons3p_200/lumRM080_H3K4me1.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_H3K4me1_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM080_H3K4me1.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
lum_H3K4me1_5p <- read.table("~/REMC/epiProfile/exons5p_200/lumRM080_H3K4me1.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_H3K4me1_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM080_H3K4me1.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
H3K4me1_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K4me1 = -1)
H3K4me1_3p$group <- interaction(H3K4me1_3p$Cell_type, H3K4me1_3p$Expression)
H3K4me1_3p[H3K4me1_3p$group == "lum.lum-specific", "H3K4me1"] <- colMeans(lum_H3K4me1_3p[lum_specific,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "lum.myo-specific", "H3K4me1"] <- colMeans(lum_H3K4me1_3p[myo_specific,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "lum.expressed_in_both", "H3K4me1"] <- colMeans(lum_H3K4me1_3p[both,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "lum.not_expressed", "H3K4me1"] <- colMeans(lum_H3K4me1_3p[neither,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "myo.lum-specific", "H3K4me1"] <- colMeans(myo_H3K4me1_3p[lum_specific,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "myo.myo-specific", "H3K4me1"] <- colMeans(myo_H3K4me1_3p[myo_specific,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "myo.expressed_in_both", "H3K4me1"] <- colMeans(myo_H3K4me1_3p[both,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "myo.not_expressed", "H3K4me1"] <- colMeans(myo_H3K4me1_3p[neither,], na.rm = T)
H3K4me1_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K4me1 = -1)
H3K4me1_5p$group <- interaction(H3K4me1_5p$Cell_type, H3K4me1_5p$Expression)
H3K4me1_5p[H3K4me1_5p$group == "lum.lum-specific", "H3K4me1"] <- colMeans(lum_H3K4me1_5p[lum_specific,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "lum.myo-specific", "H3K4me1"] <- colMeans(lum_H3K4me1_5p[myo_specific,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "lum.expressed_in_both", "H3K4me1"] <- colMeans(lum_H3K4me1_5p[both,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "lum.not_expressed", "H3K4me1"] <- colMeans(lum_H3K4me1_5p[neither,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "myo.lum-specific", "H3K4me1"] <- colMeans(myo_H3K4me1_5p[lum_specific,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "myo.myo-specific", "H3K4me1"] <- colMeans(myo_H3K4me1_5p[myo_specific,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "myo.expressed_in_both", "H3K4me1"] <- colMeans(myo_H3K4me1_5p[both,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "myo.not_expressed", "H3K4me1"] <- colMeans(myo_H3K4me1_5p[neither,], na.rm = T)
H3K4me1_boundaries <- data.frame(rbind(H3K4me1_3p, H3K4me1_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K4me1_3p)), levels = c("5-prime", "3-prime")))
(H3K4me1_boundaries_profile <- ggplot(H3K4me1_boundaries, aes(x = Position, y = H3K4me1, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type), size = 1.5) + 
   geom_point(aes(color = Expression, shape = Cell_type), size = 4) + 
   facet_wrap(~ End) + 
   ggtitle("H3K4me1 profile around exon boundaries") + 
   ylab("Average H3K4me1 signal") + 
   theme_bw())
ggsave(H3K4me1_boundaries_profile, file = "H3K4me1_boundaries_profile.pdf")

######################################################################################################
# H3K9me3 profile @ exon boundaries
lum_H3K9me3_3p <- read.table("~/REMC/epiProfile/exons3p_200/lumRM080_H3K9me3.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_H3K9me3_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM080_H3K9me3.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
lum_H3K9me3_5p <- read.table("~/REMC/epiProfile/exons5p_200/lumRM080_H3K9me3.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_H3K9me3_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM080_H3K9me3.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
H3K9me3_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K9me3 = -1)
H3K9me3_3p$group <- interaction(H3K9me3_3p$Cell_type, H3K9me3_3p$Expression)
H3K9me3_3p[H3K9me3_3p$group == "lum.lum-specific", "H3K9me3"] <- colMeans(lum_H3K9me3_3p[lum_specific,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "lum.myo-specific", "H3K9me3"] <- colMeans(lum_H3K9me3_3p[myo_specific,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "lum.expressed_in_both", "H3K9me3"] <- colMeans(lum_H3K9me3_3p[both,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "lum.not_expressed", "H3K9me3"] <- colMeans(lum_H3K9me3_3p[neither,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "myo.lum-specific", "H3K9me3"] <- colMeans(myo_H3K9me3_3p[lum_specific,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "myo.myo-specific", "H3K9me3"] <- colMeans(myo_H3K9me3_3p[myo_specific,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "myo.expressed_in_both", "H3K9me3"] <- colMeans(myo_H3K9me3_3p[both,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "myo.not_expressed", "H3K9me3"] <- colMeans(myo_H3K9me3_3p[neither,], na.rm = T)
H3K9me3_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K9me3 = -1)
H3K9me3_5p$group <- interaction(H3K9me3_5p$Cell_type, H3K9me3_5p$Expression)
H3K9me3_5p[H3K9me3_5p$group == "lum.lum-specific", "H3K9me3"] <- colMeans(lum_H3K9me3_5p[lum_specific,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "lum.myo-specific", "H3K9me3"] <- colMeans(lum_H3K9me3_5p[myo_specific,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "lum.expressed_in_both", "H3K9me3"] <- colMeans(lum_H3K9me3_5p[both,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "lum.not_expressed", "H3K9me3"] <- colMeans(lum_H3K9me3_5p[neither,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "myo.lum-specific", "H3K9me3"] <- colMeans(myo_H3K9me3_5p[lum_specific,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "myo.myo-specific", "H3K9me3"] <- colMeans(myo_H3K9me3_5p[myo_specific,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "myo.expressed_in_both", "H3K9me3"] <- colMeans(myo_H3K9me3_5p[both,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "myo.not_expressed", "H3K9me3"] <- colMeans(myo_H3K9me3_5p[neither,], na.rm = T)
H3K9me3_boundaries <- data.frame(rbind(H3K9me3_3p, H3K9me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K9me3_3p)), levels = c("5-prime", "3-prime")))
(H3K9me3_boundaries_profile <- ggplot(H3K9me3_boundaries, aes(x = Position, y = H3K9me3, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type), size = 1.5) + 
   geom_point(aes(color = Expression, shape = Cell_type), size = 4) + 
   facet_wrap(~ End) + 
   ggtitle("H3K9me3 profile around exon boundaries") + 
   ylab("Average H3K9me3 signal") + 
   theme_bw())
ggsave(H3K9me3_boundaries_profile, file = "H3K9me3_boundaries_profile.pdf")

######################################################################################################
# H3K27me3 profile @ exon boundaries
lum_H3K27me3_3p <- read.table("~/REMC/epiProfile/exons3p_200/lumRM080_H3K27me3.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_H3K27me3_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM080_H3K27me3.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
lum_H3K27me3_5p <- read.table("~/REMC/epiProfile/exons5p_200/lumRM080_H3K27me3.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_H3K27me3_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM080_H3K27me3.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
H3K27me3_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K27me3 = -1)
H3K27me3_3p$group <- interaction(H3K27me3_3p$Cell_type, H3K27me3_3p$Expression)
H3K27me3_3p[H3K27me3_3p$group == "lum.lum-specific", "H3K27me3"] <- colMeans(lum_H3K27me3_3p[lum_specific,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "lum.myo-specific", "H3K27me3"] <- colMeans(lum_H3K27me3_3p[myo_specific,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "lum.expressed_in_both", "H3K27me3"] <- colMeans(lum_H3K27me3_3p[both,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "lum.not_expressed", "H3K27me3"] <- colMeans(lum_H3K27me3_3p[neither,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "myo.lum-specific", "H3K27me3"] <- colMeans(myo_H3K27me3_3p[lum_specific,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "myo.myo-specific", "H3K27me3"] <- colMeans(myo_H3K27me3_3p[myo_specific,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "myo.expressed_in_both", "H3K27me3"] <- colMeans(myo_H3K27me3_3p[both,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "myo.not_expressed", "H3K27me3"] <- colMeans(myo_H3K27me3_3p[neither,], na.rm = T)
H3K27me3_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K27me3 = -1)
H3K27me3_5p$group <- interaction(H3K27me3_5p$Cell_type, H3K27me3_5p$Expression)
H3K27me3_5p[H3K27me3_5p$group == "lum.lum-specific", "H3K27me3"] <- colMeans(lum_H3K27me3_5p[lum_specific,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "lum.myo-specific", "H3K27me3"] <- colMeans(lum_H3K27me3_5p[myo_specific,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "lum.expressed_in_both", "H3K27me3"] <- colMeans(lum_H3K27me3_5p[both,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "lum.not_expressed", "H3K27me3"] <- colMeans(lum_H3K27me3_5p[neither,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "myo.lum-specific", "H3K27me3"] <- colMeans(myo_H3K27me3_5p[lum_specific,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "myo.myo-specific", "H3K27me3"] <- colMeans(myo_H3K27me3_5p[myo_specific,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "myo.expressed_in_both", "H3K27me3"] <- colMeans(myo_H3K27me3_5p[both,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "myo.not_expressed", "H3K27me3"] <- colMeans(myo_H3K27me3_5p[neither,], na.rm = T)
H3K27me3_boundaries <- data.frame(rbind(H3K27me3_3p, H3K27me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K27me3_3p)), levels = c("5-prime", "3-prime")))
(H3K27me3_boundaries_profile <- ggplot(H3K27me3_boundaries, aes(x = Position, y = H3K27me3, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type), size = 1.5) + 
   geom_point(aes(color = Expression, shape = Cell_type), size = 4) + 
   facet_wrap(~ End) + 
   ggtitle("H3K27me3 profile around exon boundaries") + 
   ylab("Average H3K27me3 signal") + 
   theme_bw())
ggsave(H3K27me3_boundaries_profile, file = "H3K27me3_boundaries_profile.pdf")

######################################################################################################
# H3K36me3 profile @ exon boundaries
lum_H3K36me3_3p <- read.table("~/REMC/epiProfile/exons3p_200/lumRM080_H3K36me3.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_H3K36me3_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM080_H3K36me3.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
lum_H3K36me3_5p <- read.table("~/REMC/epiProfile/exons5p_200/lumRM080_H3K36me3.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_H3K36me3_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM080_H3K36me3.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
H3K36me3_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K36me3 = -1)
H3K36me3_3p$group <- interaction(H3K36me3_3p$Cell_type, H3K36me3_3p$Expression)
H3K36me3_3p[H3K36me3_3p$group == "lum.lum-specific", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[lum_specific,], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.myo-specific", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[myo_specific,], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.expressed_in_both", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[both,], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.not_expressed", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[neither,], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.lum-specific", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[lum_specific,], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.myo-specific", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[myo_specific,], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.expressed_in_both", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[both,], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.not_expressed", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[neither,], na.rm = T)
H3K36me3_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K36me3 = -1)
H3K36me3_5p$group <- interaction(H3K36me3_5p$Cell_type, H3K36me3_5p$Expression)
H3K36me3_5p[H3K36me3_5p$group == "lum.lum-specific", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[lum_specific,], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.myo-specific", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[myo_specific,], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.expressed_in_both", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[both,], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.not_expressed", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[neither,], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.lum-specific", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[lum_specific,], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.myo-specific", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[myo_specific,], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.expressed_in_both", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[both,], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.not_expressed", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[neither,], na.rm = T)
H3K36me3_boundaries <- data.frame(rbind(H3K36me3_3p, H3K36me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K36me3_3p)), levels = c("5-prime", "3-prime")))
(H3K36me3_boundaries_profile <- ggplot(H3K36me3_boundaries, aes(x = Position, y = H3K36me3, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type), size = 1.5) + 
   geom_point(aes(color = Expression, shape = Cell_type), size = 4) + 
   facet_wrap(~ End) + 
   ggtitle("H3K36me3 profile around exon boundaries") + 
   ylab("Average H3K36me3 signal") + 
   theme_bw())
ggsave(H3K36me3_boundaries_profile, file = "H3K36me3_boundaries_profile.pdf")

save(WGBS_boundaries, MeDIP_boundaries, H3K4me3_boundaries, H3K4me1_boundaries, H3K9me3_boundaries, H3K27me3_boundaries, H3K36me3_boundaries, file = "exonProfile.Rdata")

######################################################################################################
# WGBS along exons
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
setwd("~/REMC/epiProfile/exons/")
collapse <- function(gene, raw){
  raw_data <- na.omit(as.numeric(raw[raw$X1 == gene, 2:ncol(raw)]))
  length <- length(raw_data)
  output <- vector(length = 20)
  for(i in 1:19){
    output[i] <- mean(raw_data[(as.integer((i-1)*length/20)+1):(as.integer(i*length/20)+1)])
  }
  output[20] <- mean(raw_data[(as.integer((i-1)*length/20)+1):as.integer(i*length/20)])
  return(output)
}

# less ~/REMC/epiProfile/exons/myoRM045_bismark.hg19v65_exons_for_genes.profile | awk -F' ' '{if(NF>max)max=NF}END{print max}'
# less ~/REMC/epiProfile/exons/lumRM066_bismark.hg19v65_exons_for_genes.profile | awk -F' ' '{if(NF>max)max=NF}END{print max}'
myo_bismark_exons_raw <- read.table("~/REMC/epiProfile/exons/myoRM045_bismark.hg19v65_exons_for_genes.profile", sep = " ", head = F, as.is = T, fill = T, col.names = as.character(1:1462))
myo_bismark_exons <- t(sapply(unique(myo_bismark_exons_raw$X1), collapse, raw = myo_bismark_exons_raw))
save(myo_bismark_exons, "myo_bismark_exons.Rdata")
lum_bismark_exons_raw <- read.table("~/REMC/epiProfile/exons/lumRM066_bismark.hg19v65_exons_for_genes.profile", sep = " ", head = F, as.is = T, fill = T, col.names = as.character(1:1462))
lum_bismark_exons <- t(sapply(unique(lum_bismark_exons_raw$X1), collapse, raw = lum_bismark_exons_raw))
save(lum_bismark_exons, "lum_bismark_exons.Rdata")
WGBS_exons <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(5, 100, by = 5), times = 8), WGBS = -1)
WGBS_exons$group <- interaction(WGBS_exons$Cell_type, WGBS_exons$Expression)
WGBS_exons[WGBS_exons$group == "lum.lum-specific", "WGBS"] <- colMeans(lum_bismark_exons[lum_specific,], na.rm = T)
WGBS_exons[WGBS_exons$group == "lum.myo-specific", "WGBS"] <- colMeans(lum_bismark_exons[myo_specific,], na.rm = T)
WGBS_exons[WGBS_exons$group == "lum.expressed_in_both", "WGBS"] <- colMeans(lum_bismark_exons[both,], na.rm = T)
WGBS_exons[WGBS_exons$group == "lum.not_expressed", "WGBS"] <- colMeans(lum_bismark_exons[neither,], na.rm = T)
WGBS_exons[WGBS_exons$group == "myo.lum-specific", "WGBS"] <- colMeans(myo_bismark_exons[lum_specific,], na.rm = T)
WGBS_exons[WGBS_exons$group == "myo.myo-specific", "WGBS"] <- colMeans(myo_bismark_exons[myo_specific,], na.rm = T)
WGBS_exons[WGBS_exons$group == "myo.expressed_in_both", "WGBS"] <- colMeans(myo_bismark_exons[both,], na.rm = T)
WGBS_exons[WGBS_exons$group == "myo.not_expressed", "WGBS"] <- colMeans(myo_bismark_exons[neither,], na.rm = T)
(WGBS_exons_profile <- ggplot(WGBS_exons, aes(x = Position, y = WGBS, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type), size = 1.5) + 
   geom_point(aes(color = Expression, shape = Cell_type), size = 4) + 
   ggtitle("WGBS profile along exons") + 
   xlab("Normalized position along exons") + 
   ylab("Average DNA methylation level") + 
   theme_bw())
ggsave(WGBS_exons_profile, file = "WGBS_exons_profile")

