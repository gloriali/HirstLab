# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# Epigenetic profile at intron boundaries
setwd("~/快盘/REMC/epiProfile/IR/")
# require(plyr, lib.loc = "/home/lli/bin/R-3.0.2")
library(plyr)
library(ggplot2)
library(grid) 
load("intronProfile.Rdata")
introns <- read.delim("/home/lli/hg19/hg19v65_introns_for_genes", head = F, as.is = T)
lum084_ir <- read.delim("/projects/epigenomics/acarles/intron/breast/4GloriaAnalysis/A17918.introns.retained.Id4Gloria", head = F, as.is = T)$V1
myo084_ir <- read.delim("/projects/epigenomics/acarles/intron/breast/4GloriaAnalysis/A17919.introns.retained.Id4Gloria", head = F, as.is = T)$V1
both <- intersect(lum084_ir, myo084_ir)
lum_specific <- setdiff(lum084_ir, both)
myo_specific <- setdiff(myo084_ir, both)
neither <- setdiff(introns$V5, union(lum084_ir, myo084_ir))
rm(introns)

######################################################################################################
# WGBS profile @ intron boundaries
lum_bismark_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/lumRM066_bismark.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_bismark_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM045_bismark.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_bismark_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/lumRM066_bismark.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_bismark_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM045_bismark.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
WGBS_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), WGBS = -1)
WGBS_3p$group <- interaction(WGBS_3p$Cell_type, WGBS_3p$Expression)
WGBS_3p[WGBS_3p$group == "lum.lum-specific", "WGBS"] <- colMeans(lum_bismark_3p[lum_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "lum.myo-specific", "WGBS"] <- colMeans(lum_bismark_3p[myo_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "lum.retained_in_both", "WGBS"] <- colMeans(lum_bismark_3p[both,], na.rm = T)
WGBS_3p[WGBS_3p$group == "lum.not_retained", "WGBS"] <- colMeans(lum_bismark_3p[neither,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.lum-specific", "WGBS"] <- colMeans(myo_bismark_3p[lum_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.myo-specific", "WGBS"] <- colMeans(myo_bismark_3p[myo_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.retained_in_both", "WGBS"] <- colMeans(myo_bismark_3p[both,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.not_retained", "WGBS"] <- colMeans(myo_bismark_3p[neither,], na.rm = T)
WGBS_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), WGBS = -1)
WGBS_5p$group <- interaction(WGBS_5p$Cell_type, WGBS_5p$Expression)
WGBS_5p[WGBS_5p$group == "lum.lum-specific", "WGBS"] <- colMeans(lum_bismark_5p[lum_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "lum.myo-specific", "WGBS"] <- colMeans(lum_bismark_5p[myo_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "lum.retained_in_both", "WGBS"] <- colMeans(lum_bismark_5p[both,], na.rm = T)
WGBS_5p[WGBS_5p$group == "lum.not_retained", "WGBS"] <- colMeans(lum_bismark_5p[neither,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.lum-specific", "WGBS"] <- colMeans(myo_bismark_5p[lum_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.myo-specific", "WGBS"] <- colMeans(myo_bismark_5p[myo_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.retained_in_both", "WGBS"] <- colMeans(myo_bismark_5p[both,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.not_retained", "WGBS"] <- colMeans(myo_bismark_5p[neither,], na.rm = T)
WGBS_boundaries <- data.frame(rbind(WGBS_3p, WGBS_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(WGBS_3p)), levels = c("5-prime", "3-prime")))
(WGBS_boundaries_profile <- ggplot(WGBS_boundaries, aes(x = Position, y = WGBS, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type)) + 
   geom_point(aes(color = Expression, shape = Cell_type)) + 
   facet_wrap(~ End) + 
   ggtitle("WGBS profile around intron boundaries") + 
   ylab("Average DNA methylation level") + 
   theme_bw())
ggsave(WGBS_boundaries_profile, file = "WGBS_boundaries_profile.pdf")

# add CpG content track 
CpG_content_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/CpG.hg19v65_introns_for_genes.3prime_200", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
CpG_content_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/CpG.hg19v65_introns_for_genes.5prime_200", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
CpG_3p <- data.frame(Expression = c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), Position = rep(seq(-190, 190, by = 20), times = 4), CpG = -1)
CpG_3p[CpG_3p$Expression == "lum-specific", "CpG"] <- colMeans(CpG_content_3p[lum_specific,], na.rm = T)
CpG_3p[CpG_3p$Expression == "myo-specific", "CpG"] <- colMeans(CpG_content_3p[myo_specific,], na.rm = T)
CpG_3p[CpG_3p$Expression == "retained_in_both", "CpG"] <- colMeans(CpG_content_3p[both,], na.rm = T)
CpG_3p[CpG_3p$Expression == "not_retained", "CpG"] <- colMeans(CpG_content_3p[neither,], na.rm = T)
CpG_5p <- data.frame(Expression = c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), Position = rep(seq(-190, 190, by = 20), times = 4), CpG = -1)
CpG_5p[CpG_5p$Expression == "lum-specific", "CpG"] <- colMeans(CpG_content_5p[lum_specific,], na.rm = T)
CpG_5p[CpG_5p$Expression == "myo-specific", "CpG"] <- colMeans(CpG_content_5p[myo_specific,], na.rm = T)
CpG_5p[CpG_5p$Expression == "retained_in_both", "CpG"] <- colMeans(CpG_content_5p[both,], na.rm = T)
CpG_5p[CpG_5p$Expression == "not_retained", "CpG"] <- colMeans(CpG_content_5p[neither,], na.rm = T)
CpG_boundaries <- data.frame(rbind(CpG_3p, CpG_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(CpG_3p)), levels = c("5-prime", "3-prime")))
(CpG_boundaries_profile <- ggplot(CpG_boundaries, aes(x = Position, y = CpG, group = Expression)) + 
   geom_line(aes(color = Expression)) + 
   geom_point(aes(color = Expression)) + 
   facet_wrap(~ End) + 
   ggtitle("CpG profile around intron boundaries") + 
   ylab("Average percentage of CpG") + 
   theme_bw())
ggsave(CpG_boundaries_profile, file = "CpG_boundaries_profile.pdf", height = 4)

WGBS_CpG <- data.frame(data = c(rep("WGBS", nrow(WGBS_boundaries)), rep("CpG", nrow(CpG_boundaries))), 
                       Cell_type = c(as.character(WGBS_boundaries$Cell_type), rep("lum", nrow(CpG_boundaries))), 
                       Expression = c(as.character(WGBS_boundaries$Expression), as.character(CpG_boundaries$Expression)), 
                       Position = c(WGBS_boundaries$Position, CpG_boundaries$Position), 
                       value = c(WGBS_boundaries$WGBS, CpG_boundaries$CpG), 
                       End = c(as.character(WGBS_boundaries$End), as.character(CpG_boundaries$End)))
WGBS_CpG$data <- factor(WGBS_CpG$data, levels = c("WGBS", "CpG"))
WGBS_CpG$End <- factor(WGBS_CpG$End, levels = c("5-prime", "3-prime"))
WGBS_CpG$group <- interaction(WGBS_CpG$Cell_type, WGBS_CpG$Expression)
(WGBS_CpG_profile <- ggplot(WGBS_CpG, aes(x = Position, y = value, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type)) + 
   geom_point(aes(color = Expression, shape = Cell_type)) + 
   facet_grid(data ~ End, scales = "free_y") + 
   ylab("Average DNA methylation") + 
   theme(axis.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 12, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
gt <- ggplot_gtable(ggplot_build(WGBS_CpG_profile)) 
# gt$layout 
gt$heights[[4]] <- unit(2.5, "null") 
pdf("WGBS_CpG_profile.pdf")
grid.draw(gt) 
dev.off()
######################################################################################################
# MeDIP profile @ intron boundaries
lum_MeDIP_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/lumRM035_MeDIP.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_MeDIP_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM035_MeDIP.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_MeDIP_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/lumRM035_MeDIP.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_MeDIP_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM035_MeDIP.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_MeDIP_3p <- sum(myo_MeDIP_3p)/sum(lum_MeDIP_3p) * lum_MeDIP_3p
lum_MeDIP_5p <- sum(myo_MeDIP_5p)/sum(lum_MeDIP_5p) * lum_MeDIP_5p
MeDIP_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), MeDIP = -1)
MeDIP_3p$group <- interaction(MeDIP_3p$Cell_type, MeDIP_3p$Expression)
MeDIP_3p[MeDIP_3p$group == "lum.lum-specific", "MeDIP"] <- colMeans(lum_MeDIP_3p[lum_specific,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "lum.myo-specific", "MeDIP"] <- colMeans(lum_MeDIP_3p[myo_specific,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "lum.retained_in_both", "MeDIP"] <- colMeans(lum_MeDIP_3p[both,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "lum.not_retained", "MeDIP"] <- colMeans(lum_MeDIP_3p[neither,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "myo.lum-specific", "MeDIP"] <- colMeans(myo_MeDIP_3p[lum_specific,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "myo.myo-specific", "MeDIP"] <- colMeans(myo_MeDIP_3p[myo_specific,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "myo.retained_in_both", "MeDIP"] <- colMeans(myo_MeDIP_3p[both,], na.rm = T)
MeDIP_3p[MeDIP_3p$group == "myo.not_retained", "MeDIP"] <- colMeans(myo_MeDIP_3p[neither,], na.rm = T)
MeDIP_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), MeDIP = -1)
MeDIP_5p$group <- interaction(MeDIP_5p$Cell_type, MeDIP_5p$Expression)
MeDIP_5p[MeDIP_5p$group == "lum.lum-specific", "MeDIP"] <- colMeans(lum_MeDIP_5p[lum_specific,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "lum.myo-specific", "MeDIP"] <- colMeans(lum_MeDIP_5p[myo_specific,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "lum.retained_in_both", "MeDIP"] <- colMeans(lum_MeDIP_5p[both,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "lum.not_retained", "MeDIP"] <- colMeans(lum_MeDIP_5p[neither,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "myo.lum-specific", "MeDIP"] <- colMeans(myo_MeDIP_5p[lum_specific,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "myo.myo-specific", "MeDIP"] <- colMeans(myo_MeDIP_5p[myo_specific,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "myo.retained_in_both", "MeDIP"] <- colMeans(myo_MeDIP_5p[both,], na.rm = T)
MeDIP_5p[MeDIP_5p$group == "myo.not_retained", "MeDIP"] <- colMeans(myo_MeDIP_5p[neither,], na.rm = T)
MeDIP_boundaries <- data.frame(rbind(MeDIP_3p, MeDIP_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(MeDIP_3p)), levels = c("5-prime", "3-prime")))
(MeDIP_boundaries_profile <- ggplot(MeDIP_boundaries, aes(x = Position, y = MeDIP, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type)) + 
   geom_point(aes(color = Expression, shape = Cell_type)) + 
   facet_wrap(~ End) + 
   ggtitle("MeDIP profile around intron boundaries") + 
   ylab("Average DNA methylation level") + 
   theme_bw())
ggsave(MeDIP_boundaries_profile, file = "MeDIP_boundaries_profile.pdf")

######################################################################################################
# H3K4me3 profile @ intron boundaries (no lum libraries available)
myo_H3K4me3_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM080_H3K4me3.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K4me3_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM080_H3K4me3.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
H3K4me3_3p <- data.frame(Expression = rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K4me3 = -1)
H3K4me3_3p[H3K4me3_3p$Expression == "lum-specific", "H3K4me3"] <- colMeans(myo_bismark_3p[lum_specific,], na.rm = T)
H3K4me3_3p[H3K4me3_3p$Expression == "myo-specific", "H3K4me3"] <- colMeans(myo_bismark_3p[myo_specific,], na.rm = T)
H3K4me3_3p[H3K4me3_3p$Expression == "retained_in_both", "H3K4me3"] <- colMeans(myo_bismark_3p[both,], na.rm = T)
H3K4me3_3p[H3K4me3_3p$Expression == "not_retained", "H3K4me3"] <- colMeans(myo_bismark_3p[neither,], na.rm = T)
H3K4me3_5p <- data.frame(Expression = rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K4me3 = -1)
H3K4me3_5p[H3K4me3_5p$Expression == "lum-specific", "H3K4me3"] <- colMeans(myo_bismark_5p[lum_specific,], na.rm = T)
H3K4me3_5p[H3K4me3_5p$Expression == "myo-specific", "H3K4me3"] <- colMeans(myo_bismark_5p[myo_specific,], na.rm = T)
H3K4me3_5p[H3K4me3_5p$Expression == "retained_in_both", "H3K4me3"] <- colMeans(myo_bismark_5p[both,], na.rm = T)
H3K4me3_5p[H3K4me3_5p$Expression == "not_retained", "H3K4me3"] <- colMeans(myo_bismark_5p[neither,], na.rm = T)
H3K4me3_boundaries <- data.frame(rbind(H3K4me3_3p, H3K4me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K4me3_3p)), levels = c("5-prime", "3-prime")))
(H3K4me3_boundaries_profile <- ggplot(H3K4me3_boundaries, aes(x = Position, y = H3K4me3, group = Expression)) + 
   geom_line(aes(color = Expression)) + 
   geom_point(aes(color = Expression)) + 
   facet_wrap(~ End) + 
   ggtitle("myo H3K4me3 profile around intron boundaries") + 
   ylab("Average H3K4me3 signal") + 
   theme_bw())
ggsave(H3K4me3_boundaries_profile, file = "H3K4me3_boundaries_profile.pdf")

######################################################################################################
# H3K4me1 profile @ intron boundaries
lum_H3K4me1_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/lumRM080_H3K4me1.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K4me1_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM080_H3K4me1.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K4me1_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/lumRM080_H3K4me1.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K4me1_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM080_H3K4me1.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K4me1_3p <- sum(myo_H3K4me1_3p)/sum(lum_H3K4me1_3p) * lum_H3K4me1_3p
lum_H3K4me1_5p <- sum(myo_H3K4me1_5p)/sum(lum_H3K4me1_5p) * lum_H3K4me1_5p
H3K4me1_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K4me1 = -1)
H3K4me1_3p$group <- interaction(H3K4me1_3p$Cell_type, H3K4me1_3p$Expression)
H3K4me1_3p[H3K4me1_3p$group == "lum.lum-specific", "H3K4me1"] <- colMeans(lum_H3K4me1_3p[lum_specific,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "lum.myo-specific", "H3K4me1"] <- colMeans(lum_H3K4me1_3p[myo_specific,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "lum.retained_in_both", "H3K4me1"] <- colMeans(lum_H3K4me1_3p[both,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "lum.not_retained", "H3K4me1"] <- colMeans(lum_H3K4me1_3p[neither,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "myo.lum-specific", "H3K4me1"] <- colMeans(myo_H3K4me1_3p[lum_specific,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "myo.myo-specific", "H3K4me1"] <- colMeans(myo_H3K4me1_3p[myo_specific,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "myo.retained_in_both", "H3K4me1"] <- colMeans(myo_H3K4me1_3p[both,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$group == "myo.not_retained", "H3K4me1"] <- colMeans(myo_H3K4me1_3p[neither,], na.rm = T)
H3K4me1_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K4me1 = -1)
H3K4me1_5p$group <- interaction(H3K4me1_5p$Cell_type, H3K4me1_5p$Expression)
H3K4me1_5p[H3K4me1_5p$group == "lum.lum-specific", "H3K4me1"] <- colMeans(lum_H3K4me1_5p[lum_specific,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "lum.myo-specific", "H3K4me1"] <- colMeans(lum_H3K4me1_5p[myo_specific,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "lum.retained_in_both", "H3K4me1"] <- colMeans(lum_H3K4me1_5p[both,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "lum.not_retained", "H3K4me1"] <- colMeans(lum_H3K4me1_5p[neither,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "myo.lum-specific", "H3K4me1"] <- colMeans(myo_H3K4me1_5p[lum_specific,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "myo.myo-specific", "H3K4me1"] <- colMeans(myo_H3K4me1_5p[myo_specific,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "myo.retained_in_both", "H3K4me1"] <- colMeans(myo_H3K4me1_5p[both,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$group == "myo.not_retained", "H3K4me1"] <- colMeans(myo_H3K4me1_5p[neither,], na.rm = T)
H3K4me1_boundaries <- data.frame(rbind(H3K4me1_3p, H3K4me1_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K4me1_3p)), levels = c("5-prime", "3-prime")))
(H3K4me1_boundaries_profile <- ggplot(H3K4me1_boundaries, aes(x = Position, y = H3K4me1, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type)) + 
   geom_point(aes(color = Expression, shape = Cell_type)) + 
   facet_wrap(~ End) + 
   ggtitle("H3K4me1 profile around intron boundaries") + 
   ylab("Average H3K4me1 signal") + 
   theme_bw())
ggsave(H3K4me1_boundaries_profile, file = "H3K4me1_boundaries_profile.pdf")

######################################################################################################
# H3K9me3 profile @ intron boundaries
lum_H3K9me3_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/lumRM080_H3K9me3.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K9me3_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM080_H3K9me3.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K9me3_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/lumRM080_H3K9me3.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K9me3_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM080_H3K9me3.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K9me3_3p <- sum(myo_H3K9me3_3p)/sum(lum_H3K9me3_3p) * lum_H3K9me3_3p
lum_H3K9me3_5p <- sum(myo_H3K9me3_5p)/sum(lum_H3K9me3_5p) * lum_H3K9me3_5p
H3K9me3_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K9me3 = -1)
H3K9me3_3p$group <- interaction(H3K9me3_3p$Cell_type, H3K9me3_3p$Expression)
H3K9me3_3p[H3K9me3_3p$group == "lum.lum-specific", "H3K9me3"] <- colMeans(lum_H3K9me3_3p[lum_specific,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "lum.myo-specific", "H3K9me3"] <- colMeans(lum_H3K9me3_3p[myo_specific,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "lum.retained_in_both", "H3K9me3"] <- colMeans(lum_H3K9me3_3p[both,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "lum.not_retained", "H3K9me3"] <- colMeans(lum_H3K9me3_3p[neither,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "myo.lum-specific", "H3K9me3"] <- colMeans(myo_H3K9me3_3p[lum_specific,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "myo.myo-specific", "H3K9me3"] <- colMeans(myo_H3K9me3_3p[myo_specific,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "myo.retained_in_both", "H3K9me3"] <- colMeans(myo_H3K9me3_3p[both,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$group == "myo.not_retained", "H3K9me3"] <- colMeans(myo_H3K9me3_3p[neither,], na.rm = T)
H3K9me3_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K9me3 = -1)
H3K9me3_5p$group <- interaction(H3K9me3_5p$Cell_type, H3K9me3_5p$Expression)
H3K9me3_5p[H3K9me3_5p$group == "lum.lum-specific", "H3K9me3"] <- colMeans(lum_H3K9me3_5p[lum_specific,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "lum.myo-specific", "H3K9me3"] <- colMeans(lum_H3K9me3_5p[myo_specific,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "lum.retained_in_both", "H3K9me3"] <- colMeans(lum_H3K9me3_5p[both,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "lum.not_retained", "H3K9me3"] <- colMeans(lum_H3K9me3_5p[neither,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "myo.lum-specific", "H3K9me3"] <- colMeans(myo_H3K9me3_5p[lum_specific,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "myo.myo-specific", "H3K9me3"] <- colMeans(myo_H3K9me3_5p[myo_specific,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "myo.retained_in_both", "H3K9me3"] <- colMeans(myo_H3K9me3_5p[both,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$group == "myo.not_retained", "H3K9me3"] <- colMeans(myo_H3K9me3_5p[neither,], na.rm = T)
H3K9me3_boundaries <- data.frame(rbind(H3K9me3_3p, H3K9me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K9me3_3p)), levels = c("5-prime", "3-prime")))
(H3K9me3_boundaries_profile <- ggplot(H3K9me3_boundaries, aes(x = Position, y = H3K9me3, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type)) + 
   geom_point(aes(color = Expression, shape = Cell_type)) + 
   facet_wrap(~ End) + 
   ggtitle("H3K9me3 profile around intron boundaries") + 
   ylab("Average H3K9me3 signal") + 
   theme_bw())
ggsave(H3K9me3_boundaries_profile, file = "H3K9me3_boundaries_profile.pdf")

######################################################################################################
# H3K27me3 profile @ intron boundaries
lum_H3K27me3_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/lumRM080_H3K27me3.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K27me3_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM080_H3K27me3.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K27me3_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/lumRM080_H3K27me3.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K27me3_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM080_H3K27me3.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K27me3_3p <- sum(myo_H3K27me3_3p)/sum(lum_H3K27me3_3p) * lum_H3K27me3_3p
lum_H3K27me3_5p <- sum(myo_H3K27me3_5p)/sum(lum_H3K27me3_5p) * lum_H3K27me3_5p
H3K27me3_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K27me3 = -1)
H3K27me3_3p$group <- interaction(H3K27me3_3p$Cell_type, H3K27me3_3p$Expression)
H3K27me3_3p[H3K27me3_3p$group == "lum.lum-specific", "H3K27me3"] <- colMeans(lum_H3K27me3_3p[lum_specific,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "lum.myo-specific", "H3K27me3"] <- colMeans(lum_H3K27me3_3p[myo_specific,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "lum.retained_in_both", "H3K27me3"] <- colMeans(lum_H3K27me3_3p[both,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "lum.not_retained", "H3K27me3"] <- colMeans(lum_H3K27me3_3p[neither,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "myo.lum-specific", "H3K27me3"] <- colMeans(myo_H3K27me3_3p[lum_specific,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "myo.myo-specific", "H3K27me3"] <- colMeans(myo_H3K27me3_3p[myo_specific,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "myo.retained_in_both", "H3K27me3"] <- colMeans(myo_H3K27me3_3p[both,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$group == "myo.not_retained", "H3K27me3"] <- colMeans(myo_H3K27me3_3p[neither,], na.rm = T)
H3K27me3_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "retained_in_both", "not_retained"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K27me3 = -1)
H3K27me3_5p$group <- interaction(H3K27me3_5p$Cell_type, H3K27me3_5p$Expression)
H3K27me3_5p[H3K27me3_5p$group == "lum.lum-specific", "H3K27me3"] <- colMeans(lum_H3K27me3_5p[lum_specific,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "lum.myo-specific", "H3K27me3"] <- colMeans(lum_H3K27me3_5p[myo_specific,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "lum.retained_in_both", "H3K27me3"] <- colMeans(lum_H3K27me3_5p[both,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "lum.not_retained", "H3K27me3"] <- colMeans(lum_H3K27me3_5p[neither,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "myo.lum-specific", "H3K27me3"] <- colMeans(myo_H3K27me3_5p[lum_specific,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "myo.myo-specific", "H3K27me3"] <- colMeans(myo_H3K27me3_5p[myo_specific,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "myo.retained_in_both", "H3K27me3"] <- colMeans(myo_H3K27me3_5p[both,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$group == "myo.not_retained", "H3K27me3"] <- colMeans(myo_H3K27me3_5p[neither,], na.rm = T)
H3K27me3_boundaries <- data.frame(rbind(H3K27me3_3p, H3K27me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K27me3_3p)), levels = c("5-prime", "3-prime")))
(H3K27me3_boundaries_profile <- ggplot(H3K27me3_boundaries, aes(x = Position, y = H3K27me3, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type)) + 
   geom_point(aes(color = Expression, shape = Cell_type)) + 
   facet_wrap(~ End) + 
   ggtitle("H3K27me3 profile around intron boundaries") + 
   ylab("Average H3K27me3 signal") + 
   theme_bw())
ggsave(H3K27me3_boundaries_profile, file = "H3K27me3_boundaries_profile.pdf")

######################################################################################################
# H3K36me3 signal for introns
# separate on gene RPKM for H3K36me3: <1; 1-10; 10-100; >100 
lum084gene <- read.delim("/projects/epigenomics/ep50/internal/jqc.1.7.6/A17918/coverage/A17918.G.A.rpkm.pc", head = F, as.is = T)
myo084gene <- read.delim("/projects/epigenomics/ep50/internal/jqc.1.7.6/A17919/coverage/A17919.G.A.rpkm.pc", head = F, as.is = T)
geneRPKM <- data.frame(gene = lum084gene$V1, RPKM = (lum084gene$V3 + myo084gene$V3)/2)
rownames(geneRPKM) <- geneRPKM$gene
introns_geneRPKM <- data.frame(intron = c(both, lum_specific, myo_specific, neither))
introns_geneRPKM$gene <- gsub("_chr[0-9XYM]+:[0-9_<-]+", "", introns_geneRPKM$intron)
introns_geneRPKM$geneRPKM <- geneRPKM[introns_geneRPKM$gene, "RPKM"]
introns_geneRPKM <- na.omit(introns_geneRPKM)
introns_1 <- as.character(introns_geneRPKM[introns_geneRPKM$geneRPKM <= 1, "intron"])
introns_1_10 <- as.character(introns_geneRPKM[introns_geneRPKM$geneRPKM > 1 & introns_geneRPKM$geneRPKM <= 10, "intron"])
introns_10_100 <- as.character(introns_geneRPKM[introns_geneRPKM$geneRPKM > 10 & introns_geneRPKM$geneRPKM <= 100, "intron"])
introns_100 <- as.character(introns_geneRPKM[introns_geneRPKM$geneRPKM > 100, "intron"])
rm(lum084gene, myo084gene)

lum_H3K36me3_introns <- read.delim("~/REMC/epiProfile/IR/introns/hg19v65_introns_for_genes.lumRM080_H3K36me3.coverage", head = F, as.is = T)
myo_H3K36me3_introns <- read.delim("~/REMC/epiProfile/IR/introns/hg19v65_introns_for_genes.myoRM080_H3K36me3.coverage", head = F, as.is = T)
# normalize signal
(norm <- sum(myo_H3K36me3_introns$V6)/sum(lum_H3K36me3_introns$V6))
H3K36me3_introns <- data.frame(id = c(lum_H3K36me3_introns$V4, myo_H3K36me3_introns$V4), Cell_type = c(rep("lum", nrow(lum_H3K36me3_introns)), rep("myo", nrow(myo_H3K36me3_introns))), geneRPKM = NA, Expression = NA, H3K36me3 = c(lum_H3K36me3_introns$V6 * norm, myo_H3K36me3_introns$V6))
H3K36me3_introns[H3K36me3_introns$id %in% both, "Expression"] <- "retained_in_both"
H3K36me3_introns[H3K36me3_introns$id %in% neither, "Expression"] <- "not_retained"
H3K36me3_introns[H3K36me3_introns$id %in% lum_specific, "Expression"] <- "lum-specific"
H3K36me3_introns[H3K36me3_introns$id %in% myo_specific, "Expression"] <- "myo-specific"
H3K36me3_introns$Expression <- factor(H3K36me3_introns$Expression)
H3K36me3_introns[H3K36me3_introns$id %in% introns_1, "geneRPKM"] <- "gene RPKM < 1"
H3K36me3_introns[H3K36me3_introns$id %in% introns_1_10, "geneRPKM"] <- "gene RPKM 1-10"
H3K36me3_introns[H3K36me3_introns$id %in% introns_10_100, "geneRPKM"] <- "gene RPKM 10-100"
H3K36me3_introns[H3K36me3_introns$id %in% introns_100, "geneRPKM"] <- "gene RPKM > 100"
H3K36me3_introns$geneRPKM <- factor(H3K36me3_introns$geneRPKM)
H3K36me3_introns$utilize <- interaction(H3K36me3_introns$Cell_type, H3K36me3_introns$Expression)
# fold enrichment between utilized/un-utilized isoform introns
mean(H3K36me3_introns[H3K36me3_introns$utilize == "lum.lum-specific" | H3K36me3_introns$utilize == "myo.myo-specific", "H3K36me3"])/mean(H3K36me3_introns[H3K36me3_introns$utilize == "myo.lum-specific" | H3K36me3_introns$utilize == "lum.myo-specific", "H3K36me3"])
H3K36me3_introns <- droplevels(na.omit(H3K36me3_introns[H3K36me3_introns$geneRPKM != "gene RPKM > 100", ]))
H3K36me3_introns$group <- interaction(interaction(H3K36me3_introns$Cell_type, H3K36me3_introns$Expression), H3K36me3_introns$geneRPKM)
H3K36me3_introns_stat <- ddply(H3K36me3_introns, ~ group, summarize, Cell_type = Cell_type[1], Expression = Expression[1], geneRPKM = geneRPKM[1], ymin = boxplot.stats(H3K36me3)$stats[1], lower = boxplot.stats(H3K36me3)$stats[2], middle = mean(H3K36me3), upper = boxplot.stats(H3K36me3)$stats[4], ymax = boxplot.stats(H3K36me3)$stats[5])

(H3K36me3_introns_profile <- ggplot(H3K36me3_introns_stat, aes(x = Expression, group = group)) + 
   geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax, fill = Cell_type), stat = "identity", position = "dodge", outlier.shape = NA, width = 0.8) + 
   facet_grid(geneRPKM ~ ., scales = "free") + 
   # ggtitle("H3K36me3 signal for introns") + 
   xlab("intron group") + 
   ylab("Average H3K36me3 signal") + 
   theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 15, color = "black"), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 15, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
ggsave(H3K36me3_introns_profile, file = "H3K36me3_introns_profile.pdf", width = 9, height = 8)

# t.test(H3K36me3_introns[H3K36me3_introns$group == "lum.lum-specific.gene RPKM < 1", "H3K36me3"], H3K36me3_introns[H3K36me3_introns$group == "myo.lum-specific.gene RPKM < 1", "H3K36me3"])
# t.test(H3K36me3_introns[H3K36me3_introns$group == "lum.myo-specific.gene RPKM < 1", "H3K36me3"], H3K36me3_introns[H3K36me3_introns$group == "myo.myo-specific.gene RPKM < 1", "H3K36me3"])
# t.test(H3K36me3_introns[H3K36me3_introns$group == "lum.lum-specific.gene RPKM 1-10", "H3K36me3"], H3K36me3_introns[H3K36me3_introns$group == "myo.lum-specific.gene RPKM 1-10", "H3K36me3"])
# t.test(H3K36me3_introns[H3K36me3_introns$group == "lum.myo-specific.gene RPKM 1-10", "H3K36me3"], H3K36me3_introns[H3K36me3_introns$group == "myo.myo-specific.gene RPKM 1-10", "H3K36me3"])
# t.test(H3K36me3_introns[H3K36me3_introns$group == "lum.retained_in_both.gene RPKM < 1", "H3K36me3"], H3K36me3_introns[H3K36me3_introns$group == "myo.retained_in_both.gene RPKM < 1", "H3K36me3"])
# t.test(H3K36me3_introns[H3K36me3_introns$group == "lum.not_retained.gene RPKM < 1", "H3K36me3"], H3K36me3_introns[H3K36me3_introns$group == "myo.not_retained.gene RPKM < 1", "H3K36me3"])
# t.test(H3K36me3_introns[H3K36me3_introns$group == "lum.retained_in_both.gene RPKM 1-10", "H3K36me3"], H3K36me3_introns[H3K36me3_introns$group == "myo.retained_in_both.gene RPKM 1-10", "H3K36me3"])
# t.test(H3K36me3_introns[H3K36me3_introns$group == "lum.not_retained.gene RPKM 1-10", "H3K36me3"], H3K36me3_introns[H3K36me3_introns$group == "myo.not_retained.gene RPKM 1-10", "H3K36me3"])

save(both, neither, lum_specific, myo_specific, introns_1, introns_1_10, introns_10_100, introns_100, H3K36me3_introns_stat, WGBS_CpG, 
     WGBS_boundaries, CpG_boundaries, MeDIP_boundaries, H3K4me3_boundaries, H3K4me1_boundaries, H3K9me3_boundaries, H3K27me3_boundaries, H3K36me3_introns, file = "intronProfile.Rdata")

