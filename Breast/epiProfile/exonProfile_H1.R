# Epigenetic profile at exon boundaries
setwd("~/快盘/REMC/epiProfile/")
library(ggplot2)
#load("exonProfile.Rdata")
# get H1-specific, myo-specific, expressed in both cell types exons  
isoform <- read.delim("myo_84_H1_r1a_isoform.txt", head = T, as.is = T)
isoform$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform$V1)
isoform$coord  <- gsub("chr", "", isoform$coord)
isoform$coord  <- gsub(":", "_", isoform$coord)
isoform$coord  <- gsub("-", "_", isoform$coord)
isoform$coord  <- paste0(isoform$id, "_", isoform$coord)
H1_specific <- isoform[isoform$V2 < isoform$V3, "coord"]
myo_specific <- isoform[isoform$V2 > isoform$V3, "coord"]
H1 <- read.delim("~/REMC/epiProfile/H1_r1a.G.exn.A.rpkm", head = F, as.is = T)
H1$coord <- gsub("<[-]*1", "", H1$V1)
H1$coord  <- gsub("chr", "", H1$coord)
H1$coord  <- gsub(":", "_", H1$coord)
H1$coord  <- gsub("-", "_", H1$coord)
H1$coord  <- paste0(H1$V2, "_", H1$coord)
myo084 <- read.delim("~/REMC/epiProfile/A17919.G.exn.A.rpkm", head = F, as.is = T)
myo084$coord <- gsub("<[-]*1", "", myo084$V1)
myo084$coord  <- gsub("chr", "", myo084$coord)
myo084$coord  <- gsub(":", "_", myo084$coord)
myo084$coord  <- gsub("-", "_", myo084$coord)
myo084$coord  <- paste0(myo084$V2, "_", myo084$coord)
both <- intersect(H1[H1$V4 > 0.1, "coord"], myo084[myo084$V4 > 0.1, "coord"])
neither <- setdiff(H1$coord, c(H1_specific, myo_specific, both))
rm(H1, myo084, isoform)

######################################################################################################
# WGBS profile @ exon boundaries
H1_bismark_3p <- read.table("~/REMC/epiProfile/exons3p_200/H1_bismark.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_bismark_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM045_bismark.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
H1_bismark_5p <- read.table("~/REMC/epiProfile/exons5p_200/H1_bismark.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_bismark_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM045_bismark.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
H1_bismark_3p <- H1_bismark_3p/10
H1_bismark_5p <- H1_bismark_5p/10
WGBS_3p <- data.frame(Cell_type = rep(c("H1", "myo"), each = 20*4), Expression = rep(c(rep(c("H1-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), WGBS = -1)
WGBS_3p$group <- interaction(WGBS_3p$Cell_type, WGBS_3p$Expression)
WGBS_3p[WGBS_3p$group == "H1.H1-specific", "WGBS"] <- colMeans(H1_bismark_3p[H1_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "H1.myo-specific", "WGBS"] <- colMeans(H1_bismark_3p[myo_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "H1.expressed_in_both", "WGBS"] <- colMeans(H1_bismark_3p[both,], na.rm = T)
WGBS_3p[WGBS_3p$group == "H1.not_expressed", "WGBS"] <- colMeans(H1_bismark_3p[neither,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.H1-specific", "WGBS"] <- colMeans(myo_bismark_3p[H1_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.myo-specific", "WGBS"] <- colMeans(myo_bismark_3p[myo_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.expressed_in_both", "WGBS"] <- colMeans(myo_bismark_3p[both,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.not_expressed", "WGBS"] <- colMeans(myo_bismark_3p[neither,], na.rm = T)
WGBS_5p <- data.frame(Cell_type = rep(c("H1", "myo"), each = 20*4), Expression = rep(c(rep(c("H1-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), WGBS = -1)
WGBS_5p$group <- interaction(WGBS_5p$Cell_type, WGBS_5p$Expression)
WGBS_5p[WGBS_5p$group == "H1.H1-specific", "WGBS"] <- colMeans(H1_bismark_5p[H1_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "H1.myo-specific", "WGBS"] <- colMeans(H1_bismark_5p[myo_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "H1.expressed_in_both", "WGBS"] <- colMeans(H1_bismark_5p[both,], na.rm = T)
WGBS_5p[WGBS_5p$group == "H1.not_expressed", "WGBS"] <- colMeans(H1_bismark_5p[neither,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.H1-specific", "WGBS"] <- colMeans(myo_bismark_5p[H1_specific,], na.rm = T)
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
ggsave(WGBS_boundaries_profile, file = "WGBS_boundaries_profile_H1.pdf")

# add CpG content track 
CpG_content_3p <- read.table("~/REMC/epiProfile/exons3p_200/CpG.hg19v65_exons_for_genes.3prime_200.unique", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
CpG_content_5p <- read.table("~/REMC/epiProfile/exons5p_200/CpG.hg19v65_exons_for_genes.5prime_200.unique", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
CpG_3p <- data.frame(Expression = c(rep(c("H1-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), Position = rep(seq(-190, 190, by = 20), times = 4), CpG = -1)
CpG_3p[CpG_3p$Expression == "H1-specific", "CpG"] <- colMeans(CpG_content_3p[H1_specific,], na.rm = T)
CpG_3p[CpG_3p$Expression == "myo-specific", "CpG"] <- colMeans(CpG_content_3p[myo_specific,], na.rm = T)
CpG_3p[CpG_3p$Expression == "expressed_in_both", "CpG"] <- colMeans(CpG_content_3p[both,], na.rm = T)
CpG_3p[CpG_3p$Expression == "not_expressed", "CpG"] <- colMeans(CpG_content_3p[neither,], na.rm = T)
CpG_5p <- data.frame(Expression = c(rep(c("H1-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), Position = rep(seq(-190, 190, by = 20), times = 4), CpG = -1)
CpG_5p[CpG_5p$Expression == "H1-specific", "CpG"] <- colMeans(CpG_content_5p[H1_specific,], na.rm = T)
CpG_5p[CpG_5p$Expression == "myo-specific", "CpG"] <- colMeans(CpG_content_5p[myo_specific,], na.rm = T)
CpG_5p[CpG_5p$Expression == "expressed_in_both", "CpG"] <- colMeans(CpG_content_5p[both,], na.rm = T)
CpG_5p[CpG_5p$Expression == "not_expressed", "CpG"] <- colMeans(CpG_content_5p[neither,], na.rm = T)
CpG_boundaries <- data.frame(rbind(CpG_3p, CpG_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(CpG_3p)), levels = c("5-prime", "3-prime")))
(CpG_boundaries_profile <- ggplot(CpG_boundaries, aes(x = Position, y = CpG, group = Expression)) + 
   geom_line(aes(color = Expression), size = 1.5) + 
   geom_point(aes(color = Expression), size = 4) + 
   facet_wrap(~ End) + 
   ggtitle("CpG profile around exon boundaries") + 
   ylab("Average percentage of CpG") + 
   theme_bw())
ggsave(CpG_boundaries_profile, file = "CpG_boundaries_profile.pdf", height = 4)

WGBS_CpG <- data.frame(data = c(rep("WGBS", nrow(WGBS_boundaries)), rep("CpG", nrow(CpG_boundaries))), 
                       Cell_type = c(as.character(WGBS_boundaries$Cell_type), rep("H1", nrow(CpG_boundaries))), 
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
library(grid) 
gt <- ggplot_gtable(ggplot_build(WGBS_CpG_profile)) 
# gt$layout 
gt$heights[[4]] <- unit(2.5, "null") 
pdf("WGBS_CpG_profile_H1.pdf")
grid.draw(gt) 
dev.off()
