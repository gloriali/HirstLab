# Epigenetic profile at exon boundaries
setwd("~/快盘/REMC/epiProfile/")
source("~/HirstLab/Pipeline/epiProfile.R")
library(ggplot2)
load("exonProfile.Rdata")
lum084 <- read.delim("~/REMC/epiProfile/A17918.G.exn.A.rpkm", head = F, as.is = T)
lum084$coord <- gsub("<[-]*1", "", lum084$V1)
lum084$coord  <- gsub("chr", "", lum084$coord)
lum084$coord  <- gsub(":", "_", lum084$coord)
lum084$coord  <- gsub("-", "_", lum084$coord)
lum084$coord  <- paste0(lum084$V2, "_", lum084$coord)
rownames(lum084) <- lum084$coord
myo084 <- read.delim("~/REMC/epiProfile/A17919.G.exn.A.rpkm", head = F, as.is = T)
myo084$coord <- gsub("<[-]*1", "", myo084$V1)
myo084$coord  <- gsub("chr", "", myo084$coord)
myo084$coord  <- gsub(":", "_", myo084$coord)
myo084$coord  <- gsub("-", "_", myo084$coord)
myo084$coord  <- paste0(myo084$V2, "_", myo084$coord)
rownames(myo084) <- myo084$coord
lum084gene <- read.delim("~/REMC/gene/A17918.G.A.rpkm.pc", head = F, as.is = T)
rownames(lum084gene) <- lum084gene$V1
myo084gene <- read.delim("~/REMC/gene/A17919.G.A.rpkm.pc", head = F, as.is = T)
rownames(myo084gene) <- myo084gene$V1
geneRPKM <- data.frame(gene = lum084gene$V1, RPKM = (lum084gene$V3 + myo084gene$V3)/2)
rownames(geneRPKM) <- geneRPKM$gene
exons_geneRPKM <- data.frame(exon = lum084$coord)
exons_geneRPKM$gene <- gsub("_[0-9XY]+_[0-9_]+", "", exons_geneRPKM$exon)
exons_geneRPKM$geneRPKM <- geneRPKM[exons_geneRPKM$gene, "RPKM"]
exons_geneRPKM$lum084gene <- lum084gene[exons_geneRPKM$gene, "V3"]
exons_geneRPKM$myo084gene <- myo084gene[exons_geneRPKM$gene, "V3"]
exons_geneRPKM$lum084exon <- lum084[as.character(exons_geneRPKM$exon), "V4"]
exons_geneRPKM$myo084exon <- myo084[as.character(exons_geneRPKM$exon), "V4"]
exons_geneRPKM <- na.omit(exons_geneRPKM)
exons_geneRPKM$exon <- as.character(exons_geneRPKM$exon)
lum_expressed <- as.character(exons_geneRPKM[lum084gene > 0.1, "exon"])
myo_expressed <- as.character(exons_geneRPKM[myo084gene > 0.1, "exon"])

# get lum-specific, myo-specific, expressed in both cell types exons  
isoform <- read.delim("lum084_myo084_isoform_all.txt", head = T, as.is = T)
isoform$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform$V1)
isoform$coord  <- gsub("chr", "", isoform$coord)
isoform$coord  <- gsub(":", "_", isoform$coord)
isoform$coord  <- gsub("-", "_", isoform$coord)
isoform$coord  <- paste0(isoform$id, "_", isoform$coord)
lum_specific <- isoform[isoform$V2 < isoform$V3, "coord"]  # isoform exons
myo_specific <- isoform[isoform$V2 > isoform$V3, "coord"]  # isoform exons
lum_specific <- intersect(lum_specific, intersect(lum_expressed, myo_expressed))
myo_specific <- intersect(myo_specific, intersect(lum_expressed, myo_expressed))
both <- setdiff(exons_geneRPKM[exons_geneRPKM$lum084exon > 0.1 & exons_geneRPKM$myo084exon > 0.1, "exon"], c(lum_specific, myo_specific))
neither <- setdiff(lum084$coord, c(lum_specific, myo_specific, both))
lum_exon_0 <- exons_geneRPKM[exons_geneRPKM$lum084exon == 0 & exons_geneRPKM$lum084gene > 0.1, "exon"]
myo_exon_0 <- exons_geneRPKM[exons_geneRPKM$myo084exon == 0 & exons_geneRPKM$myo084gene > 0.1, "exon"]
lum_gene_0 <- exons_geneRPKM[exons_geneRPKM$lum084gene == 0, "exon"]
myo_gene_0 <- exons_geneRPKM[exons_geneRPKM$myo084gene == 0, "exon"]
# exons_geneRPKM[exons_geneRPKM$lum084gene == 0 & exons_geneRPKM$myo084gene == 0, "express"] <- "gene_not_expressed"
# both <- intersect(lum084[lum084$V4 > 0.1, "coord"], myo084[myo084$V4 > 0.1, "coord"])     # both exon RPKM > 0.1
# neither <- intersect(lum084[lum084$V4 < 0.1, "coord"], myo084[myo084$V4 < 0.1, "coord"])  # both exon RPKM < 0.1
# rm(lum084, myo084, isoform)


# more stringent cell-type specific exons
# lum_specific <- intersect(lum_specific, myo084[myo084$V4 == 0, "coord"])
# myo_specific <- intersect(myo_specific, lum084[lum084$V4 == 0, "coord"])

# # exons with restricted gene RPKM to isoform exons gene RPKM range
# min <- min(exons_geneRPKM[exons_geneRPKM$express == "lum-specific" | exons_geneRPKM$express == "myo-specific", "geneRPKM"])
# max <- max(exons_geneRPKM[exons_geneRPKM$express == "lum-specific" | exons_geneRPKM$express == "myo-specific", "geneRPKM"])
# both <- as.character(exons_geneRPKM[exons_geneRPKM$express == "expressed_in_both" & exons_geneRPKM$geneRPKM > min & exons_geneRPKM$geneRPKM < max, "exon"])
# neither <- as.character(exons_geneRPKM[exons_geneRPKM$express == "exon_not_expressed" & exons_geneRPKM$geneRPKM > min & exons_geneRPKM$geneRPKM < max, "exon"])

# restrict exon RPKM to isoform exon RPKM range for both
# max <- max(c(lum084[c(lum_specific, myo_specific), "V4"], myo084[c(lum_specific, myo_specific), "V4"]))
# both <- intersect(both, intersect(lum084[lum084$V4 < max, "coord"], myo084[myo084$V4 < max, "coord"]))

# divide exons based on gene RPKM groups
exons_1 <- as.character(exons_geneRPKM[exons_geneRPKM$geneRPKM <= 1, "exon"])
exons_1_10 <- as.character(exons_geneRPKM[exons_geneRPKM$geneRPKM > 1 & exons_geneRPKM$geneRPKM <= 10, "exon"])
exons_10_100 <- as.character(exons_geneRPKM[exons_geneRPKM$geneRPKM > 10 & exons_geneRPKM$geneRPKM <= 100, "exon"])
exons_100 <- as.character(exons_geneRPKM[exons_geneRPKM$geneRPKM > 100, "exon"])

######################################################################################################
# WGBS profile @ exon boundaries
lum_WGBS_3p <- read.table("~/REMC/epiProfile/exons3p_200/lumRM066_WGBS.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_WGBS_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM045_WGBS.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
lum_WGBS_5p <- read.table("~/REMC/epiProfile/exons5p_200/lumRM066_WGBS.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_WGBS_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM045_WGBS.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
WGBS_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), WGBS = -1)
WGBS_3p$group <- interaction(WGBS_3p$Cell_type, WGBS_3p$Expression)
WGBS_3p[WGBS_3p$group == "lum.lum-specific", "WGBS"] <- colMeans(lum_WGBS_3p[lum_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "lum.myo-specific", "WGBS"] <- colMeans(lum_WGBS_3p[myo_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "lum.expressed_in_both", "WGBS"] <- colMeans(lum_WGBS_3p[both,], na.rm = T)
WGBS_3p[WGBS_3p$group == "lum.not_expressed", "WGBS"] <- colMeans(lum_WGBS_3p[neither,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.lum-specific", "WGBS"] <- colMeans(myo_WGBS_3p[lum_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.myo-specific", "WGBS"] <- colMeans(myo_WGBS_3p[myo_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.expressed_in_both", "WGBS"] <- colMeans(myo_WGBS_3p[both,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.not_expressed", "WGBS"] <- colMeans(myo_WGBS_3p[neither,], na.rm = T)
WGBS_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), WGBS = -1)
WGBS_5p$group <- interaction(WGBS_5p$Cell_type, WGBS_5p$Expression)
WGBS_5p[WGBS_5p$group == "lum.lum-specific", "WGBS"] <- colMeans(lum_WGBS_5p[lum_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "lum.myo-specific", "WGBS"] <- colMeans(lum_WGBS_5p[myo_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "lum.expressed_in_both", "WGBS"] <- colMeans(lum_WGBS_5p[both,], na.rm = T)
WGBS_5p[WGBS_5p$group == "lum.not_expressed", "WGBS"] <- colMeans(lum_WGBS_5p[neither,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.lum-specific", "WGBS"] <- colMeans(myo_WGBS_5p[lum_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.myo-specific", "WGBS"] <- colMeans(myo_WGBS_5p[myo_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.expressed_in_both", "WGBS"] <- colMeans(myo_WGBS_5p[both,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.not_expressed", "WGBS"] <- colMeans(myo_WGBS_5p[neither,], na.rm = T)
WGBS_boundaries <- data.frame(rbind(WGBS_3p, WGBS_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(WGBS_3p)), levels = c("5-prime", "3-prime")))
(WGBS_boundaries_profile <- ggplot(WGBS_boundaries, aes(x = Position, y = WGBS, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type), size = 1.5) + 
   geom_point(aes(color = Expression, shape = Cell_type), size = 4) + 
   facet_wrap(~ End) + 
   ggtitle("WGBS profile around exon boundaries") + 
   ylab("Average DNA methylation level") + 
   theme_bw())
ggsave(WGBS_boundaries_profile, file = "WGBS_boundaries_profile.pdf")

# expressed genes only
WGBS_boundaries <- epiProfile(mark = "WGBS", cell1 = "lum", cell2 = "myo", donor1 = "RM066", donor2 = "RM045", dirIn = "~/REMC/epiProfile/", both = both, neither = neither, cell1_specific = lum_specific, cell2_specific = myo_specific)

# add GC content track 
GC_content_3p <- read.table("~/REMC/epiProfile/exons3p_200/GC.hg19v65_exons_for_genes.3prime_200.unique", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
GC_content_5p <- read.table("~/REMC/epiProfile/exons5p_200/GC.hg19v65_exons_for_genes.5prime_200.unique", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
GC_3p <- data.frame(Expression = c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), Position = rep(seq(-190, 190, by = 20), times = 4), GC = -1)
GC_3p[GC_3p$Expression == "lum-specific", "GC"] <- colMeans(GC_content_3p[lum_specific,], na.rm = T)
GC_3p[GC_3p$Expression == "myo-specific", "GC"] <- colMeans(GC_content_3p[myo_specific,], na.rm = T)
GC_3p[GC_3p$Expression == "expressed_in_both", "GC"] <- colMeans(GC_content_3p[both,], na.rm = T)
GC_3p[GC_3p$Expression == "not_expressed", "GC"] <- colMeans(GC_content_3p[neither,], na.rm = T)
GC_5p <- data.frame(Expression = c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), Position = rep(seq(-190, 190, by = 20), times = 4), GC = -1)
GC_5p[GC_5p$Expression == "lum-specific", "GC"] <- colMeans(GC_content_5p[lum_specific,], na.rm = T)
GC_5p[GC_5p$Expression == "myo-specific", "GC"] <- colMeans(GC_content_5p[myo_specific,], na.rm = T)
GC_5p[GC_5p$Expression == "expressed_in_both", "GC"] <- colMeans(GC_content_5p[both,], na.rm = T)
GC_5p[GC_5p$Expression == "not_expressed", "GC"] <- colMeans(GC_content_5p[neither,], na.rm = T)
GC_boundaries <- data.frame(rbind(GC_3p, GC_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(GC_3p)), levels = c("5-prime", "3-prime")))
(GC_boundaries_profile <- ggplot(GC_boundaries, aes(x = Position, y = GC, group = Expression)) + 
   geom_line(aes(color = Expression), size = 1.5) + 
   geom_point(aes(color = Expression), size = 4) + 
   facet_wrap(~ End) + 
   ggtitle("GC content profile around exon boundaries") + 
   ylab("Average percentage of GC") + 
   theme_bw())
ggsave(GC_boundaries_profile, file = "GC_boundaries_profile.pdf", height = 4)

WGBS_GC <- data.frame(data = c(rep("WGBS", nrow(WGBS_boundaries)), rep("GC content", nrow(GC_boundaries))), 
                      Cell_type = c(as.character(WGBS_boundaries$Cell_type), rep("lum", nrow(GC_boundaries))), 
                      Expression = c(as.character(WGBS_boundaries$Expression), as.character(GC_boundaries$Expression)), 
                      Position = c(WGBS_boundaries$Position, GC_boundaries$Position), 
                      value = c(WGBS_boundaries$WGBS, GC_boundaries$GC), 
                      End = c(as.character(WGBS_boundaries$End), as.character(GC_boundaries$End)))
WGBS_GC$data <- factor(WGBS_GC$data, levels = c("WGBS", "GC content"))
WGBS_GC$End <- factor(WGBS_GC$End, levels = c("5-prime", "3-prime"))
WGBS_GC$group <- interaction(WGBS_GC$Cell_type, WGBS_GC$Expression)
(WGBS_GC_profile <- ggplot(WGBS_GC, aes(x = Position, y = value, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type)) + 
   geom_point(aes(color = Expression, shape = Cell_type)) + 
   facet_grid(data ~ End, scales = "free_y") + 
   ylab("Average DNA methylation / GC content") + 
   theme(axis.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 12, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
library(grid) 
gt <- ggplot_gtable(ggplot_build(WGBS_GC_profile)) 
# gt$layout 
gt$heights[[4]] <- unit(2, "null") 
gt$width[[3]] <- unit(1.1, "null") 
pdf("WGBS_GC_profile.pdf")
grid.newpage()
grid.draw(gt) 
dev.off()

# add CpG content track 
CpG_content_3p <- read.table("~/REMC/epiProfile/exons3p_200/CpG.hg19v65_exons_for_genes.3prime_200.unique", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
CpG_content_5p <- read.table("~/REMC/epiProfile/exons5p_200/CpG.hg19v65_exons_for_genes.5prime_200.unique", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
CpG_3p <- data.frame(Expression = c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), Position = rep(seq(-190, 190, by = 20), times = 4), CpG = -1)
CpG_3p[CpG_3p$Expression == "lum-specific", "CpG"] <- colMeans(CpG_content_3p[lum_specific,], na.rm = T)
CpG_3p[CpG_3p$Expression == "myo-specific", "CpG"] <- colMeans(CpG_content_3p[myo_specific,], na.rm = T)
CpG_3p[CpG_3p$Expression == "expressed_in_both", "CpG"] <- colMeans(CpG_content_3p[both,], na.rm = T)
CpG_3p[CpG_3p$Expression == "not_expressed", "CpG"] <- colMeans(CpG_content_3p[neither,], na.rm = T)
CpG_5p <- data.frame(Expression = c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), Position = rep(seq(-190, 190, by = 20), times = 4), CpG = -1)
CpG_5p[CpG_5p$Expression == "lum-specific", "CpG"] <- colMeans(CpG_content_5p[lum_specific,], na.rm = T)
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
   ylab("") + 
   scale_color_manual(values = c("lum-specific" = rgb(200,50,0, maxColorValue = 255), "myo-specific" = rgb(50,200,50, maxColorValue = 255), "expressed_in_both" = "purple", "not_expressed" = "blue")) + 
   theme(panel.margin = unit(0.75, "lines"), axis.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 12, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
library(grid) 
gt <- ggplot_gtable(ggplot_build(WGBS_CpG_profile)) 
# gt$layout 
gt$heights[[4]] <- unit(3, "null") 
# grid.newpage()
pdf("WGBS_CpG_profile.pdf", width = 9)
grid.draw(gt) 
grid.text("Average DNA methylation", x = unit(0.015, "npc"), y = unit(0.65, "npc"), rot = 90)
grid.text("No. of CpGs", x = unit(0.015, "npc"), y = unit(0.18, "npc"), rot = 90)
dev.off()
######################################################################################################
# MeDIP profile @ exon boundaries
lum_MeDIP_3p <- read.table("~/REMC/epiProfile/exons3p_200/lumRM035_MeDIP.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_MeDIP_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM035_MeDIP.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
lum_MeDIP_5p <- read.table("~/REMC/epiProfile/exons5p_200/lumRM035_MeDIP.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_MeDIP_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM035_MeDIP.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
lum_MeDIP_3p <- sum(myo_MeDIP_3p)/sum(lum_MeDIP_3p) * lum_MeDIP_3p
lum_MeDIP_5p <- sum(myo_MeDIP_5p)/sum(lum_MeDIP_5p) * lum_MeDIP_5p
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

# expressed genes only
MeDIP_boundaries <- epiProfile(mark = "MeDIP", cell1 = "lum", cell2 = "myo", donor1 = "RM035", donor2 = "RM035", dirIn = "~/REMC/epiProfile/", both = both, neither = neither, cell1_specific = lum_specific, cell2_specific = myo_specific)

######################################################################################################
# H3K4me3 profile @ exon boundaries (no lum libraries available)
myo_H3K4me3_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM080_H3K4me3.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_H3K4me3_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM080_H3K4me3.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
H3K4me3_3p <- data.frame(Expression = rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K4me3 = -1)
H3K4me3_3p[H3K4me3_3p$Expression == "lum-specific", "H3K4me3"] <- colMeans(myo_WGBS_3p[lum_specific,], na.rm = T)
H3K4me3_3p[H3K4me3_3p$Expression == "myo-specific", "H3K4me3"] <- colMeans(myo_WGBS_3p[myo_specific,], na.rm = T)
H3K4me3_3p[H3K4me3_3p$Expression == "expressed_in_both", "H3K4me3"] <- colMeans(myo_WGBS_3p[both,], na.rm = T)
H3K4me3_3p[H3K4me3_3p$Expression == "not_expressed", "H3K4me3"] <- colMeans(myo_WGBS_3p[neither,], na.rm = T)
H3K4me3_5p <- data.frame(Expression = rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K4me3 = -1)
H3K4me3_5p[H3K4me3_5p$Expression == "lum-specific", "H3K4me3"] <- colMeans(myo_WGBS_5p[lum_specific,], na.rm = T)
H3K4me3_5p[H3K4me3_5p$Expression == "myo-specific", "H3K4me3"] <- colMeans(myo_WGBS_5p[myo_specific,], na.rm = T)
H3K4me3_5p[H3K4me3_5p$Expression == "expressed_in_both", "H3K4me3"] <- colMeans(myo_WGBS_5p[both,], na.rm = T)
H3K4me3_5p[H3K4me3_5p$Expression == "not_expressed", "H3K4me3"] <- colMeans(myo_WGBS_5p[neither,], na.rm = T)
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
lum_H3K4me1_3p <- sum(myo_H3K4me1_3p)/sum(lum_H3K4me1_3p) * lum_H3K4me1_3p
lum_H3K4me1_5p <- sum(myo_H3K4me1_5p)/sum(lum_H3K4me1_5p) * lum_H3K4me1_5p
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
lum_H3K9me3_3p <- sum(myo_H3K9me3_3p)/sum(lum_H3K9me3_3p) * lum_H3K9me3_3p
lum_H3K9me3_5p <- sum(myo_H3K9me3_5p)/sum(lum_H3K9me3_5p) * lum_H3K9me3_5p
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
lum_H3K27me3_3p <- sum(myo_H3K27me3_3p)/sum(lum_H3K27me3_3p) * lum_H3K27me3_3p
lum_H3K27me3_5p <- sum(myo_H3K27me3_5p)/sum(lum_H3K27me3_5p) * lum_H3K27me3_5p
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
lum_H3K36me3_3p <- sum(myo_H3K36me3_3p)/sum(lum_H3K36me3_3p) * lum_H3K36me3_3p
lum_H3K36me3_5p <- sum(myo_H3K36me3_5p)/sum(lum_H3K36me3_5p) * lum_H3K36me3_5p
H3K36me3_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K36me3 = -1)
H3K36me3_3p <- data.frame(H3K36me3_3p, geneRPKM = rep(c("<1", "1-10", "10-100", ">100"), each = nrow(H3K36me3_3p)))
H3K36me3_3p$group <- interaction(interaction(H3K36me3_3p$Cell_type, H3K36me3_3p$Expression), H3K36me3_3p$geneRPKM)
H3K36me3_3p[H3K36me3_3p$group == "lum.lum-specific.<1", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(lum_specific, exons_1),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.myo-specific.<1", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(myo_specific, exons_1),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.expressed_in_both.<1", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(both, exons_1),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.not_expressed.<1", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(neither, exons_1),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.lum-specific.<1", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(lum_specific, exons_1),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.myo-specific.<1", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(myo_specific, exons_1),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.expressed_in_both.<1", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(both, exons_1),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.not_expressed.<1", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(neither, exons_1),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.lum-specific.1-10", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(lum_specific, exons_1_10),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.myo-specific.1-10", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(myo_specific, exons_1_10),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.expressed_in_both.1-10", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(both, exons_1_10),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.not_expressed.1-10", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(neither, exons_1_10),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.lum-specific.1-10", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(lum_specific, exons_1_10),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.myo-specific.1-10", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(myo_specific, exons_1_10),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.expressed_in_both.1-10", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(both, exons_1_10),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.not_expressed.1-10", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(neither, exons_1_10),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.lum-specific.10-100", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(lum_specific, exons_10_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.myo-specific.10-100", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(myo_specific, exons_10_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.expressed_in_both.10-100", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(both, exons_10_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.not_expressed.10-100", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(neither, exons_10_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.lum-specific.10-100", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(lum_specific, exons_10_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.myo-specific.10-100", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(myo_specific, exons_10_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.expressed_in_both.10-100", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(both, exons_10_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.not_expressed.10-100", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(neither, exons_10_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.lum-specific.>100", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(lum_specific, exons_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.myo-specific.>100", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(myo_specific, exons_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.expressed_in_both.>100", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(both, exons_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "lum.not_expressed.>100", "H3K36me3"] <- colMeans(lum_H3K36me3_3p[intersect(neither, exons_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.lum-specific.>100", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(lum_specific, exons_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.myo-specific.>100", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(myo_specific, exons_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.expressed_in_both.>100", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(both, exons_100),], na.rm = T)
H3K36me3_3p[H3K36me3_3p$group == "myo.not_expressed.>100", "H3K36me3"] <- colMeans(myo_H3K36me3_3p[intersect(neither, exons_100),], na.rm = T)

H3K36me3_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), H3K36me3 = -1)
H3K36me3_5p <- data.frame(H3K36me3_5p, geneRPKM = rep(c("<1", "1-10", "10-100", ">100"), each = nrow(H3K36me3_5p)))
H3K36me3_5p$group <- interaction(interaction(H3K36me3_5p$Cell_type, H3K36me3_5p$Expression), H3K36me3_5p$geneRPKM)
H3K36me3_5p[H3K36me3_5p$group == "lum.lum-specific.<1", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(lum_specific, exons_1),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.myo-specific.<1", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(myo_specific, exons_1),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.expressed_in_both.<1", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(both, exons_1),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.not_expressed.<1", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(neither, exons_1),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.lum-specific.<1", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(lum_specific, exons_1),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.myo-specific.<1", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(myo_specific, exons_1),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.expressed_in_both.<1", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(both, exons_1),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.not_expressed.<1", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(neither, exons_1),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.lum-specific.1-10", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(lum_specific, exons_1_10),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.myo-specific.1-10", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(myo_specific, exons_1_10),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.expressed_in_both.1-10", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(both, exons_1_10),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.not_expressed.1-10", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(neither, exons_1_10),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.lum-specific.1-10", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(lum_specific, exons_1_10),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.myo-specific.1-10", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(myo_specific, exons_1_10),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.expressed_in_both.1-10", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(both, exons_1_10),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.not_expressed.1-10", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(neither, exons_1_10),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.lum-specific.10-100", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(lum_specific, exons_10_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.myo-specific.10-100", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(myo_specific, exons_10_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.expressed_in_both.10-100", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(both, exons_10_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.not_expressed.10-100", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(neither, exons_10_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.lum-specific.10-100", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(lum_specific, exons_10_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.myo-specific.10-100", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(myo_specific, exons_10_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.expressed_in_both.10-100", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(both, exons_10_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.not_expressed.10-100", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(neither, exons_10_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.lum-specific.>100", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(lum_specific, exons_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.myo-specific.>100", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(myo_specific, exons_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.expressed_in_both.>100", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(both, exons_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "lum.not_expressed.>100", "H3K36me3"] <- colMeans(lum_H3K36me3_5p[intersect(neither, exons_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.lum-specific.>100", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(lum_specific, exons_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.myo-specific.>100", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(myo_specific, exons_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.expressed_in_both.>100", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(both, exons_100),], na.rm = T)
H3K36me3_5p[H3K36me3_5p$group == "myo.not_expressed.>100", "H3K36me3"] <- colMeans(myo_H3K36me3_5p[intersect(neither, exons_100),], na.rm = T)

H3K36me3_boundaries <- data.frame(rbind(H3K36me3_3p, H3K36me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K36me3_3p)), levels = c("5-prime", "3-prime")))
H3K36me3_boundaries <- droplevels(H3K36me3_boundaries[H3K36me3_boundaries$geneRPKM != ">100", ])
(H3K36me3_boundaries_profile <- ggplot(H3K36me3_boundaries, aes(x = Position, y = H3K36me3, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type), size = 1.5) + 
   geom_point(aes(color = Expression, shape = Cell_type), size = 4) + 
   facet_grid(geneRPKM ~ End, scales = "free_y") + 
   ggtitle("H3K36me3 profile around exon boundaries") + 
   ylab("Average H3K36me3 signal") + 
   theme_bw())
ggsave(H3K36me3_boundaries_profile, file = "H3K36me3_boundaries_profile.pdf")

######################################################################################################
# H3K36me3 signal for exons
lum_H3K36me3_exons <- read.delim("~/REMC/epiProfile/exons/hg19v65_exons_for_genes.lumRM080_H3K36me3.coverage", head = F, as.is = T)
myo_H3K36me3_exons <- read.delim("~/REMC/epiProfile/exons/hg19v65_exons_for_genes.myoRM080_H3K36me3.coverage", head = F, as.is = T)
# normalize signal
(norm <- sum(myo_H3K36me3_exons$V6)/sum(lum_H3K36me3_exons$V6))
H3K36me3_exons <- data.frame(id = c(lum_H3K36me3_exons$V4, myo_H3K36me3_exons$V4), Cell_type = c(rep("lum", nrow(lum_H3K36me3_exons)), rep("myo", nrow(myo_H3K36me3_exons))), Expression = NA, H3K36me3 = c(lum_H3K36me3_exons$V6 * norm, myo_H3K36me3_exons$V6))
H3K36me3_exons[H3K36me3_exons$id %in% both, "Expression"] <- "expressed_in_both"
# H3K36me3_exons[H3K36me3_exons$id %in% neither, "Expression"] <- "exon_not_expressed"
H3K36me3_exons[H3K36me3_exons$id %in% lum_exon_0 & as.character(H3K36me3_exons$Cell_type) == "lum", "Expression"] <- "not_expressed"
H3K36me3_exons[H3K36me3_exons$id %in% myo_exon_0 & as.character(H3K36me3_exons$Cell_type) == "myo", "Expression"] <- "not_expressed"
H3K36me3_exons[H3K36me3_exons$id %in% lum_specific, "Expression"] <- "lum-specific"
H3K36me3_exons[H3K36me3_exons$id %in% myo_specific, "Expression"] <- "myo-specific"
# H3K36me3_exons[H3K36me3_exons$id %in% not_expressed, "Expression"] <- "gene_not_expressed"
H3K36me3_exons[H3K36me3_exons$id %in% lum_gene_0 & as.character(H3K36me3_exons$Cell_type) == "lum", "Expression"] <- "gene_not_expressed"
H3K36me3_exons[H3K36me3_exons$id %in% myo_gene_0 & as.character(H3K36me3_exons$Cell_type) == "myo", "Expression"] <- "gene_not_expressed"
H3K36me3_exons <- droplevels(H3K36me3_exons[H3K36me3_exons$Expression != "gene_not_expressed", ])
H3K36me3_exons$Expression <- factor(H3K36me3_exons$Expression, levels = c("expressed_in_both", "lum-specific", "myo-specific", "not_expressed"))
# H3K36me3_exons[H3K36me3_exons$id %in% exons_1, "geneRPKM"] <- "gene RPKM < 1"
# H3K36me3_exons[H3K36me3_exons$id %in% exons_1_10, "geneRPKM"] <- "gene RPKM 1-10"
# H3K36me3_exons[H3K36me3_exons$id %in% exons_10_100, "geneRPKM"] <- "gene RPKM 10-100"
# H3K36me3_exons[H3K36me3_exons$id %in% exons_100, "geneRPKM"] <- "gene RPKM > 100"
# H3K36me3_exons$geneRPKM <- factor(H3K36me3_exons$geneRPKM)
H3K36me3_exons$utilize <- interaction(H3K36me3_exons$Cell_type, H3K36me3_exons$Expression)
H3K36me3_exons <- na.omit(H3K36me3_exons)
# H3K36me3_exons <- droplevels(na.omit(H3K36me3_exons[H3K36me3_exons$geneRPKM != "gene RPKM > 100" & H3K36me3_exons$geneRPKM != "gene RPKM 10-100", ]))
# H3K36me3_exons$group <- interaction(interaction(H3K36me3_exons$Cell_type, H3K36me3_exons$Expression), H3K36me3_exons$geneRPKM)
summary(H3K36me3_exons$utilize)
library(plyr)
H3K36me3_exons_stat <- ddply(H3K36me3_exons, ~ utilize, summarize, Cell_type = Cell_type[1], Expression = Expression[1], ymin = boxplot.stats(H3K36me3)$stats[1], lower = boxplot.stats(H3K36me3)$stats[2], middle = mean(H3K36me3), upper = boxplot.stats(H3K36me3)$stats[4], ymax = boxplot.stats(H3K36me3)$stats[5])

(H3K36me3_exons_profile <- ggplot(H3K36me3_exons_stat, aes(x = Expression, group = utilize)) + 
   geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax, fill = Cell_type), stat = "identity", position = "dodge", width = 0.8) + 
   # facet_grid(geneRPKM ~ ., scales = "free") + 
   # ggtitle("H3K36me3 signal for exons") + 
   xlab("Exon group") + 
   ylab("Average H3K36me3 signal") + 
   scale_fill_manual(values = c("lum" = rgb(200,50,0, maxColorValue = 255), "myo" = rgb(50,200,50, maxColorValue = 255))) + 
   theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 15, color = "black"), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 20, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
ggsave(H3K36me3_exons_profile, file = "H3K36me3_exons_profile.pdf", width = 9, height = 5)
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "lum.lum-specific" | H3K36me3_exons$utilize == "myo.myo-specific", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "myo.lum-specific" | H3K36me3_exons$utilize == "lum.myo-specific", "H3K36me3"], alternative = "greater")$p.value
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "lum.lum-specific", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "myo.lum-specific", "H3K36me3"])$p.value
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "myo.myo-specific", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "lum.myo-specific", "H3K36me3"])$p.value
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "lum.expressed_in_both", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "myo.expressed_in_both", "H3K36me3"])$p.value
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "lum.exon_not_expressed", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "myo.exon_not_expressed", "H3K36me3"])$p.value
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "lum.gene_not_expressed", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "myo.gene_not_expressed", "H3K36me3"])$p.value
# fold enrichment between utilized/un-utilized isoform exons
mean(H3K36me3_exons[H3K36me3_exons$utilize == "lum.lum-specific" | H3K36me3_exons$utilize == "myo.myo-specific", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.lum-specific" | H3K36me3_exons$utilize == "lum.myo-specific", "H3K36me3"])
mean(H3K36me3_exons[H3K36me3_exons$utilize == "lum.lum-specific", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.lum-specific", "H3K36me3"])
mean(H3K36me3_exons[H3K36me3_exons$utilize == "lum.myo-specific", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.myo-specific", "H3K36me3"])
mean(H3K36me3_exons[H3K36me3_exons$utilize == "lum.expressed_in_both", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.expressed_in_both", "H3K36me3"])
mean(H3K36me3_exons[H3K36me3_exons$utilize == "lum.exon_not_expressed", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.exon_not_expressed", "H3K36me3"])
mean(H3K36me3_exons[H3K36me3_exons$utilize == "lum.gene_not_expressed", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.gene_not_expressed", "H3K36me3"])

t.test(H3K36me3_exons[H3K36me3_exons$group == "lum.lum-specific.gene RPKM < 1", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$group == "myo.lum-specific.gene RPKM < 1", "H3K36me3"])
t.test(H3K36me3_exons[H3K36me3_exons$group == "lum.myo-specific.gene RPKM < 1", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$group == "myo.myo-specific.gene RPKM < 1", "H3K36me3"])
t.test(H3K36me3_exons[H3K36me3_exons$group == "lum.lum-specific.gene RPKM 1-10", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$group == "myo.lum-specific.gene RPKM 1-10", "H3K36me3"])
t.test(H3K36me3_exons[H3K36me3_exons$group == "lum.myo-specific.gene RPKM 1-10", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$group == "myo.myo-specific.gene RPKM 1-10", "H3K36me3"])
t.test(H3K36me3_exons[H3K36me3_exons$group == "lum.expressed_in_both.gene RPKM < 1", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$group == "myo.expressed_in_both.gene RPKM < 1", "H3K36me3"])
t.test(H3K36me3_exons[H3K36me3_exons$group == "lum.not_expressed.gene RPKM < 1", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$group == "myo.not_expressed.gene RPKM < 1", "H3K36me3"])
t.test(H3K36me3_exons[H3K36me3_exons$group == "lum.expressed_in_both.gene RPKM 1-10", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$group == "myo.expressed_in_both.gene RPKM 1-10", "H3K36me3"])
t.test(H3K36me3_exons[H3K36me3_exons$group == "lum.not_expressed.gene RPKM 1-10", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$group == "myo.not_expressed.gene RPKM 1-10", "H3K36me3"])

save(both, neither, not_expressed,lum_specific, myo_specific, exons_1, exons_1_10, exons_10_100, exons_100, H3K36me3_exons_stat, WGBS_CpG, CpG_boundaries, 
     WGBS_boundaries, MeDIP_boundaries, H3K4me3_boundaries, H3K4me1_boundaries, H3K9me3_boundaries, H3K27me3_boundaries, H3K36me3_boundaries, GC_boundaries, H3K36me3_exons, WGBS_GC, file = "exonProfile.Rdata")

# ################## Unable to complete: computational time too long ################################
# # WGBS along exons
# # /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# setwd("~/REMC/epiProfile/exons/")
# collapse <- function(gene, raw){
#   raw_data <- na.omit(as.numeric(raw[raw$X1 == gene, 2:ncol(raw)]))
#   length <- length(raw_data)
#   output <- vector(length = 20)
#   for(i in 1:19){
#     output[i] <- mean(raw_data[(as.integer((i-1)*length/20)+1):(as.integer(i*length/20)+1)])
#   }
#   output[20] <- mean(raw_data[(as.integer((i-1)*length/20)+1):as.integer(i*length/20)])
#   return(output)
# }
# 
# # less ~/REMC/epiProfile/exons/myoRM045_WGBS.hg19v65_exons_for_genes.profile | awk -F' ' '{if(NF>max)max=NF}END{print max}'
# # less ~/REMC/epiProfile/exons/lumRM066_WGBS.hg19v65_exons_for_genes.profile | awk -F' ' '{if(NF>max)max=NF}END{print max}'
# myo_WGBS_exons_raw <- read.table("~/REMC/epiProfile/exons/myoRM045_WGBS.hg19v65_exons_for_genes.profile", sep = " ", head = F, as.is = T, fill = T, col.names = as.character(1:1462))
# myo_WGBS_exons <- t(sapply(unique(myo_WGBS_exons_raw$X1), collapse, raw = myo_WGBS_exons_raw))
# save(myo_WGBS_exons, "myo_WGBS_exons.Rdata")
# lum_WGBS_exons_raw <- read.table("~/REMC/epiProfile/exons/lumRM066_WGBS.hg19v65_exons_for_genes.profile", sep = " ", head = F, as.is = T, fill = T, col.names = as.character(1:1462))
# lum_WGBS_exons <- t(sapply(unique(lum_WGBS_exons_raw$X1), collapse, raw = lum_WGBS_exons_raw))
# save(lum_WGBS_exons, "lum_WGBS_exons.Rdata")
# WGBS_exons <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(5, 100, by = 5), times = 8), WGBS = -1)
# WGBS_exons$group <- interaction(WGBS_exons$Cell_type, WGBS_exons$Expression)
# WGBS_exons[WGBS_exons$group == "lum.lum-specific", "WGBS"] <- colMeans(lum_WGBS_exons[lum_specific,], na.rm = T)
# WGBS_exons[WGBS_exons$group == "lum.myo-specific", "WGBS"] <- colMeans(lum_WGBS_exons[myo_specific,], na.rm = T)
# WGBS_exons[WGBS_exons$group == "lum.expressed_in_both", "WGBS"] <- colMeans(lum_WGBS_exons[both,], na.rm = T)
# WGBS_exons[WGBS_exons$group == "lum.not_expressed", "WGBS"] <- colMeans(lum_WGBS_exons[neither,], na.rm = T)
# WGBS_exons[WGBS_exons$group == "myo.lum-specific", "WGBS"] <- colMeans(myo_WGBS_exons[lum_specific,], na.rm = T)
# WGBS_exons[WGBS_exons$group == "myo.myo-specific", "WGBS"] <- colMeans(myo_WGBS_exons[myo_specific,], na.rm = T)
# WGBS_exons[WGBS_exons$group == "myo.expressed_in_both", "WGBS"] <- colMeans(myo_WGBS_exons[both,], na.rm = T)
# WGBS_exons[WGBS_exons$group == "myo.not_expressed", "WGBS"] <- colMeans(myo_WGBS_exons[neither,], na.rm = T)
# (WGBS_exons_profile <- ggplot(WGBS_exons, aes(x = Position, y = WGBS, group = group)) + 
#    geom_line(aes(color = Expression, linetype = Cell_type), size = 1.5) + 
#    geom_point(aes(color = Expression, shape = Cell_type), size = 4) + 
#    ggtitle("WGBS profile along exons") + 
#    xlab("Normalized position along exons") + 
#    ylab("Average DNA methylation level") + 
#    theme_bw())
# ggsave(WGBS_exons_profile, file = "WGBS_exons_profile")
# 
