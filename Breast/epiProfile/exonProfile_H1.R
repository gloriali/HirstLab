# Epigenetic profile at exon boundaries
setwd("~/快盘/REMC/epiProfile/")
library(ggplot2)
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
rownames(H1) <- H1$coord
myo084 <- read.delim("~/REMC/epiProfile/A17919.G.exn.A.rpkm", head = F, as.is = T)
myo084$coord <- gsub("<[-]*1", "", myo084$V1)
myo084$coord  <- gsub("chr", "", myo084$coord)
myo084$coord  <- gsub(":", "_", myo084$coord)
myo084$coord  <- gsub("-", "_", myo084$coord)
myo084$coord  <- paste0(myo084$V2, "_", myo084$coord)
rownames(myo084) <- myo084$coord
both <- intersect(H1[H1$V4 > 0.1, "coord"], myo084[myo084$V4 > 0.1, "coord"])
neither <- setdiff(H1$coord, c(H1_specific, myo_specific, both))

H1gene <- read.delim("~/REMC/gene/H1_r1a.G.A.rpkm.pc", head = F, as.is = T)
rownames(H1gene) <- H1gene$V1
myo084gene <- read.delim("~/REMC/gene/A17919.G.A.rpkm.pc", head = F, as.is = T)
rownames(myo084gene) <- myo084gene$V1
geneRPKM <- data.frame(gene = H1gene$V1, RPKM = (H1gene$V3 + myo084gene$V3)/2)
rownames(geneRPKM) <- geneRPKM$gene
exons_geneRPKM <- data.frame(exon = H1$coord)
exons_geneRPKM$gene <- gsub("_[0-9XY]+_[0-9_]+", "", exons_geneRPKM$exon)
exons_geneRPKM$geneRPKM <- geneRPKM[exons_geneRPKM$gene, "RPKM"]
exons_geneRPKM$H1gene <- H1gene[exons_geneRPKM$gene, "V3"]
exons_geneRPKM$myo084gene <- myo084gene[exons_geneRPKM$gene, "V3"]
exons_geneRPKM$H1exon <- H1[as.character(exons_geneRPKM$exon), "V4"]
exons_geneRPKM$myo084exon <- myo084[as.character(exons_geneRPKM$exon), "V4"]
exons_geneRPKM <- na.omit(exons_geneRPKM)
exons_geneRPKM$exon <- as.character(exons_geneRPKM$exon)
H1_expressed <- as.character(exons_geneRPKM[H1gene > 0.1, "exon"])
myo_expressed <- as.character(exons_geneRPKM[myo084gene > 0.1, "exon"])
H1_specific <- intersect(H1_specific, intersect(H1_expressed, myo_expressed))
myo_specific <- intersect(myo_specific, intersect(H1_expressed, myo_expressed))
both <- setdiff(exons_geneRPKM[exons_geneRPKM$H1exon > 0.1 & exons_geneRPKM$myo084exon > 0.1, "exon"], c(H1_specific, myo_specific))
neither <- setdiff(H1$coord, c(H1_specific, myo_specific, both))
H1_exon_0 <- exons_geneRPKM[exons_geneRPKM$H1exon == 0 & exons_geneRPKM$H1gene > 0.1, "exon"]
myo_exon_0 <- exons_geneRPKM[exons_geneRPKM$myo084exon == 0 & exons_geneRPKM$myo084gene > 0.1, "exon"]
H1_gene_0 <- exons_geneRPKM[exons_geneRPKM$H1gene == 0, "exon"]
myo_gene_0 <- exons_geneRPKM[exons_geneRPKM$myo084gene == 0, "exon"]

######################################################################################################
# WGBS profile @ exon boundaries
H1_WGBS_3p <- read.table("~/REMC/epiProfile/exons3p_200/H1_bismark.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_WGBS_3p <- read.table("~/REMC/epiProfile/exons3p_200/myoRM045_WGBS.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
H1_WGBS_5p <- read.table("~/REMC/epiProfile/exons5p_200/H1_bismark.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_WGBS_5p <- read.table("~/REMC/epiProfile/exons5p_200/myoRM045_WGBS.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
H1_WGBS_3p <- H1_WGBS_3p/10
H1_WGBS_5p <- H1_WGBS_5p/10
WGBS_3p <- data.frame(Cell_type = rep(c("H1", "myo"), each = 20*4), Expression = rep(c(rep(c("H1-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), WGBS = -1)
WGBS_3p$group <- interaction(WGBS_3p$Cell_type, WGBS_3p$Expression)
WGBS_3p[WGBS_3p$group == "H1.H1-specific", "WGBS"] <- colMeans(H1_WGBS_3p[H1_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "H1.myo-specific", "WGBS"] <- colMeans(H1_WGBS_3p[myo_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "H1.expressed_in_both", "WGBS"] <- colMeans(H1_WGBS_3p[both,], na.rm = T)
WGBS_3p[WGBS_3p$group == "H1.not_expressed", "WGBS"] <- colMeans(H1_WGBS_3p[neither,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.H1-specific", "WGBS"] <- colMeans(myo_WGBS_3p[H1_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.myo-specific", "WGBS"] <- colMeans(myo_WGBS_3p[myo_specific,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.expressed_in_both", "WGBS"] <- colMeans(myo_WGBS_3p[both,], na.rm = T)
WGBS_3p[WGBS_3p$group == "myo.not_expressed", "WGBS"] <- colMeans(myo_WGBS_3p[neither,], na.rm = T)
WGBS_5p <- data.frame(Cell_type = rep(c("H1", "myo"), each = 20*4), Expression = rep(c(rep(c("H1-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), WGBS = -1)
WGBS_5p$group <- interaction(WGBS_5p$Cell_type, WGBS_5p$Expression)
WGBS_5p[WGBS_5p$group == "H1.H1-specific", "WGBS"] <- colMeans(H1_WGBS_5p[H1_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "H1.myo-specific", "WGBS"] <- colMeans(H1_WGBS_5p[myo_specific,], na.rm = T)
WGBS_5p[WGBS_5p$group == "H1.expressed_in_both", "WGBS"] <- colMeans(H1_WGBS_5p[both,], na.rm = T)
WGBS_5p[WGBS_5p$group == "H1.not_expressed", "WGBS"] <- colMeans(H1_WGBS_5p[neither,], na.rm = T)
WGBS_5p[WGBS_5p$group == "myo.H1-specific", "WGBS"] <- colMeans(myo_WGBS_5p[H1_specific,], na.rm = T)
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

# add CpG content track 
CpG_content_3p <- read.table("~/hg19/CpG.hg19v65_exons_for_genes.3prime_200.unique", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
CpG_content_5p <- read.table("~/hg19/CpG.hg19v65_exons_for_genes.5prime_200.unique", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
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

WGBS_CpG <- data.frame(data = c(rep("WGBS", nrow(WGBS_boundaries)), rep("CpG", ncol(CpG_content_3p)*2)), 
                       Cell_type = c(as.character(WGBS_boundaries$Cell_type), rep("H1", ncol(CpG_content_3p)*2)), 
                       Expression = c(as.character(WGBS_boundaries$Expression), rep("CpG density", ncol(CpG_content_3p)*2)), 
                       Position = c(WGBS_boundaries$Position, rep(seq(-190, 190, by = ncol(CpG_content_3p)), times = 2)), 
                       value = c(WGBS_boundaries$WGBS, colMeans(CpG_content_3p), colMeans(CpG_content_5p)), 
                       End = c(as.character(WGBS_boundaries$End), rep(c("3-prime", "5-prime"), each = ncol(CpG_content_3p))))
WGBS_CpG$data <- factor(WGBS_CpG$data, levels = c("WGBS", "CpG"))
WGBS_CpG$End <- factor(WGBS_CpG$End, levels = c("5-prime", "3-prime"))
WGBS_CpG$Expression <- factor(WGBS_CpG$Expression, levels = c("H1-specific", "myo-specific", "expressed_in_both", "not_expressed", "CpG density"))
WGBS_CpG$group <- interaction(WGBS_CpG$Cell_type, WGBS_CpG$Expression)
(WGBS_CpG_profile <- ggplot(WGBS_CpG, aes(x = Position, y = value, group = group)) + 
   geom_line(aes(color = Expression, linetype = Cell_type)) + 
   geom_point(aes(color = Expression, shape = Cell_type)) + 
   facet_grid(data ~ End, scales = "free_y") + 
   ylab("") + 
   scale_color_manual(values = c("CpG density" = "black", "H1-specific" = rgb(200,50,0, maxColorValue = 255), "myo-specific" = rgb(50,200,50, maxColorValue = 255), "expressed_in_both" = "purple", "not_expressed" = "blue")) + 
   theme(panel.border = element_rect(linetype = "solid", fill = "transparent"), panel.margin = unit(0.75, "lines"), axis.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 12, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
library(grid) 
gt <- ggplot_gtable(ggplot_build(WGBS_CpG_profile)) 
# gt$layout 
gt$heights[[4]] <- unit(3.5, "null") 
# grid.newpage()
pdf("WGBS_CpG_profile_H1.pdf", width = 9)
grid.draw(gt) 
grid.text("Average DNA methylation", x = unit(0.015, "npc"), y = unit(0.65, "npc"), rot = 90)
grid.text("CpG density", x = unit(0.015, "npc"), y = unit(0.15, "npc"), rot = 90)
dev.off()

# statistical test 
H1_WGBS_3p_peak <- data.frame(id = rownames(H1_WGBS_3p), max = apply(H1_WGBS_3p, 1, max))
rownames(H1_WGBS_3p_peak) <- rownames(H1_WGBS_3p)
t.test(na.omit(H1_WGBS_3p_peak[both, "max"]), na.omit(H1_WGBS_3p_peak[neither, "max"]))$p.value
t.test(na.omit(H1_WGBS_3p_peak[H1_specific, "max"]), na.omit(H1_WGBS_3p_peak[neither, "max"]))$p.value
t.test(na.omit(H1_WGBS_3p_peak[myo_specific, "max"]), na.omit(H1_WGBS_3p_peak[neither, "max"]))$p.value
H1_WGBS_5p_peak <- data.frame(id = rownames(H1_WGBS_5p), max = apply(H1_WGBS_5p, 1, max))
rownames(H1_WGBS_5p_peak) <- rownames(H1_WGBS_5p)
t.test(na.omit(H1_WGBS_5p_peak[both, "max"]), na.omit(H1_WGBS_5p_peak[neither, "max"]))$p.value
t.test(na.omit(H1_WGBS_5p_peak[H1_specific, "max"]), na.omit(H1_WGBS_5p_peak[neither, "max"]))$p.value
t.test(na.omit(H1_WGBS_5p_peak[myo_specific, "max"]), na.omit(H1_WGBS_5p_peak[neither, "max"]))$p.value
myo_WGBS_3p_peak <- data.frame(id = rownames(myo_WGBS_3p), max = apply(myo_WGBS_3p, 1, max))
rownames(myo_WGBS_3p_peak) <- rownames(myo_WGBS_3p)
t.test(na.omit(myo_WGBS_3p_peak[both, "max"]), na.omit(myo_WGBS_3p_peak[neither, "max"]))$p.value
t.test(na.omit(myo_WGBS_3p_peak[H1_specific, "max"]), na.omit(myo_WGBS_3p_peak[neither, "max"]))$p.value
t.test(na.omit(myo_WGBS_3p_peak[myo_specific, "max"]), na.omit(myo_WGBS_3p_peak[neither, "max"]))$p.value
myo_WGBS_5p_peak <- data.frame(id = rownames(myo_WGBS_5p), max = apply(myo_WGBS_5p, 1, max))
rownames(myo_WGBS_5p_peak) <- rownames(myo_WGBS_5p)
t.test(na.omit(myo_WGBS_5p_peak[both, "max"]), na.omit(myo_WGBS_5p_peak[neither, "max"]))$p.value
t.test(na.omit(myo_WGBS_5p_peak[H1_specific, "max"]), na.omit(myo_WGBS_5p_peak[neither, "max"]))$p.value
t.test(na.omit(myo_WGBS_5p_peak[myo_specific, "max"]), na.omit(myo_WGBS_5p_peak[neither, "max"]))$p.value

t.test(na.omit(H1_WGBS_3p_peak[H1_specific, "max"]), na.omit(myo_WGBS_3p_peak[H1_specific, "max"]))$p.value
t.test(na.omit(H1_WGBS_5p_peak[H1_specific, "max"]), na.omit(myo_WGBS_5p_peak[H1_specific, "max"]))$p.value
t.test(na.omit(H1_WGBS_3p_peak[myo_specific, "max"]), na.omit(myo_WGBS_3p_peak[myo_specific, "max"]))$p.value
t.test(na.omit(H1_WGBS_5p_peak[myo_specific, "max"]), na.omit(myo_WGBS_5p_peak[myo_specific, "max"]))$p.value

######################################################################################################
# H3K36me3 signal for exons
H1_H3K36me3_exons <- read.delim("~/REMC/epiProfile/exons/hg19v65_exons_for_genes.H1_H3K36me3.coverage", head = F, as.is = T)
myo_H3K36me3_exons <- read.delim("~/REMC/epiProfile/exons/hg19v65_exons_for_genes.myoRM080_H3K36me3.coverage", head = F, as.is = T)
# normalize signal
(norm <- sum(myo_H3K36me3_exons$V6)/sum(H1_H3K36me3_exons$V6))
H3K36me3_exons <- data.frame(id = c(H1_H3K36me3_exons$V4, myo_H3K36me3_exons$V4), Cell_type = c(rep("H1", nrow(H1_H3K36me3_exons)), rep("myo", nrow(myo_H3K36me3_exons))), Expression = NA, H3K36me3 = c(H1_H3K36me3_exons$V6 * norm, myo_H3K36me3_exons$V6))
H3K36me3_exons[H3K36me3_exons$id %in% both, "Expression"] <- "expressed_in_both"
H3K36me3_exons[H3K36me3_exons$id %in% H1_exon_0 & as.character(H3K36me3_exons$Cell_type) == "H1", "Expression"] <- "not_expressed"
H3K36me3_exons[H3K36me3_exons$id %in% myo_exon_0 & as.character(H3K36me3_exons$Cell_type) == "myo", "Expression"] <- "not_expressed"
H3K36me3_exons[H3K36me3_exons$id %in% H1_specific, "Expression"] <- "H1-specific"
H3K36me3_exons[H3K36me3_exons$id %in% myo_specific, "Expression"] <- "myo-specific"
H3K36me3_exons[H3K36me3_exons$id %in% H1_gene_0 & as.character(H3K36me3_exons$Cell_type) == "H1", "Expression"] <- "gene_not_expressed"
H3K36me3_exons[H3K36me3_exons$id %in% myo_gene_0 & as.character(H3K36me3_exons$Cell_type) == "myo", "Expression"] <- "gene_not_expressed"
H3K36me3_exons <- droplevels(H3K36me3_exons[H3K36me3_exons$Expression != "gene_not_expressed", ])
H3K36me3_exons$Expression <- factor(H3K36me3_exons$Expression, levels = c("expressed_in_both", "H1-specific", "myo-specific", "not_expressed"))
H3K36me3_exons$utilize <- interaction(H3K36me3_exons$Cell_type, H3K36me3_exons$Expression)
H3K36me3_exons <- na.omit(H3K36me3_exons)
summary(H3K36me3_exons$utilize)
library(plyr)
H3K36me3_exons_stat <- ddply(H3K36me3_exons, ~ utilize, summarize, Cell_type = Cell_type[1], Expression = Expression[1], ymin = boxplot.stats(H3K36me3)$stats[1], lower = boxplot.stats(H3K36me3)$stats[2], middle = mean(H3K36me3), upper = boxplot.stats(H3K36me3)$stats[4], ymax = boxplot.stats(H3K36me3)$stats[5])

(H3K36me3_exons_profile <- ggplot(H3K36me3_exons_stat, aes(x = Expression, group = utilize)) + 
   geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax, fill = Cell_type), stat = "identity", position = "dodge", width = 0.8) + 
   # facet_grid(geneRPKM ~ ., scales = "free") + 
   # ggtitle("H3K36me3 signal for exons") + 
   xlab("Exon group") + 
   ylab("Average H3K36me3 signal") + 
   scale_fill_manual(values = c("H1" = rgb(200,50,0, maxColorValue = 255), "myo" = rgb(50,200,50, maxColorValue = 255))) + 
   theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 15, color = "black"), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 20, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
ggsave(H3K36me3_exons_profile, file = "H3K36me3_exons_profile.pdf", width = 9, height = 5)
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "H1.H1-specific" | H3K36me3_exons$utilize == "myo.myo-specific", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "myo.H1-specific" | H3K36me3_exons$utilize == "H1.myo-specific", "H3K36me3"], alternative = "greater")$p.value
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "H1.H1-specific", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "myo.H1-specific", "H3K36me3"])$p.value
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "myo.myo-specific", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "H1.myo-specific", "H3K36me3"])$p.value
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "H1.expressed_in_both", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "myo.expressed_in_both", "H3K36me3"])$p.value
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "H1.not_expressed", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "myo.not_expressed", "H3K36me3"])$p.value
t.test(H3K36me3_exons[H3K36me3_exons$utilize == "H1.gene_not_expressed", "H3K36me3"], H3K36me3_exons[H3K36me3_exons$utilize == "myo.gene_not_expressed", "H3K36me3"])$p.value
# fold enrichment between utilized/un-utilized isoform exons
mean(H3K36me3_exons[H3K36me3_exons$utilize == "H1.H1-specific" | H3K36me3_exons$utilize == "myo.myo-specific", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.H1-specific" | H3K36me3_exons$utilize == "H1.myo-specific", "H3K36me3"])
mean(H3K36me3_exons[H3K36me3_exons$utilize == "H1.H1-specific", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.H1-specific", "H3K36me3"])
mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.myo-specific", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "H1.myo-specific", "H3K36me3"])
mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.expressed_in_both", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "H1.expressed_in_both", "H3K36me3"])
mean(H3K36me3_exons[H3K36me3_exons$utilize == "myo.not_expressed", "H3K36me3"])/mean(H3K36me3_exons[H3K36me3_exons$utilize == "H1.not_expressed", "H3K36me3"])

# exon usage in transcripts
library(plyr)
exons <- read.delim("~/快盘/hg19/hg19v65_exons", head = F, as.is = T)
exons$id <- paste0(exons$V1, "_", exons$V3, "_", exons$V4, "_", exons$V5)
Ntranscript <- ddply(exons, ~ V1, summarize, Ntranscript = length(unique(V2)))
rownames(Ntranscript) <- Ntranscript$V1
exon <- ddply(exons, ~ id, summarize, transcript = length(unique(V2)))
exon$gene <- gsub("_[0-9XY]+_[0-9_]+", "", exon$id)
exons <- exons[!duplicated(exons$id),]
rownames(exons) <- exons$id
exon$Ntranscript <- Ntranscript[exon$gene, "Ntranscript"]
exon$percentTrans <- exon$transcript / exon$Ntranscript
rownames(exon) <- exon$id
exon_usage <- na.omit(data.frame(id = H3K36me3_exons$id, Cell_type = H3K36me3_exons$Cell_type, Expression = H3K36me3_exons$Expression, utilize = H3K36me3_exons$utilize, usage = exon[as.character(H3K36me3_exons$id), "percentTrans"]))
exon_usage_stat <- ddply(exon_usage, ~ utilize, summarize, Cell_type = Cell_type[1], Expression = Expression[1], ymin = boxplot.stats(usage)$stats[1], lower = boxplot.stats(usage)$stats[2], middle = mean(usage), upper = boxplot.stats(usage)$stats[4], ymax = boxplot.stats(usage)$stats[5])
H3K36me3_exons_stat <- data.frame(rbind(H3K36me3_exons_stat, exon_usage_stat), type = c(rep("H3K36me3", nrow(H3K36me3_exons_stat)), rep("Exon usage", nrow(exon_usage_stat))))
H3K36me3_exons_stat$type <- factor(H3K36me3_exons_stat$type, levels = c("H3K36me3", "Exon usage"))
(H3K36me3_exons_profile <- ggplot(H3K36me3_exons_stat, aes(x = Expression, group = utilize)) + 
   geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax, fill = Cell_type), stat = "identity", position = "dodge", width = 0.8) + 
   facet_grid(type ~ ., scales = "free") + 
   # ggtitle("H3K36me3 signal for exons") + 
   xlab("Exon group") + 
   ylab("") + 
   scale_fill_manual(values = c("H1" = rgb(200,50,0, maxColorValue = 255), "myo" = rgb(50,200,50, maxColorValue = 255))) + 
   theme(panel.border = element_rect(linetype = "solid", fill = "transparent"), axis.title = element_text(size = 20), axis.text.x = element_text(size = 15, color = "black"), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 20, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
library(grid) 
gt <- ggplot_gtable(ggplot_build(H3K36me3_exons_profile)) 
# gt$layout 
gt$heights[[3]] <- unit(2, "null") 
# grid.newpage()
pdf("H3K36me3_exons_profile_H1.pdf", width = 10, height = 7)
grid.draw(gt) 
grid.text("Average H3K36me3 signal", x = unit(0.015, "npc"), y = unit(0.65, "npc"), rot = 90)
grid.text("Exon usage", x = unit(0.015, "npc"), y = unit(0.25, "npc"), rot = 90)
grid.text("myo/H1=1.01", x = unit(0.18, "npc"), y = unit(0.955, "npc"))
grid.text("H1/myo=1.44", x = unit(0.34, "npc"), y = unit(0.6, "npc"))
grid.text("myo/H1=1.55", x = unit(0.51, "npc"), y = unit(0.6, "npc"))
grid.text("myo/H1=1.10", x = unit(0.68, "npc"), y = unit(0.6, "npc"))
dev.off()

