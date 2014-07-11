# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# Epigenetic profile at intron boundaries
setwd("~/快盘/REMC/epiProfile/IR/v2")
library(plyr)
library(ggplot2)
library(grid) 
load("intronProfile_v2.Rdata")
introns <- read.delim("~/hg19/hg19v65_introns_for_genes", head = F, as.is = T)
lum084_ir <- read.delim("../A17918.introns.retained.Id4Gloria", head = F, as.is = T)$V1
myo084_ir <- read.delim("../A17919.introns.retained.Id4Gloria", head = F, as.is = T)$V1
lum084_other <- setdiff(introns$V5, lum084_ir)
myo084_other <- setdiff(introns$V5, myo084_ir)
rm(introns)

######################################################################################################
# WGBS profile @ intron boundaries
lum_bismark_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/lumRM066_bismark.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_bismark_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM045_bismark.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_bismark_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/lumRM066_bismark.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_bismark_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM045_bismark.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
WGBS_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), WGBS = -1)
WGBS_3p[WGBS_3p$Expression == "lum_IR", "WGBS"] <- colMeans(lum_bismark_3p[lum084_ir,], na.rm = T)
WGBS_3p[WGBS_3p$Expression == "lum_not-retained", "WGBS"] <- colMeans(lum_bismark_3p[lum084_other,], na.rm = T)
WGBS_3p[WGBS_3p$Expression == "myo_IR", "WGBS"] <- colMeans(myo_bismark_3p[myo084_ir,], na.rm = T)
WGBS_3p[WGBS_3p$Expression == "myo_not-retained", "WGBS"] <- colMeans(myo_bismark_3p[myo084_other,], na.rm = T)
WGBS_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), WGBS = -1)
WGBS_5p[WGBS_5p$Expression == "lum_IR", "WGBS"] <- colMeans(lum_bismark_5p[lum084_ir,], na.rm = T)
WGBS_5p[WGBS_5p$Expression == "lum_not-retained", "WGBS"] <- colMeans(lum_bismark_5p[lum084_other,], na.rm = T)
WGBS_5p[WGBS_5p$Expression == "myo_IR", "WGBS"] <- colMeans(myo_bismark_5p[myo084_ir,], na.rm = T)
WGBS_5p[WGBS_5p$Expression == "myo_not-retained", "WGBS"] <- colMeans(myo_bismark_5p[myo084_other,], na.rm = T)
WGBS_boundaries_v2 <- data.frame(rbind(WGBS_3p, WGBS_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(WGBS_3p)), levels = c("5-prime", "3-prime")))
(WGBS_boundaries_v2_profile <- ggplot(WGBS_boundaries_v2, aes(x = Position, y = WGBS, group = Expression)) + 
   geom_line(aes(color = IR)) + 
   geom_point(aes(color = IR)) + 
   facet_grid(Cell_type ~ End) + 
   ggtitle("WGBS profile around intron boundaries") + 
   ylab("Average DNA methylation level") + 
   theme_bw())
ggsave(WGBS_boundaries_v2_profile, file = "WGBS_boundaries_v2_profile.pdf")

# add CpG content track 
CpG_content_3p <- read.table("~/hg19/CpG.hg19v65_introns_for_genes.3prime_200", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
CpG_content_5p <- read.table("~/hg19/CpG.hg19v65_introns_for_genes.5prime_200", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
CpG_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), CpG = -1)
CpG_3p[CpG_3p$Expression == "lum_IR", "CpG"] <- colMeans(CpG_content_3p[lum084_ir,], na.rm = T)
CpG_3p[CpG_3p$Expression == "lum_not-retained", "CpG"] <- colMeans(CpG_content_3p[lum084_other,], na.rm = T)
CpG_3p[CpG_3p$Expression == "myo_IR", "CpG"] <- colMeans(CpG_content_3p[myo084_ir,], na.rm = T)
CpG_3p[CpG_3p$Expression == "myo_not-retained", "CpG"] <- colMeans(CpG_content_3p[myo084_other,], na.rm = T)
CpG_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), CpG = -1)
CpG_5p[CpG_5p$Expression == "lum_IR", "CpG"] <- colMeans(CpG_content_5p[lum084_ir,], na.rm = T)
CpG_5p[CpG_5p$Expression == "lum_not-retained", "CpG"] <- colMeans(CpG_content_5p[lum084_other,], na.rm = T)
CpG_5p[CpG_5p$Expression == "myo_IR", "CpG"] <- colMeans(CpG_content_5p[myo084_ir,], na.rm = T)
CpG_5p[CpG_5p$Expression == "myo_not-retained", "CpG"] <- colMeans(CpG_content_5p[myo084_other,], na.rm = T)
CpG_boundaries_v2 <- data.frame(rbind(CpG_3p, CpG_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(CpG_3p)), levels = c("5-prime", "3-prime")))
(CpG_boundaries_v2_profile <- ggplot(CpG_boundaries_v2, aes(x = Position, y = CpG, group = Expression)) + 
   geom_line(aes(color = IR)) + 
   geom_point(aes(color = IR)) + 
   facet_grid(Cell_type ~ End) + 
   ggtitle("CpG profile around intron boundaries") + 
   ylab("Average percentage of CpG") + 
   theme_bw())
ggsave(CpG_boundaries_v2_profile, file = "CpG_boundaries_v2_profile.pdf", height = 4)

WGBS_CpG_v2 <- data.frame(data = c(rep("WGBS", nrow(WGBS_boundaries_v2)), rep("CpG", nrow(CpG_boundaries_v2))), 
                       Cell_type = c(as.character(WGBS_boundaries_v2$Cell_type), as.character(CpG_boundaries_v2$Cell_type)), 
                       IR = c(as.character(WGBS_boundaries_v2$IR), as.character(CpG_boundaries_v2$IR)), 
                       Expression = c(as.character(WGBS_boundaries_v2$Expression), as.character(CpG_boundaries_v2$Expression)), 
                       Position = c(WGBS_boundaries_v2$Position, CpG_boundaries_v2$Position), 
                       value = c(WGBS_boundaries_v2$WGBS, CpG_boundaries_v2$CpG), 
                       End = c(as.character(WGBS_boundaries_v2$End), as.character(CpG_boundaries_v2$End)))
WGBS_CpG_v2$IR <- gsub("non-IR", "not retained", WGBS_CpG_v2$IR)
WGBS_CpG_v2$IR <- gsub("IR", "retained introns", WGBS_CpG_v2$IR)
WGBS_CpG_v2$IR <- factor(WGBS_CpG_v2$IR, levels = c("retained introns", "not retained"))
WGBS_CpG_v2$color <- as.character(interaction(WGBS_CpG_v2$Cell_type, WGBS_CpG_v2$IR))
WGBS_CpG_v2$color <- gsub("lum.not retained", "not retained", WGBS_CpG_v2$color)
WGBS_CpG_v2$color <- gsub("myo.not retained", "not retained", WGBS_CpG_v2$color)
WGBS_CpG_v2$color <- factor(WGBS_CpG_v2$color, levels = c("lum.retained introns", "myo.retained introns", "not retained"))
WGBS_CpG_v2$group <- factor(interaction(WGBS_CpG_v2$data, WGBS_CpG_v2$Cell_type), levels = c("WGBS.lum", "CpG.lum", "WGBS.myo", "CpG.myo"))
WGBS_CpG_v2$End <- factor(WGBS_CpG_v2$End, levels = c("5-prime", "3-prime"))
(WGBS_CpG_v2_profile <- ggplot(WGBS_CpG_v2, aes(x = Position, y = value, group = Expression)) + 
   geom_line(aes(color = color)) + 
   geom_point(aes(color = color)) + 
   facet_grid(group ~ End, scales = "free_y") + 
   ylab("") + 
   scale_color_manual(name = "IR", values = c("lum.retained introns" = rgb(200,50,0, maxColorValue = 255), "myo.retained introns" = rgb(50,200,50, maxColorValue = 255), "not retained" = "blue")) + 
   theme(panel.border = element_rect(linetype = "solid", fill = "transparent"), panel.margin = unit(0.75, "lines"), axis.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 12, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
gt <- ggplot_gtable(ggplot_build(WGBS_CpG_v2_profile)) 
# gt$layout 
gt$heights[[4]] <- unit(2.5, "null") 
gt$heights[[8]] <- unit(2.5, "null") 
pdf("WGBS_CpG_v2_profile.pdf", width = 10, height = 10)
# grid.newpage()
grid.draw(gt) 
grid.text("Average DNA methylation", x = unit(0.015, "npc"), y = unit(0.82, "npc"), rot = 90)
grid.text("Average DNA methylation", x = unit(0.015, "npc"), y = unit(0.35, "npc"), rot = 90)
grid.text("No. of CpGs", x = unit(0.015, "npc"), y = unit(0.12, "npc"), rot = 90)
grid.text("No. of CpGs", x = unit(0.015, "npc"), y = unit(0.58, "npc"), rot = 90)
dev.off()
######################################################################################################
# delta WGBS vs intron/exon correlation
e <- 1e-6
lum_intron_exon <- read.delim("../A17918.introns_cov.fankingExons", , head = F, as.is = T)
rownames(lum_intron_exon) <- paste0(lum_intron_exon$V5, "_", lum_intron_exon$V1, ":", lum_intron_exon$V2, "-", lum_intron_exon$V3, "<", lum_intron_exon$V4)
lum_intron_exon$ratio3p <- (lum_intron_exon$V6 + e) / (lum_intron_exon$V8 + e)
lum_intron_exon$ratio5p <- (lum_intron_exon$V6 + e) / (lum_intron_exon$V7 + e)
lum_correlate3p <- data.frame(id = rownames(lum_intron_exon), 
                              ratio = lum_intron_exon[, "ratio3p"], 
                              delta = apply(lum_bismark_3p[, 1:10], 1, max) - apply(lum_bismark_3p[, 11:20], 1, min))
plot(x=c(0.3, 1), y=c(0, 0.2), type = "n", xlab = "intron/3p-exon", ylab = "3p-delta")
points(lum_correlate3p$ratio, lum_correlate3p$delta, pch = 19, cex = 0.5)
cor(lum_correlate3p$ratio, lum_correlate3p$delta)
lum_correlate5p <- data.frame(id = rownames(lum_intron_exon), 
                              ratio = lum_intron_exon[, "ratio5p"], 
                              delta = apply(lum_bismark_5p[, 1:10], 1, max) - apply(lum_bismark_5p[, 11:20], 1, min))
plot(x=c(0.3, 1), y=c(0, 0.2), type = "n", xlab = "intron/5p-exon", ylab = "5p-delta")
points(lum_correlate5p$ratio, lum_correlate5p$delta, pch = 19, cex = 0.5)
cor(lum_correlate5p$ratio, lum_correlate5p$delta)

myo_intron_exon <- read.delim("../A17919.introns_cov.fankingExons", , head = F, as.is = T)
rownames(myo_intron_exon) <- paste0(myo_intron_exon$V5, "_", myo_intron_exon$V1, ":", myo_intron_exon$V2, "-", myo_intron_exon$V3, "<", myo_intron_exon$V4)
myo_intron_exon$ratio3p <- (myo_intron_exon$V6 + e) / (myo_intron_exon$V8 + e)
myo_intron_exon$ratio5p <- (myo_intron_exon$V6 + e) / (myo_intron_exon$V7 + e)
myo_correlate3p <- data.frame(id = rownames(myo_intron_exon), 
                              ratio = myo_intron_exon[, "ratio3p"], 
                              delta = apply(myo_bismark_3p[, 1:10], 1, max) - apply(myo_bismark_3p[, 11:20], 1, min))
plot(x=c(0.3, 1), y=c(0, 0.2), type = "n", xlab = "intron/3p-exon", ylab = "3p-delta")
points(myo_correlate3p$ratio, myo_correlate3p$delta, pch = 19, cex = 0.5)
cor(myo_correlate3p$ratio, myo_correlate3p$delta)
myo_correlate5p <- data.frame(id = rownames(myo_intron_exon), 
                              ratio = myo_intron_exon[, "ratio5p"], 
                              delta = apply(myo_bismark_5p[, 1:10], 1, max) - apply(myo_bismark_5p[, 11:20], 1, min))
plot(x=c(0.3, 1), y=c(0, 0.2), type = "n", xlab = "intron/5p-exon", ylab = "5p-delta")
points(myo_correlate5p$ratio, myo_correlate5p$delta, pch = 19, cex = 0.5)
cor(myo_correlate5p$ratio, myo_correlate5p$delta)

# intron RPKM ecdf
pdf("Ecdf_intronRPKM.pdf")
plot(x=c(0,10), y=c(0,1), xlab = "intron RPKM", ylab = "Ecdf", type="n")
lines(ecdf(lum_intron_exon$V6), col = "red")
lines(ecdf(myo_intron_exon$V6), col = "green")
legend("bottomright", c("lum", "myo"), col = c("red", "green"), lwd = 5)
dev.off()

######################################################################################################
# MeDIP profile @ intron boundaries
lum_MeDIP_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/lumRM035_MeDIP.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_MeDIP_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM035_MeDIP.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_MeDIP_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/lumRM035_MeDIP.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_MeDIP_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM035_MeDIP.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_MeDIP_3p <- sum(myo_MeDIP_3p)/sum(lum_MeDIP_3p) * lum_MeDIP_3p
lum_MeDIP_5p <- sum(myo_MeDIP_5p)/sum(lum_MeDIP_5p) * lum_MeDIP_5p
MeDIP_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), MeDIP = -1)
MeDIP_3p[MeDIP_3p$Expression == "lum_IR", "MeDIP"] <- colMeans(lum_MeDIP_3p[lum084_ir,], na.rm = T)
MeDIP_3p[MeDIP_3p$Expression == "lum_not-retained", "MeDIP"] <- colMeans(lum_MeDIP_3p[lum084_other,], na.rm = T)
MeDIP_3p[MeDIP_3p$Expression == "myo_IR", "MeDIP"] <- colMeans(myo_MeDIP_3p[myo084_ir,], na.rm = T)
MeDIP_3p[MeDIP_3p$Expression == "myo_not-retained", "MeDIP"] <- colMeans(myo_MeDIP_3p[myo084_other,], na.rm = T)
MeDIP_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), MeDIP = -1)
MeDIP_5p[MeDIP_5p$Expression == "lum_IR", "MeDIP"] <- colMeans(lum_MeDIP_5p[lum084_ir,], na.rm = T)
MeDIP_5p[MeDIP_5p$Expression == "lum_not-retained", "MeDIP"] <- colMeans(lum_MeDIP_5p[lum084_other,], na.rm = T)
MeDIP_5p[MeDIP_5p$Expression == "myo_IR", "MeDIP"] <- colMeans(myo_MeDIP_5p[myo084_ir,], na.rm = T)
MeDIP_5p[MeDIP_5p$Expression == "myo_not-retained", "MeDIP"] <- colMeans(myo_MeDIP_5p[myo084_other,], na.rm = T)
MeDIP_boundaries_v2 <- data.frame(rbind(MeDIP_3p, MeDIP_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(MeDIP_3p)), levels = c("5-prime", "3-prime")))
MeDIP_boundaries_v2$IR <- gsub("non-IR", "not retained", MeDIP_boundaries_v2$IR)
MeDIP_boundaries_v2$IR <- gsub("IR", "retained introns", MeDIP_boundaries_v2$IR)
MeDIP_boundaries_v2$IR <- factor(MeDIP_boundaries_v2$IR, levels = c("retained introns", "not retained"))
MeDIP_boundaries_v2$color <- as.character(interaction(MeDIP_boundaries_v2$Cell_type, MeDIP_boundaries_v2$IR))
MeDIP_boundaries_v2$color <- gsub("lum.not retained", "not retained", MeDIP_boundaries_v2$color)
MeDIP_boundaries_v2$color <- gsub("myo.not retained", "not retained", MeDIP_boundaries_v2$color)
MeDIP_boundaries_v2$color <- factor(MeDIP_boundaries_v2$color, levels = c("lum.retained introns", "myo.retained introns", "not retained"))
(MeDIP_boundaries_v2_profile <- ggplot(MeDIP_boundaries_v2, aes(x = Position, y = MeDIP, group = Expression)) + 
   geom_line(aes(color = color)) + 
   geom_point(aes(color = color)) + 
   facet_grid(Cell_type ~ End) + 
   ylab("Average MeDIP signal") + 
   # ggtitle("MeDIP profile around intron boundaries") + 
   scale_color_manual(name = "IR", values = c("lum.retained introns" = rgb(200,50,0, maxColorValue = 255), "myo.retained introns" = rgb(50,200,50, maxColorValue = 255), "not retained" = "blue")) + 
   theme(panel.border = element_rect(linetype = "solid", fill = "transparent"), panel.margin = unit(0.75, "lines"), axis.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 12, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
ggsave(MeDIP_boundaries_v2_profile, file = "MeDIP_boundaries_v2_profile.pdf")

######################################################################################################
# H3K4me3 profile @ intron boundaries (no lum libraries available)
myo_H3K4me3_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM080_H3K4me3.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K4me3_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM080_H3K4me3.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
H3K4me3_3p <- data.frame(Cell_type = rep("myo", 20), IR = rep(c("IR", "non-IR"), each = 20), Expression = rep(c("myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 2), H3K4me3 = -1)
H3K4me3_3p[H3K4me3_3p$Expression == "myo_IR", "H3K4me3"] <- colMeans(myo_H3K4me3_3p[myo084_ir,], na.rm = T)
H3K4me3_3p[H3K4me3_3p$Expression == "myo_not-retained", "H3K4me3"] <- colMeans(myo_H3K4me3_3p[myo084_other,], na.rm = T)
H3K4me3_5p <- data.frame(Cell_type = rep("myo", 20), IR = rep(c("IR", "non-IR"), each = 20), Expression = rep(c("myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 2), H3K4me3 = -1)
H3K4me3_5p[H3K4me3_5p$Expression == "myo_IR", "H3K4me3"] <- colMeans(myo_H3K4me3_5p[myo084_ir,], na.rm = T)
H3K4me3_5p[H3K4me3_5p$Expression == "myo_not-retained", "H3K4me3"] <- colMeans(myo_H3K4me3_5p[myo084_other,], na.rm = T)
H3K4me3_boundaries_v2 <- data.frame(rbind(H3K4me3_3p, H3K4me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K4me3_3p)), levels = c("5-prime", "3-prime")))
(H3K4me3_boundaries_v2_profile <- ggplot(H3K4me3_boundaries_v2, aes(x = Position, y = H3K4me3, group = Expression)) + 
   geom_line(aes(color = IR)) + 
   geom_point(aes(color = IR)) + 
   facet_grid(. ~ End) + 
   ggtitle("myo H3K4me3 profile around intron boundaries") + 
   ylab("Average H3K4me3 signal in myo") + 
   theme_bw())
ggsave(H3K4me3_boundaries_v2_profile, file = "H3K4me3_boundaries_v2_profile.pdf")

######################################################################################################
# H3K4me1 profile @ intron boundaries
lum_H3K4me1_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/lumRM080_H3K4me1.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K4me1_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM080_H3K4me1.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K4me1_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/lumRM080_H3K4me1.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K4me1_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM080_H3K4me1.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K4me1_3p <- sum(myo_H3K4me1_3p)/sum(lum_H3K4me1_3p) * lum_H3K4me1_3p
lum_H3K4me1_5p <- sum(myo_H3K4me1_5p)/sum(lum_H3K4me1_5p) * lum_H3K4me1_5p
H3K4me1_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K4me1 = -1)
H3K4me1_3p[H3K4me1_3p$Expression == "lum_IR", "H3K4me1"] <- colMeans(lum_H3K4me1_3p[lum084_ir,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$Expression == "lum_not-retained", "H3K4me1"] <- colMeans(lum_H3K4me1_3p[lum084_other,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$Expression == "myo_IR", "H3K4me1"] <- colMeans(myo_H3K4me1_3p[myo084_ir,], na.rm = T)
H3K4me1_3p[H3K4me1_3p$Expression == "myo_not-retained", "H3K4me1"] <- colMeans(myo_H3K4me1_3p[myo084_other,], na.rm = T)
H3K4me1_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K4me1 = -1)
H3K4me1_5p[H3K4me1_5p$Expression == "lum_IR", "H3K4me1"] <- colMeans(lum_H3K4me1_5p[lum084_ir,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$Expression == "lum_not-retained", "H3K4me1"] <- colMeans(lum_H3K4me1_5p[lum084_other,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$Expression == "myo_IR", "H3K4me1"] <- colMeans(myo_H3K4me1_5p[myo084_ir,], na.rm = T)
H3K4me1_5p[H3K4me1_5p$Expression == "myo_not-retained", "H3K4me1"] <- colMeans(myo_H3K4me1_5p[myo084_other,], na.rm = T)
H3K4me1_boundaries_v2 <- data.frame(rbind(H3K4me1_3p, H3K4me1_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K4me1_3p)), levels = c("5-prime", "3-prime")))
(H3K4me1_boundaries_v2_profile <- ggplot(H3K4me1_boundaries_v2, aes(x = Position, y = H3K4me1, group = Expression)) + 
   geom_line(aes(color = IR)) + 
   geom_point(aes(color = IR)) + 
   facet_grid(Cell_type ~ End) + 
   ggtitle("H3K4me1 profile around intron boundaries") + 
   ylab("Average H3K4me1 signal") + 
   theme_bw())
ggsave(H3K4me1_boundaries_v2_profile, file = "H3K4me1_boundaries_v2_profile.pdf")

######################################################################################################
# H3K9me3 profile @ intron boundaries
lum_H3K9me3_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/lumRM080_H3K9me3.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K9me3_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM080_H3K9me3.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K9me3_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/lumRM080_H3K9me3.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K9me3_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM080_H3K9me3.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K9me3_3p <- sum(myo_H3K9me3_3p)/sum(lum_H3K9me3_3p) * lum_H3K9me3_3p
lum_H3K9me3_5p <- sum(myo_H3K9me3_5p)/sum(lum_H3K9me3_5p) * lum_H3K9me3_5p
H3K9me3_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K9me3 = -1)
H3K9me3_3p[H3K9me3_3p$Expression == "lum_IR", "H3K9me3"] <- colMeans(lum_H3K9me3_3p[lum084_ir,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$Expression == "lum_not-retained", "H3K9me3"] <- colMeans(lum_H3K9me3_3p[lum084_other,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$Expression == "myo_IR", "H3K9me3"] <- colMeans(myo_H3K9me3_3p[myo084_ir,], na.rm = T)
H3K9me3_3p[H3K9me3_3p$Expression == "myo_not-retained", "H3K9me3"] <- colMeans(myo_H3K9me3_3p[myo084_other,], na.rm = T)
H3K9me3_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K9me3 = -1)
H3K9me3_5p[H3K9me3_5p$Expression == "lum_IR", "H3K9me3"] <- colMeans(lum_H3K9me3_5p[lum084_ir,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$Expression == "lum_not-retained", "H3K9me3"] <- colMeans(lum_H3K9me3_5p[lum084_other,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$Expression == "myo_IR", "H3K9me3"] <- colMeans(myo_H3K9me3_5p[myo084_ir,], na.rm = T)
H3K9me3_5p[H3K9me3_5p$Expression == "myo_not-retained", "H3K9me3"] <- colMeans(myo_H3K9me3_5p[myo084_other,], na.rm = T)
H3K9me3_boundaries_v2 <- data.frame(rbind(H3K9me3_3p, H3K9me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K9me3_3p)), levels = c("5-prime", "3-prime")))
(H3K9me3_boundaries_v2_profile <- ggplot(H3K9me3_boundaries_v2, aes(x = Position, y = H3K9me3, group = Expression)) + 
   geom_line(aes(color = IR)) + 
   geom_point(aes(color = IR)) + 
   facet_grid(Cell_type ~ End) + 
   ggtitle("H3K9me3 profile around intron boundaries") + 
   ylab("Average H3K9me3 signal") + 
   theme_bw())
ggsave(H3K9me3_boundaries_v2_profile, file = "H3K9me3_boundaries_v2_profile.pdf")

######################################################################################################
# H3K27me3 profile @ intron boundaries
lum_H3K27me3_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/lumRM080_H3K27me3.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K27me3_3p <- read.table("~/REMC/epiProfile/IR/introns3p_200/myoRM080_H3K27me3.hg19v65_introns_for_genes.3prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K27me3_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/lumRM080_H3K27me3.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
myo_H3K27me3_5p <- read.table("~/REMC/epiProfile/IR/introns5p_200/myoRM080_H3K27me3.hg19v65_introns_for_genes.5prime_200.profile", sep = " ", head = F, as.is = T, row.names = 1)[, -21]
lum_H3K27me3_3p <- sum(myo_H3K27me3_3p)/sum(lum_H3K27me3_3p) * lum_H3K27me3_3p
lum_H3K27me3_5p <- sum(myo_H3K27me3_5p)/sum(lum_H3K27me3_5p) * lum_H3K27me3_5p
H3K27me3_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K27me3 = -1)
H3K27me3_3p[H3K27me3_3p$Expression == "lum_IR", "H3K27me3"] <- colMeans(lum_H3K27me3_3p[lum084_ir,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$Expression == "lum_not-retained", "H3K27me3"] <- colMeans(lum_H3K27me3_3p[lum084_other,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$Expression == "myo_IR", "H3K27me3"] <- colMeans(myo_H3K27me3_3p[myo084_ir,], na.rm = T)
H3K27me3_3p[H3K27me3_3p$Expression == "myo_not-retained", "H3K27me3"] <- colMeans(myo_H3K27me3_3p[myo084_other,], na.rm = T)
H3K27me3_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*2), IR = rep(rep(c("IR", "non-IR"), each = 20), times = 2), Expression = rep(c("lum_IR", "lum_not-retained", "myo_IR", "myo_not-retained"), each = 20), Position = rep(seq(-190, 190, by = 20), times = 4), H3K27me3 = -1)
H3K27me3_5p[H3K27me3_5p$Expression == "lum_IR", "H3K27me3"] <- colMeans(lum_H3K27me3_5p[lum084_ir,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$Expression == "lum_not-retained", "H3K27me3"] <- colMeans(lum_H3K27me3_5p[lum084_other,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$Expression == "myo_IR", "H3K27me3"] <- colMeans(myo_H3K27me3_5p[myo084_ir,], na.rm = T)
H3K27me3_5p[H3K27me3_5p$Expression == "myo_not-retained", "H3K27me3"] <- colMeans(myo_H3K27me3_5p[myo084_other,], na.rm = T)
H3K27me3_boundaries_v2 <- data.frame(rbind(H3K27me3_3p, H3K27me3_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(H3K27me3_3p)), levels = c("5-prime", "3-prime")))
(H3K27me3_boundaries_v2_profile <- ggplot(H3K27me3_boundaries_v2, aes(x = Position, y = H3K27me3, group = Expression)) + 
   geom_line(aes(color = IR)) + 
   geom_point(aes(color = IR)) + 
   facet_grid(Cell_type ~ End) + 
   ggtitle("H3K27me3 profile around intron boundaries") + 
   ylab("Average H3K27me3 signal") + 
   theme_bw())
ggsave(H3K27me3_boundaries_v2_profile, file = "H3K27me3_boundaries_v2_profile.pdf")

######################################################################################################
# H3K36me3 signal for introns
# separate on gene RPKM for H3K36me3: <1; 1-10; 10-100; >100 
lum084gene <- read.delim("~/REMC/gene/A17918.G.A.rpkm.pc", head = F, as.is = T)
myo084gene <- read.delim("~/REMC/gene/A17919.G.A.rpkm.pc", head = F, as.is = T)
geneRPKM <- data.frame(gene = lum084gene$V1, RPKM = (lum084gene$V3 + myo084gene$V3)/2)
rownames(geneRPKM) <- geneRPKM$gene
introns_geneRPKM <- data.frame(intron = c(lum084_other, lum084_ir, myo084_ir, myo084_other))
introns_geneRPKM$gene <- gsub("_chr[0-9XYM]+:[0-9_<-]+", "", introns_geneRPKM$intron)
introns_geneRPKM$geneRPKM <- geneRPKM[introns_geneRPKM$gene, "RPKM"]
introns_geneRPKM <- na.omit(introns_geneRPKM)
introns_1 <- as.character(introns_geneRPKM[introns_geneRPKM$geneRPKM <= 1, "intron"])
introns_1_10 <- as.character(introns_geneRPKM[introns_geneRPKM$geneRPKM > 1 & introns_geneRPKM$geneRPKM <= 10, "intron"])
introns_10_100 <- as.character(introns_geneRPKM[introns_geneRPKM$geneRPKM > 10 & introns_geneRPKM$geneRPKM <= 100, "intron"])
introns_100 <- as.character(introns_geneRPKM[introns_geneRPKM$geneRPKM > 100, "intron"])
rm(lum084gene, myo084gene)

lum_H3K36me3_introns_v2 <- read.delim("~/REMC/epiProfile/IR/introns/hg19v65_introns_for_genes.lumRM080_H3K36me3.coverage", head = F, as.is = T)
myo_H3K36me3_introns_v2 <- read.delim("~/REMC/epiProfile/IR/introns/hg19v65_introns_for_genes.myoRM080_H3K36me3.coverage", head = F, as.is = T)
# normalize signal
(norm <- sum(myo_H3K36me3_introns_v2$V6)/sum(lum_H3K36me3_introns_v2$V6))
H3K36me3_introns_v2 <- data.frame(id = c(lum_H3K36me3_introns_v2$V4, myo_H3K36me3_introns_v2$V4), Cell_type = c(rep("lum", nrow(lum_H3K36me3_introns_v2)), rep("myo", nrow(myo_H3K36me3_introns_v2))), geneRPKM = NA, Expression = NA, H3K36me3 = c(lum_H3K36me3_introns_v2$V6 * norm, myo_H3K36me3_introns_v2$V6))
H3K36me3_introns_v2[H3K36me3_introns_v2$id %in% lum084_other & H3K36me3_introns_v2$Cell_type == "lum", "Expression"] <- "not-retained"
H3K36me3_introns_v2[H3K36me3_introns_v2$id %in% myo084_other & H3K36me3_introns_v2$Cell_type == "myo", "Expression"] <- "not-retained"
H3K36me3_introns_v2[H3K36me3_introns_v2$id %in% lum084_ir & H3K36me3_introns_v2$Cell_type == "lum", "Expression"] <- "IR"
H3K36me3_introns_v2[H3K36me3_introns_v2$id %in% myo084_ir & H3K36me3_introns_v2$Cell_type == "myo", "Expression"] <- "IR"
H3K36me3_introns_v2$Expression <- factor(H3K36me3_introns_v2$Expression)
H3K36me3_introns_v2[H3K36me3_introns_v2$id %in% introns_1, "geneRPKM"] <- "gene RPKM < 1"
H3K36me3_introns_v2[H3K36me3_introns_v2$id %in% introns_1_10, "geneRPKM"] <- "gene RPKM 1-10"
H3K36me3_introns_v2[H3K36me3_introns_v2$id %in% introns_10_100, "geneRPKM"] <- "gene RPKM 10-100"
H3K36me3_introns_v2[H3K36me3_introns_v2$id %in% introns_100, "geneRPKM"] <- "gene RPKM > 100"
H3K36me3_introns_v2$geneRPKM <- factor(H3K36me3_introns_v2$geneRPKM)
# # fold enrichment between utilized/un-utilized isoform introns
# mean(H3K36me3_introns_v2[H3K36me3_introns_v2$utilize == "lum.lum_IR" | H3K36me3_introns_v2$utilize == "myo.myo_IR", "H3K36me3"])/mean(H3K36me3_introns_v2[H3K36me3_introns_v2$utilize == "myo.lum_IR" | H3K36me3_introns_v2$utilize == "lum.myo_IR", "H3K36me3"])
H3K36me3_introns_v2 <- droplevels(na.omit(H3K36me3_introns_v2[H3K36me3_introns_v2$geneRPKM != "gene RPKM > 100", ]))
H3K36me3_introns_v2$group <- interaction(interaction(H3K36me3_introns_v2$Cell_type, H3K36me3_introns_v2$Expression), H3K36me3_introns_v2$geneRPKM)
H3K36me3_introns_v2_stat_v2 <- ddply(H3K36me3_introns_v2, ~ group, summarize, Cell_type = Cell_type[1], Expression = Expression[1], geneRPKM = geneRPKM[1], ymin = boxplot.stats(H3K36me3)$stats[1], lower = boxplot.stats(H3K36me3)$stats[2], middle = mean(H3K36me3), upper = boxplot.stats(H3K36me3)$stats[4], ymax = boxplot.stats(H3K36me3)$stats[5])
H3K36me3_introns_v2_stat_v2$Expression <- gsub("not-retained", "not retained", H3K36me3_introns_v2_stat_v2$Expression)
H3K36me3_introns_v2_stat_v2$Expression <- gsub("IR", "retained introns", H3K36me3_introns_v2_stat_v2$Expression)
H3K36me3_introns_v2_stat_v2$Expression <- factor(H3K36me3_introns_v2_stat_v2$Expression, levels = c("retained introns", "not retained"))

(H3K36me3_introns_v2_profile <- ggplot(H3K36me3_introns_v2_stat_v2, aes(x = Expression, group = group)) + 
   geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax, fill = Expression), color = "grey", stat = "identity", position = "dodge", width = 0.8) + 
   facet_grid(geneRPKM ~ Cell_type, scales = "free") + 
   # ggtitle("H3K36me3 signal for introns") + 
   xlab("intron retention") + 
   ylab("Average H3K36me3 signal") + 
   scale_fill_manual(name = "IR", values = c("retained introns" = "red", "not retained" = "blue")) + 
   theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 15, color = "black"), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 15, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
ggsave(H3K36me3_introns_v2_profile, file = "H3K36me3_introns_v2_profile.pdf", width = 10, height = 8)

t.test(H3K36me3_introns_v2[H3K36me3_introns_v2$group == "lum.IR.gene RPKM < 1", "H3K36me3"], H3K36me3_introns_v2[H3K36me3_introns_v2$group == "lum.not-retained.gene RPKM < 1", "H3K36me3"])
t.test(H3K36me3_introns_v2[H3K36me3_introns_v2$group == "myo.IR.gene RPKM < 1", "H3K36me3"], H3K36me3_introns_v2[H3K36me3_introns_v2$group == "myo.not-retained.gene RPKM < 1", "H3K36me3"])
t.test(H3K36me3_introns_v2[H3K36me3_introns_v2$group == "lum.IR.gene RPKM 1-10", "H3K36me3"], H3K36me3_introns_v2[H3K36me3_introns_v2$group == "lum.not-retained.gene RPKM 1-10", "H3K36me3"])
t.test(H3K36me3_introns_v2[H3K36me3_introns_v2$group == "myo.IR.gene RPKM 1-10", "H3K36me3"], H3K36me3_introns_v2[H3K36me3_introns_v2$group == "myo.not-retained.gene RPKM 1-10", "H3K36me3"])
t.test(H3K36me3_introns_v2[H3K36me3_introns_v2$group == "lum.IR.gene RPKM 10-100", "H3K36me3"], H3K36me3_introns_v2[H3K36me3_introns_v2$group == "lum.not-retained.gene RPKM 10-100", "H3K36me3"])
t.test(H3K36me3_introns_v2[H3K36me3_introns_v2$group == "myo.IR.gene RPKM 10-100", "H3K36me3"], H3K36me3_introns_v2[H3K36me3_introns_v2$group == "myo.not-retained.gene RPKM 10-100", "H3K36me3"])

save(lum084_other, myo084_other, lum084_ir, myo084_ir, introns_1, introns_1_10, introns_10_100, introns_100, H3K36me3_introns_v2_stat_v2, WGBS_CpG_v2, lum_correlate3p, lum_correlate5p, myo_correlate3p, myo_correlate5p, 
     WGBS_boundaries_v2, CpG_boundaries_v2, MeDIP_boundaries_v2, H3K4me3_boundaries_v2, H3K4me1_boundaries_v2, H3K9me3_boundaries_v2, H3K27me3_boundaries_v2, H3K36me3_introns_v2, file = "intronProfile_v2.Rdata")

