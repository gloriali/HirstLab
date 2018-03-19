# NHA VitC - comparing wt - mut - vitc

library(ggplot2)
library(plyr)
library(VennDiagram)
library(grid)
library(gridExtra)
library(gplots)
library(dendextend)
library(reshape2)
library(dplyr)
library(RCircos)
library(stringr)
library(scales)
source('~/HirstLab/Pipeline/R/enrich.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/wt_mut_vitc/")
load("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/wt_mut_vitc/wt_mut_vitc.Rdata")
RPKM <- read.delim("../RNAseq/RPKM/vitc.RPKM", as.is = T)
rownames(RPKM) <- RPKM$ENSG

## ========= 5mC modifiers RPKM =========
DNAme_regulators_RPKM <- read.delim("DNAme_regulators.RPKM", as.is = T) %>% select(-Name) %>% melt(id = c("ENSG", "gene")) 
(DNAme_regulators_RPKM_figure <- ggplot(DNAme_regulators_RPKM, aes(gene, log10(value), color = variable)) + 
		geom_point(position = position_jitter(width = 0.2)) + 
		coord_flip() + 
		guides(color = guide_legend(title = NULL)) + 
		xlab("") + 
		ylab("log10 RPKM") + 
		theme_bw())
ggsave(DNAme_regulators_RPKM_figure, file = "DNAme_regulators_RPKM_figure.pdf", height = 5, width = 6)

## ============= H3K27ac ==========
e <- 1e-4
H3K27ac_ER_union_RPKM <- read.delim("./H3K27ac/H3K27ac.ER.union.RPKM", as.is = T)
pdf("./H3K27ac/H3K27ac_ER_union_RPKM_figure.pdf", height = 5, width = 5)
plot(ecdf(H3K27ac_ER_union_RPKM$X.NHA_control.), xlim = c(0, 100), main = "Ecdf of H3K27ac union ERs RPKM", xlab = "RPKM", ylab = "ecdf")
lines(ecdf(H3K27ac_ER_union_RPKM$X.NHAR_control.), col = "red")
lines(ecdf(H3K27ac_ER_union_RPKM$X.NHAR_vitc.), col = "blue")
legend("bottomright", c("wt", "mut", "vitc"), col = c("black", "red", "blue"), lwd = 3, cex = 0.8)
plot(density(log2((H3K27ac_ER_union_RPKM$X.NHAR_vitc.+e)/(H3K27ac_ER_union_RPKM$X.NHAR_control.+e))), xlim = c(-2, 2), main = "Density of H3K27ac RPKM FC", xlab = "log2 RPKM FC", ylab = "density")
lines(density(log2((H3K27ac_ER_union_RPKM$X.NHAR_control.+e)/(H3K27ac_ER_union_RPKM$X.NHA_control.+e))), col = "red")
abline(v = 0)
legend("topright", c("vitc/mut", "mut/wt"), col = c("black", "red"), lwd = 3, cex = 0.8)
dev.off()
H3K27ac_DMR_summary <- read.delim("./H3K27ac/DMR.summary.stats", as.is = T)

## ============= DMR ============
DMR_DhMR_summary <- read.delim("./methylation/DMR.DhMR.summary.stats", as.is = T)
PBAL_delta_hMeDIP_FC <- read.delim("./methylation/DMR.vitc-mut.cor", as.is = T, head = F, col.names = c("id", "FC", "delta"))
(PBAL_delta_hMeDIP_FC_figure <- ggplot(PBAL_delta_hMeDIP_FC, aes(log2(FC), delta)) + 
		geom_point() + 
		coord_cartesian(xlim = c(-8, 8)) + 
		ggtitle("PBAL delta vs hMeDIP FC - vitC vs mut") + 
		theme_bw())
ggsave(PBAL_delta_hMeDIP_FC_figure, file = "./methylation/PBAL_delta_hMeDIP_FC_figure.pdf")
PBAL_450K_DMR <- read.delim("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/Chan/DMR.PBAL_mut-wt.bed", head = F, as.is = T, col.names = c("chr", "start", "end", "id", "450K", "PBAL"))
(PBAL_450K_DMR_figure <- ggplot(PBAL_450K_DMR, aes(X450K, PBAL)) + 
		geom_point(size = 0.5, alpha = 0.2) + 
		xlab("450K") + 
		ylab("PBAL") + 
		ggtitle("delta 5mC mut-wt") + 
		theme_bw())
ggsave(PBAL_450K_DMR_figure, file = "/projects/epigenomics3/epigenomics3_results/users/lli/NHA/Chan/DMR.PBAL_mut-wt.pdf")

## ============= save ===========
save(list = c(ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "enrich"), ls(pattern = "DAVID"), ls(pattern = "GREAT"), ls(pattern = "venn")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/NHA/wt_mut_vitc/wt_mut_vitc.Rdata")
