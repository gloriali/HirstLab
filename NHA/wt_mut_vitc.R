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



## ============= save ===========
save(list = c(ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "enrich"), ls(pattern = "DAVID"), ls(pattern = "GREAT"), ls(pattern = "venn")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/NHA/wt_mut_vitc/wt_mut_vitc.Rdata")
