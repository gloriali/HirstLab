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

## ============= hMeDIP =========
### -------- enrichment in genomic regions ------
genomic_breakdown <- read.delim("./hMeDIP/intersect/genomic.breakdown.summary", as.is = T) %>% 
	mutate(group = gsub("group.wt-mut-vitc.", "", Name), NCpG = NULL, Name = NULL)
genomic_breakdown_tall <- melt(genomic_breakdown, id = "group") 
(genomic_breakdown_figure <- ggplot(genomic_breakdown_tall, aes(variable, log2(value), fill = group)) + 
		geom_bar(stat = "identity", position = position_dodge(), width = 0.8) + 
		geom_hline(yintercept = c(-1, 1)) + 
		facet_wrap(~ group) + 
		guides(fill = "none") + 
		xlab("") + 
		ylab("Fold enrichment") + 
		coord_flip() + 
		theme_bw())
ggsave(genomic_breakdown_figure, file = "./hMeDIP/genomic_breakdown_figure.pdf", height = 7, width = 7)
### ------------ intersect with enhancer ---------------
(DN_UP_enhancer_enrich <- enrich_GREAT("DN-UP_enhancer", "DN-UP_enhancer", top = 10, dirIn = "./hMeDIP/enrich/", dirOut = "./hMeDIP/enrich/", height = 6, width = 7))
TF_ENSG <- read.delim("/home/lli/hg19/TF.ENSG", as.is = T) %>% merge(RPKM)
homer_DN_UP_enhancer <- read.delim("./hMeDIP/homer/DN-UP_enhancer/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 20, q.value..Benjamini. <= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) %>% merge(TF_ENSG)
(homer_DN_UP_enhancer_figure <- ggplot(homer_DN_UP_enhancer, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("5hmC DN-UP enhancere") + 
		theme_bw())
ggsave(homer_DN_UP_enhancer_figure, file = "./hMeDIP/homer/homer_DN_UP_enhancer_figure.pdf", height = 5, width = 5)
homer_DN_UP_enhancer_RPKM <- homer_DN_UP_enhancer %>% select(TF, MGG119_control:NHA_vitc) %>% melt(id = "TF") %>% mutate(TF = factor(TF, levels = levels(homer_DN_UP_enhancer$TF)))
(homer_DN_UP_enhancer_RPKM_figure <- ggplot(homer_DN_UP_enhancer_RPKM, aes(TF, log10(value), color = variable)) + 
		geom_point(size = 2, position = position_jitter(width = 0.1)) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("NHARvitc_NHARcontrol NHARcontrol unique") + 
		guides(color = guide_legend(title = "")) + 
		theme_bw() + 
		theme(axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(homer_DN_UP_enhancer_RPKM_figure, file = "./hMeDIP/homer/homer_DN_UP_enhancer_RPKM_figure.pdf", height = 5, width = 5)
pdf("./hMeDIP/homer/homer_DN_UP_enhancer_figure.pdf", height = 5, width = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(homer_DN_UP_enhancer_figure + theme(plot.margin = unit(c(1,0,1,1), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(homer_DN_UP_enhancer_RPKM_figure + ggtitle("") + theme(plot.margin = unit(c(1,1,1,-0.5), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()


## ============= save ===========
save(list = c(ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "enrich"), ls(pattern = "DAVID"), ls(pattern = "GREAT"), ls(pattern = "venn")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/NHA/wt_mut_vitc/wt_mut_vitc.Rdata")
