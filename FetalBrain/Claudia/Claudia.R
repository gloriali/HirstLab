# FetalBrain - data from Claudia

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
library(preprocessCore)
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/")
load("/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/Claudia.Rdata")

## ------- 450K array --------
DNAme <- read.csv("DNAme/180212-FB-betas-@lli.csv", row.names = 1) 
DNAme_norm <- DNAme %>% as.matrix() %>% normalize.quantiles() %>% as.data.frame()
colnames(DNAme_norm) <- colnames(DNAme)
DNAme_norm$ID <- row.names(DNAme)
DNAme_norm$diff.M291 <- DNAme_norm$GW17.M594 - DNAme_norm$GW13.M291
DNAme_norm$diff.M318 <- DNAme_norm$GW17.M594 - DNAme_norm$GW13.M318
DNAme_norm$diff.M446 <- DNAme_norm$GW17.M594 - DNAme_norm$GW13.M446
DNAme_norm$diff.M498 <- DNAme_norm$GW17.M594 - DNAme_norm$GW13.M498
cut <- 0.15 
DM_sum <- data.frame(GW13 = c("M291", "M318", "M446", "M498"), 
										 hyper = c(sum(DNAme_norm$diff.M291 > cut), sum(DNAme_norm$diff.M318 > cut), sum(DNAme_norm$diff.M446 > cut), sum(DNAme_norm$diff.M498 > cut)), 
										 hypo = c(sum(DNAme_norm$diff.M291 < -cut), sum(DNAme_norm$diff.M318 < -cut), sum(DNAme_norm$diff.M446 < -cut), sum(DNAme_norm$diff.M498 < -cut)))
(DM_sum_figure <- ggplot(DM_sum %>% melt(id = "GW13"), aes(GW13, value, fill = variable)) + 
		geom_bar(stat = "identity", position = position_dodge()) + 
		scale_fill_manual(values = c("hyper" = "red", "hypo" = "blue"), name = "DM in GW17") + 
		ylab("No. of probes") + 
		theme_bw())
ggsave(DM_sum_figure, file = "DM_sum_figure.pdf", height = 4, width = 4)
DNAme_ave <- data.frame(ID = DNAme_norm$ID)
for(GW in c("GW8", "GW9", "GW10", "GW11", "GW12", "GW13", "GW14", "GW17")){
	DNAme_ave[, GW] <- rowMeans(DNAme_norm %>% select(starts_with(GW)))
}
DNAme_ave$diff <- DNAme_ave$GW17 - DNAme_ave$GW13
DNAme_long <- DNAme_ave %>% select(ID, starts_with("GW")) %>% melt(id = "ID") %>% mutate(GW = factor(variable, levels = c("GW8", "GW9", "GW10", "GW11", "GW12", "GW13", "GW14", "GW17")))
Olig2_450K <- read.delim("DNAme/Olig2.promoter.450K.bed", head = F, as.is = T)
(Olig2_450K_figure <- ggplot(DNAme_long %>% filter(ID %in% Olig2_450K$V4), aes(GW, as.numeric(value), color = variable)) + 
		geom_point() + 
		facet_wrap(~ ID) + 
		guides(color = "none") + 
		xlab("") + 
		ylab("beta value") + 
		ggtitle("Olig2 promoter") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(Olig2_450K_figure, file = "Olig2_450K_figure.pdf", height = 3, width = 6)
GW_hyper_450K <- read.delim("DNAme/DMR.HuFNSC02_HuFNSC04.450K.hyper.bed", head = F, as.is = T)
(GW_hyper_450K_figure <- ggplot(DNAme_long %>% filter(ID %in% GW_hyper_450K$V4), aes(GW, value, color = variable)) + 
		geom_point() + 
		facet_wrap(~ ID) + 
		guides(color = "none") + 
		xlab("") + 
		ylab("beta value") + 
		ggtitle("GW17 hyper") + 
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)))
ggsave(GW_hyper_450K_figure, file = "GW_hyper_450K_figure.pdf", height = 8, width = 8)
GW_hypo_450K <- read.delim("DNAme/DMR.HuFNSC02_HuFNSC04.450K.hypo.bed", head = F, as.is = T)
cut <- 0.05
GW_hypo_450K_validate <- GW_hypo_450K %>% filter(V4 %in% DNAme_norm[DNAme_norm$diff.M291 < -cut | DNAme_norm$diff.M318 < -cut | DNAme_norm$diff.M446 < -cut | DNAme_norm$diff.M498 < -cut, "ID"])
write.table(GW_hypo_450K_validate, file = "GW_hypo_450K_validate.bed", sep = "\t", row.names = F, col.names = F, quote = F)
(GW_hypo_450K_validate_enrich <- enrich_GREAT("GW_hypo_450K_validate", "GW_hypo_450K_validate", dirIn = "./", dirOut = "./", width = 6))
(GW_hypo_450K_figure <- ggplot(DNAme_long %>% filter(ID %in% intersect(GW_hypo_450K$V4, DNAme_ave[DNAme_ave$diff < -0.01, "ID"])), aes(GW, value, color = ID)) + 
		geom_smooth(aes(group = ID), se = F) + 
		#facet_wrap(~ID) +
		guides(color = "none") + 
		xlab("") + 
		ylab("beta value") + 
		ggtitle("GW17 hypo") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(GW_hypo_450K_figure, file = "GW_hypo_450K_figure.pdf", height = 20, width = 20)
DNAme$sd <- apply(DNAme[, 2:ncol(DNAme)], 1, sd)
(GW_450K_figure <- ggplot(DNAme_long %>% filter(ID %in% DNAme[DNAme$sd>quantile(DNAme$sd, 0.8), "ID"]), aes(GW, value, color = variable)) + 
		geom_boxplot() + 
		guides(color = "none") + 
		xlab("") + 
		ylab("beta value") + 
		ggtitle("Top 20% most variable probes") + 
		theme_bw())
ggsave(GW_450K_figure, file = "GW_450K_figure.pdf", height = 4, width = 4)

## --------- RNAseq ----------
RPKM <- read.delim("RNAseq/RPKM/Claudia.RPKM", as.is = T) 
RPKM_long <- RPKM %>% melt(id = "ENSG") %>% mutate(GW = factor(gsub("\\..*", "", variable), levels = c("GW8", "GW10", "GW11", "GW12", "GW13", "GW14", "GW17")))
(Olig2_RPKM_figure <- ggplot(RPKM_long %>% filter(ENSG == "ENSG00000233863"), aes(GW, value, color = variable)) + 
		geom_point() + 
		guides(color = "none") + 
		xlab("") + 
		ylab("RPKM") + 
		ggtitle("Olig2") + 
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)))
ggsave(Olig2_RPKM_figure, file = "Olig2_RPKM_figure.pdf", height = 4, width = 4)
GW_17_13_DN <- read.delim("/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_DN_duplicated.txt", as.is = T, head = F)
(GW_17_13_DN_figure <- ggplot(RPKM_long %>% filter(ENSG %in% intersect(GW_17_13_DN$V1, RPKM[RPKM$GW17.M594-RPKM$GW13.M318 < -5, "ENSG"])), aes(GW, log10(value), color = variable)) + 
		geom_point() + 
		facet_wrap(~ENSG) + 
		guides(color = "none") + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("GW17 DN") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(GW_17_13_DN_figure, file = "GW_17_13_DN_figure.pdf", height = 8, width = 8)
GW_17_13_UP <- read.delim("/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_UP_duplicated.txt", as.is = T, head = F)
(GW_17_13_UP_figure <- ggplot(RPKM_long %>% filter(ENSG %in% intersect(GW_17_13_UP$V1, RPKM[RPKM$GW17.M594-RPKM$GW13.M318>1, "ENSG"])), aes(GW, log10(value), color = variable)) + 
		geom_boxplot() + 
		guides(color = "none") + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("GW17 UP") + 
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)))
ggsave(GW_17_13_UP_figure, file = "GW_17_13_UP_figure.pdf", height = 4, width = 4)
RPKM$sd <- apply(RPKM[, 2:ncol(RPKM)], 1, sd)
(GW_RPKM_figure <- ggplot(RPKM_long %>% filter(ENSG %in% RPKM[RPKM$sd>quantile(RPKM$sd, 0.8), "ENSG"]), aes(GW, value, color = variable)) + 
		geom_boxplot() + 
		coord_cartesian(ylim = c(0, 100)) +
		guides(color = "none") + 
		xlab("") + 
		ylab("RPKM") + 
		ggtitle("Top 20% most variable genes") + 
		theme_bw())
ggsave(GW_RPKM_figure, file = "GW_RPKM_figure.pdf", height = 4, width = 4)
BrainSpan <- read.csv("./BrainSpan/expression_matrix.csv", head = F, row.names = 1)
row_metadata <- read.csv("./BrainSpan/rows_metadata.csv")
rownames(BrainSpan) <- row_metadata$ensembl_gene_id
col_metadata <- read.csv("./BrainSpan/columns_metadata.csv") %>% mutate(ID = paste0(donor_name, "_", age, "_", structure_acronym))
colnames(BrainSpan) <- col_metadata$ID
BrainSpan$ID <- rownames(BrainSpan)
BrainSpan_long <- BrainSpan %>% select(ID, contains("pcw")) %>% melt(id = "ID") %>% mutate(donor = gsub("_.*", "", variable), structure = gsub(".*_", "", variable), GW = as.numeric(gsub(".*_", "", gsub(" pcw.*", "", variable))) + 2) %>%
	filter(GW < 20, !(structure %in% c("CB", "CBC", "CGE", "LGE", "MGE", "DTH", "M1C", "M1C-S1C", "STR", "MD", "HIP", "Ocx", "PCx", "S1C", "TCx", "URL")))
BrainSpan$GW18mean <- BrainSpan %>% select(contains("16 pcw")) %>% rowMeans()
BrainSpan$GW14mean <- BrainSpan %>% select(contains("12 pcw")) %>% rowMeans()
(BrainSpan_GW_17_13_DN_figure <- ggplot(BrainSpan_long %>% filter(ID %in% intersect(GW_17_13_DN$V1, BrainSpan[BrainSpan$GW14mean - BrainSpan$GW18mean > 1, "ID"])), aes(as.character(GW), log10(value), color = variable)) + 
		geom_point() + 
		facet_grid(ID ~ structure) + 
		guides(color = "none") + 
		xlab("GW") + 
		ylab("log10 RPKM") + 
		ggtitle("GW17 DN") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(BrainSpan_GW_17_13_DN_figure, file = "BrainSpan_GW_17_13_DN_figure.pdf", height = 20, width = 20)
(BrainSpan_GW_17_13_UP_figure <- ggplot(BrainSpan_long %>% filter(ID %in% intersect(GW_17_13_UP$V1, BrainSpan[BrainSpan$GW18mean - BrainSpan$GW14mean > 5, "ID"])), aes(as.character(GW), log10(value), color = variable)) + 
		geom_boxplot() + 
		facet_wrap(~ structure) + 
		guides(color = "none") + 
		xlab("GW") + 
		ylab("log10 RPKM") + 
		ggtitle("GW17 UP") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(BrainSpan_GW_17_13_UP_figure, file = "BrainSpan_GW_17_13_UP_figure.pdf", height = 6, width = 7)


save(list = ls(pattern = "figure"), file = "/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/Claudia.Rdata")
