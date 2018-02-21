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
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/")
load("/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/Claudia.Rdata")

## ------- 450K array --------
DNAme <- read.csv("DNAme/180212-FB-betas-@lli.csv") 
DNAme_long <- DNAme %>% select(-starts_with("NA")) %>% melt(id = "ID") %>% mutate(GW = factor(gsub("\\..*", "", variable), levels = c("GW8", "GW9", "GW10", "GW11", "GW12", "GW13", "GW14", "GW17")))
Olig2_450K <- read.delim("DNAme/Olig2.promoter.450K.bed", head = F, as.is = T)
(Olig2_450K_figure <- ggplot(DNAme_long %>% filter(ID %in% Olig2_450K$V4), aes(GW, value, color = variable)) + 
		geom_boxplot() + 
		guides(color = "none") + 
		xlab("") + 
		ylab("beta value") + 
		ggtitle("Olig2 promoter") + 
		theme_bw())
ggsave(Olig2_450K_figure, file = "Olig2_450K_figure.pdf", height = 4, width = 4)
GW_hyper_450K <- read.delim("DNAme/DMR.HuFNSC02_HuFNSC04.450K.hyper.bed", head = F, as.is = T)
(GW_hyper_450K_figure <- ggplot(DNAme_long %>% filter(ID %in% intersect(GW_hyper_450K$V4, DNAme[DNAme$GW17.M594-DNAme$GW13.M291>0.001, "ID"])), aes(GW, value, color = variable)) + 
		geom_boxplot() + 
		guides(color = "none") + 
		xlab("") + 
		ylab("beta value") + 
		ggtitle("GW17 hyper") + 
		theme_bw())
ggsave(GW_hyper_450K_figure, file = "GW_hyper_450K_figure.pdf", height = 4, width = 4)
GW_hypo_450K <- read.delim("DNAme/DMR.HuFNSC02_HuFNSC04.450K.hypo.bed", head = F, as.is = T)
(GW_hypo_450K_figure <- ggplot(DNAme_long %>% filter(ID %in% intersect(GW_hypo_450K$V4, DNAme[DNAme$GW17.M594-DNAme$GW13.M291<0.001, "ID"])), aes(GW, value, color = variable)) + 
		geom_boxplot() + 
		guides(color = "none") + 
		xlab("") + 
		ylab("beta value") + 
		ggtitle("GW17 hypo") + 
		theme_bw())
ggsave(GW_hypo_450K_figure, file = "GW_hypo_450K_figure.pdf", height = 4, width = 4)
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
		theme_bw())
ggsave(Olig2_RPKM_figure, file = "Olig2_RPKM_figure.pdf", height = 4, width = 4)
GW_17_13_DN <- read.delim("/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_DN_duplicated.txt", as.is = T, head = F)
(GW_17_13_DN_figure <- ggplot(RPKM_long %>% filter(ENSG %in% GW_17_13_DN$V1), aes(GW, log10(value), color = variable)) + 
		geom_boxplot() + 
		guides(color = "none") + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("GW17 DN") + 
		theme_bw())
ggsave(GW_17_13_DN_figure, file = "GW_17_13_DN_figure.pdf", height = 4, width = 4)
GW_17_13_UP <- read.delim("/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_UP_duplicated.txt", as.is = T, head = F)
(GW_17_13_UP_figure <- ggplot(RPKM_long %>% filter(ENSG %in% intersect(GW_17_13_UP$V1, RPKM[RPKM$GW17.M594-RPKM$GW13.M318>1, "ENSG"])), aes(GW, log10(value), color = variable)) + 
		geom_boxplot() + 
		guides(color = "none") + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("GW17 UP") + 
		theme_bw())
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

save(list = ls(pattern = "figure"), file = "/projects/epigenomics3/epigenomics3_results/users/lli/Claudia/Claudia.Rdata")
