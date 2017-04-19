# CTCF in glioma

library(ggplot2)
library(plyr)
library(VennDiagram)
library(gridExtra)
library(gplots)
library(dendextend)
library(reshape2)
library(wq)
library(dplyr)
library(RCircos)
library(stringr)
source('~/HirstLab/Pipeline/R/enrich.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
load("/projects/epigenomics2/users/lli/glioma/CTCF/CTCF.Rdata")
setwd("/projects/epigenomics2/users/lli/glioma/CTCF")
libs <- c("CEMT_19", "CEMT_21", "CEMT_22", "CEMT_23", "CEMT_47")

## ------- 5mC at CTCF loss regions -------
CTCF_loss_5mC <- read.delim("./WGBS/CTCF_IDHwt_unique.bed.5mC", as.is = T) %>% mutate(type = ifelse(grepl("NPC", sample), "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = gsub("NPC.", "", sample))
(CTCF_loss_5mC_figure <- ggplot(CTCF_loss_5mC, aes(sample, fractional, fill = type)) + 
		geom_boxplot() + 
		xlab("") + 
		ylab("Fractional methylation") + 
		ggtitle("CFCT loss") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CTCF_loss_5mC_figure, file = "./WGBS/CTCF_loss_5mC_figure.pdf", height = 4, width = 5)
CTCF_loss_5mC_wide <- CTCF_loss_5mC %>% select(ID, fractional, sample) %>% dcast(ID ~ sample, value.var = "fractional") %>% 
	mutate(chr = gsub(":.*", "", ID), start = gsub(".*:", "", gsub("-.*", "", ID)), end = gsub(".*-", "", ID), mut = (CEMT_19 + CEMT_22 + CEMT_47)/3, wt = CEMT_23, diff = mut - wt)
responder <- CTCF_loss_5mC_wide %>% filter(diff >= 0.2) %>% select(chr, start, end, ID, mut, wt, diff)
nonresponder <- CTCF_loss_5mC_wide %>% filter(diff <= 0) %>% select(chr, start, end, ID, mut, wt, diff)
write.table(responder, file = "./WGBS/responder.bed", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(nonresponder, file = "./WGBS/nonresponder.bed", sep = "\t", quote = F, row.names = F, col.names = F)
CTCF_gain_5mC <- read.delim("./WGBS/CTCF_IDHmut_unique.bed.5mC", as.is = T) %>% mutate(type = ifelse(grepl("NPC", sample), "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = gsub("NPC.", "", sample))
(CTCF_gain_5mC_figure <- ggplot(CTCF_gain_5mC, aes(sample, fractional, fill = type)) + 
		geom_boxplot() + 
		xlab("") + 
		ylab("Fractional methylation") + 
		ggtitle("CFCT gain") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CTCF_gain_5mC_figure, file = "./WGBS/CTCF_gain_5mC_figure.pdf", height = 4, width = 5)

## ------- H3K36me3 at CTCF loss regions -------
e <- 1e-4
CTCF_nonresponder_K36 <- read.delim("./H3K36me3/resonpder.H3K36me3.bed", as.is = T, col.names = c("chr", "start", "end", "ID", "mut", "wt", "diff", "N", "sample", "RPKM")) %>% 
	mutate(type = ifelse(grepl("NPC", sample), "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = gsub("NPC_", "", sample))
(CTCF_nonresponder_K36_figure <- ggplot(CTCF_nonresponder_K36, aes(sample, log10(RPKM), fill = type)) + 
		geom_boxplot() + 
		xlab("") + 
		ylab("log10 H3K36me3 RPKM") + 
		ggtitle("non-responders") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CTCF_nonresponder_K36_figure, file = "./H3K36me3/CTCF_nonresponder_K36_figure.pdf", height = 4, width = 5)
CTCF_nonresponder_K36_wide <- rbind(CTCF_nonresponder_K36 %>% filter(sample %in% c("CEMT_19", "CEMT_23")) %>%  dcast(ID + diff ~ sample, value.var = "RPKM") %>% 
																 	mutate(FC = (CEMT_19+e)/(CEMT_23+e), sample = "CEMT_19") %>% select(ID, diff, FC, sample), 
																 CTCF_nonresponder_K36 %>% filter(sample %in% c("CEMT_22", "CEMT_23")) %>%  dcast(ID + diff ~ sample, value.var = "RPKM") %>%
																 	mutate(FC = (CEMT_22+e)/(CEMT_23+e), sample = "CEMT_22") %>% select(ID, diff, FC, sample), 
																 CTCF_nonresponder_K36 %>% filter(sample %in% c("CEMT_47", "CEMT_23")) %>%  dcast(ID + diff ~ sample, value.var = "RPKM") %>%
																 	mutate(FC = (CEMT_47+e)/(CEMT_23+e), sample = "CEMT_47") %>% select(ID, diff, FC, sample))
(CTCF_nonresponder_K36_2D_figure <- ggplot(CTCF_nonresponder_K36_wide, aes(diff, log2(FC))) + 
		geom_point(size = 0.5) + 
		facet_wrap(~ sample, nrow = 3) + 
		xlab("5mC IDHmut-IDHwt") + 
		ylab("H3K36me3 RPKM log2 IDHmut/IDHwt") + 
		ggtitle("non-responders") + 
		theme_bw())
ggsave(CTCF_nonresponder_K36_2D_figure, file = "./H3K36me3/CTCF_nonresponder_K36_2D_figure.pdf", height = 9, width = 7)
CTCF_responder_K36 <- read.delim("./H3K36me3/resonpder.H3K36me3.bed", as.is = T, col.names = c("chr", "start", "end", "ID", "mut", "wt", "diff", "N", "sample", "RPKM")) %>% 
	mutate(type = ifelse(grepl("NPC", sample), "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = gsub("NPC_", "", sample))
(CTCF_responder_K36_figure <- ggplot(CTCF_responder_K36, aes(sample, log10(RPKM), fill = type)) + 
		geom_boxplot() + 
		xlab("") + 
		ylab("log10 H3K36me3 RPKM") + 
		ggtitle("responders") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CTCF_responder_K36_figure, file = "./H3K36me3/CTCF_responder_K36_figure.pdf", height = 4, width = 5)
CTCF_responder_K36_wide <- rbind(CTCF_responder_K36 %>% filter(sample %in% c("CEMT_19", "CEMT_23")) %>%  dcast(ID + diff ~ sample, value.var = "RPKM") %>% 
																 	mutate(FC = (CEMT_19+e)/(CEMT_23+e), sample = "CEMT_19") %>% select(ID, diff, FC, sample), 
	CTCF_responder_K36 %>% filter(sample %in% c("CEMT_22", "CEMT_23")) %>%  dcast(ID + diff ~ sample, value.var = "RPKM") %>%
		mutate(FC = (CEMT_22+e)/(CEMT_23+e), sample = "CEMT_22") %>% select(ID, diff, FC, sample), 
	CTCF_responder_K36 %>% filter(sample %in% c("CEMT_47", "CEMT_23")) %>%  dcast(ID + diff ~ sample, value.var = "RPKM") %>%
		mutate(FC = (CEMT_47+e)/(CEMT_23+e), sample = "CEMT_47") %>% select(ID, diff, FC, sample))
(CTCF_responder_K36_2D_figure <- ggplot(CTCF_responder_K36_wide, aes(diff, log2(FC))) + 
		geom_point(size = 0.5) + 
		facet_wrap(~ sample, nrow = 3) + 
		xlab("5mC IDHmut-IDHwt") + 
		ylab("H3K36me3 RPKM log2 IDHmut/IDHwt") + 
		ggtitle("responders") + 
		theme_bw())
ggsave(CTCF_responder_K36_2D_figure, file = "./H3K36me3/CTCF_responder_K36_2D_figure.pdf", height = 9, width = 7)

## ------- H3K27me3 at CTCF loss regions -------
CTCF_nonresponder_K27 <- read.delim("./H3K27me3/resonpder.H3K27me3.bed", as.is = T, col.names = c("chr", "start", "end", "ID", "mut", "wt", "diff", "N", "sample", "RPKM")) %>% 
	mutate(type = ifelse(grepl("NPC", sample), "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = gsub("NPC_", "", sample))
(CTCF_nonresponder_K27_figure <- ggplot(CTCF_nonresponder_K27, aes(sample, log10(RPKM), fill = type)) + 
		geom_boxplot() + 
		xlab("") + 
		ylab("log10 H3K27me3 RPKM") + 
		ggtitle("non-responders") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CTCF_nonresponder_K27_figure, file = "./H3K27me3/CTCF_nonresponder_K27_figure.pdf", height = 4, width = 5)
CTCF_responder_K27 <- read.delim("./H3K27me3/resonpder.H3K27me3.bed", as.is = T, col.names = c("chr", "start", "end", "ID", "mut", "wt", "diff", "N", "sample", "RPKM")) %>% 
	mutate(type = ifelse(grepl("NPC", sample), "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = gsub("NPC_", "", sample))
(CTCF_responder_K27_figure <- ggplot(CTCF_responder_K27, aes(sample, log10(RPKM), fill = type)) + 
		geom_boxplot() + 
		xlab("") + 
		ylab("log10 H3K27me3 RPKM") + 
		ggtitle("responders") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CTCF_responder_K27_figure, file = "./H3K27me3/CTCF_responder_K27_figure.pdf", height = 4, width = 5)

## ---------- H3K27ac and RNA-seq ----------
e <- 1e-4
gene <- read.delim("/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.EnsID.HUGO_sorted", head = F, as.is = T, row.names = 1)
H3K27ac_TSS_CTCF_loss <- read.delim("./H3K27ac/H3K27ac.TSS.CTCF_IDHwt_unique.bed.RPKM", as.is = T) %>% 
	filter(sample %in% c("CEMT_19", "CEMT_22", "CEMT_47")) %>% arrange(ENSG, d) %>% distinct(ENSG, CTCF, sample) %>%
	mutate(IDHmut = ifelse(sample == "CEMT_19", CEMT_19, ifelse(sample == "CEMT_22", CEMT_22, CEMT_47)), logFC = log2((IDHmut + e)/(CEMT_23 + e)), Gene = gene[ENSG, "V2"]) %>% 
	filter(logFC >= 1) %>% select(ENSG, Gene, CTCF, enhancer, d, IDHmut, CEMT_23, logFC) %>% distinct(ENSG, CTCF) %>% filter(d <= 50000) %>% arrange(desc(logFC))
H3K27ac_TSS_CTCF_gain <- read.delim("./H3K27ac/H3K27ac.TSS.CTCF_IDHmut_unique.bed.RPKM", as.is = T) %>% 
	filter(sample == "CEMT_23") %>% arrange(ENSG, d) %>% distinct(ENSG, CTCF, sample) %>% 
	mutate(IDHmut = (CEMT_19 + CEMT_22 + CEMT_47)/3, logFC = log2((IDHmut + e)/(CEMT_23 + e)), Gene = gene[ENSG, "V2"]) %>% filter(logFC <= -1) %>% 
	select(ENSG, Gene, CTCF, enhancer, d, IDHmut, CEMT_23, logFC) %>% filter(d <= 50000) %>% arrange(logFC)


save(list = c(ls(pattern = "figure"), "CTCF_loss_5mC_wide", "responder", "nonresponder", "H3K27ac_TSS_CTCF_gain", "H3K27ac_TSS_CTCF_loss"), 
		 file = "/projects/epigenomics2/users/lli/glioma/CTCF/CTCF.Rdata")

