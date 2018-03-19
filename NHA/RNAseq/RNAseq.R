# NHA VitC - ChIP-seq analysis

library(ggplot2)
library(plyr)
library(VennDiagram)
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
load("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/RNAseq.Rdata")
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/")
RPKM <- read.delim("./RPKM/vitc.RPKM", as.is = T)
# name <- read.delim("/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.EnsID_sorted.HUGO", head = F, as.is = T, col.names = c("ENSG", "Name"))
# RPKM <- merge(name, RPKM)
# write.table(RPKM, file = "./RPKM/vitc.RPKM", sep = "\t", col.names = T, row.names = F, quote = F)

## ============= QC =============
QC_summary <- read.delim("./bam/QC_summary.txt", as.is = T) %>% mutate(Sample = gsub("MGG_", "MGG119_", Library))
QC_summary_reads <- QC_summary %>% select(Sample, Total_Number_Of_Reads, Number_Reads_Aligned, Number_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = "Sample")
(QC_summary_reads_figure <- ggplot(QC_summary_reads, aes(Sample, value/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ .) + 
		guides(fill = "none") + 
		xlab("") + 
		ylab("No. of million reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_reads_figure, file = "./bam/QC_summary_reads_figure.pdf", height = 8, width = 3)
QC_summary_percent <- QC_summary %>% select(Sample, Mapping_Efficiency, Percent_of_dups, Percent_of_Paired_Alignments, Percent_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = "Sample")
(QC_summary_percent_figure <- ggplot(QC_summary_percent, aes(Sample, value, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ ., scales = "free_y") + 
		guides(fill = "none") + 
		xlab("") + 
		ylab("Percent of total number of reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_percent_figure, file = "./bam/QC_summary_percent_figure.pdf", height = 10, width = 3)

## ============= DE =============
MGG119vitc_MGG119control_UP <- read.delim("./DEfine/UP.MGG119_vitc_MGG119_control.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = c("ENSG", "RPKM1", "RPKM2", "FDR"))
MGG119vitc_MGG119control_DN <- read.delim("./DEfine/DN.MGG119_vitc_MGG119_control.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = c("ENSG", "RPKM1", "RPKM2", "FDR"))
NHARvitc_NHARcontrol_UP <- read.delim("./DEfine/UP.NHAR_vitc_NHAR_control.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = c("ENSG", "RPKM1", "RPKM2", "FDR"))
NHARvitc_NHARcontrol_DN <- read.delim("./DEfine/DN.NHAR_vitc_NHAR_control.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = c("ENSG", "RPKM1", "RPKM2", "FDR"))
NHARcontrol_NHAcontrol_UP <- read.delim("./DEfine/UP.NHAR_control_NHA_control.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = c("ENSG", "RPKM1", "RPKM2", "FDR"))
NHARcontrol_NHAcontrol_DN <- read.delim("./DEfine/DN.NHAR_control_NHA_control.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = c("ENSG", "RPKM1", "RPKM2", "FDR"))
### intersect
vitc_control_UP <- list(MGG119 = MGG119vitc_MGG119control_UP$ENSG, NHAR = NHARvitc_NHARcontrol_UP$ENSG)
vitc_control_UP_venn <- venn.diagram(vitc_control_UP, filename = NULL, fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./DEfine/vitc_control_UP_venn.pdf", height = 3, width = 4)
grid.draw(vitc_control_UP_venn)
dev.off()
vitc_control_DN <- list(MGG119 = MGG119vitc_MGG119control_DN$ENSG, NHAR = NHARvitc_NHARcontrol_DN$ENSG)
vitc_control_DN_venn <- venn.diagram(vitc_control_DN, filename = NULL, fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./DEfine/vitc_control_DN_venn.pdf", height = 3, width = 4)
grid.draw(vitc_control_DN_venn)
dev.off()
vitc_UP_mut_DN <- list(vitc_UP = NHARvitc_NHARcontrol_UP$ENSG, mut_DN = NHARcontrol_NHAcontrol_DN$ENSG)
vitc_UP_mut_DN_venn <- venn.diagram(vitc_UP_mut_DN, filename = NULL, category.names = c("NHAR vitc_control UP", "NHAR_NHA control DN"), margin = 0.1, fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05)
pdf("./DEfine/vitc_UP_mut_DN_venn.pdf", height = 3, width = 4)
grid.draw(vitc_UP_mut_DN_venn)
dev.off()
vitc_DN_mut_UP <- list(vitc_DN = NHARvitc_NHARcontrol_DN$ENSG, mut_UP = NHARcontrol_NHAcontrol_UP$ENSG)
vitc_DN_mut_UP_venn <- venn.diagram(vitc_DN_mut_UP, filename = NULL, category.names = c("NHAR vitc_control DN", "NHAR_NHA control UP"), margin = 0.1, fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05)
pdf("./DEfine/vitc_DN_mut_UP_venn.pdf", height = 3, width = 4)
grid.draw(vitc_DN_mut_UP_venn)
dev.off()
write.table(merge(NHARvitc_NHARcontrol_UP, NHARcontrol_NHAcontrol_DN, by = "ENSG"), file = "./DEfine/vitc_UP_mut_DN.txt", sep = "\t", col.names = T, row.names = F, quote = F)
### DAVID
(enrich_NHAR_NHA_control_UP <- enrich(name = "UP.NHAR_control_NHA_control", dirIn = "./DEfine/enrich/", dirOut = "./DEfine/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 7, width = 8))
(enrich_NHAR_NHA_control_DN <- enrich(name = "DN.NHAR_control_NHA_control", dirIn = "./DEfine/enrich/", dirOut = "./DEfine/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 10, width = 8))
(enrich_NHAR_vitc_control_UP <- enrich(name = "UP.NHAR_vitc_NHAR_control", dirIn = "./DEfine/enrich/", dirOut = "./DEfine/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 7, width = 8))
(enrich_NHAR_vitc_control_DN <- enrich(name = "DN.NHAR_vitc_NHAR_control", dirIn = "./DEfine/enrich/", dirOut = "./DEfine/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 4, width = 8))
(enrich_MGG119_vitc_control_UP <- enrich(name = "UP.MGG119_vitc_MGG119_control", dirIn = "./DEfine/enrich/", dirOut = "./DEfine/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 7, width = 8))
(enrich_MGG119_vitc_control_DN <- enrich(name = "DN.MGG119_vitc_MGG119_control", dirIn = "./DEfine/enrich/", dirOut = "./DEfine/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 2, width = 8))
(enrich_vitc_UP_mut_DN <- enrich(name = "vitc_UP_mut_DN", dirIn = "./DEfine/enrich/", dirOut = "./DEfine/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 7, width = 8))

## ============= save ===========
save(list = c("RPKM", ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "enrich"), ls(pattern = "GREAT"), ls(pattern = "venn")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/RNAseq.Rdata")

