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
load("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/ChIPseq.Rdata")
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/")
marks <- c("H3K27me3", "H3K9me3", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3")

## ============= QC =============
QC_summary <- read.delim("./bam/QC_summary.txt", as.is = T) %>% mutate(Mark = gsub("_.*", "", Library), Sample = gsub(".*_NHA", "NHA", Library))
QC_summary_reads <- QC_summary %>% select(Mark, Sample, Total_Number_Of_Reads, Number_Reads_Aligned, Number_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = c("Mark", "Sample"))
(QC_summary_reads_figure <- ggplot(QC_summary_reads, aes(Sample, value/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ Mark) + 
		xlab("") + 
		ylab("No. of million reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_reads_figure, file = "./bam/QC_summary_reads_figure.pdf", height = 8, width = 8)
QC_summary_percent <- QC_summary %>% select(Mark, Sample, Mapping_Efficiency, Percent_of_dups, Percent_of_Paired_Alignments, Percent_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = c("Mark", "Sample"))
(QC_summary_percent_figure <- ggplot(QC_summary_percent, aes(Sample, value, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ Mark, scales = "free_y") + 
		xlab("") + 
		ylab("Percent of total number of reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_percent_figure, file = "./bam/QC_summary_percent_figure.pdf", height = 10, width = 8)

## ============= save ===========
save(list = c(ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "DAVID"), ls(pattern = "GREAT"), ls(pattern = "venn")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/ChIPseq.Rdata")

