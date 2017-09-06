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
QC_target_region <- read.delim("./bam/QC.target_regions.txt", as.is = T)
QC_summary_domain <- merge(QC_summary, QC_target_region) %>% mutate(Percent_Domain_Reads = Number_Target_Regions*100/Number_Reads_Aligned)  %>% filter(!grepl("input", Mark))
(QC_summary_domain_figure <- ggplot(QC_summary_domain, aes(Sample, Percent_Domain_Reads, fill = Sample)) + 
		geom_bar(stat = "identity", width = 0.5) + 
		facet_grid(Mark ~ ., scales = "free_y") + 
		xlab("") + 
		ylab("Domain reads (% of mapped)") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_domain_figure, file = "./bam/QC_summary_domain_figure.pdf", height = 7, width = 4)
NHAR_H3K9me3_H3K4me3 <- read.delim("./bam/NHAR.H3K9me3.H3K4me3.bed", as.is = T, head = F, col.names = c("chr", "start", "end", "ID", "N", "Sample", "RPKM")) %>% select(ID, Sample, RPKM) %>%
	reshape(direction = "wide", idvar = "ID", timevar = "Sample")
(NHAR_H3K9me3_H3K4me3_figure <- ggplot(NHAR_H3K9me3_H3K4me3, aes(log10(RPKM.H3K9me3_NHAR_vitc), log10(RPKM.H3K4me3_NHAR_vitc))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.1, 0.2, 0.3, 0.04, 0.5, 0.8, 1), guide = "none") +	
		coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) + 
		ggtitle(cor(NHAR_H3K9me3_H3K4me3$RPKM.H3K9me3_NHAR_vitc, NHAR_H3K9me3_H3K4me3$RPKM.H3K4me3_NHAR_vitc)) + 
		theme_bw())
ggsave(NHAR_H3K9me3_H3K4me3_figure, file = "./bam/NHAR_H3K9me3_H3K4me3_figure.pdf", width = 5, height = 5)
(H3K9me3_NHAR_vitC_NTC_figure <- ggplot(NHAR_H3K9me3_H3K4me3, aes(log10(RPKM.H3K9me3_NHAR_vitc), log10(RPKM.H3K9me3_NHAR_control))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.1, 0.2, 0.3, 0.04, 0.5, 0.8, 1), guide = "none") +	
		coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) + 
		ggtitle(cor(NHAR_H3K9me3_H3K4me3$RPKM.H3K9me3_NHAR_vitc, NHAR_H3K9me3_H3K4me3$RPKM.H3K9me3_NHAR_control)) + 
		theme_bw())
ggsave(H3K9me3_NHAR_vitC_NTC_figure, file = "./bam/H3K9me3_NHAR_vitC_NTC_figure.pdf", width = 5, height = 5)
NHAR_H3K9me3_H3K4me3_venn <- draw.pairwise.venn(area1 = 44496864, area2 = 54317631, cross.area = 41973315, category = c("H3K9me3_NHAR_vitc", "H3K4me3_NHAR_vitc"), fill = c("red", "blue"))
H3K9me3_NHAR_vitC_NTC_venn <- draw.pairwise.venn(area1 = 44496864, area2 = 121812813, cross.area = 41856862, category = c("H3K9me3_NHAR_vitc", "H3K9me3_NHAR_control"), fill = c("red", "blue"))

## ============= ER =============
ER_summary_MACS2 <- read.delim("./MACS2/ER_summary.txt", as.is = T)
(ER_summary_MACS2_figure <- ggplot(ER_summary_MACS2, aes(Sample, Total_length/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(Mark ~., scales = "free_y") + 
		xlab("") + 
		ylab("Total length (MB)") + 
		ggtitle("MACS2 ER") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(ER_summary_MACS2_figure, file = "./MACS2/ER_summary_MACS2_figure.pdf", height = 6, width = 5)
ER_summary_FindER <- read.delim("./FindER/ER_summary.txt", as.is = T)
(ER_summary_FindER_figure <- ggplot(ER_summary_FindER, aes(Sample, Total_length/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(Mark ~., scales = "free_y") + 
		xlab("") + 
		ylab("Total length (MB)") + 
		ggtitle("FindER ER") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(ER_summary_FindER_figure, file = "./FindER/ER_summary_FindER_figure.pdf", height = 6, width = 5)

## ============= save ===========
save(list = c(ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "DAVID"), ls(pattern = "GREAT"), ls(pattern = "venn")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/ChIPseq.Rdata")

