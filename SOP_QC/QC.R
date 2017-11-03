# QC for optimizing MeDIP and hMeDIP SOP

## QC for MeDIP positive and negative qPCR primers
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
setwd("/projects/epigenomics2/users/lli/")
meDIP_qPCR <- read.delim("meDIP_qPCR_primers_hg19.5mC.bed", as.is = T) %>%
  mutate(type = gsub("[0-9]+", "", ID), CpG_density = n/(end-start)*100)
(meDIP_qPCR_fractional_figure <- ggplot(meDIP_qPCR, aes(ID, fractional, fill = type)) + 
    geom_boxplot() + 
    scale_fill_manual(values = c("red", "blue")) + 
    xlab("") + 
    ylab("Fractional methylation") +
    theme_bw())
ggsave(meDIP_qPCR_fractional_figure, file = "meDIP_qPCR_fractional_figure.pdf", height = 4, width = 5)
meDIP_qPCR_CGdensity <- read.delim("meDIP_qPCR_primers_hg19.CpGdensity.bed", as.is = T, head = F, col.names = c("chr", "start", "end", "ID", "n", "density")) %>%
  mutate(type = gsub("[0-9]+", "", ID))
(meDIP_qPCR_CGdensity_figure <- ggplot(meDIP_qPCR_CGdensity, aes(ID, density, color = type)) + 
    geom_point(size = 3) + 
    scale_color_manual(values = c("red", "blue")) + 
    xlab("") + 
    ylab("CpG density - #CpG per 100bp") +
    theme_bw())
ggsave(meDIP_qPCR_CGdensity_figure, file = "meDIP_qPCR_CGdensity_figure.pdf", height = 4, width = 5)

## QC for SOP
setwd("/projects/epigenomics3/bams/hg19/UBC_meDIPQC_01NOV2017/")
QC_summary <- read.delim("SOP_QC_summary.txt", as.is = T) %>% mutate(Sample = gsub(".*-", "", gsub("_S.*", "", Library)), Type = gsub("IP-", "", gsub("-\\dX.*", "", Library)))
QC_summary_reads <- QC_summary %>% select(Sample, Type, Total_Number_Of_Reads, Number_Reads_Aligned, Number_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = c("Type", "Sample"))
(QC_summary_reads_figure <- ggplot(QC_summary_reads, aes(Sample, value/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ Type) + 
		guides(fill = "none") + 
		xlab("") + 
		ylab("No. of million reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_reads_figure, file = "QC_summary_reads_figure.pdf", height = 8, width = 5)
QC_summary_percent <- QC_summary %>% select(Sample, Type, Mapping_Efficiency, Percent_of_dups, Percent_of_Paired_Alignments, Percent_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = c("Type", "Sample"))
(QC_summary_percent_figure <- ggplot(QC_summary_percent, aes(Sample, value, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ Type, scales = "free_y") + 
		guides(fill = "none") + 
		xlab("") + 
		ylab("Percent of total number of reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_percent_figure, file = "QC_summary_percent_figure.pdf", height = 10, width = 5)

## coverage vs GC content
GC_coverage_summary <- read.delim("/projects/epigenomics3/bams/hg19/UBC_meDIPQC_01NOV2017/GCcontent.coverage.summary", as.is = T) %>% mutate(antibody = gsub(".*-", "", gsub("_S.*", "", sample)), type = gsub("IP-", "", gsub("-\\dX.*", "", sample)))
(GC_coverage_summary_figure <- ggplot(GC_coverage_summary, aes(GC, average, color = antibody)) + 
		geom_line() + 
		facet_wrap(~type, nrow = 2, scales = "free_y") + 
		xlab("GC content") + 
		ylab("Average coverage") + 
		theme_bw())
ggsave(GC_coverage_summary_figure, file = "/projects/epigenomics3/bams/hg19/UBC_meDIPQC_01NOV2017/GC_coverage_summary_figure.pdf")

## spike-in coverage
spikein_coverage <- read.delim("spike_in.coverage", as.is = T) %>% mutate(Sample = gsub(".*-", "", gsub("_S.*", "", Library)), Type = gsub("IP-", "", gsub("-\\dX.*", "", Library)), Genome = revalue(Genome, c("NC_001416.1" = "Lambda", "NC_001604.1" = "T7", "NC_003287.2" = "M13")))
(spikein_coverage_5mC_figure <- ggplot(spikein_coverage %>% filter(Type == "5mC"), aes(Sample, Coverage, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		guides(fill = "none") + 
		facet_grid(Genome ~ Type, scales = "free") + 
		xlab("") + 
		ylab("No. of uniquely aligned reads") + 
		theme_bw())
(spikein_coverage_5hmC_figure <- ggplot(spikein_coverage %>% filter(Type == "5hmC"), aes(Sample, Coverage, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		guides(fill = "none") + 
		facet_grid(Genome ~ Type, scales = "free") + 
		xlab("") + 
		ylab("") + 
		theme_bw())
pdf("spikein_coverage_figure.pdf", height = 6, width = 5)
grid.arrange(spikein_coverage_5mC_figure, spikein_coverage_5hmC_figure, nrow = 1)
dev.off()
