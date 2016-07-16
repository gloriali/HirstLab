# ChIPseq analysis for PB - Nada Jabado's sample

library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
load("/projects/epigenomics2/Brain/NJabado/NJabado.Rdata")
setwd("/projects/epigenomics2/Brain/NJabado/")

## -------- QC --------
ChIPseq_QC <- read.delim("./ChIPseq/QC/summary.txt", as.is = T) %>%
	mutate(Sample = gsub("PC.", "", Sample), Mark = str_split_fixed(Sample, "_", 3)[, 2], Sample = gsub("_.*", "", Sample)) 
ChIPseq_QC_target_region <- read.delim("./ChIPseq/QC/QC.target_regions.txt", as.is = T)
ChIPseq_QC <- merge(ChIPseq_QC, ChIPseq_QC_target_region) %>% mutate(Percent_Domain_Reads = Number_Target_Regions*100/Number_Reads_Aligned) 
ChIPseq_QC_long <- ChIPseq_QC %>% 
	select(Sample, Mark, Number_Reads_Aligned, Mapping_Efficiency, Percent_of_dups, Percent_of_Paired_Alignments, Number_Uniquely_Aligned_Reads_without_Dups_and_Q_.._10, Percent_Uniquely_Aligned_Reads_without_Dups_and_Q_.._10, Av_Insert_Size) %>% 
	melt(id = c("Sample", "Mark")) %>% 
	mutate(variable = gsub("Uniquely_Aligned_Reads_without_Dups_and_Q_.._10", "after_filter", variable), variable = gsub("Percent_of_Paired_Alignments", "Percent_of_Paired", variable), variable = factor(variable, levels = c("Number_Reads_Aligned", "Number_after_filter", "Mapping_Efficiency", "Percent_of_dups", "Percent_of_Paired", "Percent_after_filter", "Percent_Domain_Reads", "Av_Insert_Size")))
(ChIPseq_QC_figure <- ggplot(ChIPseq_QC_long, aes(Sample, value, fill = Mark)) + 
	geom_bar(stat = "identity", position = position_dodge(), width = 0.8) + 
	facet_grid(variable ~ ., scales = "free_y") + 
	xlab("") + 
	ylab("") + 
	theme_bw())
ggsave(ChIPseq_QC_figure, file = "./ChIPseq/QC/ChIPseq_QC_figure.pdf", height = 12, width = 7)
(ChIPseq_QC_Domain_figure <- ggplot(ChIPseq_QC %>% filter(!grepl("Input", Mark)), aes(Sample, Percent_Domain_Reads)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	facet_grid(Mark ~ ., scales = "free_y") + 
	xlab("") + 
	ylab("Domain reads (% of mapped)") + 
	theme_bw())
ggsave(ChIPseq_QC_Domain_figure, file = "./ChIPseq/QC/ChIPseq_QC_Domain_figure.pdf", height = 7, width = 4)

save(list = c("ChIPseq_QC", "ChIPseq_QC_figure", "ChIPseq_QC_Domain_figure"), 
		 file = "/projects/epigenomics2/Brain/NJabado/NJabado.Rdata")
