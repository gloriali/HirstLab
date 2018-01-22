# NHA VitC - hMeDIP analysis

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
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/")
load("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/hMeDIP")
RPKM <- read.delim("../RNAseq/RPKM/vitc.RPKM", as.is = T)
rownames(RPKM) <- RPKM$ENSG

## ============= QC =============
QC_summary <- read.delim("./bam/summary.txt", as.is = T) %>% filter(grepl("trim", Library)) %>% mutate(Sample = gsub(".trim", "", Library))
QC_summary_reads <- QC_summary %>% select(Sample, Total_Number_Of_Reads, Number_Reads_Aligned, Number_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = "Sample")
(QC_summary_reads_figure <- ggplot(QC_summary_reads, aes(Sample, value/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ ., scales = "free_x") + 
		xlab("") + 
		ylab("No. of million reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_reads_figure, file = "./bam/QC_summary_reads_figure.pdf", height = 8, width = 5)
QC_summary_percent <- QC_summary %>% select(Sample, Mapping_Efficiency, Percent_of_dups, Percent_of_Paired_Alignments, Percent_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = "Sample")
(QC_summary_percent_figure <- ggplot(QC_summary_percent, aes(Sample, value, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ ., scales = "free") + 
		xlab("") + 
		ylab("Percent of total number of reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_percent_figure, file = "./bam/QC_summary_percent_figure.pdf", height = 10, width = 5)
### correlation between replicates
cov_MGG_control <- read.delim("./wig/hg19.chrlen.autoXY.1KB.bed.MGG_control.coverage", head = F, as.is = T, col.names = c("ID", "MGG_control1", "MGG_control2"))
(cov_MGG_control_figure <- ggplot(cov_MGG_control, aes(log10(MGG_control1), log10(MGG_control2))) + 
		stat_density_2d(aes(fill = ..level..), geom = "polygon") + 
		theme_bw())
ggsave(cov_MGG_control_figure, file = "./wig/cov_MGG_control_figure.pdf", height = 5, width = 5)
cov_MGG_vitc <- read.delim("./wig/hg19.chrlen.autoXY.1KB.bed.MGG_vitc.coverage", head = F, as.is = T, col.names = c("ID", "MGG_vitc1", "MGG_vitc2"))
(cov_MGG_vitc_figure <- ggplot(cov_MGG_vitc, aes(log10(MGG_vitc1), log10(MGG_vitc2))) + 
		stat_density_2d(aes(fill = ..level..), geom = "polygon") + 
		theme_bw())
ggsave(cov_MGG_vitc_figure, file = "./wig/cov_MGG_vitc_figure.pdf", height = 5, width = 5)
cov_MGG_control_vitc <- read.delim("./wig/hg19.chrlen.autoXY.1KB.bed.MGG_control1_vitc1.coverage", head = F, as.is = T, col.names = c("ID", "MGG_control1", "MGG_vitc1"))
(cov_MGG_control_vitc_figure <- ggplot(cov_MGG_control_vitc, aes(log10(MGG_control1), log10(MGG_vitc1))) + 
		stat_density_2d(aes(fill = ..level..), geom = "polygon") + 
		theme_bw())
ggsave(cov_MGG_control_vitc_figure, file = "./wig/cov_MGG_control_vitc_figure.pdf", height = 5, width = 5)

## ============= MACS2 ==========
ER_summary <- read.delim("./MACS2/ER_summary.txt", as.is = T) %>% mutate(Sample = gsub(".trim", "", Sample))
(ER_summary_figure <- ggplot(ER_summary, aes(Sample, Total_length/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		xlab("") + 
		ylab("Total length (MB)") + 
		ggtitle("MACS2 ER") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(ER_summary_figure, file = "./MACS2/ER_summary_figure.pdf", height = 6, width = 5)

### unique ERs
ER_unique_summary <- read.delim("./unique/ER_unique_summary.txt", as.is = T) 
#### enrichment in genomic features
ER_unique_genomic_breakdown <- read.delim("./unique/intersect/genomic.breakdown.summary", as.is = T) %>% mutate(Comparison = gsub("\\..*", "", Name), Name = str_replace(Name, "[A-Za-z_]+\\.", ""), NCpG = NULL)
ER_unique_genomic_breakdown_tall <- melt(ER_unique_genomic_breakdown, id = c("Comparison", "Name")) 
(ER_unique_genomic_breakdown_figure <- ggplot(ER_unique_genomic_breakdown_tall, aes(variable, value, fill = Name)) + 
		geom_bar(position = position_dodge(), stat = "identity", width = 0.5) + 
		facet_wrap(~ Comparison, scales = "free_x") + 
		xlab("") + 
		ylab("Fold enrichment") + 
		coord_flip() + 
		theme_bw())
ggsave(ER_unique_genomic_breakdown_figure, file = "./unique/ER_unique_genomic_breakdown_figure.pdf", height = 6, width = 10)
#### vitc reversed
NHAR_wt_unique <- as.numeric(system("less ./unique/NHAR_control_NHA_control.NHA_control.unique.bed | awk '{s=s+$3-$2}END{print s}'", intern = T))
NHAR_vitc_unique <- as.numeric(system("less ./unique/NHAR_vitc_NHAR_control.NHAR_vitc.unique.bed | awk '{s=s+$3-$2}END{print s}'", intern = T))
NHAR_wt_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique/NHAR_control_NHA_control.NHA_control.unique.bed -b ./unique/NHAR_vitc_NHAR_control.NHAR_vitc.unique.bed | awk '{s=s+$3-$2}END{print s}'", intern = T))
NHAR_wt_vitc_venn <- draw.pairwise.venn(area1 = NHAR_wt_unique, area2 = NHAR_vitc_unique, cross.area = NHAR_wt_vitc, category = c("NHAR mut loss", "NHAR vitc gain"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique/NHAR_wt_vitc_venn.pdf", height = 3, width = 4)
grid.draw(NHAR_wt_vitc_venn)
dev.off()
MGG_vitc_unique <- as.numeric(system("less ./unique/MGG_vitc_MGG_control.MGG_vitc.unique.bed | awk '{s=s+$3-$2}END{print s}'", intern = T))
NHAR_vitc_unique <- as.numeric(system("less ./unique/NHAR_vitc_NHAR_control.NHAR_vitc.unique.bed | awk '{s=s+$3-$2}END{print s}'", intern = T))
NHAR_MGG_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique/MGG_vitc_MGG_control.MGG_vitc.unique.bed -b ./unique/NHAR_vitc_NHAR_control.NHAR_vitc.unique.bed | awk '{s=s+$3-$2}END{print s}'", intern = T))
NHAR_MGG_vitc_venn <- draw.pairwise.venn(area1 = MGG_vitc_unique, area2 = NHAR_vitc_unique, cross.area = NHAR_MGG_vitc, category = c("MGG vitc gain", "NHAR vitc gain"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique/NHAR_MGG_vitc_venn.pdf", height = 3, width = 4)
grid.draw(NHAR_MGG_vitc_venn)
dev.off()
#### homer
homer_unique_MGG119vitc_MGG119control <- read.delim("./unique/homer/MGG_vitc_MGG_control.MGG_vitc/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) 
(homer_unique_MGG119vitc_MGG119control_figure <- ggplot(homer_unique_MGG119vitc_MGG119control, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("MGG119vitc_MGG119control MGG119vitc unique") + 
		theme_bw())
ggsave(homer_unique_MGG119vitc_MGG119control_figure, file = "./unique/homer/homer_unique_MGG119vitc_MGG119control_figure.pdf", height = 5, width = 5)
homer_unique_NHARvitc_NHARcontrol <- read.delim("./unique/homer/NHAR_vitc_NHAR_control.NHAR_vitc/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) 
(homer_unique_NHARvitc_NHARcontrol_figure <- ggplot(homer_unique_NHARvitc_NHARcontrol, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("NHARvitc_NHARcontrol NHARvitc unique") + 
		theme_bw())
ggsave(homer_unique_NHARvitc_NHARcontrol_figure, file = "./unique/homer/homer_unique_NHARvitc_NHARcontrol_figure.pdf", height = 5, width = 5)
homer_unique_NHAcontrol_NHARcontrol <- read.delim("./unique/homer/NHAR_control_NHA_control.NHA_control/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) 
(homer_unique_NHAcontrol_NHARcontrol_figure <- ggplot(homer_unique_NHAcontrol_NHARcontrol, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("NHARvitc_NHARcontrol NHARvitc unique") + 
		theme_bw())
ggsave(homer_unique_NHAcontrol_NHARcontrol_figure, file = "./unique/homer/homer_unique_NHAcontrol_NHARcontrol_figure.pdf", height = 5, width = 5)


## ============= save ===========
save(list = c(ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "enrich"), ls(pattern = "DAVID"), ls(pattern = "GREAT"), ls(pattern = "venn")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/hMeDIP.Rdata")
