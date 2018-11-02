# MGG VitC 

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
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/MGG/")
load("/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MGG.Rdata")

## ========== ChIPseq ==========
ER_nearest_dis <- read.delim("./ChIPseq/FindER2/ER.nearest.dis.bed", as.is = T)
(ER_nearest_dis_figure <- ggplot(ER_nearest_dis, aes(dis, color = sample)) + 
		geom_density(adjust = 0.1) + 
		facet_grid(category ~ mark, scales = "free_y") + 
		coord_cartesian(xlim = c(0, 50000)) + 
		theme_bw())
ggsave(ER_nearest_dis_figure, file = "./ChIPseq/FindER2/ER_nearest_dis_figure.pdf", height = 4, width = 6)
(ER_length_figure <- ggplot(ER_nearest_dis, aes(len, color = sample)) + 
		geom_density(adjust = 0.001) + 
		facet_grid(category ~ mark, scales = "free_y") + 
		coord_cartesian(xlim = c(0, 3000)) + 
		theme_bw())
ggsave(ER_length_figure, file = "./ChIPseq/FindER2/ER_length_figure.pdf", height = 4, width = 6)

## ============= WGBS ===========
### enhancer 5mC
enhancer_5mC <- read.delim("./WGBS/enhancer/promoter.enhancer.5mC", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample")) %>%
	mutate(group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median"))) 
enhancer_5mC_summary <- merge(enhancer_5mC %>% group_by(sample, group) %>% summarize(N=n()), enhancer_5mC %>% group_by(sample) %>% summarize(Total=n())) %>% 
	mutate(percent = N/Total, signal = 2.2, fractional = ifelse(group == "hypo", 0.15, ifelse(group == "hyper", 0.85, 0.5)), labels = paste0(round(percent*100, 0), "%"))
(enhancer_5mC_figure <- ggplot(enhancer_5mC, aes(fractional, log10(signal))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.01, 0.02, 0.04, 0.08, 0.3, 0.5, 1), guide = "none") +		
		geom_vline(xintercept = 0.3, color = "white") + 
		geom_vline(xintercept = 0.7, color = "white") + 
		geom_text(data = enhancer_5mC_summary, aes(fractional, signal, label = labels), color = "white", size = 3) + 
		facet_wrap(~ sample) + 
		coord_cartesian(ylim = c(0.3, 2.3)) + 
		xlab("Fractional methylation") + 
		ylab("log10 H3K27ac normalized signal") + 
		theme_bw())
ggsave(enhancer_5mC_figure, file = "./WGBS/enhancer_5mC_figure.pdf", height = 4, width = 6)
### enrichment in genomic features
DMR_genomic_breakdown <- read.delim("./WGBS/DMR/intersect/genomic.breakdown.summary", as.is = T) %>% mutate(Name = gsub(".*\\.", "", Name), NCpG = NULL)
DMR_genomic_breakdown_tall <- melt(DMR_genomic_breakdown, id = c("Name")) 
(DMR_genomic_breakdown_figure <- ggplot(DMR_genomic_breakdown_tall, aes(variable, log2(value), fill = Name)) + 
		geom_bar(position = position_dodge(), stat = "identity", width = 0.5) + 
		xlab("") + 
		ylab("log2 Fold enrichment") + 
		coord_flip() + 
		theme_bw())
ggsave(DMR_genomic_breakdown_figure, file = "./WGBS/DMR/DMR_genomic_breakdown_figure.pdf", height = 4, width = 5)


## ============= hMeDIP ==========
### enrichment in genomic features
ER_unique_genomic_breakdown <- read.delim("./hMeDIP/FindER2/intersect/genomic.breakdown.summary", as.is = T) %>% mutate(NCpG = NULL)
ER_unique_genomic_breakdown_tall <- melt(ER_unique_genomic_breakdown, id = c("Name")) 
(ER_unique_genomic_breakdown_figure <- ggplot(ER_unique_genomic_breakdown_tall, aes(variable, log2(value), fill = Name)) + 
		geom_bar(position = position_dodge(), stat = "identity", width = 0.5) + 
		xlab("") + 
		ylab("log2 Fold enrichment") + 
		coord_flip() + 
		theme_bw())
ggsave(ER_unique_genomic_breakdown_figure, file = "./hMeDIP/FindER2/ER_unique_genomic_breakdown_figure.pdf", height = 4, width = 5)
### homer
homer_unique_MGGvitc <- read.delim("./hMeDIP/FindER2/homer/MGG_vitc.MGG_control/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 20, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) 
(homer_unique_MGGvitc_figure <- ggplot(homer_unique_MGGvitc, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of region with motif") + 
		ggtitle("MGGvitc unique") + 
		theme_bw())
ggsave(homer_unique_MGGvitc_figure, file = "./hMeDIP/FindER2/homer_unique_MGGvitc_figure.pdf", height = 5, width = 5)
homer_unique_MGGcontrol <- read.delim("./hMeDIP/FindER2/homer/MGG_control.MGG_vitc/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 20, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) 
(homer_unique_MGGcontrol_figure <- ggplot(homer_unique_MGGcontrol, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of region with motif") + 
		ggtitle("MGGcontrol unique") + 
		theme_bw())
ggsave(homer_unique_MGGcontrol_figure, file = "./hMeDIP/FindER2/homer_unique_MGGcontrol_figure.pdf", height = 5, width = 5)
homer_unique_enhaner_MGGvitc <- read.delim("./hMeDIP/FindER2/enhancer/homer/MGG_vitc.MGG_control/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 20, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) 
(homer_unique_enhaner_MGGvitc_figure <- ggplot(homer_unique_enhaner_MGGvitc, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of region with motif") + 
		ggtitle("MGGvitc unique enhancer") + 
		theme_bw())
ggsave(homer_unique_enhaner_MGGvitc_figure, file = "./hMeDIP/FindER2/homer_unique_enhaner_MGGvitc_figure.pdf", height = 5, width = 5)
homer_unique_enhaner_MGGcontrol <- read.delim("./hMeDIP/FindER2/enhancer/homer/MGG_control.MGG_vitc/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 20, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) 
(homer_unique_enhaner_MGGcontrol_figure <- ggplot(homer_unique_enhaner_MGGcontrol, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of region with motif") + 
		ggtitle("MGGcontrol unique enhacer") + 
		theme_bw())
ggsave(homer_unique_enhaner_MGGcontrol_figure, file = "./hMeDIP/FindER2/homer_unique_enhaner_MGGcontrol_figure.pdf", height = 3, width = 5)


## ============= save ===========
save(list = c(ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "enrich"), ls(pattern = "DAVID"), ls(pattern = "GREAT"), ls(pattern = "venn")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MGG.Rdata")
