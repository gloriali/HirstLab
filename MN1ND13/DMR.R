# MN1 ND13 DP vs DN DMR analysis

library(ggplot2)
library(plyr)
library(VennDiagram)
library(gridExtra)
library(gplots)
library(reshape2)
library(wq)
library(dplyr)
library(RCircos)
library(stringr)
source('~/HirstLab/Pipeline/R/DMR.figures.R')
source('~/HirstLab/Pipeline/R/enrich.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
load("/projects/epigenomics2/PING/DMR/DMR.Rdata")
setwd("/projects/epigenomics2/PING/DMR/")

## ========== global QC ==========
DP_coverage <- read.delim("A04_MN1ND13_DP.bed.CpG.txt.gz.coverage.txt", head = F, as.is = T, col.names = c("coverage", "N")) %>% mutate(sample = "DP")
DN_coverage <- read.delim("A05_MN1ND13_DN.bed.CpG.txt.gz.coverage.txt", head = F, as.is = T, col.names = c("coverage", "N")) %>% mutate(sample = "DN")
CD34_coverage <- read.delim("CEMT_32_CD34.5mC.CpG.coverage.txt", head = F, as.is = T, col.names = c("coverage", "N")) %>% mutate(sample = "CD34")
CpG_coverage <- rbind(DP_coverage, DN_coverage, CD34_coverage)
(CpG_coverage_figure <- ggplot(data = CpG_coverage, aes(coverage, N/1e6, color = sample)) + 
	geom_point() + 
	geom_line() + 
	geom_vline(xintercept = 3) + 
	coord_cartesian(xlim = c(0, 50)) + 
	ylab("No. of million CpGs") +
	theme_bw())
ggsave(CpG_coverage_figure, file = "CpG_coverage_figure.pdf", height = 5, width = 5)
quantile_5mC <- read.delim("qc.5mC.quantile", as.is = T) %>% mutate(sample = gsub(".*_", "", sample))
(quantile_5mC_figure <- ggplot(quantile_5mC, aes(x = type, fill = sample, lower = lower, middle = median, upper = upper, ymin = ymin, ymax = ymax)) + 
	geom_boxplot(stat = "identity") + 
	xlab("") + 
	ylab("Fractional methylation") + 
	scale_fill_manual(values = c("grey", "coral", "royal blue")) + 
	theme_bw())
ggsave(quantile_5mC_figure, file = "quantile_5mC_figure.pdf", height = 4, width = 3)
CGI_5mC <- read.delim("CGI.5mC", as.is = T, col.names = c("ID", "fractional", "sample")) 
(CGI_5mC_figure <- ggplot(CGI_5mC, aes(fractional, color = sample)) + 
		geom_density(size = 1.2) + 
		xlab("Fractional methylation") + 
		ylab("") + 
		scale_color_manual(values = c("darkgrey", "coral", "royal blue")) + 
		theme_bw())
ggsave(CGI_5mC_figure, file = "CGI_5mC_figure.pdf", height = 4, width = 5)
CGI_5mC_wide <- CGI_5mC %>% dcast(ID ~ sample, value.var = "fractional")
t.test(CGI_5mC_wide$CD34, CGI_5mC_wide$DN)$p.value
t.test(CGI_5mC_wide$CD34, CGI_5mC_wide$DP)$p.value
t.test(CGI_5mC_wide$DP, CGI_5mC_wide$DN)$p.value
delta_5mC_DP_DN <- read.delim("DP_DN.5mC.delta.summary", head = F, as.is = T, col.names = c("delta", "N")) %>% filter(delta != 1) %>% mutate(samples = "DP-DN")
delta_5mC_DP_CD34 <- read.delim("DP_CD34.5mC.delta.summary", head = F, as.is = T, col.names = c("delta", "N")) %>% filter(delta != 1) %>% mutate(samples = "DP-CD34")
delta_5mC_DN_CD34 <- read.delim("DN_CD34.5mC.delta.summary", head = F, as.is = T, col.names = c("delta", "N")) %>% filter(delta != 1) %>% mutate(samples = "DN-CD34")
delta_5mC <- rbind(delta_5mC_DP_DN, delta_5mC_DP_CD34, delta_5mC_DN_CD34)
(delta_5mC_figure <- ggplot(delta_5mC, aes(delta + 0.05, N/1e6, color = samples)) + 
	geom_point() + 
	geom_line() + 
	geom_vline(xintercept = 0) + 
	xlab("Delta fractional methylation") + 
	ylab("No. of million CpGs") +
	theme_bw())
ggsave(delta_5mC_figure, file = "delta_5mC_figure.pdf", height = 5, width = 5)

## ============ DMRs ================
colname <- c("chr", "start", "end", "ID", "DM", "count", "length")
DMR_DP_DN <- read.delim("DMR.DP_DN.s500.c3", head = F, as.is = T, col.names = colname) %>% mutate(chr = paste0("chr", chr))
DMR_DP_DN_figure <- DMR_figures(DMR_DP_DN, "DP", "DN", figures = c("length", "count", "adjacentDis", "frequency", "circos"), hist_width = 20)
(DMR_DP_DN_figure$length)
(DMR_DP_DN_figure$count)
(DMR_DP_DN_figure$dis)
(DMR_DP_DN_figure$freq <- DMR_DP_DN_figure$freq + 
		scale_fill_manual(values = c("coral", "royal blue"), name = "", labels = c("hypo in DN", "hypo in DP")) + 
		coord_flip(ylim = c(-600, 600)))
(GREAT_DMR_DP_DN_hyper_figure <- enrich_GREAT("DP_DN_hyper", "DP_DN_hyper", categories = c("GOBP", "DiseaseOntology", "MousePhenotype"), height = 10, width = 8))
DMR_DP_CD34 <- read.delim("DMR.DP_CD34.s500.c3", head = F, as.is = T, col.names = colname) %>% mutate(chr = paste0("chr", chr))
DMR_DP_CD34_figure <- DMR_figures(DMR_DP_CD34, "DP", "CD34", figures = c("length", "count", "adjacentDis", "frequency", "circos"), hist_width = 20)
(DMR_DP_CD34_figure$length)
(DMR_DP_CD34_figure$count)
(DMR_DP_CD34_figure$dis)
(DMR_DP_CD34_figure$freq)
(GREAT_DMR_DP_CD34_hyper_figure <- enrich_GREAT("DP_CD34_hyper", "DP_CD34_hyper", categories = c("GOBP", "DiseaseOntology", "MousePhenotype"), height = 10, width = 8))
(GREAT_DMR_DP_CD34_hypo_figure <- enrich_GREAT("DP_CD34_hypo", "DP_CD34_hypo", categories = c("GOBP", "GOMF", "InterPro"), height = 7, width = 8))
DMR_DN_CD34 <- read.delim("DMR.DN_CD34.s500.c3", head = F, as.is = T, col.names = colname) %>% mutate(chr = paste0("chr", chr))
DMR_DN_CD34_figure <- DMR_figures(DMR_DN_CD34, "DN", "CD34", figures = c("length", "count", "adjacentDis", "frequency", "circos"), hist_width = 20)
(DMR_DN_CD34_figure$length)
(DMR_DN_CD34_figure$count)
(DMR_DN_CD34_figure$dis)
(DMR_DN_CD34_figure$freq)
(GREAT_DMR_DN_CD34_hyper_figure <- enrich_GREAT("DN_CD34_hyper", "DN_CD34_hyper", categories = c("GOBP", "DiseaseOntology", "MousePhenotype"), height = 10, width = 8))
(GREAT_DMR_DN_CD34_hypo_figure <- enrich_GREAT("DN_CD34_hypo", "DN_CD34_hypo", categories = c("GOMF", "DiseaseOntology", "MousePhenotype"), height = 10, width = 8))
genomic_breakdown <- read.delim("./CpG/genomic.breakdown.summary", as.is = T) %>% 
	mutate(DM = gsub(".*\\.", "", Name), Sample = gsub("\\.s500.*", "", Name), NCpG = NULL, Name = NULL)
genomic_breakdown_tall <- melt(genomic_breakdown, id = c("DM", "Sample")) 
(genomic_breakdown_figure <- ggplot(genomic_breakdown_tall, aes(variable, log2(value), fill = DM)) + 
	geom_bar(position = position_dodge(), stat = "identity", width = 0.5) + 	
	facet_wrap(~ Sample) + 
	xlab("") + 
	ylab("log2 fold enrichment") + 
	scale_fill_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw())
ggsave(genomic_breakdown_figure, file = "genomic_breakdown_figure.pdf", height = 5, width = 8)

save(list = c(ls(pattern = ".*figure")),
		 file = "/projects/epigenomics2/PING/DMR/DMR.Rdata")
