# NHA VitC - PBAL analysis

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
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/")
load("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/PBAL.Rdata")
RPKM <- read.delim("../RNAseq/RPKM/vitc.RPKM", as.is = T)
rownames(RPKM) <- RPKM$ENSG

## ------- QC -------
QC_summary <- read.delim("./bam/summary.xls", as.is = T, row.names = 2, head = F) %>% select(-V1) %>% t() %>% as.data.frame() %>% mutate(Sample = Library)
QC_summary_reads <- QC_summary %>% select(Sample, Total_Number_Of_Reads, Number_Reads_Aligned, Number_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = "Sample")
(QC_summary_reads_figure <- ggplot(QC_summary_reads, aes(Sample, as.numeric(value)/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ ., scales = "free_x") + 
		xlab("") + 
		ylab("No. of million reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_reads_figure, file = "./bam/QC_summary_reads_figure.pdf", height = 8, width = 5)
QC_summary_percent <- QC_summary %>% select(Sample, Mapping_Efficiency, Percent_of_dups, Percent_of_Paired_Alignments, Percent_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = "Sample")
(QC_summary_percent_figure <- ggplot(QC_summary_percent, aes(Sample, as.numeric(value), fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ ., scales = "free") + 
		xlab("") + 
		ylab("Percent of total number of reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_percent_figure, file = "./bam/QC_summary_percent_figure.pdf", height = 10, width = 5)
qc_5mC_coverage1 <- read.delim("/projects/epigenomics2/users/lli/glioma/WGBS/qc_5mC_coverage.txt", as.is = T) %>% mutate(category = gsub("_.*", "", sample), batch = gsub("[0-9]+", "", gsub("_.*", "", gsub("IDH.*t_", "", sample))))
qc_5mC_coverage <- read.delim("qc_5mC_coverage.txt", as.is = T) %>% mutate(category = sample, batch = "PBAL") %>% rbind(qc_5mC_coverage1)
(qc_5mC_coverage_figure <- ggplot(qc_5mC_coverage, aes(coverage, N/1e6, color = batch, group = sample)) + 
		geom_line() + 
		guides(color = guide_legend(title = NULL)) + 
		ylab("No. of million CpGs") + 
		coord_cartesian(xlim = c(0, 50)) + 
		theme_bw())
ggsave(qc_5mC_coverage_figure, file = "qc_5mC_coverage_figure.pdf", height = 5, width = 5)

## ------- clustering ---------
genome_5mC <- read.table("matrix_genome.5mC", sep = " ", row.names = 1, head = T)
spearman_genome_5mC <- cor(genome_5mC, method = "spearman")
pearson_genome_5mC <- cor(genome_5mC, method = "pearson")
write.table(spearman_genome_5mC, file = "cluster_spearman_genome_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(pearson_genome_5mC, file = "cluster_pearson_genome_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
CGI_5mC <- read.table("matrix_CGI.5mC", sep = " ", row.names = 1, head = T)
spearman_CGI_5mC <- cor(CGI_5mC, method = "spearman")
pearson_CGI_5mC <- cor(CGI_5mC, method = "pearson")
write.table(spearman_CGI_5mC, file = "cluster_spearman_CGI_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(pearson_CGI_5mC, file = "cluster_pearson_CGI_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
op <- par(mar = c(5, 4, 4, 8))
spearman_genome_5mC <- read.delim("cluster_spearman_genome_5mC.cor", as.is = T)
spearman_genome_5mC_dend <-  hclust(1 - spearman_genome_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(spearman_genome_5mC_dend, main = "genome-wide spearman", horiz = TRUE)
pearson_genome_5mC <- read.delim("cluster_pearson_genome_5mC.cor", as.is = T)
pearson_genome_5mC_dend <- hclust(1 - pearson_genome_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(pearson_genome_5mC_dend, main = "genom-wide pearson", horiz = TRUE)
spearman_CGI_5mC <- read.delim("cluster_spearman_CGI_5mC.cor", as.is = T)
spearman_CGI_5mC_dend <- hclust(1 - spearman_CGI_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(spearman_CGI_5mC_dend, main = "CGI spearman", horiz = TRUE)
pearson_CGI_5mC <- read.delim("cluster_pearson_CGI_5mC.cor", as.is = T)
pearson_CGI_5mC_dend <-  hclust(1 - pearson_CGI_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(pearson_CGI_5mC_dend, main = "CGI pearson", horiz = TRUE)
par(op)
pdf("cluster.pdf", height = 5, width = 5)
op <- par(mar = c(5, 4, 4, 8))
plot(spearman_genome_5mC_dend, main = "genome-wide spearman", horiz = TRUE)
plot(pearson_genome_5mC_dend, main = "genome-wide pearson", horiz = TRUE)
plot(spearman_CGI_5mC_dend, main = "CGI spearman", horiz = TRUE)
plot(pearson_CGI_5mC_dend, main = "CGI pearson", horiz = TRUE)
par(op)
dev.off()
plot(ecdf(genome_5mC$MGG_vitc))
lines(ecdf(genome_5mC$NHAR_control), col = "red")
lines(ecdf(genome_5mC$NHAR_vitc), col = "blue")
plot(density(genome_5mC$MGG_vitc - genome_5mC$MGG_control), xlim = c(-.3, .3))
lines(density(genome_5mC$NHAR_control - genome_5mC$NHA_control), col = "red")
lines(density(genome_5mC$NHAR_vitc - genome_5mC$NHAR_control), col = "blue")

## ------- 5mC distribution --------
qc_5mC_violin <- rbind(melt(genome_5mC) %>% mutate(type = "genome"), melt(CGI_5mC) %>% mutate(type = "CGI"))
(qc_5mC_violin_figure <- ggplot(qc_5mC_violin, aes(variable, value, color = variable)) + 
		geom_violin() + 
		facet_grid(type ~ ., scales = "free_y") + 
		guides(color = guide_legend(title = NULL)) + 
		xlab("") + 
		ylab("Fractional methylation") + 
		coord_flip() + 
		theme_bw())
ggsave(qc_5mC_violin_figure, file = "qc_5mC_violin_figure.pdf", height = 7, width = 6)
qc_5mC_profile <- read.delim("qc_5mC_profile.txt", as.is = T) %>% mutate(category = gsub("_.*", "", sample), N = ifelse(type == "genome", N/1e6, N/1e3)) 
(qc_5mC_profile_figure <- ggplot(qc_5mC_profile, aes(fractional, N, color = sample)) + 
		geom_smooth(se = F, span = 0.1, size = 0.5) + 
		facet_grid(type ~ ., scales = "free_y") + 
		guides(color = guide_legend(title = NULL), override.aes = list(size = 5)) + 
		coord_cartesian(ylim = c(0, 6)) + 
		xlab("Fractional methylation") + 
		ylab("No. of million CpGs                    No. of thousand CGIs") + 
		theme_bw())
ggsave(qc_5mC_profile_figure, file = "qc_5mC_profile_figure.pdf", height = 5, width = 6)
qc_5mC_quantile <- read.delim("qc_5mC_quantile.txt", as.is = T) %>% mutate(category = gsub("_.*", "", sample)) %>% filter(sample != "IDHmut_CEMT_21")
(qc_5mC_quantile_figure <- ggplot(qc_5mC_quantile, aes(x = sample, lower = lower, middle = median, upper = upper, ymin = ymin, ymax = ymax, fill = category)) + 
		geom_boxplot(stat = "identity", color = "grey") + 
		facet_grid(type ~ .) + 
		guides(color = guide_legend(title = NULL)) + 
		xlab("") + 
		ylab("Fractional methylation") + 
		coord_flip() + 
		theme_bw())
ggsave(qc_5mC_quantile_figure, file = "qc_5mC_quantile_figure.pdf", height = 7, width = 6)

## ----------- DMRs -------------
### -------- enrichment in genomic regions ------
genomic_breakdown <- read.delim("./DMR/intersect/genomic.breakdown.summary", as.is = T) %>% 
	mutate(sample = gsub("\\..*", "", gsub("DhMR.", "", Name)), DM = gsub(".*\\.", "", Name), NCpG = NULL, Name = NULL)
genomic_breakdown_tall <- melt(genomic_breakdown, id = c("sample", "DM")) 
(genomic_breakdown_figure <- ggplot(genomic_breakdown_tall, aes(variable, log2(value), fill = DM)) + 
		geom_bar(position = position_dodge(), stat = "identity", width = 0.5) + 
		geom_hline(yintercept = c(-1, 1)) + 
		facet_wrap(~sample) + 
		xlab("") + 
		ylab("log2 Fold enrichment") + 
		scale_fill_manual(name = "", values = c("red", "blue")) + 
		coord_flip() + 
		theme_bw())
ggsave(genomic_breakdown_figure, file = "./DMR/genomic_breakdown_figure.pdf", height = 6, width = 7)


## ============= save ===========
save(list = c(ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "dend"), ls(pattern = "enrich"), ls(pattern = "DAVID"), ls(pattern = "GREAT"), ls(pattern = "venn")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/PBAL.Rdata")
