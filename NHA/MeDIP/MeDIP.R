# MGG - MeDIP analysis

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
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/")
load("/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/MeDIP.Rdata")

## ============= QC =============
QC_summary <- read.delim("./bam/summary.xls", as.is = T, row.names = 2, head = F) %>% select(-V1) %>% t() %>% as.data.frame()
colnames(QC_summary) <- gsub("without_Dups_and_Q_>=_10", "After_Filter", colnames(QC_summary))
QC_summary_reads <- QC_summary %>% select(Library, Total_Number_Of_Reads, Number_Reads_Aligned, Number_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = "Library")
(QC_summary_reads_figure <- ggplot(QC_summary_reads, aes(Library, as.numeric(value)/10^6, fill = Library)) + 
		geom_bar(stat = "identity", position = position_dodge()) + 
		facet_grid(variable ~ ., scales = "free_x") + 
		xlab("") + 
		ylab("No. of million reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_reads_figure, file = "./bam/QC_summary_reads_figure.pdf", height = 8, width = 5)
QC_summary_percent <- QC_summary %>% select(Library, Mapping_Efficiency, Percent_of_dups, Percent_of_Paired_Alignments, Percent_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = "Library")
(QC_summary_percent_figure <- ggplot(QC_summary_percent, aes(Library, as.numeric(value), fill = Library)) + 
		geom_bar(stat = "identity", position = position_dodge()) + 
		facet_grid(variable ~ ., scales = "free") + 
		xlab("") + 
		ylab("Percent of total number of reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_percent_figure, file = "./bam/QC_summary_percent_figure.pdf", height = 10, width = 5)
### fractional calls: signal/background coverage Ecdf
names <- c('MGG_control.24h','MGG_vitc.24h','MGG_vitc.48h','MGG_vitc.72h','MGG_vitc.6d')
chrs <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
for(i in 1:5){
	name <- names[i]
	sig_total <- data.frame(cov = 0:1000, total = 0)
	back_total <- data.frame(cov = 0:1000, total = 0)
	for(j in 1:24){
		chr <- chrs[j]
		print(paste(name, chr))
		sig <- read.delim(paste0("./fractional/CG_25_around_chr/", name, "/", chr, "/", chr, ".gz.", name, ".covDist"), head = F, col.names = c("cov", "n", "total")) 
		sig <- sig %>% mutate(cdf = total/sig[nrow(sig), "total"])
		back <- read.delim(paste0("./fractional/CG_empty_500_chr/", name, "/", chr, "/", chr, ".gz.", name, ".covDist"), head = F, col.names = c("cov", "n", "total")) 
		back <- back %>% mutate(cdf = total/back[nrow(back), "total"])
		pdf(paste0("./fractional/CDF_cov_plots/CDF_cov_", name, ".", chr, ".pdf"))
		plot(x = c(0, 50), y = c(0, 1), main = paste(name, chr), type = "n", xlab = "Average coverage", ylab = "Ecdf")
		lines(sig$cov, sig$cdf, col = "red")
		lines(back$cov, back$cdf, col = "black")
		abline(h = 1)
		legend("bottomright", c("signal", "background"), col = c("red", "black"), lwd = 3)
		dev.off()
		sig_total$total <- sig_total$total + sig$total
		back_total$total <- back_total$total + back$total
	}
	sig_total <- sig_total %>% mutate(cdf = total/sig_total[nrow(sig_total), "total"])
	back_total <- back_total %>% mutate(cdf = total/back_total[nrow(back_total), "total"])
	pdf(paste0("./fractional/CDF_cov_plots/CDF_cov_", name, ".pdf"))
	plot(x = c(0, 50), y = c(0, 1), main = name, type = "n", xlab = "Average coverage", ylab = "Ecdf")
	lines(sig_total$cov, sig_total$cdf, col = "red")
	lines(back_total$cov, back_total$cdf, col = "black")
	abline(h = 1)
	legend("bottomright", c("signal", "background"), col = c("red", "black"), lwd = 3)
	dev.off()
}

## ------- 5mC distribution --------
qc_5mC_profile <- read.delim("./fractional/qc_5mC_profile.txt", as.is = T) %>% mutate(N = ifelse(type == "genome", N/1e6, N/1e3)) 
(qc_5mC_profile_figure <- ggplot(qc_5mC_profile, aes(fractional, N, color = sample)) + 
		geom_smooth(se = F, span = 0.1, size = 0.5) + 
		facet_grid(type ~ ., scales = "free_y") + 
		guides(color = guide_legend(title = NULL), override.aes = list(size = 5)) + 
		coord_cartesian(ylim = c(0, 8)) + 
		xlab("Fractional methylation") + 
		ylab("No. of million CpGs                    No. of thousand CGIs") + 
		theme_bw())
ggsave(qc_5mC_profile_figure, file = "./fractional/qc_5mC_profile_figure.pdf", height = 5, width = 6)
qc_5mC_quantile <- read.delim("./fractional/qc_5mC_quantile.txt", as.is = T) 
(qc_5mC_quantile_figure <- ggplot(qc_5mC_quantile, aes(x = sample, lower = lower, middle = median, upper = upper, ymin = ymin, ymax = ymax, fill = sample)) + 
		geom_boxplot(stat = "identity", color = "grey") + 
		facet_grid(type ~ .) + 
		guides(color = guide_legend(title = NULL)) + 
		xlab("") + 
		ylab("Fractional methylation") + 
		coord_flip() + 
		theme_bw())
ggsave(qc_5mC_quantile_figure, file = "./fractional/qc_5mC_quantile_figure.pdf", height = 5, width = 8)

## ------- clustering ---------
genome_5mC <- read.table("./fractional/matrix_genome.5mC", sep = " ", head = T, row.names = 1) 
spearman_genome_5mC <- cor(genome_5mC, method = "spearman")
pearson_genome_5mC <- cor(genome_5mC, method = "pearson")
write.table(spearman_genome_5mC, file = "./fractional/cluster_spearman_genome_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(pearson_genome_5mC, file = "./fractional/cluster_pearson_genome_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
CGI_5mC <- read.table("./fractional/matrix_CGI.5mC", sep = " ", head = T, row.names = 1) 
spearman_CGI_5mC <- cor(CGI_5mC, method = "spearman")
pearson_CGI_5mC <- cor(CGI_5mC, method = "pearson")
write.table(spearman_CGI_5mC, file = "./fractional/cluster_spearman_CGI_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(pearson_CGI_5mC, file = "./fractional/cluster_pearson_CGI_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
op <- par(mar = c(5, 4, 4, 10))
spearman_genome_5mC <- read.delim("./fractional/cluster_spearman_genome_5mC.cor", as.is = T)
spearman_genome_5mC_dend <-  hclust(1 - spearman_genome_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(spearman_genome_5mC_dend, main = "genome-wide spearman", horiz = TRUE)
pearson_genome_5mC <- read.delim("./fractional/cluster_pearson_genome_5mC.cor", as.is = T)
pearson_genome_5mC_dend <- hclust(1 - pearson_genome_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(pearson_genome_5mC_dend, main = "genom-wide pearson", horiz = TRUE)
spearman_CGI_5mC <- read.delim("./fractional/cluster_spearman_CGI_5mC.cor", as.is = T)
spearman_CGI_5mC_dend <- hclust(1 - spearman_CGI_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(spearman_CGI_5mC_dend, main = "CGI spearman", horiz = TRUE)
pearson_CGI_5mC <- read.delim("./fractional/cluster_pearson_CGI_5mC.cor", as.is = T)
pearson_CGI_5mC_dend <-  hclust(1 - pearson_CGI_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(pearson_CGI_5mC_dend, main = "CGI pearson", horiz = TRUE)
par(op)
pdf("./fractional/cluster.pdf", height = 5, width = 8)
op <- par(mar = c(5, 4, 4, 10))
plot(spearman_genome_5mC_dend, main = "genome-wide spearman", horiz = TRUE)
plot(pearson_genome_5mC_dend, main = "genome-wide pearson", horiz = TRUE)
plot(spearman_CGI_5mC_dend, main = "CGI spearman", horiz = TRUE)
plot(pearson_CGI_5mC_dend, main = "CGI pearson", horiz = TRUE)
par(op)
dev.off()

