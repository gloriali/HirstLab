# Glioma - WGBS analysis

library(ggplot2)
library(plyr)
library(VennDiagram)
library(gridExtra)
library(gplots)
library(dendextend)
library(reshape2)
library(wq)
library(dplyr)
library(RCircos)
library(stringr)
source('~/HirstLab/Pipeline/R/DMR.figures.R')
source('~/HirstLab/Pipeline/R/enrich.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
setwd("/projects/epigenomics2/users/lli/glioma/WGBS/")
load("/projects/epigenomics2/users/lli/glioma/WGBS/WGBS.Rdata")
libs <- c("IDHmut_CEMT_19", "IDHmut_CEMT_21", "IDHmut_CEMT_22", "IDHmut_CEMT_47", "IDHmut_TCGA060128", "IDHmut_TCGA161460", "IDHmut_TCGA191788", "IDHwt_CEMT_23", "IDHwt_TCGA141401", "IDHwt_TCGA141454", "IDHwt_TCGA143477")

## ------- 5mC modifiers RPKM --------
DNAme_regulators_RPKM <- read.delim("DNAme_regulators.RPKM", as.is = T) %>% select(-contains("Brain")) %>% melt(id = c("ENSG", "gene")) %>% 
	mutate(type = ifelse(variable == "CEMT_23", "IDHwt", ifelse(grepl("CEMT", variable), "IDHmut", "NPC")))
(DNAme_regulators_RPKM_figure <- ggplot(DNAme_regulators_RPKM, aes(gene, value, color = type)) + 
		geom_point(position = position_jitter(width = 0.2)) + 
		coord_flip() + 
		guides(color = guide_legend(title = NULL)) + 
		xlab("") + 
		ylab("RPKM") + 
		theme_bw())
ggsave(DNAme_regulators_RPKM_figure, file = "DNAme_regulators_RPKM_figure.pdf", height = 5, width = 6)

## ------- CpG coverage ---------
qc_5mC_coverage <- read.delim("qc_5mC_coverage.txt", as.is = T) %>% mutate(category = gsub("_.*", "", sample), batch = gsub("[0-9]+", "", gsub("_.*", "", gsub("IDH.*t_", "", sample))))
(qc_5mC_coverage_figure <- ggplot(qc_5mC_coverage, aes(coverage, N/1e6, color = batch)) + 
		geom_line(aes(group = sample)) + 
		guides(color = guide_legend(title = NULL)) + 
		ylab("No. of million CpGs") + 
		coord_cartesian(xlim = c(0, 50)) + 
		theme_bw())
ggsave(qc_5mC_coverage_figure, file = "qc_5mC_coverage_figure.pdf", height = 5, width = 5)

## ------- 5mC distribution --------
qc_5mC_profile <- read.delim("qc_5mC_profile.txt", as.is = T) %>% mutate(category = gsub("_.*", "", sample), N = ifelse(type == "genome", N/1e6, N/1e3)) %>% filter(sample != "IDHmut_CEMT_21")
(qc_5mC_profile_figure <- ggplot(qc_5mC_profile, aes(fractional, N, color = category, group = sample)) + 
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

## ------- clustering ---------
setwd("/projects/epigenomics2/users/lli/glioma/WGBS/")
genome_5mC <- read.table("matrix_genome.5mC", sep = " ", head = F, row.names = 1, col.names = c("ID", "IDHmut_CEMT_19", "IDHmut_CEMT_21", "IDHmut_CEMT_22", "IDHmut_CEMT_47", "IDHmut_TCGA060128", "IDHmut_TCGA161460", "IDHmut_TCGA191788", "IDHwt_CEMT_23", "IDHwt_TCGA141401", "IDHwt_TCGA141454", "IDHwt_TCGA143477", "NPC_Cortex02", "NPC_Cortex04", "NPC_GE02", "NPC_GE04"))
spearman_genome_5mC <- cor(genome_5mC, method = "spearman")
pearson_genome_5mC <- cor(genome_5mC, method = "pearson")
write.table(spearman_genome_5mC, file = "cluster_spearman_genome_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(pearson_genome_5mC, file = "cluster_pearson_genome_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
CGI_5mC <- read.table("matrix_CGI.5mC", sep = " ", head = F, row.names = 1, col.names = c("ID", "IDHmut_CEMT_19", "IDHmut_CEMT_21", "IDHmut_CEMT_22", "IDHmut_CEMT_47", "IDHmut_TCGA060128", "IDHmut_TCGA161460", "IDHmut_TCGA191788", "IDHwt_CEMT_23", "IDHwt_TCGA141401", "IDHwt_TCGA141454", "IDHwt_TCGA143477", "NPC_Cortex02", "NPC_Cortex04", "NPC_GE02", "NPC_GE04"))
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
pdf("cluster.pdf", height = 5, width = 8)
op <- par(mar = c(5, 4, 4, 8))
plot(spearman_genome_5mC_dend, main = "genome-wide spearman", horiz = TRUE)
plot(pearson_genome_5mC_dend, main = "genome-wide pearson", horiz = TRUE)
plot(spearman_CGI_5mC_dend, main = "CGI spearman", horiz = TRUE)
plot(pearson_CGI_5mC_dend, main = "CGI pearson", horiz = TRUE)
par(op)
dev.off()

## ------- changes at CGI edges -------
CGI_edge <- read.delim("./CGI_edge/CGI.edge.profile") %>% mutate(category = gsub("_.*", "", sample), edge = revalue(edge, c("L" = "5-prime", "R" = "3-prime"))) %>% filter(sample != "IDHmut_CEMT_21")
(CGI_edge_figure <- ggplot(CGI_edge %>% filter(sample %in% libs), aes(-distance, fractional, color = category, group = sample)) + 
	geom_smooth(se = F, span = 0.2, size = 0.5) + 
	facet_wrap(~ edge, scales = "free_x") + 
	guides(color = guide_legend(title = NULL)) + 
	xlab("Distance to CGI edge (bp)") + 
	ylab("Fractional methylation") + 
	theme_bw())
ggsave(CGI_edge_figure, file = "./CGI_edge/CGI_edge_figure.pdf", width = 8, height = 5)
CGI_edge_delta <- read.delim("./CGI_edge/CGI.edge.delta.profile") %>% mutate(category = gsub("_.*", "", sample1), sample1 = gsub("(IDH.*t)_", "\\1\n", sample1), edge = revalue(edge, c("L" = "5-prime", "R" = "3-prime"))) %>% filter(!grepl("CEMT_21", sample1))
(CGI_edge_delta_figure <- ggplot(CGI_edge_delta %>% filter(sample1 %in% libs), aes(-distance, delta, color = category, linetype = sample2, group = interaction(sample1, sample2))) + 
		geom_smooth(se = F, span = 0.4) + 
		geom_hline(yintercept = 0) + 
		facet_grid(. ~ edge, scales = "free_x") + 
		scale_color_manual(values = c("IDHmut" = "red", "IDHwt" = "blue")) + 
		guides(color = guide_legend(title = NULL), linetype = guide_legend(title = NULL, color = "black")) + 
		xlab("Distance to CGI edge (bp)") + 
		ylab("Difference in fractional methylation\nglioma - NPC") + 
		theme_bw())
ggsave(CGI_edge_delta_figure, file = "./CGI_edge/CGI_edge_delta_figure.pdf", width = 8, height = 5)

## ------- 5mC at enhancer regions --------
enhancer_5mC <- read.delim("./H3K27ac/enhancer.5mC", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample")) %>% filter(sample != "CEMT_21") %>% 
	mutate(category = gsub(".*_", "", ID), type = ifelse(sample == "NPC_GE04", "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = paste0(type, "\n", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median")))
enhancer_5mC_summary <- merge(enhancer_5mC %>% group_by(sample, category, group) %>% summarize(N=n()), enhancer_5mC %>% group_by(sample, category) %>% summarize(Total=n())) %>% 
	mutate(percent = N/Total, signal = 2.1, fractional = ifelse(group == "hypo", 0.15, ifelse(group == "hyper", 0.85, 0.5)), labels = paste0(round(percent*100, 0), "%"))
(enhancer_5mC_figure <- ggplot(enhancer_5mC, aes(fractional, log10(signal))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.01, 0.02, 0.04, 0.08, 0.3, 0.5, 1), guide = "none") +		
		geom_vline(xintercept = 0.3, color = "white") + 
		geom_vline(xintercept = 0.7, color = "white") + 
		geom_text(data = enhancer_5mC_summary, aes(fractional, signal, label = labels), color = "white", size = 3) + 
		facet_grid(sample ~ category) + 
		coord_cartesian(ylim = c(0.5, 2.2)) + 
		xlab("Fractional methylation") + 
		ylab("log10 H3K27ac normalized signal") + 
		theme_bw())
ggsave(enhancer_5mC_figure, file = "./H3K27ac/enhancer_5mC_figure.pdf", height = 6, width = 6)
enhancer_5mC_RPKM <- read.delim("./H3K27ac/enhancer.5mC.RPKM", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "category", "signal", "fractional", "sample", "ENSG", "RPKM")) %>% 
	mutate(RPKM = RPKM + 1e-5, type = ifelse(sample == "NPC_GE04", "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = paste0(type, "\n", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median")), group = factor(group, levels = c("hypo", "median", "hyper")))
(enhancer_5mC_RPKM_figure <- ggplot(enhancer_5mC_RPKM, aes(group, log10(RPKM), fill = type)) + 
		geom_violin() + 
		coord_cartesian(ylim = c(-1.5, 2.5)) + 
		facet_grid(sample ~ category) + 
		xlab("") + 
		theme_bw())
ggsave(enhancer_5mC_RPKM_figure, file = "./H3K27ac/enhancer_5mC_RPKM_figure.pdf", height = 6, width = 6)
RPKM_p <- data.frame(sample = c(), category = c(), p = c())
for(lib in c("CEMT_19", "CEMT_22", "CEMT_23", "CEMT_47", "GE04")){
	for(c in c("promoter", "regular", "super")){
		RPKM_p <- rbind(RPKM_p, data.frame(sample = lib, category = c, p = t.test((enhancer_5mC_RPKM %>% filter(category == c, grepl(lib, sample), group == "hyper"))$RPKM, (enhancer_5mC_RPKM %>% filter(category == c, grepl(lib, sample), group == "hypo"))$RPKM)$p.value))
	}
}
enhancer_5mC_CpG <- read.delim("./H3K27ac/enhancer.5mC.CpG", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample", "NCpG", "CpGdensity")) %>% filter(sample != "CEMT_21") %>% 
	mutate(category = gsub(".*_", "", ID), type = ifelse(sample == "NPC_GE04", "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = paste0(type, "\n", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median")))
(enhancer_5mC_CpG_figure <- ggplot(enhancer_5mC_CpG, aes(CpGdensity, color = group)) + 
		geom_density(adjust = 2) +
		facet_grid(sample ~ category, scale = "free_y") + 
		coord_cartesian(xlim = c(0, 100)) + 
		xlab("No. of CpGs per kB") + 
		ylab("density") + 
		theme_bw())
ggsave(enhancer_5mC_CpG_figure, file = "./H3K27ac/enhancer_5mC_CpG_figure.pdf", height = 6, width = 6)
enhancer_5mC_CGI <- read.delim("./H3K27ac/enhancer.5mC.CGI", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample")) %>% filter(sample != "CEMT_21") %>% 
	mutate(category = gsub(".*_", "", ID), type = ifelse(sample == "NPC_GE04", "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = paste0(type, "\n", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median")))
enhancer_5mC_CGI_summary <- merge(enhancer_5mC_CGI %>% group_by(sample, category, group) %>% summarize(N=n()), enhancer_5mC_CGI %>% group_by(sample, category) %>% summarize(Total=n())) %>% 
	mutate(percent = N/Total, signal = 2.1, fractional = ifelse(group == "hypo", 0.15, ifelse(group == "hyper", 0.85, 0.5)), labels = paste0(round(percent*100, 0), "%"))
(enhancer_5mC_CGI_figure <- ggplot(enhancer_5mC_CGI, aes(fractional, log10(signal))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.01, 0.02, 0.04, 0.08, 0.3, 0.5, 1), guide = "none") +		
		geom_vline(xintercept = 0.3, color = "white") + 
		geom_vline(xintercept = 0.7, color = "white") + 
		geom_text(data = enhancer_5mC_CGI_summary, aes(fractional, signal, label = labels), color = "white", size = 3) + 
		facet_grid(sample ~ category) + 
		coord_cartesian(ylim = c(0.5, 2.2)) + 
		xlab("Fractional methylation") + 
		ylab("log10 H3K27ac normalized signal") + 
		ggtitle("CGI enhancers") + 
		theme_bw())
ggsave(enhancer_5mC_CGI_figure, file = "./H3K27ac/enhancer_5mC_CGI_figure.pdf", height = 6, width = 6)
enhancer_5mC_nonCGI <- read.delim("./H3K27ac/enhancer.5mC.nonCGI", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample")) %>% filter(sample != "CEMT_21") %>% 
	mutate(category = gsub(".*_", "", ID), type = ifelse(sample == "NPC_GE04", "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = paste0(type, "\n", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median")))
enhancer_5mC_nonCGI_summary <- merge(enhancer_5mC_nonCGI %>% group_by(sample, category, group) %>% summarize(N=n()), enhancer_5mC_nonCGI %>% group_by(sample, category) %>% summarize(Total=n())) %>% 
	mutate(percent = N/Total, signal = 2.1, fractional = ifelse(group == "hypo", 0.15, ifelse(group == "hyper", 0.85, 0.5)), labels = paste0(round(percent*100, 0), "%"))
(enhancer_5mC_nonCGI_figure <- ggplot(enhancer_5mC_nonCGI, aes(fractional, log10(signal))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.01, 0.02, 0.04, 0.08, 0.3, 0.5, 1), guide = "none") +		
		geom_vline(xintercept = 0.3, color = "white") + 
		geom_vline(xintercept = 0.7, color = "white") + 
		geom_text(data = enhancer_5mC_nonCGI_summary, aes(fractional, signal, label = labels), color = "white", size = 3) + 
		facet_grid(sample ~ category) + 
		coord_cartesian(ylim = c(0.5, 2.2)) + 
		xlab("Fractional methylation") + 
		ylab("log10 H3K27ac normalized signal") + 
		ggtitle("nonCGI enhancers") + 
		theme_bw())
ggsave(enhancer_5mC_nonCGI_figure, file = "./H3K27ac/enhancer_5mC_nonCGI_figure.pdf", height = 6, width = 6)
enhancer_5mC_TCGA <- read.delim("./H3K27ac/TCGA.enhancer.coverage.5mC", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample")) %>%
	mutate(category = gsub(".*_", "", ID), type = gsub("_.*", "", sample), sample = gsub("_", "\n", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median")))
enhancer_5mC_TCGA_summary <- merge(enhancer_5mC_TCGA %>% group_by(sample, category, group) %>% summarize(N=n()), enhancer_5mC_TCGA %>% group_by(sample, category) %>% summarize(Total=n())) %>% 
	mutate(percent = N/Total, signal = 2.1, fractional = ifelse(group == "hypo", 0.15, ifelse(group == "hyper", 0.85, 0.5)), labels = paste0(round(percent*100, 0), "%"))
(enhancer_5mC_TCGA_figure <- ggplot(enhancer_5mC_TCGA, aes(fractional, log10(signal))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.01, 0.02, 0.04, 0.08, 0.3, 0.5, 1), guide = "none") +		
		geom_vline(xintercept = 0.3, color = "white") + 
		geom_vline(xintercept = 0.7, color = "white") + 
		geom_text(data = enhancer_5mC_TCGA_summary, aes(fractional, signal, label = labels), color = "white", size = 3) + 
		facet_grid(sample ~ category) + 
		coord_cartesian(ylim = c(0.5, 2.2)) + 
		xlab("Fractional methylation") + 
		ylab("log10 H3K27ac normalized signal") + 
		theme_bw())
ggsave(enhancer_5mC_TCGA_figure, file = "./H3K27ac/enhancer_5mC_TCGA_figure.pdf", height = 7, width = 8)
enhancer_5mC_CpG_TCGA <- read.delim("./H3K27ac/TCGA.enhancer.coverage.5mC.CpG", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample", "NCpG", "CpGdensity")) %>% 
	mutate(category = gsub(".*_", "", ID), sample = gsub("_", "\n", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median")))
(enhancer_5mC_CpG_TCGA_figure <- ggplot(enhancer_5mC_CpG_TCGA, aes(CpGdensity, color = group)) + 
		geom_density(adjust = 2) +
		facet_grid(sample ~ category, scale = "free_y") + 
		coord_cartesian(xlim = c(0, 100)) + 
		xlab("No. of CpGs per kB") + 
		ylab("density") + 
		theme_bw())
ggsave(enhancer_5mC_CpG_TCGA_figure, file = "./H3K27ac/enhancer_5mC_CpG_TCGA_figure.pdf", height = 6, width = 6)
CEMT <- read.delim("../CEMT.txt", head = F, as.is = T) %>% mutate(V1 = paste0("CEMT_", V1))
row.names(CEMT) <- CEMT$V1
enhancer_5mC_other <- read.delim("./H3K27ac/others/promoter.enhancer.5mC", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample")) %>%
	mutate(sample = paste0(CEMT[sample, "V2"], "_", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median"))) %>% filter(!grepl("Brain", sample))
enhancer_5mC_other_summary <- merge(enhancer_5mC_other %>% group_by(sample, group) %>% summarize(N=n()), enhancer_5mC_other %>% group_by(sample) %>% summarize(Total=n())) %>% 
	mutate(percent = N/Total, signal = 2.2, fractional = ifelse(group == "hypo", 0.15, ifelse(group == "hyper", 0.85, 0.5)), labels = paste0(round(percent*100, 0), "%"))
(enhancer_5mC_other_figure <- ggplot(enhancer_5mC_other, aes(fractional, log10(signal))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.01, 0.02, 0.04, 0.08, 0.3, 0.5, 1), guide = "none") +		
		geom_vline(xintercept = 0.3, color = "white") + 
		geom_vline(xintercept = 0.7, color = "white") + 
		geom_text(data = enhancer_5mC_other_summary, aes(fractional, signal, label = labels), color = "white", size = 3) + 
		facet_wrap(~ sample) + 
		coord_cartesian(ylim = c(0.3, 2.3)) + 
		xlab("Fractional methylation") + 
		ylab("log10 H3K27ac normalized signal") + 
		theme_bw())
ggsave(enhancer_5mC_other_figure, file = "./H3K27ac/enhancer_5mC_other_figure.pdf", height = 6, width = 6)
enhancer_5mC_CpG_other <- read.delim("./H3K27ac/others/promoter.enhancer.5mC.CpG", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample", "NCpG", "CpGdensity")) %>% 
	mutate(sample = paste0(CEMT[sample, "V2"], "_", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median"))) %>% filter(!grepl("Brain", sample))
(enhancer_5mC_CpG_other_figure <- ggplot(enhancer_5mC_CpG_other, aes(CpGdensity, color = group)) + 
		geom_density(adjust = 2) +
		facet_wrap(~ sample) + 
		coord_cartesian(xlim = c(0, 100)) + 
		xlab("No. of CpGs per kB") + 
		ylab("density") + 
		theme_bw())
ggsave(enhancer_5mC_CpG_other_figure, file = "./H3K27ac/enhancer_5mC_CpG_other_figure.pdf", height = 6, width = 6)
for(g in c("hypo", "median", "hyper")){
	for(c in c("promoter", "regular", "super")){
		assign(paste0("IDHmut_gene_", g, "_", c), list(CEMT_19 = (enhancer_5mC_RPKM %>% filter(as.character(group) == g, as.character(category) == c, grepl("CEMT_19", sample)))[, "ENSG"], 
																						CEMT_22 = (enhancer_5mC_RPKM %>% filter(as.character(group) == g, as.character(category) == c, grepl("CEMT_22", sample)))[, "ENSG"], 
																						CEMT_47 = (enhancer_5mC_RPKM %>% filter(as.character(group) == g, as.character(category) == c, grepl("CEMT_47", sample)))[, "ENSG"]))
		assign(paste0("IDHmut_venn_", g, "_", c), venn.diagram(get(paste0("IDHmut_gene_", g, "_", c)), filename = NULL, fill = c("green", "red", "blue"), cat.cex = 0.5))
	}
}
grid.arrange(gTree(children = IDHmut_venn_hyper_promoter), gTree(children = IDHmut_venn_hyper_regular), gTree(children = IDHmut_venn_hyper_super), 
						 gTree(children = IDHmut_venn_median_promoter), gTree(children = IDHmut_venn_median_regular), gTree(children = IDHmut_venn_median_super), 
						 gTree(children = IDHmut_venn_hypo_promoter), gTree(children = IDHmut_venn_hypo_regular), gTree(children = IDHmut_venn_hypo_super), 
						 nrow = 3)
pdf("./H3K27ac/venn_IDHmut_gene.pdf", height = 7, width = 7)
grid.arrange(gTree(children = IDHmut_venn_hyper_promoter), gTree(children = IDHmut_venn_hyper_regular), gTree(children = IDHmut_venn_hyper_super), 
						 gTree(children = IDHmut_venn_median_promoter), gTree(children = IDHmut_venn_median_regular), gTree(children = IDHmut_venn_median_super), 
						 gTree(children = IDHmut_venn_hypo_promoter), gTree(children = IDHmut_venn_hypo_regular), gTree(children = IDHmut_venn_hypo_super), 
						 nrow = 3)
dev.off()
for(c in c("promoter", "regular", "super")){
	assign(paste0("enhancer_5mC_homer_known_", c), read.delim(paste0("./H3K27ac/homer/homer.knownResults.summary.", c), as.is = T) %>% 
				 	mutate(significant = ifelse(q <= 0.05, TRUE, FALSE), percent_with_motif = ifelse(group == "hypo", -percent_with_motif, percent_with_motif), type = ifelse(sample == "NPC_GE04", "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut"))))
	tf <- (get(paste0("enhancer_5mC_homer_known_", c)) %>% filter(abs(percent_with_motif) >= 20, q <= 0.05) %>% arrange(percent_with_motif) %>% distinct(motif))$TF
	assign(paste0("enhancer_5mC_homer_known_", c), get(paste0("enhancer_5mC_homer_known_", c)) %>% filter(TF %in% tf) %>% mutate(TF = factor(TF, levels = tf)))
	assign(paste0("enhancer_5mC_homer_known_", c, "_figure"), ggplot(get(paste0("enhancer_5mC_homer_known_", c)), aes(TF, percent_with_motif, fill = type, alpha = significant)) + 
				 	geom_bar(stat = "identity", position = position_dodge()) + 
				 	geom_hline(yintercept = 0) + 
				 	scale_alpha_manual(values = c(0.2, 1), guide = "none") +
				 	scale_y_continuous(breaks = c(-60, -30, 0 , 30, 60), labels = c(60, 30, 0, 30, 60)) + 
				 	xlab("") + 
				 	ylab("Percent of enhancers with motif") + 
				 	coord_flip() + 
				 	theme_bw())
}
ggsave(enhancer_5mC_homer_known_promoter_figure, file = "./H3K27ac/homer/enhancer_5mC_homer_known_promoter_figure.pdf", height = 5, width = 6)
ggsave(enhancer_5mC_homer_known_regular_figure, file = "./H3K27ac/homer/enhancer_5mC_homer_known_regular_figure.pdf", height = 7, width = 6)
ggsave(enhancer_5mC_homer_known_super_figure, file = "./H3K27ac/homer/enhancer_5mC_homer_known_super_figure.pdf", height = 6, width = 6)

## ------- DMR glioma vs NPC -------
### -------- summary ---------
DMR_summary <- read.delim("./DMR/intermediate/DMR.summary.stats", head = T, as.is = T)
DMR_summary_tall <- data.frame(glioma = rep(gsub("_NPC.*", "", DMR_summary$sample), 2), NPC = rep(gsub(".*_NPC_", "NPC_", DMR_summary$sample), 2), DM = rep(c("hyper", "hypo"), each = nrow(DMR_summary)), length = c(DMR_summary$hyper, -DMR_summary$hypo)/10^6)
#DMR_summary_tall <- DMR_summary_tall %>% filter(glioma != "CEMT_21") %>% mutate(glioma = ifelse(glioma == "CEMT_23", paste0("IDHwt\n", glioma), paste0("IDHmut\n", glioma)))
(DMR_summary_figure <- ggplot(DMR_summary_tall, aes(x = glioma, y = length, color = DM, shape = NPC)) + 
	geom_point(position = position_jitter(width = 0.1), size = 3) + 
	geom_hline(yintercept = 0) + 
	xlab("") + 
	ylab("Total DMR length (Mb)") + 
	scale_color_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw() + 
	theme(axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15)))
ggsave(DMR_summary_figure, file = "./DMR/DMR_summary_figure.pdf", height = 5, width = 6)

### ------- visualization ----- 
colname <- c("chr", "start", "end", "ID", "DM", "length")
for(lib in libs){
	print(lib)
	assign(paste0("DMR_", lib), read.delim(paste0("./DMR/DMR.", lib, "_NPC"), head = F, as.is = T, col.names = colname))
	assign(paste0("DMR_", lib, "_figure"), DMR_figures(get(paste0("DMR_", lib)), lib, "NPCs", dirOut = "./DMR/", figures = c("length", "frequency", "circos"), colname = colname, hist_width = 3))
}

### -------- enrichment in genomic regions ------
genomic_breakdown <- read.delim("./DMR/intersect/genomic.breakdown.summary", as.is = T) %>% 
	mutate(sample = gsub("_NPC.*", "", Name), DM = gsub(".*NPC\\.", "", Name), NCpG = NULL, Name = NULL)
genomic_breakdown_tall <- melt(genomic_breakdown, id = c("sample", "DM")) %>% 
	mutate(value = ifelse(DM == "hyper", value, -value))
(genomic_breakdown_figure <- ggplot(genomic_breakdown_tall %>% filter(sample %in% libs), aes(variable, value, fill = DM)) + 
	geom_bar(position = "identity", stat = "identity", width = 0.5) + 
	geom_hline(yintercept = c(-2, 2)) + 
	facet_wrap(~sample) + 
	xlab("") + 
	ylab("Fold enrichment") + 
	scale_fill_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw())
ggsave(genomic_breakdown_figure, file = "./DMR/genomic_breakdown_figure.pdf", height = 6, width = 7)

### -------- hyper CGI with K36me3 ---------
CGI_DMR_hyper_summary <- read.delim("./DMR/CGI/CGI.DMR.hyper.summary", as.is = T)
CGI_hyper_summary <- read.delim("./DMR/CGI/CGI.hyper.H3K36me3.summary", as.is = T)

### -------- % of hyper CpGs in hyper CGIs -----------
CGI_DMR_hyper_DM_all <- read.delim("./DMR/CGI/CGI.DMR.hyper.DM.all", as.is = T)
(CGI_DMR_hyper_DM_figure <- ggplot(CGI_DMR_hyper_DM_all, aes(glioma, percent, fill = NPC)) + 
		geom_boxplot(position = position_dodge()) + 
		xlab("")+
		ylab("Percent of hyper CpGs in each hyper CGIs") + 
		facet_grid(Type ~.) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CGI_DMR_hyper_DM_figure, file = "./DMR/CGI/CGI_DMR_hyper_DM_figure.pdf", height = 6, width = 6)

### -------- distance to closest CGI --------
colname <- c("chr", "start", "end", "ID", "CGI_chr", "CGI_start", "CGI_end", "CGI_ID", "distance", "norm_dis")
for(lib in libs){
	print(lib)
	assign(paste0("DMR_", lib, "_CGI_dis"), rbind((read.delim(paste0("./DMR/CGI_dis/DMR.", lib, "_NPC.hyper.CGI.dis"), head = F, as.is = T, col.names = colname) %>% mutate(DM = "hyper")), 
																								(read.delim(paste0("./DMR/CGI_dis/DMR.", lib, "_NPC.hypo.CGI.dis"), head = F, as.is = T, col.names = colname) %>% mutate(DM = "hypo"))))
	assign(paste0("DMR_", lib, "_CGI_dis_figure"), ggplot(get(paste0("DMR_", lib, "_CGI_dis")) %>% filter(abs(norm_dis) <= 3), aes(norm_dis, color = DM)) + 
				 	geom_density() + 
				 	geom_vline(xintercept = c(-1, 1)) + 
				 	scale_color_manual(name = "", values = c("red", "blue")) + 
				 	xlab("Normalized distance to CGI") + 
				 	ylab("density") + 
				 	ggtitle(lib) + 
				 	theme_bw())
	ggsave(get(paste0("DMR_", lib, "_CGI_dis_figure")), file = paste0("./DMR/CGI_dis/DMR_", lib, "_CGI_dis_figure.pdf"), height = 5, width = 6)
}

### -------- GREAT --------
(GREAT_DMR_CEMT_19_hyper_figure <- enrich_GREAT("CEMT_19_hyper", "CEMT_19_hyper", dirIn = "./DMR/enrich/", dirOut = "./DMR/enrich/", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 8, width = 7))
(GREAT_DMR_CEMT_19_hypo_figure <- enrich_GREAT("CEMT_19_hypo", "CEMT_19_hypo", dirIn = "./DMR/enrich/", dirOut = "./DMR/enrich/", categories = c("GOBP"), height = 3, width = 7))
(GREAT_DMR_CEMT_21_hyper_figure <- enrich_GREAT("CEMT_21_hyper", "CEMT_21_hyper", dirIn = "./DMR/enrich/", dirOut = "./DMR/enrich/", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 8, width = 7))
(GREAT_DMR_CEMT_21_hypo_figure <- enrich_GREAT("CEMT_21_hypo", "CEMT_21_hypo", dirIn = "./DMR/enrich/", dirOut = "./DMR/enrich/", categories = c("MSigPerturbation"), height = 4, width = 7))
(GREAT_DMR_CEMT_22_hyper_figure <- enrich_GREAT("CEMT_22_hyper", "CEMT_22_hyper", dirIn = "./DMR/enrich/", dirOut = "./DMR/enrich/", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 8, width = 7))
# GREAT_DMR_CEMT_22_hypo_figure : no enrichment
(GREAT_DMR_CEMT_23_hyper_figure <- enrich_GREAT("CEMT_23_hyper", "CEMT_23_hyper", dirIn = "./DMR/enrich/", dirOut = "./DMR/enrich/", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 9, width = 7))
(GREAT_DMR_CEMT_23_hypo_figure <- enrich_GREAT("CEMT_23_hypo", "CEMT_23_hypo", dirIn = "./DMR/enrich/", dirOut = "./DMR/enrich/", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))
(GREAT_DMR_CEMT_47_hyper_figure <- enrich_GREAT("CEMT_47_hyper", "CEMT_47_hyper", dirIn = "./DMR/enrich/", dirOut = "./DMR/enrich/", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 9, width = 7))
(GREAT_DMR_CEMT_47_hypo_figure <- enrich_GREAT("CEMT_47_hypo", "CEMT_47_hypo", dirIn = "./DMR/enrich/", dirOut = "./DMR/enrich/", categories = c("MSigPerturbation"), height = 2, width = 7))

### -------- enrichment in chromatin states --------
DMR_ChromHMM_summary <- read.delim("./DMR/intersect/DMR.chromHMM.enrich.summary", as.is = T) 
(DMR_ChromHMM_summary_figure <- ggplot(DMR_ChromHMM_summary, aes(Name, Enrichment, fill = Sample)) + 
	geom_bar(stat = "identity", position = position_dodge()) + 
	facet_grid(. ~ DM) + 
	coord_flip() + 
	xlab("") + 
	ylab("Fold enrichment") +
	theme_bw())
ggsave(DMR_ChromHMM_summary_figure, file = "./DMR/DMR_ChromHMM_summary_figure.pdf")

### -------- enrichment in differentially marked histone mods -------
DMR_DHM_enrich <- read.delim("./DMR/DMR.uniqueHM.enrichment.summary", as.is = T) %>% 
	mutate(HM = gsub("CEMT.*.unique", "Histone modification gain", HM), HM = gsub("NPC.*.unique", "Histone modification loss", HM), sig = ifelse(p.value <= 0.01, "p-value <= 0.01", "p-value > 0.01"))
(DMR_DHM_enrich_figure <- ggplot(DMR_DHM_enrich, aes(Sample, log2(Fold), fill = Mark, alpha = sig)) + 
	geom_bar(stat = "identity", position = position_dodge(), width = 0.8) + 
	scale_alpha_manual(guide = "none", values = c(1, 0.2)) + 
	facet_grid(HM ~ DMR) + 
	coord_flip() + 
	xlab("") + 
	ylab("log2 Fold Enrichemnt") + 
	theme_bw())
ggsave(DMR_DHM_enrich_figure, file = "./DMR/DMR_DHM_enrich_figure.pdf")
DMR_enhancer_enrich <- DMR_DHM_enrich %>% filter(Sample != "CEMT_21", Mark == "H3K27ac")
pdf("./DMR/DMR_enhancer_enrich_venn.pdf")
DMR_enhancer_enrich_venn <- draw.pairwise.venn(area1 = 9128, area2 = 10940, cross.area = 826, category = c("Hypermethylation", "Loss of H3K27ac"), fill = c("orange", "blue"), cat.pos = c(200, 160), cat.cex = 1.5)
dev.off()

### -------- associated with DE genes ----------
DMR_DE <- read.delim("./DMR/DE/DMR.DE.summary", as.is = T) %>% mutate(Significant = p_Fisher < 0.05) %>% filter(Sample != "CEMT_21")
(DMR_DE_figure <- ggplot(DMR_DE, aes(DM, Percent_intersect, fill = DE, alpha = Significant)) + 
	geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
	scale_fill_manual(name = "", values = c("blue", "red")) + 
	scale_alpha_manual(values = c(0.3, 1), guide = "none") + 
	facet_wrap(~ Sample) + 
	xlab("") + 
	ylab("Fraction of DE genes") + 
	theme_bw())
ggsave(DMR_DE_figure, file = "./DMR/DMR_DE_figure.pdf")
DMR_DE_HM_summary <- read.delim("./DMR/DE/DMR.DE.HM.summary", as.is = T)
hyper_UP_K27ac_summary <- read.delim("./DMR/DE/hyper.UP_2FC.H2K27ac.summary", as.is = T)

save(list = c(ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "venn"), ls(pattern = "dend")),
		 file = "/projects/epigenomics2/users/lli/glioma/WGBS/WGBS.Rdata")

