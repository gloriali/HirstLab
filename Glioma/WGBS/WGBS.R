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
load("/projects/epigenomics2/users/lli/glioma/WGBS/WGBS.Rdata")
setwd("/projects/epigenomics2/users/lli/glioma/WGBS/DMR")
libs <- c("CEMT_19", "CEMT_21", "CEMT_22", "CEMT_23", "CEMT_47", "TCGA060128_IDHmut", "TCGA141401_IDHwt", "TCGA141454_IDHwt", "TCGA143477_IDHwt", "TCGA161460_IDHmut", "TCGA191788_IDHmut")

## ------- 5mC modifiers RPKM --------
modifier_RPKM <- read.delim("../DNAme_regulators.RPKM", as.is = T) %>% select(-contains("Brain")) %>% melt(id = c("ENSG", "gene")) %>% 
	mutate(type = ifelse(variable == "CEMT_23", "IDHwt", ifelse(grepl("CEMT", variable), "IDHmut", "NPC")))
(modifier_RPKM_figure <- ggplot(modifier_RPKM, aes(gene, value, color = type)) + 
		geom_point(position = position_jitter(width = 0.2)) + 
		coord_flip() + 
		xlab("") + 
		ylab("RPKM") + 
		theme_bw())
ggsave(modifier_RPKM_figure, file = "../modifier_RPKM_figure.pdf", height = 5, width = 6)

## ------- CpG coverage ---------
coverage <- read.delim("../coverage.txt", as.is = T) %>% mutate(type = ifelse(grepl("NPC", sample), "NPC", ifelse(grepl("IDHwt|CEMT_23", sample), "IDHwt", "IDHmut")), category = gsub("[0-9]*_.*", "", gsub("\\..*", "", sample)))
(coverage_figure <- ggplot(coverage, aes(coverage, N/1e6, color = category)) + 
		geom_line(aes(group = sample)) + 
		ylab("No. of million CpGs") + 
		coord_cartesian(xlim = c(0, 50)) + 
		theme_bw())
ggsave(coverage_figure, file = "../coverage_figure.pdf", height = 5, width = 5)

## ------- 5mC distribution------
quantile_5mC <- read.delim("../qc.5mC.quantile", as.is = T) %>% mutate(type = ifelse(grepl("CGI", sample), "CGI", "genome"), category = ifelse(grepl("NPC", sample), "NPC", ifelse(grepl("IDHwt|CEMT_23", sample), "IDHwt", "IDHmut")), sample = gsub("_[gC].*", "", gsub("_IDH.*", "", sample))) %>%
	filter(sample != "CEMT_21")
(quantile_5mC_figure <- ggplot(quantile_5mC, aes(x = sample, lower = lower, middle = median, upper = upper, ymin = ymin, ymax = ymax, fill = category)) + 
	geom_boxplot(stat = "identity", outlier.shape = NA, width = 0.5, color = "grey") + 
	facet_grid(type ~ .) + 
	xlab("") + 
	ylab("Fractional methylation") + 
	coord_flip() + 
	theme_bw())
ggsave(quantile_5mC_figure, file = "../quantile_5mC_figure.pdf", height = 7, width = 6)

## ------- clustering -----------------
all_5mC <- read.table("/projects/epigenomics2/users/lli/glioma/WGBS/all.combine.5mC", sep = " ", head = F, row.names = 1, col.names = c("ID", "IDHmut_CEMT_19", "IDHmut_CEMT_21", "IDHmut_CEMT_22", "IDHwt_CEMT_23", "IDHmut_CEMT_47", "NPC_Cortex02", "NPC_Cortex04", "NPC_GE02", "NPC_GE04", "IDHmut_TCGA060128", "IDHwt_TCGA141401", "IDHwt_TCGA141454", "IDHwt_TCGA143477", "IDHmut_TCGA161460", "IDHmut_TCGA191788"))
spearman_all_5mC <- cor(all_5mC, method = "spearman")
pearson_all_5mC <- cor(all_5mC, method = "pearson")
write.table(spearman_all_5mC, file = "/projects/epigenomics2/users/lli/glioma/WGBS/spearman_all_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(pearson_all_5mC, file = "/projects/epigenomics2/users/lli/glioma/WGBS/pearson_all_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
CGI_5mC <- read.table("/projects/epigenomics2/users/lli/glioma/WGBS/CGI.combine.5mC", sep = " ", head = F, row.names = 1, col.names = c("ID", "IDHmut_CEMT_19", "IDHmut_CEMT_21", "IDHmut_CEMT_22", "IDHwt_CEMT_23", "IDHmut_CEMT_47", "NPC_Cortex02", "NPC_Cortex04", "NPC_GE02", "NPC_GE04", "IDHmut_TCGA060128", "IDHwt_TCGA141401", "IDHwt_TCGA141454", "IDHwt_TCGA143477", "IDHmut_TCGA161460", "IDHmut_TCGA191788"))
spearman_CGI_5mC <- cor(CGI_5mC, method = "spearman")
pearson_CGI_5mC <- cor(CGI_5mC, method = "pearson")
write.table(spearman_CGI_5mC, file = "/projects/epigenomics2/users/lli/glioma/WGBS/spearman_CGI_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(pearson_CGI_5mC, file = "/projects/epigenomics2/users/lli/glioma/WGBS/pearson_CGI_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
op <- par(mar = c(5, 4, 4, 8))
spearman_all_5mC <- read.delim("../spearman_all_5mC.cor", as.is = T)
spearman_all_5mC_dend <-  hclust(1 - spearman_all_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(spearman_all_5mC_dend, main = "all CpGs spearman", horiz = TRUE)
pearson_all_5mC <- read.delim("../pearson_all_5mC.cor", as.is = T)
pearson_all_5mC_dend <- hclust(1 - pearson_all_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(pearson_all_5mC_dend, main = "all CpGs pearson", horiz = TRUE)
spearman_CGI_5mC <- read.delim("../spearman_CGI_5mC.cor", as.is = T)
spearman_CGI_5mC_dend <- hclust(1 - spearman_CGI_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(spearman_CGI_5mC_dend, main = "CGI spearman", horiz = TRUE)
pearson_CGI_5mC <- read.delim("../pearson_CGI_5mC.cor", as.is = T)
pearson_CGI_5mC_dend <-  hclust(1 - pearson_CGI_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
plot(pearson_CGI_5mC_dend, main = "CGI pearson", horiz = TRUE)
par(op)
pdf("../cluster.pdf", height = 5, width = 8)
op <- par(mar = c(5, 4, 4, 8))
plot(spearman_all_5mC_dend, main = "genome-wide spearman", horiz = TRUE)
plot(pearson_all_5mC_dend, main = "genome-wide pearson", horiz = TRUE)
plot(spearman_CGI_5mC_dend, main = "CGI spearman", horiz = TRUE)
plot(pearson_CGI_5mC_dend, main = "CGI pearson", horiz = TRUE)
par(op)
dev.off()

## ------- changes at CGI edges -------
CGI_edge <- read.delim("../CGI_edge/CGI.edge.all", head = F, col.names = c("chr", "start", "end", "ID", "mC", "CGI", "edge", "dis", "sample")) %>% 
	mutate(edge = revalue(edge, c("L" = "5-prime", "R" = "3-prime")))
CGI_edge_CpG_ID <- CGI_edge %>% group_by(ID) %>% summarise(n = n()) %>% filter(n == 9)
CGI_edge <- CGI_edge %>% filter(ID %in% CGI_edge_CpG_ID$ID)
(CGI_edge_figure <- ggplot(CGI_edge, aes(-dis, mC, color = sample)) + 
	geom_smooth() + 
	facet_wrap(~ edge, scales = "free_x") + 
	xlab("Distance to CGI edge (bp)") + 
	ylab("Fractional methylation") + 
	theme_bw())
ggsave(CGI_edge_figure, file = "../CGI_edge/CGI_edge_figure.pdf", width = 8, height = 5)
CGI_edge_delta <- read.delim("../CGI_edge/CGI.edge.delta.all", head = F, col.names = c("ID", "CGI", "edge", "dis", "delta", "sample", "NPC")) %>% 
	mutate(NPC = gsub(".+-", "", NPC), edge = revalue(edge, c("L" = "5-prime", "R" = "3-prime")))
(CGI_edge_delta_figure <- ggplot(CGI_edge_delta, aes(-dis, delta, color = NPC)) + 
		geom_smooth() + 
		geom_hline(yintercept = 0) + 
		facet_grid(sample ~ edge, scales = "free_x") + 
		xlab("Distance to CGI edge (bp)") + 
		ylab("Difference in fractional methylation\nglioma - NPC") + 
		theme_bw())
ggsave(CGI_edge_delta_figure, file = "../CGI_edge/CGI_edge_delta_figure.pdf", width = 8, height = 8)

## ------- 5mC at enhancer regions --------
enhancer_5mC <- read.delim("../H3K27ac/enhancer.5mC", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample")) %>% filter(sample != "CEMT_21") %>% 
	mutate(category = gsub(".*_", "", ID), type = ifelse(sample == "NPC_GE04", "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = paste0(type, "\n", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median")))
enhancer_5mC_summary <- merge(enhancer_5mC %>% group_by(sample, category, group) %>% summarize(N=n()), enhancer_5mC %>% group_by(sample, category) %>% summarize(Total=n())) %>% 
	mutate(percent = N/Total, signal = 3, fractional = ifelse(group == "hypo", 0.15, ifelse(group == "hyper", 0.85, 0.5)), labels = paste0(round(percent*100, 0), "%"))
(enhancer_5mC_figure <- ggplot(enhancer_5mC, aes(fractional, log10(signal))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.01, 0.02, 0.04, 0.08, 0.3, 0.5, 1), guide = "none") +		
		geom_vline(xintercept = 0.3, color = "white") + 
		geom_vline(xintercept = 0.7, color = "white") + 
		geom_text(data = enhancer_5mC_summary, aes(fractional, signal, label = labels), color = "white", size = 3) + 
		facet_grid(sample ~ category) + 
		xlab("Fractional methylation") + 
		ylab("log10 H3K27ac normalized signal") + 
		theme_bw())
ggsave(enhancer_5mC_figure, file = "../H3K27ac/enhancer_5mC_figure.pdf", height = 6, width = 6)
enhancer_5mC_RPKM <- read.delim("../H3K27ac/enhancer.5mC.RPKM", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "category", "signal", "fractional", "sample", "ENSG", "RPKM")) %>% 
	mutate(RPKM = RPKM + 1e-5, type = ifelse(sample == "NPC_GE04", "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = paste0(type, "\n", sample), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median")), group = factor(group, levels = c("hypo", "median", "hyper")))
(enhancer_5mC_RPKM_figure <- ggplot(enhancer_5mC_RPKM, aes(group, log10(RPKM), fill = type)) + 
		geom_violin() + 
		coord_cartesian(ylim = c(-1.5, 2.5)) + 
		facet_grid(sample ~ category) + 
		xlab("") + 
		theme_bw())
ggsave(enhancer_5mC_RPKM_figure, file = "../H3K27ac/enhancer_5mC_RPKM_figure.pdf", height = 6, width = 6)
RPKM_p <- data.frame(sample = c(), category = c(), p = c())
for(lib in c("CEMT_19", "CEMT_22", "CEMT_23", "CEMT_47", "GE04")){
	for(c in c("promoter", "regular", "super")){
		RPKM_p <- rbind(RPKM_p, data.frame(sample = lib, category = c, p = t.test((enhancer_5mC_RPKM %>% filter(category == c, grepl(lib, sample), group == "hyper"))$RPKM, (enhancer_5mC_RPKM %>% filter(category == c, grepl(lib, sample), group == "hypo"))$RPKM)$p.value))
	}
}
enhancer_5mC_TCGA <- read.delim("../H3K27ac/TCGA.enhancer.coverage.5mC", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "signal", "fractional", "sample")) %>%
	mutate(category = gsub(".*_", "", ID), type = gsub(".*_", "", sample), sample = paste0(type, "\n", gsub("_.*", "", sample)), group = ifelse(fractional < 0.3, "hypo", ifelse(fractional > 0.7, "hyper", "median")))
enhancer_5mC_TCGA_summary <- merge(enhancer_5mC_TCGA %>% group_by(sample, category, group) %>% summarize(N=n()), enhancer_5mC_TCGA %>% group_by(sample, category) %>% summarize(Total=n())) %>% 
	mutate(percent = N/Total, signal = 2.8, fractional = ifelse(group == "hypo", 0.15, ifelse(group == "hyper", 0.85, 0.5)), labels = paste0(round(percent*100, 0), "%"))
(enhancer_5mC_TCGA_figure <- ggplot(enhancer_5mC_TCGA, aes(fractional, log10(signal))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.01, 0.02, 0.04, 0.08, 0.3, 0.5, 1), guide = "none") +		
		geom_vline(xintercept = 0.3, color = "white") + 
		geom_vline(xintercept = 0.7, color = "white") + 
		geom_text(data = enhancer_5mC_TCGA_summary, aes(fractional, signal, label = labels), color = "white", size = 3) + 
		facet_grid(sample ~ category) + 
		xlab("Fractional methylation") + 
		ylab("log10 H3K27ac normalized signal") + 
		theme_bw())
ggsave(enhancer_5mC_TCGA_figure, file = "../H3K27ac/enhancer_5mC_TCGA_figure.pdf", height = 6, width = 6)
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
pdf("../H3K27ac/venn_IDHmut_gene.pdf", height = 7, width = 7)
grid.arrange(gTree(children = IDHmut_venn_hyper_promoter), gTree(children = IDHmut_venn_hyper_regular), gTree(children = IDHmut_venn_hyper_super), 
						 gTree(children = IDHmut_venn_median_promoter), gTree(children = IDHmut_venn_median_regular), gTree(children = IDHmut_venn_median_super), 
						 gTree(children = IDHmut_venn_hypo_promoter), gTree(children = IDHmut_venn_hypo_regular), gTree(children = IDHmut_venn_hypo_super), 
						 nrow = 3)
dev.off()
for(c in c("promoter", "regular", "super")){
	assign(paste0("enhancer_5mC_homer_known_", c), read.delim(paste0("../H3K27ac/homer/homer.knownResults.summary.", c), as.is = T) %>% 
				 	mutate(significant = ifelse(q <= 0.05, TRUE, FALSE), percent_with_motif = ifelse(group == "hypo", -percent_with_motif, percent_with_motif), type = ifelse(sample == "NPC_GE04", "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut"))))
	tf <- (get(paste0("enhancer_5mC_homer_known_", c)) %>% filter(abs(percent_with_motif) >= 20, q <= 0.05) %>% arrange(percent_with_motif) %>% distinct(motif))$TF
	assign(paste0("enhancer_5mC_homer_known_", c), get(paste0("enhancer_5mC_homer_known_", c)) %>% filter(TF %in% tf) %>% mutate(TF = factor(TF, levels = tf)))
	assign(paste0("enhancer_5mC_homer_known_", c, "_figure"), ggplot(get(paste0("enhancer_5mC_homer_known_", c)), aes(TF, percent_with_motif, fill = type, color = significant)) + 
				 	geom_bar(stat = "identity", position = position_dodge()) + 
				 	geom_hline(yintercept = 0) + 
				 	scale_color_manual(values = c("grey", "black")) +
				 	scale_y_continuous(breaks = c(-60, -30, 0 , 30, 60), labels = c(60, 30, 0, 30, 60)) + 
				 	xlab("") + 
				 	ylab("Percent of enhancers with motif") + 
				 	coord_flip() + 
				 	theme_bw())
}
ggsave(enhancer_5mC_homer_known_promoter_figure, file = "../H3K27ac/homer/enhancer_5mC_homer_known_promoter_figure.pdf", height = 5, width = 6)
ggsave(enhancer_5mC_homer_known_regular_figure, file = "../H3K27ac/homer/enhancer_5mC_homer_known_regular_figure.pdf", height = 7, width = 6)
ggsave(enhancer_5mC_homer_known_super_figure, file = "../H3K27ac/homer/enhancer_5mC_homer_known_super_figure.pdf", height = 6, width = 6)

## ------- 5mC at CTCF loss regions -------
CTCF_loss_5mC <- read.delim("../CTCF/CTCF.loss.5mC", as.is = T) %>% mutate(type = ifelse(grepl("NPC", sample), "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = gsub("NPC.", "", sample))
(CTCF_loss_5mC_figure <- ggplot(CTCF_loss_5mC, aes(sample, fractional, fill = type)) + 
		geom_boxplot() + 
		xlab("") + 
		ylab("Fractional methylation") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CTCF_loss_5mC_figure, file = "../CTCF/CTCF_loss_5mC_figure.pdf", height = 4, width = 5)

## ------- DMR glioma vs NPC -------
### -------- summary ---------
DMR_summary <- read.delim("./intermediate/DMR.summary.stats", head = T, as.is = T)
DMR_summary_tall <- data.frame(glioma = rep(gsub("_NPC.*", "", DMR_summary$sample), 2), NPC = rep(gsub(".*_NPC\\.", "", DMR_summary$sample), 2), DM = rep(c("hyper", "hypo"), each = nrow(DMR_summary)), length = c(DMR_summary$hyper, -DMR_summary$hypo)/10^6)
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
ggsave(DMR_summary_figure, file = "DMR_summary_figure.pdf", height = 5, width = 6)

### ------- visualization ----- 
colname <- c("chr", "start", "end", "ID", "DM", "length")
for(lib in libs){
	print(lib)
	assign(paste0("DMR_", lib), read.delim(paste0("DMR.", lib, "_NPC"), head = F, as.is = T, col.names = colname))
	assign(paste0("DMR_", lib, "_figure"), DMR_figures(get(paste0("DMR_", lib)), lib, "NPCs", figures = c("length", "frequency", "circos"), colname = colname, hist_width = 3))
}

### -------- enrichment in genomic regions ------
genomic_breakdown <- read.delim("./intersect/genomic.breakdown.summary", as.is = T) %>% 
	mutate(sample = gsub("_NPC.*", "", Name), DM = gsub(".*NPC\\.", "", Name), NCpG = NULL, Name = NULL)
genomic_breakdown_tall <- melt(genomic_breakdown, id = c("sample", "DM")) %>% 
	mutate(value = ifelse(DM == "hyper", value, -value))
genomic_breakdown_figure <- ggplot(genomic_breakdown_tall, aes(variable, value, fill = DM)) + 
	geom_bar(position = "identity", stat = "identity", width = 0.5) + 
	geom_hline(yintercept = c(-2, 2)) + 
	facet_wrap(~sample) + 
	xlab("") + 
	ylab("Fold enrichment") + 
	scale_fill_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw()
ggsave(genomic_breakdown_figure, file = "genomic_breakdown_figure.pdf", height = 6, width = 7)

### -------- hyper CGI with K36me3 ---------
CGI_DMR_hyper_summary <- read.delim("./CGI/CGI.DMR.hyper.summary", as.is = T)
CGI_hyper_summary <- read.delim("./CGI/CGI.hyper.H3K36me3.summary", as.is = T)

### -------- % of hyper CpGs in hyper CGIs -----------
CGI_DMR_hyper_DM_all <- read.delim("./CGI/CGI.DMR.hyper.DM.all", as.is = T)
(CGI_DMR_hyper_DM_figure <- ggplot(CGI_DMR_hyper_DM_all, aes(glioma, percent, fill = NPC)) + 
		geom_boxplot(position = position_dodge()) + 
		xlab("")+
		ylab("Percent of hyper CpGs in each hyper CGIs") + 
		facet_grid(Type ~.) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CGI_DMR_hyper_DM_figure, file = "./CGI/CGI_DMR_hyper_DM_figure.pdf", height = 6, width = 6)

### -------- distance to closest CGI --------
colname <- c("chr", "start", "end", "ID", "CGI_chr", "CGI_start", "CGI_end", "CGI_ID", "distance", "norm_dis")
for(lib in libs){
	print(lib)
	assign(paste0("DMR_", lib, "_CGI_dis"), rbind((read.delim(paste0("./CGI_dis/DMR.", lib, "_NPC.hyper.CGI.dis"), head = F, as.is = T, col.names = colname) %>% mutate(DM = "hyper")), 
																								(read.delim(paste0("./CGI_dis/DMR.", lib, "_NPC.hypo.CGI.dis"), head = F, as.is = T, col.names = colname) %>% mutate(DM = "hypo"))))
	assign(paste0("DMR_", lib, "_CGI_dis_figure"), ggplot(get(paste0("DMR_", lib, "_CGI_dis")) %>% filter(abs(norm_dis) <= 3), aes(norm_dis, color = DM)) + 
				 	geom_density() + 
				 	geom_vline(xintercept = c(-1, 1)) + 
				 	scale_color_manual(name = "", values = c("red", "blue")) + 
				 	xlab("Normalized distance to CGI") + 
				 	ylab("density") + 
				 	ggtitle(lib) + 
				 	theme_bw())
	ggsave(get(paste0("DMR_", lib, "_CGI_dis_figure")), file = paste0("DMR_", lib, "_CGI_dis_figure.pdf"), height = 5, width = 6)
}

### -------- GREAT --------
(GREAT_DMR_CEMT_19_hyper <- enrich_GREAT("CEMT_19_hyper", "CEMT_19_hyper", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 8, width = 7))
(GREAT_DMR_CEMT_19_hypo <- enrich_GREAT("CEMT_19_hypo", "CEMT_19_hypo", categories = c("GOBP"), height = 3, width = 7))
(GREAT_DMR_CEMT_21_hyper <- enrich_GREAT("CEMT_21_hyper", "CEMT_21_hyper", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 8, width = 7))
(GREAT_DMR_CEMT_21_hypo <- enrich_GREAT("CEMT_21_hypo", "CEMT_21_hypo", categories = c("MSigPerturbation"), height = 4, width = 7))
(GREAT_DMR_CEMT_22_hyper <- enrich_GREAT("CEMT_22_hyper", "CEMT_22_hyper", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 8, width = 7))
# GREAT_DMR_CEMT_22_hypo : no enrichment
(GREAT_DMR_CEMT_23_hyper <- enrich_GREAT("CEMT_23_hyper", "CEMT_23_hyper", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 9, width = 7))
(GREAT_DMR_CEMT_23_hypo <- enrich_GREAT("CEMT_23_hypo", "CEMT_23_hypo", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))
(GREAT_DMR_CEMT_47_hyper <- enrich_GREAT("CEMT_47_hyper", "CEMT_47_hyper", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 9, width = 7))
(GREAT_DMR_CEMT_47_hypo <- enrich_GREAT("CEMT_47_hypo", "CEMT_47_hypo", categories = c("MSigPerturbation"), height = 2, width = 7))

### -------- Intersect --------
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
DMR_intersect <- data.frame(sample1 = NA, sample2 = NA, hyper1 = NA, hyper2 = NA, hyper_intersect = NA, hyper_percent1 = NA, hyper_percent2 = NA, hyper_jaccard = NA, hypo1 = NA, hypo2 = NA, hypo_intersect = NA, hypo_percent1 = NA, hypo_percent2 = NA, hypo_jaccard = NA)
for(i1 in 1:(length(libs)-1)){
	for(i2 in (i1+1):length(libs)){
		print(c(i1, i2, libs[i1], libs[i2]))
		a1_hyper <- as.numeric(system(paste0("less DMR.", libs[i1], "_NPC.hyper.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		a2_hyper <- as.numeric(system(paste0("less DMR.", libs[i2], "_NPC.hyper.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		cross_hyper <- as.numeric(system(paste0(BEDTOOLS, "/intersectBed -a  DMR.", libs[i1], "_NPC.hyper.bed -b DMR.", libs[i2], "_NPC.hyper.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		jaccard_hyper <- as.numeric(str_split_fixed(system(paste0(BEDTOOLS, "bedtools jaccard -a DMR.", libs[i1], "_NPC.hyper.bed -b DMR.", libs[i2], "_NPC.hyper.bed"), intern = T)[2], pattern = "\\t", n = 4)[1,3])
		assign(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hyper"), 
					 draw.pairwise.venn(a1_hyper, a2_hyper, cross_hyper, category = c(libs[i1], libs[i2]), cat.pos = 0, ext.pos = 180))
		a1_hypo <- as.numeric(system(paste0("less DMR.", libs[i1], "_NPC.hypo.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		a2_hypo <- as.numeric(system(paste0("less DMR.", libs[i2], "_NPC.hypo.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		cross_hypo <- as.numeric(system(paste0(BEDTOOLS, "/intersectBed -a  DMR.", libs[i1], "_NPC.hypo.bed -b DMR.", libs[i2], "_NPC.hypo.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		jaccard_hypo <- as.numeric(str_split_fixed(system(paste0(BEDTOOLS, "bedtools jaccard -a DMR.", libs[i1], "_NPC.hypo.bed -b DMR.", libs[i2], "_NPC.hypo.bed"), intern = T)[2], pattern = "\\t", n = 4)[1,3])
		assign(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hypo"), 
					 draw.pairwise.venn(a1_hypo, a2_hypo, cross_hypo, category = c(libs[i1], libs[i2]), cat.pos = 0, ext.pos = 180))
		DMR_intersect <- rbind(DMR_intersect, data.frame(sample1 = libs[i1], sample2 = libs[i2], hyper1 = a1_hyper, hyper2 = a2_hyper, hyper_intersect = cross_hyper, hyper_percent1 = cross_hyper/a1_hyper, hyper_percent2 = cross_hyper/a2_hyper, hyper_jaccard = jaccard_hyper, hypo1 = a1_hypo, hypo2 = a2_hypo, hypo_intersect = cross_hypo, hypo_percent1 = cross_hypo/a1_hypo, hypo_percent2 = cross_hypo/a2_hypo, hypo_jaccard = jaccard_hypo))
	}
}
DMR_intersect <- na.omit(DMR_intersect)
DMR_jaccard_hyper <- data.frame(sample1 = c(DMR_intersect$sample1, DMR_intersect$sample2, libs), sample2 = c(DMR_intersect$sample2, DMR_intersect$sample1, libs), jaccard = c(rep(DMR_intersect$hyper_jaccard, 2), rep(1, length(libs)))) 
DMR_jaccard_hyper_figure <- ggplot(DMR_jaccard_hyper, aes(sample1, sample2, fill = jaccard)) + 
	geom_tile() + 
	xlab("") + 
	ylab("") + 
	ggtitle("Jaccard similarity - hyper DMRs") + 
	theme_bw()
ggsave(DMR_jaccard_hyper_figure, file = "DMR_jaccard_hyper_figure.pdf", height = 5, width = 6)
DMR_jaccard_hypo <- data.frame(sample1 = c(DMR_intersect$sample1, DMR_intersect$sample2, libs), sample2 = c(DMR_intersect$sample2, DMR_intersect$sample1, libs), jaccard = c(rep(DMR_intersect$hypo_jaccard, 2), rep(1, length(libs)))) 
DMR_jaccard_hypo_figure <- ggplot(DMR_jaccard_hypo, aes(sample1, sample2, fill = jaccard)) + 
	geom_tile() + 
	xlab("") + 
	ylab("") + 
	ggtitle("Jaccard similarity - hypo DMRs") + 
	theme_bw()
ggsave(DMR_jaccard_hypo_figure, file = "DMR_jaccard_hypo_figure.pdf", height = 5, width = 6)
pdf("Venn_DMR.pdf")
for(i1 in 1:(length(libs)-1)){
	for(i2 in (i1+1):length(libs)){
		grid.arrange(gTree(children = get(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hyper"))), gTree(children = get(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hypo"))), nrow = 2)
		grid.text("hyper", x = unit(0.1, "npc"), y = unit(0.98, "npc"))
		grid.text("hypo", x = unit(0.1, "npc"), y = unit(0.48, "npc"))
	}
}
dev.off()

### -------- enrichment in chromatin states --------
DMR_ChromHMM_summary <- read.delim("./CpG/DMR.chromHMM.enrich.summary", as.is = T) 
(DMR_ChromHMM_summary_figure <- ggplot(DMR_ChromHMM_summary, aes(Name, Enrichment, fill = Sample)) + 
	geom_bar(stat = "identity", position = position_dodge()) + 
	facet_grid(. ~ DM) + 
	coord_flip() + 
	xlab("") + 
	ylab("Fold enrichment") +
	theme_bw())
ggsave(DMR_ChromHMM_summary_figure, file = "DMR_ChromHMM_summary_figure.pdf")

### -------- enrichment in differentially marked histone mods -------
DMR_DHM_enrich <- read.delim("./DMR.uniqueHM.enrichment.summary", as.is = T) %>% 
	mutate(HM = gsub("CEMT.*.unique", "Histone modification gain", HM), HM = gsub("NPC.*.unique", "Histone modification loss", HM), sig = ifelse(p.value <= 0.01, "p-value <= 0.01", "p-value > 0.01"))
(DMR_DHM_enrich_figure <- ggplot(DMR_DHM_enrich, aes(Sample, log2(Fold), fill = Mark, color = sig)) + 
	geom_bar(stat = "identity", position = position_dodge(), width = 0.8) + 
	scale_color_manual(name = "", values = c("red", "transparent")) + 
	facet_grid(HM ~ DMR) + 
	coord_flip() + 
	xlab("") + 
	ylab("log2 Fold Enrichemnt") + 
	theme_bw())
ggsave(DMR_DHM_enrich_figure, file = "DMR_DHM_enrich_figure.pdf")
DMR_enhancer_enrich <- DMR_DHM_enrich %>% filter(Sample != "CEMT_21", Mark == "H3K27ac")
pdf("DMR_enhancer_enrich_venn.pdf")
DMR_enhancer_enrich_venn <- draw.pairwise.venn(area1 = 9128, area2 = 10940, cross.area = 826, category = c("Hypermethylation", "Loss of H3K27ac"), fill = c("orange", "blue"), cat.pos = c(200, 160), cat.cex = 1.5)
dev.off()

### -------- associated with DE genes ----------
DMR_DE <- read.delim("./DE/DMR.DE.summary", as.is = T) %>% mutate(Significant = p_Fisher < 0.01)
(DMR_DE_figure <- ggplot(DMR_DE, aes(DM, Percent_intersect, fill = DE, color = Significant)) + 
	geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
	scale_fill_manual(name = "", values = c("blue", "red")) + 
	scale_color_manual(values = c("transparent", "green")) + 
	facet_wrap(~ Sample) + 
	xlab("") + 
	ylab("Fraction of DE genes") + 
	theme_bw())
ggsave(DMR_DE_figure, file = "DMR_DE_figure.pdf")
DMR_DE_HM_summary <- read.delim("./DE/DMR.DE.HM.summary", as.is = T)
hyper_UP_K27ac_summary <- read.delim("./DE/hyper.UP_2FC.H2K27ac.summary", as.is = T)

save(list = c("DMR_intersect",  
							ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "venn"), ls(pattern = "dend"),  
							ls(pattern = "GREAT_DMR_*"), ls(pattern = "Venn_DMR_*")),
		 file = "/projects/epigenomics2/users/lli/glioma/WGBS/WGBS.Rdata")

