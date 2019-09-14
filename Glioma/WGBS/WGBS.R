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
library(car)
library(limma)
library(pheatmap)
source('~/HirstLab/Pipeline/R/DMR.figures.R')
source('~/HirstLab/Pipeline/R/enrich.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/glioma/WGBS/")
load("/projects/epigenomics3/epigenomics3_results/users/lli/glioma/WGBS/WGBS.Rdata")

## ------- 5mC modifiers RPKM --------
DNAme_regulators_RPKM <- read.delim("DNAme_regulators.RPKM", as.is = T) %>% melt(id = c("ENSG", "gene")) %>% 
	mutate(type = gsub("\\..*", "", variable))
(DNAme_regulators_RPKM_figure <- ggplot(DNAme_regulators_RPKM, aes(gene, log10(value), color = type)) + 
		geom_point(position = position_jitter(width = 0.2)) + 
		scale_color_manual(values = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[5], l = 65, c = 100))) + 
		coord_flip() + 
		guides(color = guide_legend(title = NULL)) + 
		xlab("") + 
		ylab("log10 RPKM") + 
		theme_bw())
ggsave(DNAme_regulators_RPKM_figure, file = "DNAme_regulators_RPKM_figure.pdf", height = 5, width = 6)

## ------- CpG coverage ---------
qc_5mC_coverage <- read.delim("qc_5mC_coverage.txt", as.is = T) %>% mutate(category = gsub("\\..*", "", sample))
(qc_5mC_coverage_figure <- ggplot(qc_5mC_coverage, aes(coverage, N/1e6, color = category)) + 
		geom_line(aes(group = sample)) + 
		guides(color = guide_legend(title = NULL)) + 
		ylab("No. of million CpGs") + 
		coord_cartesian(xlim = c(0, 50)) + 
		theme_bw())
ggsave(qc_5mC_coverage_figure, file = "qc_5mC_coverage_figure.pdf", height = 5, width = 7)

## ------- 5mC distribution --------
qc_5mC_profile <- read.delim("qc_5mC_profile.txt", as.is = T) %>% mutate(category = gsub("\\..*", "", sample), N = ifelse(type == "genome", N/1e6, N/1e3)) 
(qc_5mC_profile_figure <- ggplot(qc_5mC_profile, aes(fractional, N, color = category, group = sample)) + 
		geom_smooth(se = F, span = 0.1, size = 0.5) + 
		facet_grid(type ~ ., scales = "free_y") + 
		guides(color = guide_legend(title = NULL), override.aes = list(size = 5)) + 
		coord_cartesian(ylim = c(0, 8)) + 
		xlab("Fractional methylation") + 
		ylab("No. of million CpGs                    No. of thousand CGIs") + 
		theme_bw())
ggsave(qc_5mC_profile_figure, file = "qc_5mC_profile_figure.pdf", height = 5, width = 6)
qc_5mC_quantile <- read.delim("qc_5mC_quantile.txt", as.is = T) %>% mutate(category = gsub("\\..*", "", sample)) 
(qc_5mC_quantile_figure <- ggplot(qc_5mC_quantile %>% filter(!(category %in% c("NHAR", "Normal", "MGG"))), aes(x = sample, lower = lower, middle = median, upper = upper, ymin = ymin, ymax = ymax, fill = category)) + 
		geom_boxplot(stat = "identity", color = "grey") + 
		scale_fill_manual(values = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[5], l = 65, c = 100))) + 
		facet_grid(type ~ .) + 
		guides(color = guide_legend(title = NULL)) + 
		xlab("") + 
		ylab("Fractional methylation") + 
		coord_flip() + 
		theme_bw())
ggsave(qc_5mC_quantile_figure, file = "qc_5mC_quantile_figure.pdf", height = 7, width = 8)

## ------- cluster ---------
genome_5mC <- read.table("matrix_genome.5mC", sep = " ", head = T, row.names = 1) %>% select(-contains("Normal."), -contains("NHAR"), -contains("MGG"))
spearman_genome_5mC <- cor(genome_5mC, method = "spearman")
pearson_genome_5mC <- cor(genome_5mC, method = "pearson")
write.table(spearman_genome_5mC, file = "cluster_spearman_genome_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(pearson_genome_5mC, file = "cluster_pearson_genome_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
CGI_5mC <- read.table("matrix_CGI.5mC", sep = " ", head = T, row.names = 1) %>% select(-contains("Normal."), -contains("NHAR"), -contains("MGG"))
spearman_CGI_5mC <- cor(CGI_5mC, method = "spearman")
pearson_CGI_5mC <- cor(CGI_5mC, method = "pearson")
write.table(spearman_CGI_5mC, file = "cluster_spearman_CGI_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(pearson_CGI_5mC, file = "cluster_pearson_CGI_5mC.cor", sep = "\t", quote = F, row.names = T, col.names = T)
spearman_genome_5mC_dend <-  hclust(1 - spearman_genome_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
pearson_genome_5mC_dend <- hclust(1 - pearson_genome_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
spearman_CGI_5mC_dend <- hclust(1 - spearman_CGI_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
pearson_CGI_5mC_dend <-  hclust(1 - pearson_CGI_5mC %>% as.dist, method = "ward.D2") %>% as.dendrogram
pdf("cluster.pdf", height = 5, width = 8)
op <- par(mar = c(5, 4, 4, 10))
plot(spearman_genome_5mC_dend, main = "genome-wide spearman", horiz = TRUE)
plot(pearson_genome_5mC_dend, main = "genome-wide pearson", horiz = TRUE)
plot(spearman_CGI_5mC_dend, main = "CGI spearman", horiz = TRUE)
plot(pearson_CGI_5mC_dend, main = "CGI pearson", horiz = TRUE)
par(op)
dev.off()
ann <- data.frame(category = gsub("\\..*", "", row.names(pearson_genome_5mC))) 
rownames(ann) <- row.names(pearson_genome_5mC)
ann_color <- list(category = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), 
															 hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100), 
															 hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100), 
															 hcl(h = seq(15, 375, length = 5 + 1)[5], l = 65, c = 100)))
names(ann_color$category) <- c("IDHmut", "IDHwt", "NormalAdjacent", "NPC")
heatmap_pearson_genome_5mC <- pheatmap(pearson_genome_5mC, main = "genome-wide 5mC pearson", filename = "./heatmap_pearson_genome_5mC.pdf", height = 7, width = 6, color = colorRampPalette(c("forest green ","white","purple"))(100), clustering_method = "ward.D2", cutree_cols = 4, cutree_rows = 4, annotation_row = ann, annotation_col = ann, annotation_colors = ann_color, show_rownames = F)
CGI_5mC_variable <- CGI_5mC %>% mutate(id = rownames(.), sd = apply(., 1, sd)) %>% filter(sd > quantile(sd, 0.8))
rownames(CGI_5mC_variable) <- CGI_5mC_variable$id
heatmap_pearson_CGI_5mC <- pheatmap(CGI_5mC_variable %>% select(-id, -sd), main = "20% variable CGI 5mC pearson", filename = "./heatmap_pearson_CGI_5mC.pdf", height = 8, width = 6, color = colorRampPalette(c("forest green ","white","purple"))(100), clustering_distance_cols = "correlation", clustering_method = "ward.D2", cutree_rows = 4, cutree_cols = 4, annotation_col = ann, annotation_colors = ann_color, show_rownames = F)

#plot(density(genome_5mC$IDHmut_CEMT_19 - genome_5mC$IDHwt_CEMT_23), xlim = c(-.6, .6))
#lines(density(genome_5mC$IDHmut_TCGA161460 - genome_5mC$IDHwt_TCGA141454), col = "red")
#lines(density(genome_5mC$IDHmut_CEMT_19 - genome_5mC$NPC_GE04), col = "blue")
#lines(density(genome_5mC$IDHmut_CEMT_19 - genome_5mC$IDHmut_CEMT_22), col = "green")

## ------- DMR limma ------------------
cutoff <- 0.05
DM_logit <- genome_5mC %>% logit(adjust = 0.0001)
DM_design <- data.frame(Group = relevel(factor(gsub("\\..*", "", colnames(DM_logit))), ref = "NPC"), row.names = colnames(DM_logit))
DM_DesMat <- model.matrix(~ Group, DM_design)
DM_fit <- lmFit(DM_logit, DM_DesMat)
DM_fitEb <- eBayes(DM_fit)
IDHmut_NPC_DM <- topTable(DM_fitEb, coef = 'GroupIDHmut', number = Inf, adjust.method="BH", p.value = cutoff) 
write.table(IDHmut_NPC_DM, file = "./limma/DM.IDHmut_NPC_limma.txt", sep = "\t", quote = F, row.names = T, col.names = T)
IDHwt_NPC_DM <- topTable(DM_fitEb, coef = 'GroupIDHwt', number = Inf, adjust.method="BH", p.value = cutoff) 
write.table(IDHwt_NPC_DM, file = "./limma/DM.IDHwt_NPC_limma.txt", sep = "\t", quote = F, row.names = T, col.names = T)
NormalAdjacent_NPC_DM <- topTable(DM_fitEb, coef = 'GroupNormalAdjacent', number = Inf, adjust.method="BH", p.value = cutoff) 
write.table(NormalAdjacent_NPC_DM, file = "./limma/DM.NormalAdjacent_NPC_limma.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# DM_DesMat <- model.matrix(~ 0 + Group, DM_design)
# colnames(DM_DesMat) <- gsub("Group", "", colnames(DM_DesMat))
# DM_fit <- lmFit(DM_logit, DM_DesMat)
# contrast.matrix <- makeContrasts(IDHmut-NPC, IDHwt-NPC, NormalAdjacent-NPC, levels = DM_DesMat)
# DM_fit <- contrasts.fit(DM_fit, contrast.matrix)
# DM_fitEb <- eBayes(DM_fit)
# IDHmut_NPC_DM <- topTable(DM_fitEb, coef = 'IDHmut - NPC', number = Inf, adjust.method="BH", p.value = cutoff) 

### ------- DMR visualization ----- 
DMR_summary <- read.delim("./limma/DMR.summary.stats", as.is = T) %>% select(Name, hyper, hypo) %>% melt(id = "Name")
(DMR_summary_figure <- ggplot(DMR_summary, aes(Name, value/1e6, fill = variable)) + 
		geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
		xlab("") + 
		ylab("Total DMR length (Mb)") + 
		scale_fill_manual(name = "", values = c("red", "blue")) + 
		theme_bw())
ggsave(DMR_summary_figure, file = "./limma/DMR_summary_figure.pdf", height = 4, width = 5)
DMR_limma_IDHmut_NPC <- read.delim("./limma/DMR.IDHmut_NPC.s500.c3", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "DM", "count", "length")) %>% mutate(chr = paste0("chr", chr))
DMR_limma_IDHmut_NPC_figure <- DMR_figures(DMR_limma_IDHmut_NPC, "IDHmut", "NPC", dirOut = "./limma/", figures = c("count", "length", "frequency", "circos"), hist_width = 3)
DMR_limma_IDHwt_NPC <- read.delim("./limma/DMR.IDHwt_NPC.s500.c3", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "DM", "count", "length")) %>% mutate(chr = paste0("chr", chr))
DMR_limma_IDHwt_NPC_figure <- DMR_figures(DMR_limma_IDHwt_NPC, "IDHwt", "NPC", dirOut = "./limma/", figures = c("count", "length", "frequency", "circos"), hist_width = 3)
DMR_limma_NormalAdjacent_NPC <- read.delim("./limma/DMR.NormalAdjacent_NPC.s500.c3", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "DM", "count", "length")) %>% mutate(chr = paste0("chr", chr))
DMR_limma_NormalAdjacent_NPC_figure <- DMR_figures(DMR_limma_NormalAdjacent_NPC, "NormalAdjacent", "NPC", dirOut = "./limma/", figures = c("count", "length", "frequency", "circos"), hist_width = 3)

### -------- enrichment in genomic regions ------
genomic_enrichment <- read.delim("./limma/intersect/genomic.breakdown.summary", as.is = T) %>% 
	mutate(DM = gsub(".*\\.", "", Name), Name = gsub("\\..*", "", Name), NCpG = NULL)
genomic_enrichment_tall <- melt(genomic_enrichment, id = c("DM", "Name")) %>% 
	mutate(value = ifelse(DM == "hyper", value, -value))
(genomic_enrichment_figure <- ggplot(genomic_enrichment_tall, aes(variable, value, fill = DM)) + 
	geom_bar(position = "identity", stat = "identity", width = 0.5) + 
	#geom_hline(yintercept = c(-2, 2)) + 
	facet_wrap(~ Name) + 
	xlab("") + 
	ylab("Fold enrichment") + 
	scale_fill_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw())
ggsave(genomic_enrichment_figure, file = "./limma/genomic_enrichment_figure.pdf", height = 6, width = 8)
DMR_CGI_category <- read.delim("./limma/intersect/DMR.CGI.category", as.is = T)
(DMR_CGI_category_figure <- ggplot(DMR_CGI_category, aes(sample, fill = category)) + 
		geom_bar(position = "fill") + 
		xlab("") + 
		ylab("") + 
		facet_grid(DM ~ .) + 
		coord_flip() + 
		theme_bw())
ggsave(DMR_CGI_category_figure, file = "./limma/DMR_CGI_category_figure.pdf", height = 3, width = 6)
DMR_enhancer_category <- read.delim("./limma/intersect/DMR.enhancer.category", as.is = T)
(DMR_enhancer_category_figure <- ggplot(DMR_enhancer_category, aes(sample, fill = category)) + 
		geom_bar(position = "fill") + 
		scale_fill_manual(values = c("black", "light grey")) + 
		xlab("") + 
		ylab("") + 
		facet_grid(DM ~ .) + 
		coord_flip() + 
		theme_bw())
ggsave(DMR_enhancer_category_figure, file = "./limma/DMR_enhancer_category_figure.pdf", height = 3, width = 6)

### -------- associated with DE genes ----------
DMR_DE <- read.delim("./limma/DE/DMR.DE.summary", as.is = T) 
(DMR_DE_figure <- ggplot(DMR_DE, aes(DM, Percent_intersect, fill = DE)) + 
		geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
		scale_fill_manual(name = "", values = c("blue", "red")) + 
		facet_grid(region ~ Sample) + 
		xlab("") + 
		ylab("Fraction of DE genes") + 
		theme_bw())
ggsave(DMR_DE_figure, file = "./limma/DMR_DE_figure.pdf", height = 5, width = 6)

### ------- enhancer DMR homer ------------
RPKM <- read.delim("/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/RPKM/RPKM.matrix", as.is = T)
HUGO <- read.delim("/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.EnsID_sorted.HUGO", as.is = T, head = F, col.names = c("ENSG", "HUGO")) %>% distinct(HUGO) %>% mutate(HUGO = gsub("-", ".", HUGO))
rownames(HUGO) <- HUGO$HUGO
e <- 1e-6
homer_enhancer_IDHmut_NPC_hyper <- read.delim("./limma/homer/IDHmut_NPC.hyper/knownResults.txt", as.is = T) %>%
	mutate(TF = toupper(gsub("\\(.*", "", Motif.Name)), ENSG = HUGO[TF, "ENSG"], Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	filter(Percent_with_motif >= 20, q.value..Benjamini. < 0.01) %>% select(TF, ENSG, Percent_with_motif) %>% merge(RPKM) %>% arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_enhancer_IDHmut_NPC_hyper_figure <- ggplot(homer_enhancer_IDHmut_NPC_hyper, aes(TF, Percent_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "red") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percentage of enhancers with motif") + 
		ggtitle("IDHmut vs NPC hyper") + 
		theme_bw() + 
		theme(plot.margin = unit(c(0.5,0.2,0.5,0.5), "cm")))
homer_enhancer_IDHmut_NPC_hyper_RPKM <- homer_enhancer_IDHmut_NPC_hyper %>% melt(id = c("TF", "ENSG", "Percent_with_motif")) %>% mutate(type = gsub("\\..*", "", variable))
(homer_enhancer_IDHmut_NPC_hyper_RPKM_figure <- ggplot(homer_enhancer_IDHmut_NPC_hyper_RPKM, aes(TF, log10(value+e), color = type)) + 
		geom_point(position = position_jitter(width = 0.1)) + 
		scale_color_manual(values = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[5], l = 65, c = 100))) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("IDHmut vs NPC hyper") + 
		theme_bw() + 
		theme(plot.title = element_text(color = "white"), axis.text.y = element_text(size = 0), plot.margin = unit(c(0.5,0.5,0.5,-0.5), "cm")))
pdf("./limma/homer_enhancer_IDHmut_NPC_hyper_figure.pdf", height = 4, width = 7)
layOut(list(homer_enhancer_IDHmut_NPC_hyper_figure , 1, 1),
			 list(homer_enhancer_IDHmut_NPC_hyper_RPKM_figure, 1, 2)) 
dev.off()
homer_enhancer_IDHmut_NPC_hypo <- read.delim("./limma/homer/IDHmut_NPC.hypo/knownResults.txt", as.is = T) %>%
	mutate(TF = toupper(gsub("\\(.*", "", Motif.Name)), ENSG = HUGO[TF, "ENSG"], Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	filter(Percent_with_motif >= 20, q.value..Benjamini. < 0.01) %>% select(TF, ENSG, Percent_with_motif) %>% merge(RPKM) %>% arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_enhancer_IDHmut_NPC_hypo_figure <- ggplot(homer_enhancer_IDHmut_NPC_hypo, aes(TF, Percent_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percentage of enhancers with motif") + 
		ggtitle("IDHmut vs NPC hypo") + 
		theme_bw() + 
		theme(plot.margin = unit(c(0.5,0.2,0.5,0.5), "cm")))
homer_enhancer_IDHmut_NPC_hypo_RPKM <- homer_enhancer_IDHmut_NPC_hypo %>% melt(id = c("TF", "ENSG", "Percent_with_motif")) %>% mutate(type = gsub("\\..*", "", variable))
(homer_enhancer_IDHmut_NPC_hypo_RPKM_figure <- ggplot(homer_enhancer_IDHmut_NPC_hypo_RPKM, aes(TF, log10(value+e), color = type)) + 
		geom_point(position = position_jitter(width = 0.1)) + 
		scale_color_manual(values = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[5], l = 65, c = 100))) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("IDHmut vs NPC hypo") + 
		theme_bw() + 
		theme(plot.title = element_text(color = "white"), axis.text.y = element_text(size = 0), plot.margin = unit(c(0.5,0.5,0.5,-0.5), "cm")))
pdf("./limma/homer_enhancer_IDHmut_NPC_hypo_figure.pdf", height = 4, width = 7)
layOut(list(homer_enhancer_IDHmut_NPC_hypo_figure , 1, 1),
			 list(homer_enhancer_IDHmut_NPC_hypo_RPKM_figure, 1, 2)) 
dev.off()

homer_enhancer_IDHwt_NPC_hyper <- read.delim("./limma/homer/IDHwt_NPC.hyper/knownResults.txt", as.is = T) %>%
	mutate(TF = toupper(gsub("\\(.*", "", Motif.Name)), ENSG = HUGO[TF, "ENSG"], Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	filter(Percent_with_motif >= 20, q.value..Benjamini. < 0.01) %>% select(TF, ENSG, Percent_with_motif) %>% merge(RPKM) %>% arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_enhancer_IDHwt_NPC_hyper_figure <- ggplot(homer_enhancer_IDHwt_NPC_hyper, aes(TF, Percent_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "red") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percentage of enhancers with motif") + 
		ggtitle("IDHwt vs NPC hyper") + 
		theme_bw() + 
		theme(plot.margin = unit(c(0.5,0.2,0.5,0.5), "cm")))
homer_enhancer_IDHwt_NPC_hyper_RPKM <- homer_enhancer_IDHwt_NPC_hyper %>% melt(id = c("TF", "ENSG", "Percent_with_motif")) %>% mutate(type = gsub("\\..*", "", variable))
(homer_enhancer_IDHwt_NPC_hyper_RPKM_figure <- ggplot(homer_enhancer_IDHwt_NPC_hyper_RPKM, aes(TF, log10(value+e), color = type)) + 
		geom_point(position = position_jitter(width = 0.1)) + 
		scale_color_manual(values = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[5], l = 65, c = 100))) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("IDHwt vs NPC hyper") + 
		theme_bw() + 
		theme(plot.title = element_text(color = "white"), axis.text.y = element_text(size = 0), plot.margin = unit(c(0.5,0.5,0.5,-0.5), "cm")))
pdf("./limma/homer_enhancer_IDHwt_NPC_hyper_figure.pdf", height = 4, width = 7)
layOut(list(homer_enhancer_IDHwt_NPC_hyper_figure , 1, 1),
			 list(homer_enhancer_IDHwt_NPC_hyper_RPKM_figure, 1, 2)) 
dev.off()
homer_enhancer_IDHwt_NPC_hypo <- read.delim("./limma/homer/IDHwt_NPC.hypo/knownResults.txt", as.is = T) %>%
	mutate(TF = toupper(gsub("\\(.*", "", Motif.Name)), ENSG = HUGO[TF, "ENSG"], Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	filter(Percent_with_motif >= 20, q.value..Benjamini. < 0.01) %>% select(TF, ENSG, Percent_with_motif) %>% merge(RPKM) %>% arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_enhancer_IDHwt_NPC_hypo_figure <- ggplot(homer_enhancer_IDHwt_NPC_hypo, aes(TF, Percent_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percentage of enhancers with motif") + 
		ggtitle("IDHwt vs NPC hypo") + 
		theme_bw() + 
		theme(plot.margin = unit(c(0.5,0.2,0.5,0.5), "cm")))
homer_enhancer_IDHwt_NPC_hypo_RPKM <- homer_enhancer_IDHwt_NPC_hypo %>% melt(id = c("TF", "ENSG", "Percent_with_motif")) %>% mutate(type = gsub("\\..*", "", variable))
(homer_enhancer_IDHwt_NPC_hypo_RPKM_figure <- ggplot(homer_enhancer_IDHwt_NPC_hypo_RPKM, aes(TF, log10(value+e), color = type)) + 
		geom_point(position = position_jitter(width = 0.1)) + 
		scale_color_manual(values = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[5], l = 65, c = 100))) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("IDHwt vs NPC hypo") + 
		theme_bw() + 
		theme(plot.title = element_text(color = "white"), axis.text.y = element_text(size = 0), plot.margin = unit(c(0.5,0.5,0.5,-0.5), "cm")))
pdf("./limma/homer_enhancer_IDHwt_NPC_hypo_figure.pdf", height = 4, width = 7)
layOut(list(homer_enhancer_IDHwt_NPC_hypo_figure , 1, 1),
			 list(homer_enhancer_IDHwt_NPC_hypo_RPKM_figure, 1, 2)) 
dev.off()

homer_enhancer_NormalAdjacent_NPC_hyper <- read.delim("./limma/homer/NormalAdjacent_NPC.hyper/knownResults.txt", as.is = T) %>%
	mutate(TF = toupper(gsub("\\(.*", "", Motif.Name)), ENSG = HUGO[TF, "ENSG"], Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	filter(Percent_with_motif >= 20, q.value..Benjamini. < 0.01) %>% select(TF, ENSG, Percent_with_motif) %>% merge(RPKM) %>% arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_enhancer_NormalAdjacent_NPC_hyper_figure <- ggplot(homer_enhancer_NormalAdjacent_NPC_hyper, aes(TF, Percent_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "red") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percentage of enhancers with motif") + 
		ggtitle("NormalAdjacent vs NPC hyper") + 
		theme_bw() + 
		theme(plot.margin = unit(c(0.5,0.2,0.5,0.5), "cm")))
homer_enhancer_NormalAdjacent_NPC_hyper_RPKM <- homer_enhancer_NormalAdjacent_NPC_hyper %>% melt(id = c("TF", "ENSG", "Percent_with_motif")) %>% mutate(type = gsub("\\..*", "", variable))
(homer_enhancer_NormalAdjacent_NPC_hyper_RPKM_figure <- ggplot(homer_enhancer_NormalAdjacent_NPC_hyper_RPKM, aes(TF, log10(value+e), color = type)) + 
		geom_point(position = position_jitter(width = 0.1)) + 
		scale_color_manual(values = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[5], l = 65, c = 100))) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("NormalAdjacent vs NPC hyper") + 
		theme_bw() + 
		theme(plot.title = element_text(color = "white"), axis.text.y = element_text(size = 0), plot.margin = unit(c(0.5,0.5,0.5,-0.5), "cm")))
pdf("./limma/homer_enhancer_NormalAdjacent_NPC_hyper_figure.pdf", height = 4, width = 7)
layOut(list(homer_enhancer_NormalAdjacent_NPC_hyper_figure , 1, 1),
			 list(homer_enhancer_NormalAdjacent_NPC_hyper_RPKM_figure, 1, 2)) 
dev.off()
homer_enhancer_NormalAdjacent_NPC_hypo <- read.delim("./limma/homer/NormalAdjacent_NPC.hypo/knownResults.txt", as.is = T) %>%
	mutate(TF = toupper(gsub("\\(.*", "", Motif.Name)), ENSG = HUGO[TF, "ENSG"], Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	filter(Percent_with_motif >= 20, q.value..Benjamini. < 0.01) %>% select(TF, ENSG, Percent_with_motif) %>% merge(RPKM) %>% arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_enhancer_NormalAdjacent_NPC_hypo_figure <- ggplot(homer_enhancer_NormalAdjacent_NPC_hypo, aes(TF, Percent_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percentage of enhancers with motif") + 
		ggtitle("NormalAdjacent vs NPC hypo") + 
		theme_bw() + 
		theme(plot.margin = unit(c(0.5,0.2,0.5,0.5), "cm")))
homer_enhancer_NormalAdjacent_NPC_hypo_RPKM <- homer_enhancer_NormalAdjacent_NPC_hypo %>% melt(id = c("TF", "ENSG", "Percent_with_motif")) %>% mutate(type = gsub("\\..*", "", variable))
(homer_enhancer_NormalAdjacent_NPC_hypo_RPKM_figure <- ggplot(homer_enhancer_NormalAdjacent_NPC_hypo_RPKM, aes(TF, log10(value+e), color = type)) + 
		geom_point(position = position_jitter(width = 0.1)) + 
		scale_color_manual(values = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[5], l = 65, c = 100))) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("NormalAdjacent vs NPC hypo") + 
		theme_bw() + 
		theme(plot.title = element_text(color = "white"), axis.text.y = element_text(size = 0), plot.margin = unit(c(0.5,0.5,0.5,-0.5), "cm")))
pdf("./limma/homer_enhancer_NormalAdjacent_NPC_hypo_figure.pdf", height = 4, width = 7)
layOut(list(homer_enhancer_NormalAdjacent_NPC_hypo_figure , 1, 1),
			 list(homer_enhancer_NormalAdjacent_NPC_hypo_RPKM_figure, 1, 2)) 
dev.off()

homer_enhancer <- rbind(homer_enhancer_IDHmut_NPC_hyper %>% mutate(sample = "IDHmut hyper", color = "hyper"), 
												homer_enhancer_IDHmut_NPC_hypo %>% mutate(sample = "IDHmut hypo", color = "hypo"), 
												homer_enhancer_IDHwt_NPC_hyper %>% mutate(sample = "IDHwt hyper", color = "hyper"), 
												homer_enhancer_IDHwt_NPC_hypo %>% mutate(sample = "IDHwt hypo", color = "hypo"), 
												homer_enhancer_NormalAdjacent_NPC_hyper %>% mutate(sample = "NormalAdjacent hyper", color = "hyper"), 
												homer_enhancer_NormalAdjacent_NPC_hypo %>% mutate(sample = "NormalAdjacent hypo", color = "hypo")) %>% dcast(TF ~ sample, value.var = "color", fill = "none") %>% melt(id = "TF")
(homer_enhancer_figure <- ggplot(homer_enhancer, aes(variable, TF, fill = value)) + 
		geom_tile() + 
		scale_fill_manual(guide = "none", values = c("red", "blue", "light grey")) + 
		labs(x = "", y = "") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(homer_enhancer_figure, file = "./limma/homer_enhancer_figure.pdf", height = 6, width = 2)

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
(GREAT_DMR_IDHmut_NPC_hyper_figure <- enrich_GREAT("IDHmut_NPC_hyper", "IDHmut_NPC_hyper", dirIn = "./enrich/", dirOut = "./", height = 8, width = 7, top = 10))
(GREAT_DMR_IDHmut_NPC_hypo_figure <- enrich_GREAT("IDHmut_NPC_hypo", "IDHmut_NPC_hypo", dirIn = "./enrich/", dirOut = "./", height = 8, width = 7, top = 10))
(GREAT_DMR_IDHwt_NPC_hyper_figure <- enrich_GREAT("IDHwt_NPC_hyper", "IDHwt_NPC_hyper", dirIn = "./enrich/", dirOut = "./", height = 8, width = 7, top = 10))

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

## ------- DMR pairwise -------
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
colname <- c("chr", "start", "end", "ID", "DM", "length")
for(lib in libs){
	print(lib)
	assign(paste0("DMR_", lib), read.delim(paste0("./DMR/DMR.", lib, "_NPC"), head = F, as.is = T, col.names = colname))
	assign(paste0("DMR_", lib, "_figure"), DMR_figures(get(paste0("DMR_", lib)), lib, "NPCs", dirOut = "./DMR/", figures = c("length", "frequency", "circos"), colname = colname, hist_width = 3))
}

## ------- DMR NB141 ------------------
DMR_NB141_summary <- read.delim("./NB141/DMR.summary.stats", as.is = T) %>% mutate(sample = gsub("_Normal.NB141", "", gsub("IDHmut.", "", sample))) %>% select(sample, hyper, hypo) %>% melt(id = "sample")
(DMR_NB141_summary_figure <- ggplot(DMR_NB141_summary, aes(sample, value/1e6, fill = variable)) + 
		geom_bar(stat = "identity", width = 0.5, position = position_dodge()) + 
		coord_flip() + 
		xlab("") + 
		ylab("Total length (Mb)") + 
		ggtitle("IDHmut vs NB141") + 
		theme_bw())
ggsave(DMR_NB141_summary_figure, file = "./NB141/DMR_NB141_summary_figure.pdf", height = 5, width = 5)
### genomic breakdown
genomic_breakdown <- read.delim("./NB141/intersect/genomic.breakdown.summary", as.is = T) %>% 
	mutate(DM = gsub(".*\\.", "", Name), Name = gsub("_Normal.*", "", gsub("IDHmut.", "", Name)), NCpG = NULL)
genomic_breakdown_tall <- melt(genomic_breakdown, id = c("DM", "Name")) %>% 
	mutate(value = ifelse(DM == "hyper", value, -value))
(genomic_breakdown_figure <- ggplot(genomic_breakdown_tall, aes(variable, value, fill = DM)) + 
		geom_bar(position = "identity", stat = "identity", width = 0.5) + 
		#geom_hline(yintercept = c(-2, 2)) + 
		facet_wrap(~ Name) + 
		xlab("") + 
		ylab("Fold enrichment") + 
		scale_fill_manual(name = "", values = c("red", "blue")) + 
		coord_flip() + 
		theme_bw())
ggsave(genomic_breakdown_figure, file = "./NB141/genomic_breakdown_figure.pdf", height = 7, width = 9)
### enhancer DMR homer
homer_enhancer_IDHmut_NB141_hyper <- read.delim("./NB141/homer/IDHmut.hyper/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	filter(Percent_with_motif >= 20, q.value..Benjamini. < 0.01) %>% arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_enhancer_IDHmut_NB141_hyper_figure <- ggplot(homer_enhancer_IDHmut_NB141_hyper, aes(TF, Percent_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("IDHmut vs NB141 hyper") + 
		theme_bw())
ggsave(homer_enhancer_IDHmut_NB141_hyper_figure, file = "./NB141/homer_enhancer_IDHmut_NB141_hyper_figure.pdf", height = 6, width = 5)
homer_enhancer_IDHmut_NB141_hypo <- read.delim("./NB141/homer/IDHmut.hypo/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	filter(Percent_with_motif >= 20, q.value..Benjamini. < 0.01) %>% arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_enhancer_IDHmut_NB141_hypo_figure <- ggplot(homer_enhancer_IDHmut_NB141_hypo, aes(TF, Percent_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("IDHmut vs NB141 hypo") + 
		theme_bw())
ggsave(homer_enhancer_IDHmut_NB141_hypo_figure, file = "./NB141/homer_enhancer_IDHmut_NB141_hypo_figure.pdf", height = 6, width = 5)
homer_IDHmut_NB141_hyper_vitC_hMeDIP <- read.delim("./NB141/homer/IDHmut.hyper.vitc_hMeDIP/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	filter(Percent_with_motif >= 20, q.value..Benjamini. < 0.01) %>% arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_IDHmut_NB141_hyper_vitC_hMeDIP_figure <- ggplot(homer_IDHmut_NB141_hyper_vitC_hMeDIP, aes(TF, Percent_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("IDHmut 5mC hyper intersect with vitC 5hmC gain") + 
		theme_bw())
ggsave(homer_IDHmut_NB141_hyper_vitC_hMeDIP_figure, file = "./NB141/homer_IDHmut_NB141_hyper_vitC_hMeDIP_figure.pdf", height = 6, width = 5)

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
(enhancer_5mC_bar_figure <- ggplot(enhancer_5mC_summary %>% mutate(DNAme = factor(group, levels = c("hyper", "median", "hypo")), type = gsub("\\n.*", "", sample), sample = gsub(".*\\n", "", sample)), aes(sample, percent, fill = DNAme)) + 
		geom_bar(stat = "identity") + 
		facet_grid(category ~ type, scales = "free_x", space = "free_x") + 
		xlab("") + 
		ylab("% of enhancers") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(enhancer_5mC_bar_figure, file = "./H3K27ac/enhancer_5mC_bar_figure.pdf", height = 4, width = 4)
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
IDHmut_promoter_hyper_gene <- intersect((enhancer_5mC_RPKM %>% filter(as.character(group) == "hyper", as.character(category) == "promoter", grepl("CEMT_19", sample)))[, "ENSG"], intersect((enhancer_5mC_RPKM %>% filter(as.character(group) == "hyper", as.character(category) == "promoter", grepl("CEMT_22", sample)))[, "ENSG"], (enhancer_5mC_RPKM %>% filter(as.character(group) == "hyper", as.character(category) == "promoter", grepl("CEMT_47", sample)))[, "ENSG"]))
write.table(IDHmut_promoter_hyper_gene, file = "./H3K27ac/IDHmut_promoter_hyper_gene.txt", quote = F, row.names = F, col.names = F)
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
IDHmut_promoter_hyper_states <- read.delim("./H3K27ac/IDHmut.enhancer.5mC.promoter.hyper.states", as.is = T, head = F) %>% rename(chr = V1, start = V2, end = V3, ID = V4, sample = V7, DE = V12) %>% 
	mutate(DM = paste0(V8, "-", V10), H3K27ac = paste0(gsub(".*_", "", V13), "-", gsub(".*_", "", V14)), H3K4me1 = paste0(gsub(".*_", "", V15), "-", gsub(".*_", "", V16)), H3K4me3 = paste0(gsub(".*_", "", V17), "-", gsub(".*_", "", V18)), H3K27me3 = paste0(gsub(".*_", "", V19), "-", gsub(".*_", "", V20))) %>%
	select(-starts_with("V")) 
IDHmut_promoter_hyper_states_summary <- IDHmut_promoter_hyper_states %>% group_by(sample, DM, H3K27ac, H3K4me1, H3K4me3, H3K27me3) %>% summarize(N = n()) %>% filter(N >= 10) %>% dcast(DM + H3K27ac + H3K4me1 + H3K4me3 + H3K27me3 ~ sample) %>% na.omit() 
IDHmut_promoter_hyper_DE_summary <- IDHmut_promoter_hyper_states %>% group_by(sample, DE, DM, H3K27ac) %>% summarize(N = n()) %>%  dcast(DE + DM + H3K27ac ~ sample) %>% na.omit() 
(enrich_IDHmut.enhancer.5mC.promoter.hyper.states.UP.DM <- enrich(name = "IDHmut.enhancer.5mC.promoter.hyper.states.UP.DM", dirIn = "./H3K27ac/", dirOut = "./H3K27ac/", fdr = 0.01, p = "PValue", erminej = F, height = 5, width = 8))
TF.ENSG <- read.delim("./H3K27ac/homer/TF.ENSG", as.is = T) %>% merge(RPKM) %>% select(-contains("Brain"), -CEMT_21) 
for(c in c("promoter", "regular", "super")){
	assign(paste0("enhancer_5mC_homer_known_", c), read.delim(paste0("./H3K27ac/homer/homer.knownResults.summary.", c), as.is = T) %>% 
				 	mutate(TF = gsub("\\(.*", "", TF), significant = ifelse(q <= 0.05 & abs(percent_with_motif) >= 20, TRUE, FALSE), percent_with_motif = ifelse(group == "hypo", -percent_with_motif, percent_with_motif), type = ifelse(sample == "NPC_GE04", "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut"))))
	tf <- (get(paste0("enhancer_5mC_homer_known_", c)) %>% filter(abs(percent_with_motif) >= 20, q <= 0.05) %>% arrange(percent_with_motif) %>% distinct(motif))$TF
	assign(paste0("enhancer_5mC_homer_known_", c), get(paste0("enhancer_5mC_homer_known_", c)) %>% filter(TF %in% tf) %>% mutate(TF = factor(TF, levels = tf)) %>% group_by(TF, type, group) %>% summarise(n = sum(significant)) %>% mutate(significant = ifelse(n>0, TRUE, FALSE)))
	assign(paste0("enhancer_5mC_homer_known_", c, "_figure"), ggplot(get(paste0("enhancer_5mC_homer_known_", c)), aes(type, TF, fill = significant)) + 
				 	geom_tile() + 
				 	facet_wrap(~group, nrow = 1) + 
				 	scale_fill_manual(values = c("white", "red")) +
				 	xlab("") + 
				 	ylab("") + 
				 	ggtitle(paste0(c, " enhancers")) + 
				 	theme_bw() + 
				 	theme(axis.text.x = element_text(angle = 90, size = 12)))
	# assign(paste0("enhancer_5mC_homer_known_", c, "_RPKM"), get(paste0("enhancer_5mC_homer_known_", c)) %>% merge(TF.ENSG) %>% select(-type:significant) %>% melt(id = c("ENSG", "TF")) %>% mutate(type = ifelse(variable == "CEMT_23", "IDHwt", ifelse(grepl("CEMT", variable), "IDHmut", "NPC"))))
	# assign(paste0("enhancer_5mC_homer_known_", c, "_RPKM_figure"), ggplot(get(paste0("enhancer_5mC_homer_known_", c, "_RPKM")), aes(TF, log10(value), color = type)) + 
	# 			 	geom_point(position = position_jitter(width = 0.2)) + 
	# 			 	coord_flip() + 
	# 			 	guides(color = guide_legend(title = NULL)) + 
	# 			 	xlab("") + 
	# 			 	ylab("log10 RPKM") + 
	# 			 	theme_bw())
	# pdf(paste0("./H3K27ac/homer/enhancer_5mC_homer_known_", c, "_figure.pdf"), height = 5, width = 6)
	# grid.newpage()
	# pushViewport(viewport(layout = grid.layout(1, 2)))
	# print(get(paste0("enhancer_5mC_homer_known_", c, "_figure")) + theme(plot.margin = unit(c(1,0,1,1), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	# print(get(paste0("enhancer_5mC_homer_known_", c, "_RPKM_figure")) + theme(plot.margin = unit(c(1,1,1,-0.5), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	# dev.off()
}
ggsave(enhancer_5mC_homer_known_promoter_figure, file = "./H3K27ac/homer/enhancer_5mC_homer_known_promoter_figure.pdf", height = 7, width = 5)
ggsave(enhancer_5mC_homer_known_regular_figure, file = "./H3K27ac/homer/enhancer_5mC_homer_known_regular_figure.pdf", height = 7, width = 5)
ggsave(enhancer_5mC_homer_known_super_figure, file = "./H3K27ac/homer/enhancer_5mC_homer_known_super_figure.pdf", height = 7, width = 5)
e <- 1e-5
IDHmut_enhancer_5mC_promoter_hyper_DE_DM_IDHwt <- read.delim("./H3K27ac/IDHmut.enhancer.5mC.promoter.hyper.states.DE.DM.CEMT_23.5mC", as.is = T, head = F, col.names = c("ID", "wt5mC"))
IDHmut_enhancer_5mC_promoter_hyper_DE_DM <- read.delim("./H3K27ac/IDHmut.enhancer.5mC.promoter.hyper.states.DE.DM", as.is = T, head = F, col.names = c("chr", "start", "end", "ID", "mut5mC", "NPC5mC", "ENSG", "mutK27ac", "NPCK27ac")) %>% merge(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_IDHwt) %>%
	merge(RPKM) %>% select(-CEMT_21, -contains("Brain")) %>% distinct(ID) %>% arrange(mut5mC - NPC5mC) %>% mutate(ID = factor(ID, levels = ID), mutK27ac = gsub("H3K27ac_glioma_", "", mutK27ac), NPCK27ac = gsub("H3K27ac_NPC_", "", NPCK27ac), wtRPKM = CEMT_23) 
IDHmut_enhancer_5mC_promoter_hyper_DE_DM$mutRPKM <- rowMeans(select(IDHmut_enhancer_5mC_promoter_hyper_DE_DM, CEMT_19, CEMT_22, CEMT_47))
IDHmut_enhancer_5mC_promoter_hyper_DE_DM$NPCRPKM <- rowMeans(select(IDHmut_enhancer_5mC_promoter_hyper_DE_DM, contains("Cortex"), contains("GE")))
IDHmut_enhancer_5mC_promoter_hyper_DE_DM <- IDHmut_enhancer_5mC_promoter_hyper_DE_DM %>% mutate(delta_5mC = mut5mC - NPC5mC, RPKM_FC = log2(mutRPKM+e/NPCRPKM+e)) %>% filter(RPKM_FC > 0)
IDHmut_enhancer_5mC_promoter_hyper_DE_DM_5mC <- IDHmut_enhancer_5mC_promoter_hyper_DE_DM %>% select(ID, mut5mC, wt5mC, NPC5mC) %>% melt(id = "ID") %>% mutate(variable = gsub("5mC", "", variable))
(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_5mC_figure <- ggplot(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_5mC, aes(x = variable, y = ID, fill = value)) + 
		geom_tile() + 
		scale_fill_gradient(name = " Fractional\nmethylation", low = "black", high = "darkred") + 
		xlab("") + 
		ylab("") + 
		theme_bw() + 
		theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent"), legend.position = "bottom"))
ggsave(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_5mC_figure, file = "./H3K27ac/IDHmut_enhancer_5mC_promoter_hyper_DE_DM_5mC_figure.pdf", height = 5, width = 3)
IDHmut_enhancer_5mC_promoter_hyper_DE_DM_K27ac <- IDHmut_enhancer_5mC_promoter_hyper_DE_DM %>% select(ID, mutK27ac, NPCK27ac) %>% melt(id = "ID") %>% mutate(variable = gsub("K27ac", "", variable))
(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_K27ac_figure <- ggplot(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_K27ac, aes(x = variable, y = ID, fill = value)) + 
		geom_tile() + 
		scale_fill_manual(name = " H3K27ac", values = c("T" = "darkgreen", "F" = "white"), labels = c("T" = "Marked", "F" = "Unmarked")) + 
		xlab("") + 
		ylab("") + 
		theme_bw() + 
		theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent"), legend.position = "bottom"))
ggsave(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_K27ac_figure, file = "./H3K27ac/IDHmut_enhancer_5mC_promoter_hyper_DE_DM_K27ac_figure.pdf", height = 5, width = 3)
IDHmut_enhancer_5mC_promoter_hyper_DE_DM_RPKM <- IDHmut_enhancer_5mC_promoter_hyper_DE_DM %>% select(ID, mutRPKM, NPCRPKM, wtRPKM) %>% melt(id = "ID") %>% mutate(variable = gsub("RPKM", "", variable))
(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_RPKM_figure <- ggplot(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_RPKM, aes(x = variable, y = ID, fill = log10(value))) + 
		geom_tile() + 
		scale_fill_gradient(name = " log10 RPKM") + 
		xlab("") + 
		ylab("") + 
		theme_bw() + 
		theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent"), legend.position = "bottom"))
ggsave(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_RPKM_figure, file = "./H3K27ac/IDHmut_enhancer_5mC_promoter_hyper_DE_DM_RPKM_figure.pdf", height = 5, width = 3)
IDHmut_enhancer_5mC_promoter_hyper_DE_DM_5mC <- IDHmut_enhancer_5mC_promoter_hyper_DE_DM %>% select(ID, mut5mC, wt5mC, NPC5mC) %>% melt(id = "ID") %>% mutate(variable = gsub("5mC", "", variable), type = "5mC")
(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_5mC_figure <- ggplot(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_5mC, aes(x = variable, y = value, fill = variable)) + 
		geom_violin() + 
		xlab("") + 
		ylab("Fractional methylation") + 
		facet_wrap(~ type) + 
		theme_bw() + 
		theme(axis.text.x = element_text(size = 12), legend.title = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_5mC_figure, file = "./H3K27ac/IDHmut_enhancer_5mC_promoter_hyper_DE_DM_5mC_figure.pdf", height = 5, width = 4)
IDHmut_enhancer_5mC_promoter_hyper_DE_DM_RPKM <- IDHmut_enhancer_5mC_promoter_hyper_DE_DM %>% select(ID, mutRPKM, NPCRPKM, wtRPKM) %>% melt(id = "ID") %>% mutate(variable = gsub("RPKM", "", variable), type = "RPKM")
(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_RPKM_figure <- ggplot(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_RPKM, aes(x = variable, y = log10(value), fill = variable)) + 
		geom_violin() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		facet_wrap(~ type) + 
		theme_bw() + 
		theme(axis.text.x = element_text(size = 12), legend.title = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(IDHmut_enhancer_5mC_promoter_hyper_DE_DM_RPKM_figure, file = "./H3K27ac/IDHmut_enhancer_5mC_promoter_hyper_DE_DM_RPKM_figure.pdf", height = 5, width = 4)

## ----------- save ------------
save(list = c(ls(pattern = "figure"), ls(pattern = "dend"), ls(pattern = "heatmap"), "pearson_genome_5mC", "CGI_5mC_variable", ls(pattern = "ann"), "DM_fitEb", ls(pattern = "_DM$")),
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/glioma/WGBS/WGBS.Rdata")

