# Glioma RNA-seq expression analysis

library(ggplot2)
library(dplyr)
library(reshape2)
library(limma)
library(pheatmap)
library(VennDiagram)
source('~/HirstLab/Pipeline/R/enrich.R')
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/")
load("/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/RNAseq.Rdata")

## ----- cluster -----------------------
RPKM <- read.delim("./RPKM/RPKM.matrix", row.names = 1, as.is = T)
spearman_RPKM <- cor(RPKM, method = "spearman")
spearman_RPKM_dend <-  hclust(1 - spearman_RPKM %>% as.dist, method = "ward.D2") %>% as.dendrogram
pdf("./RPKM_spearman.pdf")
op <- par(mar = c(5, 4, 4, 10))
plot(spearman_RPKM_dend, main = "pc genes spearman", horiz = TRUE)
par(op)
dev.off()
ann <- data.frame(category = gsub("\\..*", "", row.names(spearman_RPKM))) 
rownames(ann) <- row.names(spearman_RPKM)
ann_color <- list(category = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), 
															 hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100), 
															 hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100), 
															 hcl(h = seq(15, 375, length = 5 + 1)[5], l = 65, c = 100)))
names(ann_color$category) <- c("IDHmut", "IDHwt", "NormalAdjacent", "NPC")
pheatmap_spearman_RPKM <- pheatmap(spearman_RPKM, main = "pc genes spearman", filename = "./heatmap_spearman_RPKM.pdf", height = 7, width = 7, color = colorRampPalette(c("forest green ","white","purple"))(100), clustering_method = "ward.D2", cutree_cols = 3, cutree_rows = 3, annotation_row = ann, annotation_col = ann, annotation_colors = ann_color, show_rownames = F)
RPKM_variable <- RPKM %>% mutate(id = rownames(.), sd = apply(., 1, sd)) %>% filter(sd > quantile(sd, 0.8))
rownames(RPKM_variable) <- RPKM_variable$id
pheatmap_pearson_RPKM_variable <- pheatmap(RPKM_variable %>% select(-id, -sd), main = "20% variable pc genes", filename = "./heatmap_pearson_RPKM_variable.pdf", height = 8, width = 6, scale = "row", color = colorRampPalette(c("forest green ","white","purple"))(100), clustering_distance_cols = "correlation", clustering_method = "ward.D2", cutree_rows = 4, cutree_cols = 3, annotation_col = ann, annotation_colors = ann_color, show_rownames = F)

## ----- DE limma ----------------------
cutoff <- 0.01
fc <- 1
e <- 1e-6
RPKM_log10 <- as.data.frame(RPKM + e) %>% log10() %>% mutate(ENSG = row.names(.), exp = rowSums(. > -1)) %>% filter(exp > 0)
rownames(RPKM_log10) <- RPKM_log10$ENSG
RPKM_log10 <- RPKM_log10 %>% select(-ENSG, -exp)
DE_design <- data.frame(Group = relevel(factor(gsub("\\..*", "", colnames(RPKM_log10))), ref = "NPC"), row.names = colnames(RPKM_log10))
DE_DesMat <- model.matrix(~ 0 + Group, DE_design)
colnames(DE_DesMat) <- gsub("Group", "", colnames(DE_DesMat))
DE_fit <- lmFit(RPKM_log10, DE_DesMat)
contrast.matrix <- makeContrasts(IDHmut-NPC, IDHwt-NPC, NormalAdjacent-NPC, IDHmut-IDHwt, levels = DE_DesMat)
DE_fit <- contrasts.fit(DE_fit, contrast.matrix)
DE_fitEb <- eBayes(DE_fit)
IDHmut_NPC_DE <- topTable(DE_fitEb, coef = 'IDHmut - NPC', number = Inf, adjust.method="BH", p.value = cutoff, lfc = fc) %>% mutate(ENSG = row.names(.), DE = ifelse(logFC > 0, "UP", "DN")) %>% select(ENSG, DE, logFC:adj.P.Val)
IDHmut_NPC_UP <- IDHmut_NPC_DE %>% filter(DE == "UP")
IDHmut_NPC_DN <- IDHmut_NPC_DE %>% filter(DE == "DN")
write.table(IDHmut_NPC_DE, file = "./limma/DE.IDHmut_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHmut_NPC_UP, file = "./limma/UP.IDHmut_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHmut_NPC_DN, file = "./limma/DN.IDHmut_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
NormalAdjacent_NPC_DE <- topTable(DE_fitEb, coef = 'NormalAdjacent - NPC', number = Inf, adjust.method="BH", p.value = cutoff, lfc = fc) %>% mutate(ENSG = row.names(.), DE = ifelse(logFC > 0, "UP", "DN")) %>% select(ENSG, DE, logFC:adj.P.Val) 
NormalAdjacent_NPC_UP <- NormalAdjacent_NPC_DE %>% filter(DE == "UP")
NormalAdjacent_NPC_DN <- NormalAdjacent_NPC_DE %>% filter(DE == "DN")
write.table(NormalAdjacent_NPC_DE, file = "./limma/DE.NormalAdjacent_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(NormalAdjacent_NPC_UP, file = "./limma/UP.NormalAdjacent_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(NormalAdjacent_NPC_DN, file = "./limma/DN.NormalAdjacent_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
IDHwt_NPC_DE <- topTable(DE_fitEb, coef = 'IDHwt - NPC', number = Inf, adjust.method="BH", p.value = cutoff, lfc = fc) %>% mutate(ENSG = row.names(.), DE = ifelse(logFC > 0, "UP", "DN")) %>% select(ENSG, DE, logFC:adj.P.Val) 
IDHwt_NPC_UP <- IDHwt_NPC_DE %>% filter(DE == "UP")
IDHwt_NPC_DN <- IDHwt_NPC_DE %>% filter(DE == "DN")
write.table(IDHwt_NPC_DE, file = "./limma/DE.IDHwt_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHwt_NPC_UP, file = "./limma/UP.IDHwt_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHwt_NPC_DN, file = "./limma/DN.IDHwt_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
IDHmut_IDHwt_DE <- topTable(DE_fitEb, coef = 'IDHmut - IDHwt', number = Inf, adjust.method="BH", p.value = cutoff, lfc = fc) %>% mutate(ENSG = row.names(.), DE = ifelse(logFC > 0, "UP", "DN")) %>% select(ENSG, DE, logFC:adj.P.Val)
IDHmut_IDHwt_UP <- IDHmut_IDHwt_DE %>% filter(DE == "UP")
IDHmut_IDHwt_DN <- IDHmut_IDHwt_DE %>% filter(DE == "DN")
write.table(IDHmut_IDHwt_DE, file = "./limma/DE.IDHmut_IDHwt_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHmut_IDHwt_UP, file = "./limma/UP.IDHmut_IDHwt_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHmut_IDHwt_DN, file = "./limma/DN.IDHmut_IDHwt_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### ----- volcano plots ----------------------
IDHmut_NPC_limma <- topTable(DE_fitEb, coef = 'IDHmut - NPC', number = Inf, adjust.method="BH") %>% mutate(ENSG = row.names(.), DE = ifelse(ENSG %in% IDHmut_NPC_UP$ENSG, "UP", ifelse(ENSG %in% IDHmut_NPC_DN$ENSG, "DN", "ST"))) 
(IDHmut_NPC_volcano <- ggplot(IDHmut_NPC_limma, aes(logFC, -log10(adj.P.Val), color = DE)) + 
		geom_point(alpha = 0.2, size = 0.5) + 
		scale_colour_manual(values = c("UP" = "red", "DN" = "blue", "ST" = "grey"), guide = "none") + 
		geom_vline(xintercept = c(-1, 1), color = "black", linetype = 2) + 
		geom_hline(yintercept = 2, color = "black", linetype = 2) + 
		annotate("text", x = -5, y = 20, label = paste0("DN=", nrow(IDHmut_NPC_DN)), color = "blue", size = 6) + 
		annotate("text", x = 5, y = 20, label = paste0("UP=", nrow(IDHmut_NPC_UP)), color = "red", size = 6) + 
		labs(title = "IDHmut vs NPC", x = "log2 FC", y = "-log10 FDR") + 
		theme_bw())
ggsave(IDHmut_NPC_volcano, file = "./limma/IDHmut_NPC_volcano.pdf", height = 5, width = 5)
NormalAdjacent_NPC_limma <- topTable(DE_fitEb, coef = 'NormalAdjacent - NPC', number = Inf, adjust.method="BH") %>% mutate(ENSG = row.names(.), DE = ifelse(ENSG %in% NormalAdjacent_NPC_UP$ENSG, "UP", ifelse(ENSG %in% NormalAdjacent_NPC_DN$ENSG, "DN", "ST"))) 
(NormalAdjacent_NPC_volcano <- ggplot(NormalAdjacent_NPC_limma, aes(logFC, -log10(adj.P.Val), color = DE)) + 
		geom_point(alpha = 0.2, size = 0.5) + 
		scale_colour_manual(values = c("UP" = "red", "DN" = "blue", "ST" = "grey"), guide = "none") + 
		geom_vline(xintercept = c(-1, 1), color = "black", linetype = 2) + 
		geom_hline(yintercept = 2, color = "black", linetype = 2) + 
		annotate("text", x = -5, y = 20, label = paste0("DN=", nrow(NormalAdjacent_NPC_DN)), color = "blue", size = 6) + 
		annotate("text", x = 5, y = 20, label = paste0("UP=", nrow(NormalAdjacent_NPC_UP)), color = "red", size = 6) + 
		labs(title = "NormalAdjacent vs NPC", x = "log2 FC", y = "-log10 FDR") + 
		theme_bw())
ggsave(NormalAdjacent_NPC_volcano, file = "./limma/NormalAdjacent_NPC_volcano.pdf", height = 5, width = 5)
IDHwt_NPC_limma <- topTable(DE_fitEb, coef = 'IDHwt - NPC', number = Inf, adjust.method="BH") %>% mutate(ENSG = row.names(.), DE = ifelse(ENSG %in% IDHwt_NPC_UP$ENSG, "UP", ifelse(ENSG %in% IDHwt_NPC_DN$ENSG, "DN", "ST"))) 
(IDHwt_NPC_volcano <- ggplot(IDHwt_NPC_limma, aes(logFC, -log10(adj.P.Val), color = DE)) + 
		geom_point(alpha = 0.2, size = 0.5) + 
		scale_colour_manual(values = c("UP" = "red", "DN" = "blue", "ST" = "grey"), guide = "none") + 
		geom_vline(xintercept = c(-1, 1), color = "black", linetype = 2) + 
		geom_hline(yintercept = 2, color = "black", linetype = 2) + 
		annotate("text", x = -5, y = 20, label = paste0("DN=", nrow(IDHwt_NPC_DN)), color = "blue", size = 6) + 
		annotate("text", x = 5, y = 20, label = paste0("UP=", nrow(IDHwt_NPC_UP)), color = "red", size = 6) + 
		labs(title = "IDHwt vs NPC", x = "log2 FC", y = "-log10 FDR") + 
		theme_bw())
ggsave(IDHwt_NPC_volcano, file = "./limma/IDHwt_NPC_volcano.pdf", height = 5, width = 5)
IDHmut_IDHwt_limma <- topTable(DE_fitEb, coef = 'IDHmut - IDHwt', number = Inf, adjust.method="BH") %>% mutate(ENSG = row.names(.), DE = ifelse(ENSG %in% IDHmut_IDHwt_UP$ENSG, "UP", ifelse(ENSG %in% IDHmut_IDHwt_DN$ENSG, "DN", "ST"))) 
(IDHmut_IDHwt_volcano <- ggplot(IDHmut_IDHwt_limma, aes(logFC, -log10(adj.P.Val), color = DE)) + 
		geom_point(alpha = 0.2, size = 0.5) + 
		scale_colour_manual(values = c("UP" = "red", "DN" = "blue", "ST" = "grey"), guide = "none") + 
		geom_vline(xintercept = c(-1, 1), color = "black", linetype = 2) + 
		geom_hline(yintercept = 2, color = "black", linetype = 2) + 
		annotate("text", x = -4, y = 11, label = paste0("DN=", nrow(IDHmut_IDHwt_DN)), color = "blue", size = 6) + 
		annotate("text", x = 4, y = 11, label = paste0("UP=", nrow(IDHmut_IDHwt_UP)), color = "red", size = 6) + 
		labs(title = "IDHmut vs IDHwt", x = "log2 FC", y = "-log10 FDR") + 
		theme_bw())
ggsave(IDHmut_IDHwt_volcano, file = "./limma/IDHmut_IDHwt_volcano.pdf", height = 5, width = 5)

### ----- intersect venn ----------------------
limma_venn_UP <- venn.diagram(list(IDHmut = IDHmut_NPC_UP$ENSG, NormalAdjacent = NormalAdjacent_NPC_UP$ENSG, IDHwt = IDHwt_NPC_UP$ENSG), filename = NULL, main = "UP compared to NPC", cat.pos = c(-30, 30, 180), fill = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100)))
pdf("./limma/limma_venn_UP.pdf", height = 4, width = 4)
grid.draw(limma_venn_UP)
dev.off()
limma_venn_DN <- venn.diagram(list(IDHmut = IDHmut_NPC_DN$ENSG, NormalAdjacent = NormalAdjacent_NPC_DN$ENSG, IDHwt = IDHwt_NPC_DN$ENSG), filename = NULL, main = "DN compared to NPC", cat.pos = c(-30, 30, 180), fill = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100)))
pdf("./limma/limma_venn_DN.pdf", height = 4, width = 4)
grid.draw(limma_venn_DN)
dev.off()
IDHmut_NPC_UP_specific <- IDHmut_NPC_UP %>% filter(!(ENSG %in% IDHwt_NPC_UP$ENSG)) %>% filter(!(ENSG %in% NormalAdjacent_NPC_UP$ENSG))
write.table(IDHmut_NPC_UP_specific, file = "./limma/UP.IDHmut_NPC_limma_specific.txt", sep = "\t", quote = F, row.names = F, col.names = T)
IDHmut_NPC_DN_specific <- IDHmut_NPC_DN %>% filter(!(ENSG %in% IDHwt_NPC_DN$ENSG)) %>% filter(!(ENSG %in% NormalAdjacent_NPC_UP$ENSG))
write.table(IDHmut_NPC_DN_specific, file = "./limma/DN.IDHmut_NPC_limma_specific.txt", sep = "\t", quote = F, row.names = F, col.names = T)
UP_common <- IDHmut_NPC_UP %>% filter(ENSG %in% IDHwt_NPC_UP$ENSG) %>% filter(ENSG %in% NormalAdjacent_NPC_UP$ENSG)
write.table(UP_common, file = "./limma/UP.limma_NPC_common.txt", sep = "\t", quote = F, row.names = F, col.names = T)
DN_common <- IDHmut_NPC_DN %>% filter(ENSG %in% IDHwt_NPC_DN$ENSG) %>% filter(ENSG %in% NormalAdjacent_NPC_DN$ENSG)
write.table(DN_common, file = "./limma/DN.limma_NPC_common.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### ----- DAVID ----------------------
(NPC_common_UP_DAVID <- enrich(name = "UP.NPC_common", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 1e-10, p = "FDR", erminej = F, height = 7, width = 6))
(NPC_common_DN_DAVID <- enrich(name = "DN.NPC_common", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 1e-2, p = "FDR", erminej = F, height = 7, width = 6))
(IDHmut_IDHwt_UP_DAVID <- enrich(name = "UP.IDHmut_IDHwt", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 1e-2, p = "FDR", erminej = F, height = 3, width = 6))
(IDHmut_IDHwt_DN_DAVID <- enrich(name = "DN.IDHmut_IDHwt", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 1e-2, p = "FDR", erminej = F, height = 6, width = 6))
(IDHmut_IDHwt_UP_TCGA_LGG_DAVID <- enrich(name = "UP.IDHmut_IDHwt.TCGA_LGG", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 1e-2, p = "FDR", erminej = F, height = 2, width = 6))
(IDHmut_IDHwt_DN_TCGA_LGG_DAVID <- enrich(name = "DN.IDHmut_IDHwt.TCGA_LGG", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 1e-2, p = "FDR", erminej = F, height = 6, width = 7))

## ----- TCGA -------------------------
load("./tcga_rna.rdata")
ann_color_tcga <- list(IDH = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100)), TET = c("black", "lightgrey"))
names(ann_color_tcga$IDH) <- c("mut", "wt")
names(ann_color_tcga$TET) <- c("mut", "wt")
ann_lgg <- data.frame(IDH = ifelse(data.lgg.label$IDH == 1, "mut", "wt")) 
rownames(ann_lgg) <- data.lgg.label$Sample.ID
spearman_RPKM_lgg <- cor(data.rna.lgg, method = "spearman")
heatmap_spearman_RPKM_lgg <- pheatmap(spearman_RPKM_lgg, main = "TCGA LGG pc genes spearman", filename = "./heatmap_spearman_RPKM_lgg.pdf", height = 5, width = 5, color = colorRampPalette(c("forest green ","white","purple"))(100), clustering_method = "ward.D2", cutree_cols = 2, cutree_rows = 2, annotation_row = ann_lgg, annotation_col = ann_lgg, annotation_colors = ann_color_tcga, show_rownames = F, show_colnames = F)
ann_aml <- data.frame(IDH = ifelse(data.aml.label$IDH == 1, "mut", "wt"), TET = ifelse(data.aml.label$TET == 1, "mut", "wt")) 
rownames(ann_aml) <- data.aml.label$Sample.ID
spearman_RPKM_aml <- cor(data.rna.aml, method = "spearman")
heatmap_spearman_RPKM_aml <- pheatmap(spearman_RPKM_aml, main = "TCGA AML pc genes spearman", filename = "./heatmap_spearman_RPKM_aml.pdf", height = 5, width = 5, color = colorRampPalette(c("forest green ","white","purple"))(100), clustering_method = "ward.D2", annotation_row = ann_aml, annotation_col = ann_aml, annotation_colors = ann_color_tcga, show_rownames = F, show_colnames = F)

cutoff <- 0.01
fc <- 1
e <- 1e-6
RPKM_log10_lgg <- as.data.frame(data.rna.lgg + e) %>% log10() %>% mutate(Entrez = row.names(.), exp = rowSums(. > -1)) %>% filter(exp > 0)
rownames(RPKM_log10_lgg) <- RPKM_log10_lgg$Entrez
RPKM_log10_lgg <- RPKM_log10_lgg %>% select(-Entrez, -exp)
DE_design_lgg <- data.frame(Group = relevel(factor(paste0("IDH", ann_lgg$IDH)), ref = "IDHwt"), row.names = rownames(ann_lgg)) 
DE_DesMat_lgg <- model.matrix(~ 0 + Group, DE_design_lgg)
colnames(DE_DesMat_lgg) <- gsub("Group", "", colnames(DE_DesMat_lgg))
DE_fit_lgg <- lmFit(RPKM_log10_lgg, DE_DesMat_lgg)
contrast.matrix_lgg <- makeContrasts(IDHmut-IDHwt, levels = DE_DesMat_lgg)
DE_fit_lgg <- contrasts.fit(DE_fit_lgg, contrast.matrix_lgg)
DE_fitEb_lgg <- eBayes(DE_fit_lgg)
IDHmut_IDHwt_DE_lgg <- topTable(DE_fitEb_lgg, coef = 'IDHmut - IDHwt', number = Inf, adjust.method="BH", p.value = cutoff, lfc = fc) %>% mutate(ENSG = row.names(.), DE = ifelse(logFC > 0, "UP", "DN")) %>% select(ENSG, DE, logFC:adj.P.Val)
IDHmut_IDHwt_UP_lgg <- IDHmut_IDHwt_DE_lgg %>% filter(DE == "UP")
IDHmut_IDHwt_DN_lgg <- IDHmut_IDHwt_DE_lgg %>% filter(DE == "DN")
write.table(IDHmut_IDHwt_DE_lgg, file = "./limma/DE.IDHmut_IDHwt_limma_lgg.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHmut_IDHwt_UP_lgg, file = "./limma/UP.IDHmut_IDHwt_limma_lgg.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHmut_IDHwt_DN_lgg, file = "./limma/DN.IDHmut_IDHwt_limma_lgg.txt", sep = "\t", quote = F, row.names = F, col.names = T)
IDHmut_IDHwt_limma_lgg <- topTable(DE_fitEb_lgg, coef = 'IDHmut - IDHwt', number = Inf, adjust.method="BH") %>% mutate(ENSG = row.names(.), DE = ifelse(ENSG %in% IDHmut_IDHwt_UP_lgg$ENSG, "UP", ifelse(ENSG %in% IDHmut_IDHwt_DN_lgg$ENSG, "DN", "ST"))) 
(IDHmut_IDHwt_volcano_lgg <- ggplot(IDHmut_IDHwt_limma_lgg, aes(logFC, -log10(adj.P.Val), color = DE)) + 
		geom_point(alpha = 0.2, size = 0.5) + 
		scale_colour_manual(values = c("UP" = "red", "DN" = "blue", "ST" = "grey"), guide = "none") + 
		geom_vline(xintercept = c(-1, 1), color = "black", linetype = 2) + 
		geom_hline(yintercept = 2, color = "black", linetype = 2) + 
		annotate("text", x = -2.5, y = 70, label = paste0("DN=", nrow(IDHmut_IDHwt_DN_lgg)), color = "blue", size = 6) + 
		annotate("text", x = 2.8, y = 70, label = paste0("UP=", nrow(IDHmut_IDHwt_UP_lgg)), color = "red", size = 6) + 
		labs(title = "TCGA LGG IDHmut vs IDHwt", x = "log2 FC", y = "-log10 FDR") + 
		theme_bw())
ggsave(IDHmut_IDHwt_volcano_lgg, file = "./limma/IDHmut_IDHwt_volcano_lgg.pdf", height = 5, width = 5)
lgg_venn_UP <- venn.diagram(list(InHouse = IDHmut_IDHwt_UP$ENSG, TCGA = IDHmut_IDHwt_UP_lgg$ENSG), filename = NULL, main = "compared to TCGA UP", cat.pos = c(-30, 30), fill = c("red", "blue"))
pdf("./limma/lgg_venn_UP.pdf", height = 4, width = 4)
grid.draw(lgg_venn_UP)
dev.off()
lgg_venn_DN <- venn.diagram(list(InHouse = IDHmut_IDHwt_DN$ENSG, TCGA = IDHmut_IDHwt_DN_lgg$ENSG), filename = NULL, main = "compared to TCGA DN", cat.pos = c(-30, 30), fill = c("red", "blue"))
pdf("./limma/lgg_venn_DN.pdf", height = 4, width = 4)
grid.draw(lgg_venn_DN)
dev.off()

RPKM_log10_aml <- as.data.frame(data.rna.aml + e) %>% log10() %>% mutate(Entrez = row.names(.), exp = rowSums(. > -1)) %>% filter(exp > 0)
rownames(RPKM_log10_aml) <- RPKM_log10_aml$Entrez
RPKM_log10_aml <- RPKM_log10_aml %>% select(-Entrez, -exp)
DE_design_aml <- data.frame(Group = relevel(factor(paste0("IDH", ann_aml$IDH)), ref = "IDHwt"), row.names = rownames(ann_aml)) 
DE_DesMat_aml <- model.matrix(~ 0 + Group, DE_design_aml)
colnames(DE_DesMat_aml) <- gsub("Group", "", colnames(DE_DesMat_aml))
DE_fit_aml <- lmFit(RPKM_log10_aml, DE_DesMat_aml)
contrast.matrix_aml <- makeContrasts(IDHmut-IDHwt, levels = DE_DesMat_aml)
DE_fit_aml <- contrasts.fit(DE_fit_aml, contrast.matrix_aml)
DE_fitEb_aml <- eBayes(DE_fit_aml)
IDHmut_IDHwt_DE_aml <- topTable(DE_fitEb_aml, coef = 'IDHmut - IDHwt', number = Inf, adjust.method="BH", p.value = cutoff, lfc = fc) %>% mutate(ENSG = row.names(.), DE = ifelse(logFC > 0, "UP", "DN")) %>% select(ENSG, DE, logFC:adj.P.Val)
IDHmut_IDHwt_UP_aml <- IDHmut_IDHwt_DE_aml %>% filter(DE == "UP")
IDHmut_IDHwt_DN_aml <- IDHmut_IDHwt_DE_aml %>% filter(DE == "DN")
write.table(IDHmut_IDHwt_DE_aml, file = "./limma/DE.IDHmut_IDHwt_limma_aml.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHmut_IDHwt_UP_aml, file = "./limma/UP.IDHmut_IDHwt_limma_aml.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHmut_IDHwt_DN_aml, file = "./limma/DN.IDHmut_IDHwt_limma_aml.txt", sep = "\t", quote = F, row.names = F, col.names = T)
IDHmut_IDHwt_limma_aml <- topTable(DE_fitEb_aml, coef = 'IDHmut - IDHwt', number = Inf, adjust.method="BH") %>% mutate(ENSG = row.names(.), DE = ifelse(ENSG %in% IDHmut_IDHwt_UP_aml$ENSG, "UP", ifelse(ENSG %in% IDHmut_IDHwt_DN_aml$ENSG, "DN", "ST"))) 
(IDHmut_IDHwt_volcano_aml <- ggplot(IDHmut_IDHwt_limma_aml, aes(logFC, -log10(adj.P.Val), color = DE)) + 
		geom_point(alpha = 0.2, size = 0.5) + 
		scale_colour_manual(values = c("UP" = "red", "DN" = "blue", "ST" = "grey"), guide = "none") + 
		geom_vline(xintercept = c(-1, 1), color = "black", linetype = 2) + 
		geom_hline(yintercept = 2, color = "black", linetype = 2) + 
		annotate("text", x = -2, y = 20, label = paste0("DN=", nrow(IDHmut_IDHwt_DN_aml)), color = "blue", size = 6) + 
		annotate("text", x = 2, y = 20, label = paste0("UP=", nrow(IDHmut_IDHwt_UP_aml)), color = "red", size = 6) + 
		labs(title = "TCGA AML IDHmut vs IDHwt", x = "log2 FC", y = "-log10 FDR") + 
		theme_bw())
ggsave(IDHmut_IDHwt_volcano_aml, file = "./limma/IDHmut_IDHwt_volcano_aml.pdf", height = 5, width = 5)

## ----- IDH1/2 expression levels ------
IDH_glioma <- read.delim("./RPKM/IDH.RPKM", head = F, col.names = c("ID", "Gene", "Sample", "RPKM"))
IDH_NPC <- read.delim("/projects/epigenomics/users/lli/FetalBrain/Tables/rpkm_pc.txt") %>% 
	filter(Ensembl %in% IDH_glioma$ID)
colnames(IDH_NPC) <- gsub("\\.HuFNSC", "", colnames(IDH_NPC))
IDH_NPC <- melt(IDH_NPC, id = "Ensembl") %>% mutate(Gene = ifelse(Ensembl == "ENSG00000138413", "IDH1", "IDH2"))
colnames(IDH_NPC) <- c("ID", "Sample", "RPKM", "Gene")
IDH_RPKM <- rbind(IDH_glioma, IDH_NPC)
(IDH_RPKM_figure <- ggplot(IDH_RPKM, aes(Sample, RPKM, color = Gene)) + 
	geom_point(size = 5) + 
	geom_line(aes(group = Gene)) + 
	scale_color_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw())
ggsave(IDH_RPKM_figure, file = "./RPKM/IDH_RPKM_figure.pdf", height = 5, width = 6)

## ----- DE glioma vs NPC ------
DE_pc_summary <- read.delim("./DEfine/DE.pc.summary", as.is = T) %>%
	melt(id = c("Sample",	"DE"), variable.name = "NPC") %>% mutate(value = ifelse(DE == "DN", -value, value))
(DE_pc_summary_figure <- ggplot(DE_pc_summary, aes(Sample, value, fill = NPC)) + 
	geom_bar(stat = "identity", position = position_dodge(), width = 0.8) + 
	geom_hline(yintercept = 0) + 
	xlab("") + 
	ylab("No. of DE genes") + 
	theme_bw())
ggsave(DE_pc_summary_figure, file = "./DEfine/DE_pc_summary_figure.pdf", height = 6, width = 6)
for(lib in libs){
	assign(paste0("DAVID_UP_", lib, "_figure"), enrich(paste0("UP.", lib, "_NPC"), dirIn = "./DEfine/enrich/", dirOut = "./DEfine/enrich/", erminej = F))
	assign(paste0("DAVID_DN_", lib, "_figure"), enrich(paste0("DN.", lib, "_NPC"), dirIn = "./DEfine/enrich/", dirOut = "./DEfine/enrich/", erminej = F))
}

save(list = c("RPKM", "RPKM_variable", "spearman_RPKM", "spearman_RPKM_dend", ls(pattern = "heatmap"), ls(pattern = "ann"), ls(pattern = "DE_fitEb"), 
							ls(pattern = "DE"), ls(pattern = "UP"), ls(pattern = "DN"), ls(pattern = "limma"), ls(pattern = "volcano"), ls(pattern = "venn"), ls(pattern = "DAVID")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/RNAseq.Rdata")
