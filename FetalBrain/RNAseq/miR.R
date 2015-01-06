# miRNA analysis on FetalBrain

# =========== data input ============
library(ggplot2)
library(dplyr)
library(reshape)
library(dendextend)
library(VennDiagram)
library(gridExtra)
library(gplots)
source('~/HirstLab/Pipeline/R/DE_FC.R')
setwd("/projects/epigenomics/users/lli/FetalBrain/miR/")
miR <- read.delim("/projects/epigenomics/ep50/miRs/expn_matrix_norm.txt", as.is = T)
miR_FetalBrain <- data.frame(Gene = miR$Gene,
                             Brain01 = miR$MX0319_ATTGGC, 
                             Brain02 = miR$MX0319_GATCTG, 
                             Cortex01 = miR$MX0318_ACATCG, 
                             Cortex02 = miR$MX0318_TGGTCA, 
                             Cortex03 = miR$MX0320_TGGTCA, 
                             Cortex04 = miR$MX0673_GGGGTT, 
                             GE01 = miR$MX0318_GCCTAA, 
                             GE02 = miR$MX0319_CACTGT, 
                             GE03 = miR$MX0673_CGTACG, 
                             GE04 = miR$MX0673_CAAGTT)
row.names(miR_FetalBrain) <- miR_FetalBrain$Gene
write.table(miR_FetalBrain, file = "miR_FetalBrain.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# ========== summary ==========
miR_summary <- data.frame(Sample = colnames(miR_FetalBrain)[-1], 
                          N_0.1 = c(filter(miR_FetalBrain, Brain01 > 0.1) %>% nrow, filter(miR_FetalBrain, Brain02 > 0.1) %>% nrow, filter(miR_FetalBrain, Cortex01 > 0.1) %>% nrow, filter(miR_FetalBrain, Cortex02 > 0.1) %>% nrow, filter(miR_FetalBrain, Cortex03 > 0.1) %>% nrow, filter(miR_FetalBrain, Cortex04 > 0.1) %>% nrow, filter(miR_FetalBrain, GE01 > 0.1) %>% nrow, filter(miR_FetalBrain, GE02 > 0.1) %>% nrow, filter(miR_FetalBrain, GE03 > 0.1) %>% nrow, filter(miR_FetalBrain, GE04 > 0.1) %>% nrow), 
                          N_1 = c(filter(miR_FetalBrain, Brain01 > 1) %>% nrow, filter(miR_FetalBrain, Brain02 > 1) %>% nrow, filter(miR_FetalBrain, Cortex01 > 1) %>% nrow, filter(miR_FetalBrain, Cortex02 > 1) %>% nrow, filter(miR_FetalBrain, Cortex03 > 1) %>% nrow, filter(miR_FetalBrain, Cortex04 > 1) %>% nrow, filter(miR_FetalBrain, GE01 > 1) %>% nrow, filter(miR_FetalBrain, GE02 > 1) %>% nrow, filter(miR_FetalBrain, GE03 > 1) %>% nrow, filter(miR_FetalBrain, GE04 > 1) %>% nrow), 
                          N_10 = c(filter(miR_FetalBrain, Brain01 > 10) %>% nrow, filter(miR_FetalBrain, Brain02 > 10) %>% nrow, filter(miR_FetalBrain, Cortex01 > 10) %>% nrow, filter(miR_FetalBrain, Cortex02 > 10) %>% nrow, filter(miR_FetalBrain, Cortex03 > 10) %>% nrow, filter(miR_FetalBrain, Cortex04 > 10) %>% nrow, filter(miR_FetalBrain, GE01 > 10) %>% nrow, filter(miR_FetalBrain, GE02 > 10) %>% nrow, filter(miR_FetalBrain, GE03 > 10) %>% nrow, filter(miR_FetalBrain, GE04 > 10) %>% nrow), 
                          N_100 = c(filter(miR_FetalBrain, Brain01 > 100) %>% nrow, filter(miR_FetalBrain, Brain02 > 100) %>% nrow, filter(miR_FetalBrain, Cortex01 > 100) %>% nrow, filter(miR_FetalBrain, Cortex02 > 100) %>% nrow, filter(miR_FetalBrain, Cortex03 > 100) %>% nrow, filter(miR_FetalBrain, Cortex04 > 100) %>% nrow, filter(miR_FetalBrain, GE01 > 100) %>% nrow, filter(miR_FetalBrain, GE02 > 100) %>% nrow, filter(miR_FetalBrain, GE03 > 100) %>% nrow, filter(miR_FetalBrain, GE04 > 100) %>% nrow))
rownames(miR_summary) <- colnames(miR_FetalBrain)[-1]
miR_summary_tall <- melt(miR_summary, id = "Sample")
(miR_summary_figure <- ggplot(miR_summary_tall, aes(x = Sample, y = value, fill = variable)) + 
   geom_bar(stat = "identity", position = "dodge") + 
   xlab("") + 
   ylab("No. of genes") + 
   scale_fill_discrete(labels=c("N > 0.1","N > 1","N > 10", "N > 100"), name = "") + 
   theme_bw())
ggsave(miR_summary_figure, file = "/projects/epigenomics/users/lli/FetalBrain/miR/miR_summary_figure.pdf")

# =========== cluster ===========
highExpr <- as.character(miR_FetalBrain[rowSums(miR_FetalBrain[, 2:11] > 100) > 0, "Gene"])
miR_dend <- (1- (miR_FetalBrain[highExpr, -1] %>% as.matrix %>% cor(method = "spearman"))) %>% as.dist %>% hclust %>% as.dendrogram %>% 
  set("by_labels_branches_col", value = c("Brain01", "Brain02"), TF_value = "green", type = "all") %>% 
  set("by_labels_branches_col", value = c("Cortex01", "Cortex02", "Cortex03", "Cortex04"), TF_value = "red", type = "all") %>% 
  set("by_labels_branches_col", value = c("GE01", "GE02", "GE03", "GE04"), TF_value = "blue", type = "all") %>% 
  set("by_labels_branches_lwd", value = colnames(miR_FetalBrain[, -1]), TF_value = 2)
pdf("/projects/epigenomics/users/lli/FetalBrain/miR/miR_dend.pdf")
op <- par(mar = c(5, 4, 4, 4))
plot(miR_dend, main = "miRNA", horiz = TRUE)
legend("topleft", c("Brain", "Cortex", "GE"), col = c("green", "red", "blue"), lwd = 8)
par(op)
dev.off()

# =========== DE ==============
## MZ twins 
Brain01_Brain02_miR_DE <- DE_FC(miR_FetalBrain[, c("Brain01", "Brain02")], log_cut = 2)
Cortex01_Cortex02_miR_DE <- DE_FC(miR_FetalBrain[, c("Cortex01", "Cortex02")], log_cut = 2)
GE01_GE02_miR_DE <- DE_FC(miR_FetalBrain[, c("GE01", "GE02")], log_cut = 2)
miR_DE_MZ_summary <- rbind(Brain01_Brain02_miR_DE$summary, Cortex01_Cortex02_miR_DE$summary, GE01_GE02_miR_DE$summary)
miR_DE_MZ_UP <- list(Brain = Brain01_Brain02_miR_DE$UP$geneID, Cortex = Cortex01_Cortex02_miR_DE$UP$geneID, GE = GE01_GE02_miR_DE$UP$geneID)
venn_miR_DE_MZ_UP <- venn.diagram(miR_DE_MZ_UP, filename = NULL, fill = c("green", "red", "blue"), main = "Venn diagram of MZ UP miRNAs")
miR_DE_MZ_DN <- list(Brain = Brain01_Brain02_miR_DE$DN$geneID, Cortex = Cortex01_Cortex02_miR_DE$DN$geneID, GE = GE01_GE02_miR_DE$DN$geneID)
venn_miR_DE_MZ_DN <- venn.diagram(miR_DE_MZ_DN, filename = NULL, fill = c("green", "red", "blue"), main = "Venn diagram of MZ DN miRNAs")
pdf("/projects/epigenomics/users/lli/FetalBrain/miR/DE/venn_miR_DE_MZ.pdf", height = 4, width = 8)
grid.arrange(gTree(children = venn_miR_DE_MZ_UP), gTree(children = venn_miR_DE_MZ_DN), nrow = 1)
dev.off()
miR_DE_MZ_heat <- miR_FetalBrain[unique(c(Brain01_Brain02_miR_DE$UP$geneID, Cortex01_Cortex02_miR_DE$UP$geneID, GE01_GE02_miR_DE$UP$geneID, Brain01_Brain02_miR_DE$DN$geneID, Cortex01_Cortex02_miR_DE$DN$geneID, GE01_GE02_miR_DE$DN$geneID)), c("Brain01", "Brain02", "Cortex01", "Cortex02", "GE01", "GE02")]
pdf("/projects/epigenomics/users/lli/FetalBrain/miR/DE/heatmap_miR_DE_MZ.pdf", height = 20)
heatmap.2(as.matrix(miR_DE_MZ_heat), scale = "row", trace = "none", margins = c(11, 6), keysize = 1, density.info = "none", col = bluered(256), key.title = "", key.xlab = "")
dev.off()
## Cortex vs GE
Cortex01_GE01_miR_DE <- DE_FC(miR_FetalBrain[, c("Cortex01", "GE01")], log_cut = 2)
Cortex02_GE02_miR_DE <- DE_FC(miR_FetalBrain[, c("Cortex02", "GE02")], log_cut = 2)
Cortex03_GE03_miR_DE <- DE_FC(miR_FetalBrain[, c("Cortex03", "GE03")], log_cut = 2)
Cortex04_GE04_miR_DE <- DE_FC(miR_FetalBrain[, c("Cortex04", "GE04")], log_cut = 2)
miR_DE_neurospheres_summary <- rbind(Cortex01_GE01_miR_DE$summary, Cortex02_GE02_miR_DE$summary, Cortex03_GE03_miR_DE$summary, Cortex04_GE04_miR_DE$summary)
miR_DE_neurospheres_UP <- list(HuFNSC01 = Cortex01_GE01_miR_DE$UP$geneID, HuFNSC02 = Cortex02_GE02_miR_DE$UP$geneID, HuFNSC03 = Cortex03_GE03_miR_DE$UP$geneID, HuFNSC04 = Cortex04_GE04_miR_DE$UP$geneID)
venn_miR_DE_neurospheres_UP <- venn.diagram(miR_DE_neurospheres_UP, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of neurospheres UP miRNAs")
miR_DE_neurospheres_DN <- list(HuFNSC01 = Cortex01_GE01_miR_DE$DN$geneID, HuFNSC02 = Cortex02_GE02_miR_DE$DN$geneID, HuFNSC03 = Cortex03_GE03_miR_DE$DN$geneID, HuFNSC04 = Cortex04_GE04_miR_DE$DN$geneID)
venn_miR_DE_neurospheres_DN <- venn.diagram(miR_DE_neurospheres_DN, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of neurospheres DN miRNAs")
pdf("/projects/epigenomics/users/lli/FetalBrain/miR/DE/venn_miR_DE_neurospheres.pdf", height = 4, width = 8)
grid.arrange(gTree(children = venn_miR_DE_neurospheres_UP), gTree(children = venn_miR_DE_neurospheres_DN), nrow = 1)
dev.off()
miR_DE_neurospheres_heat <- miR_FetalBrain[unique(c(Cortex01_GE01_miR_DE$UP$geneID, Cortex02_GE02_miR_DE$UP$geneID, Cortex03_GE03_miR_DE$UP$geneID, Cortex04_GE04_miR_DE$UP$geneID, Cortex01_GE01_miR_DE$DN$geneID, Cortex02_GE02_miR_DE$DN$geneID, Cortex03_GE03_miR_DE$DN$geneID, Cortex04_GE04_miR_DE$DN$geneID)), c("Cortex01", "Cortex02", "Cortex03", "Cortex04", "GE01", "GE02", "GE03", "GE04")]
pdf("/projects/epigenomics/users/lli/FetalBrain/miR/DE/heatmap_miR_DE_neurospheres.pdf", height = 20)
heatmap.2(as.matrix(miR_DE_neurospheres_heat), scale = "row", trace = "none", margins = c(11, 6), keysize = 1, density.info = "none", col = bluered(256), key.title = "", key.xlab = "")
dev.off()
## GW
Cortex01_Cortex03_miR_DE <- DE_FC(miR_FetalBrain[, c("Cortex01", "Cortex03")], log_cut = 2)
Cortex01_Cortex04_miR_DE <- DE_FC(miR_FetalBrain[, c("Cortex01", "Cortex04")], log_cut = 2)
Cortex02_Cortex03_miR_DE <- DE_FC(miR_FetalBrain[, c("Cortex02", "Cortex03")], log_cut = 2)
Cortex02_Cortex04_miR_DE <- DE_FC(miR_FetalBrain[, c("Cortex02", "Cortex04")], log_cut = 2)
Cortex03_Cortex04_miR_DE <- DE_FC(miR_FetalBrain[, c("Cortex03", "Cortex04")], log_cut = 2)
miR_DE_GW_Cortex_UP <- list(HuFNSC01_03 = Cortex01_Cortex03_miR_DE$UP$geneID, HuFNSC01_04 = Cortex01_Cortex04_miR_DE$UP$geneID, HuFNSC01_03 = Cortex02_Cortex03_miR_DE$UP$geneID, HuFNSC02_04 = Cortex02_Cortex04_miR_DE$UP$geneID)
venn_miR_DE_GW_Cortex_UP <- venn.diagram(miR_DE_GW_Cortex_UP, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of GW_Cortex UP miRNAs")
miR_DE_GW_Cortex_DN <- list(HuFNSC01_03 = Cortex01_Cortex03_miR_DE$DN$geneID, HuFNSC01_04 = Cortex01_Cortex04_miR_DE$DN$geneID, HuFNSC01_03 = Cortex02_Cortex03_miR_DE$DN$geneID, HuFNSC02_04 = Cortex02_Cortex04_miR_DE$DN$geneID)
venn_miR_DE_GW_Cortex_DN <- venn.diagram(miR_DE_GW_Cortex_DN, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of GW_Cortex DN miRNAs")
GE01_GE03_miR_DE <- DE_FC(miR_FetalBrain[, c("GE01", "GE03")], log_cut = 2)
GE01_GE04_miR_DE <- DE_FC(miR_FetalBrain[, c("GE01", "GE04")], log_cut = 2)
GE02_GE03_miR_DE <- DE_FC(miR_FetalBrain[, c("GE02", "GE03")], log_cut = 2)
GE02_GE04_miR_DE <- DE_FC(miR_FetalBrain[, c("GE02", "GE04")], log_cut = 2)
GE03_GE04_miR_DE <- DE_FC(miR_FetalBrain[, c("GE03", "GE04")], log_cut = 2)
miR_DE_GW_GE_UP <- list(HuFNSC01_03 = GE01_GE03_miR_DE$UP$geneID, HuFNSC01_04 = GE01_GE04_miR_DE$UP$geneID, HuFNSC01_03 = GE02_GE03_miR_DE$UP$geneID, HuFNSC02_04 = GE02_GE04_miR_DE$UP$geneID)
venn_miR_DE_GW_GE_UP <- venn.diagram(miR_DE_GW_GE_UP, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of GW_GE UP miRNAs")
miR_DE_GW_GE_DN <- list(HuFNSC01_03 = GE01_GE03_miR_DE$DN$geneID, HuFNSC01_04 = GE01_GE04_miR_DE$DN$geneID, HuFNSC01_03 = GE02_GE03_miR_DE$DN$geneID, HuFNSC02_04 = GE02_GE04_miR_DE$DN$geneID)
venn_miR_DE_GW_GE_DN <- venn.diagram(miR_DE_GW_GE_DN, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of GW_GE DN miRNAs")
miR_DE_GW_summary <- rbind(Cortex01_Cortex03_miR_DE$summary, Cortex01_Cortex04_miR_DE$summary, Cortex02_Cortex03_miR_DE$summary, Cortex02_Cortex04_miR_DE$summary, Cortex03_Cortex04_miR_DE$summary, 
                           GE01_GE03_miR_DE$summary, GE01_GE04_miR_DE$summary, GE02_GE03_miR_DE$summary, GE02_GE04_miR_DE$summary, GE03_GE04_miR_DE$summary) %>% 
  mutate(GW = rep(c("17 vs 15", "17 vs 13", "17 vs 15", "17 vs 13", "15 vs 13"), 2))
grid.arrange(gTree(children = venn_miR_DE_GW_Cortex_UP), gTree(children = venn_miR_DE_GW_Cortex_DN), gTree(children = venn_miR_DE_GW_GE_UP), gTree(children = venn_miR_DE_GW_GE_DN), nrow = 2)
miR_DE_GW_Cortex <- list(GW17_GW15_1 = c(Cortex01_Cortex03_miR_DE$UP$geneID, Cortex01_Cortex03_miR_DE$DE$geneID), GW17_GW15_2 = c(Cortex02_Cortex03_miR_DE$UP$geneID, Cortex02_Cortex03_miR_DE$DE$geneID), GW15_GW13 = c(Cortex03_Cortex04_miR_DE$UP$geneID, Cortex03_Cortex04_miR_DE$DE$geneID))
venn_miR_DE_GW_Cortex <- venn.diagram(miR_DE_GW_Cortex, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of GW Cortex miRNAs")
miR_DE_GW_GE <- list(GW17_GW15_1 = c(GE01_GE03_miR_DE$UP$geneID, GE01_GE03_miR_DE$DE$geneID), GW17_GW15_2 = c(GE02_GE03_miR_DE$UP$geneID, GE02_GE03_miR_DE$DE$geneID), GW15_GW13 = c(GE03_GE04_miR_DE$UP$geneID, GE03_GE04_miR_DE$DE$geneID))
venn_miR_DE_GW_GE <- venn.diagram(miR_DE_GW_GE, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of GW GE miRNAs")
pdf("/projects/epigenomics/users/lli/FetalBrain/miR/DE/venn_miR_DE_GW.pdf", height = 4, width = 8)
grid.arrange(gTree(children = venn_miR_DE_GW_Cortex), gTree(children = venn_miR_DE_GW_GE), nrow = 1)
dev.off()
miR_DE_GW_Cortex_UP_intersect <- merge(merge(Cortex01_Cortex03_miR_DE$UP, Cortex01_Cortex04_miR_DE$UP, by = "geneID"), merge(Cortex02_Cortex03_miR_DE$UP, Cortex02_Cortex04_miR_DE$UP, by = "geneID"), by = "geneID")
miR_DE_GW_Cortex_UP_intersect <- miR_DE_GW_Cortex_UP_intersect %>% arrange(-(miR_DE_GW_Cortex_UP_intersect %>% select(grep("logFC", colnames(miR_DE_GW_Cortex_UP_intersect))) %>% rowMeans))
miR_DE_GW_Cortex_DN_intersect <- merge(merge(Cortex01_Cortex03_miR_DE$DN, Cortex01_Cortex04_miR_DE$DN, by = "geneID"), merge(Cortex02_Cortex03_miR_DE$DN, Cortex02_Cortex04_miR_DE$DN, by = "geneID"), by = "geneID")
miR_DE_GW_Cortex_DN_intersect <- miR_DE_GW_Cortex_DN_intersect %>% arrange(miR_DE_GW_Cortex_DN_intersect %>% select(grep("logFC", colnames(miR_DE_GW_Cortex_DN_intersect))) %>% rowMeans)
miR_DE_GW_GE_UP_intersect <- merge(merge(GE01_GE03_miR_DE$UP, GE01_GE04_miR_DE$UP, by = "geneID"), merge(GE02_GE03_miR_DE$UP, GE02_GE04_miR_DE$UP, by = "geneID"), by = "geneID")
miR_DE_GW_GE_UP_intersect <- miR_DE_GW_GE_UP_intersect %>% arrange(-(miR_DE_GW_GE_UP_intersect %>% select(grep("logFC", colnames(miR_DE_GW_GE_UP_intersect))) %>% rowMeans))
miR_DE_GW_GE_DN_intersect <- merge(merge(GE01_GE03_miR_DE$DN, GE01_GE04_miR_DE$DN, by = "geneID"), merge(GE02_GE03_miR_DE$DN, GE02_GE04_miR_DE$DN, by = "geneID"), by = "geneID")
miR_DE_GW_GE_DN_intersect <- miR_DE_GW_GE_DN_intersect %>% arrange(miR_DE_GW_GE_DN_intersect %>% select(grep("logFC", colnames(miR_DE_GW_GE_DN_intersect))) %>% rowMeans)
miR_DE_GW_heat <- miR_FetalBrain[unique(c(GE01_GE03_miR_DE$UP$geneID, GE01_GE04_miR_DE$UP$geneID, GE02_GE03_miR_DE$UP$geneID, GE02_GE04_miR_DE$UP$geneID, GE03_GE04_miR_DE$UP$geneID, 
                                          GE01_GE03_miR_DE$DN$geneID, GE01_GE04_miR_DE$DN$geneID, GE02_GE03_miR_DE$DN$geneID, GE02_GE04_miR_DE$DN$geneID, GE03_GE04_miR_DE$DN$geneID, 
                                          Cortex01_Cortex03_miR_DE$UP$CortexneID, Cortex01_Cortex04_miR_DE$UP$CortexneID, Cortex02_Cortex03_miR_DE$UP$CortexneID, Cortex02_Cortex04_miR_DE$UP$CortexneID, Cortex03_Cortex04_miR_DE$UP$CortexneID, 
                                          Cortex01_Cortex03_miR_DE$DN$CortexneID, Cortex01_Cortex04_miR_DE$DN$CortexneID, Cortex02_Cortex03_miR_DE$DN$CortexneID, Cortex02_Cortex04_miR_DE$DN$CortexneID, Cortex03_Cortex04_miR_DE$DN$CortexneID)), 
                                 c("Cortex01", "Cortex02", "Cortex03", "Cortex04", "GE01", "GE02", "GE03", "GE04")]
pdf("/projects/epigenomics/users/lli/FetalBrain/miR/DE/heatmap_miR_DE_GW.pdf", height = 20)
heatmap.2(as.matrix(miR_DE_GW_heat), scale = "row", trace = "none", margins = c(11, 6), keysize = 1, density.info = "none", col = bluered(256), key.title = "", key.xlab = "")
dev.off()

save(miR_FetalBrain, miR_summary, miR_summary_figure, miR_dend, 
     Brain01_Brain02_miR_DE, Cortex01_Cortex02_miR_DE, GE01_GE02_miR_DE, miR_DE_MZ_summary, venn_miR_DE_MZ_UP, venn_miR_DE_MZ_DN, 
     Cortex01_GE01_miR_DE, Cortex02_GE02_miR_DE, Cortex03_GE03_miR_DE, Cortex04_GE04_miR_DE, miR_DE_neurospheres_summary, venn_miR_DE_neurospheres_UP, venn_miR_DE_neurospheres_DN, 
     Cortex03_Cortex04_miR_DE, Cortex01_Cortex03_miR_DE, Cortex01_Cortex04_miR_DE, Cortex02_Cortex03_miR_DE, Cortex02_Cortex04_miR_DE, venn_miR_DE_GW_Cortex_UP, venn_miR_DE_GW_Cortex_DN, 
     GE03_GE04_miR_DE, GE01_GE03_miR_DE, GE01_GE04_miR_DE, GE02_GE03_miR_DE, GE02_GE04_miR_DE, miR_DE_GW_summary, venn_miR_DE_GW_Cortex, venn_miR_DE_GW_GE, 
     miR_DE_MZ_heat, miR_DE_neurospheres_heat, miR_DE_GW_heat, 
     file = "/projects/epigenomics/users/lli/FetalBrain/miR/FetalBrain_miR.Rdata")


