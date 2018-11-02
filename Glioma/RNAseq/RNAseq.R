# Glioma RNA-seq expression analysis
library(ggplot2)
library(dplyr)
library(reshape2)
library(limma)
library(VennDiagram)
source('~/HirstLab/Pipeline/R/enrich.R')
load("/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/RNAseq.Rdata")
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/")

## ----- cluster -----------------------
RPKM <- read.delim("./RPKM/RPKM.matrix", row.names = 1, as.is = T)
spearman_RPKM <- cor(RPKM, method = "spearman")
spearman_RPKM_dend <-  hclust(1 - spearman_RPKM %>% as.dist, method = "ward.D2") %>% as.dendrogram
pdf("./RPKM_spearman.pdf")
op <- par(mar = c(5, 4, 4, 10))
plot(spearman_RPKM_dend, main = "pc genes spearman", horiz = TRUE)
par(op)
dev.off()

## ----- DE limma ----------------------
cutoff <- 0.01
e <- 1e-6
IDHmut_NPC_log10 <- as.data.frame(RPKM + e) %>% select(contains("IDHmut"), contains("NPC")) %>% log10() %>% mutate(ENSG = row.names(.), exp = rowSums(. > -1)) %>% filter(exp > 0)
rownames(IDHmut_NPC_log10) <- IDHmut_NPC_log10$ENSG
IDHmut_NPC_log10 <- IDHmut_NPC_log10 %>% select(-ENSG, -exp)
IDHmut_NPC_design <- data.frame(Group = relevel(factor(gsub("\\..*", "", colnames(IDHmut_NPC_log10))), ref = "NPC"), row.names = colnames(IDHmut_NPC_log10))
IDHmut_NPC_DesMat <- model.matrix(~ Group, IDHmut_NPC_design)
IDHmut_NPC_DEfit <- lmFit(IDHmut_NPC_log10, IDHmut_NPC_DesMat)
IDHmut_NPC_DEfitEb <- eBayes(IDHmut_NPC_DEfit)
IDHmut_NPC_DE <- topTable(IDHmut_NPC_DEfitEb, coef = 'GroupIDHmut', number = Inf, adjust.method="BH", p.value = cutoff) %>% mutate(ENSG = row.names(.), DE = ifelse(logFC > 0, "UP", "DN")) %>% select(ENSG, DE, logFC:adj.P.Val) %>% filter(abs(logFC) >= 2)
IDHmut_NPC_UP <- IDHmut_NPC_DE %>% filter(DE == "UP")
IDHmut_NPC_DN <- IDHmut_NPC_DE %>% filter(DE == "DN")
write.table(IDHmut_NPC_DE, file = "./limma/DE.IDHmut_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHmut_NPC_UP, file = "./limma/UP.IDHmut_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHmut_NPC_DN, file = "./limma/DN.IDHmut_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
IDHmut_NPC_limma <- topTable(IDHmut_NPC_DEfitEb, coef = 'GroupIDHmut', number = Inf, adjust.method="BH", p.value = 1) %>% mutate(ENSG = row.names(.), DE = ifelse(ENSG %in% IDHmut_NPC_UP$ENSG, "UP", ifelse(ENSG %in% IDHmut_NPC_DN$ENSG, "DN", "ST"))) 
(IDHmut_NPC_volcano <- ggplot(IDHmut_NPC_limma, aes(logFC, -log10(adj.P.Val), color = DE)) + 
		geom_point(alpha = 0.2, size = 0.5) + 
		scale_colour_manual(values = c("UP" = "red", "DN" = "blue", "ST" = "grey"), guide = "none") + 
		geom_vline(xintercept = c(-2, 2), color = "black", linetype = 2) + 
		geom_hline(yintercept = 2, color = "black", linetype = 2) + 
		annotate("text", x = -5, y = 20, label = paste0("DN=", nrow(IDHmut_NPC_DN)), color = "blue", size = 6) + 
		annotate("text", x = 5, y = 20, label = paste0("UP=", nrow(IDHmut_NPC_UP)), color = "red", size = 6) + 
		labs(title = "IDHmut vs NPC", x = "log2 FC", y = "-log10 BH adj.p.value") + 
		theme_bw())
ggsave(IDHmut_NPC_volcano, file = "./limma/IDHmut_NPC_volcano.pdf", height = 5, width = 5)
(IDHmut_NPC_UP_DAVID <- enrich(name = "UP.IDHmut_NPC_limma", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 1e-5, p = "FDR", erminej = F, height = 7, width = 8))
(IDHmut_NPC_DN_DAVID <- enrich(name = "DN.IDHmut_NPC_limma", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 7, width = 8))
IDHwt_NPC_log10 <- as.data.frame(RPKM + e) %>% select(contains("IDHwt"), contains("NPC")) %>% log10() %>% mutate(ENSG = row.names(.), exp = rowSums(. > -1)) %>% filter(exp > 0)
rownames(IDHwt_NPC_log10) <- IDHwt_NPC_log10$ENSG
IDHwt_NPC_log10 <- IDHwt_NPC_log10 %>% select(-ENSG, -exp)
IDHwt_NPC_design <- data.frame(Group = relevel(factor(gsub("\\..*", "", colnames(IDHwt_NPC_log10))), ref = "NPC"), row.names = colnames(IDHwt_NPC_log10))
IDHwt_NPC_DesMat <- model.matrix(~ Group, IDHwt_NPC_design)
IDHwt_NPC_DEfit <- lmFit(IDHwt_NPC_log10, IDHwt_NPC_DesMat)
IDHwt_NPC_DEfitEb <- eBayes(IDHwt_NPC_DEfit)
IDHwt_NPC_DE <- topTable(IDHwt_NPC_DEfitEb, coef = 'GroupIDHwt', number = Inf, adjust.method="BH", p.value = cutoff) %>% mutate(ENSG = row.names(.), DE = ifelse(logFC > 0, "UP", "DN")) %>% select(ENSG, DE, logFC:adj.P.Val) %>% filter(abs(logFC) >= 2)
IDHwt_NPC_UP <- IDHwt_NPC_DE %>% filter(DE == "UP")
IDHwt_NPC_DN <- IDHwt_NPC_DE %>% filter(DE == "DN")
write.table(IDHwt_NPC_DE, file = "./limma/DE.IDHwt_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHwt_NPC_UP, file = "./limma/UP.IDHwt_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(IDHwt_NPC_DN, file = "./limma/DN.IDHwt_NPC_limma.txt", sep = "\t", quote = F, row.names = F, col.names = T)
IDHwt_NPC_limma <- topTable(IDHwt_NPC_DEfitEb, coef = 'GroupIDHwt', number = Inf, adjust.method="BH", p.value = 1) %>% mutate(ENSG = row.names(.), DE = ifelse(ENSG %in% IDHwt_NPC_UP$ENSG, "UP", ifelse(ENSG %in% IDHwt_NPC_DN$ENSG, "DN", "ST"))) 
(IDHwt_NPC_volcano <- ggplot(IDHwt_NPC_limma, aes(logFC, -log10(adj.P.Val), color = DE)) + 
		geom_point(alpha = 0.2, size = 0.5) + 
		scale_colour_manual(values = c("UP" = "red", "DN" = "blue", "ST" = "grey"), guide = "none") + 
		geom_vline(xintercept = c(-2, 2), color = "black", linetype = 2) + 
		geom_hline(yintercept = 2, color = "black", linetype = 2) + 
		annotate("text", x = -5, y = 20, label = paste0("DN=", nrow(IDHwt_NPC_DN)), color = "blue", size = 6) + 
		annotate("text", x = 5, y = 20, label = paste0("UP=", nrow(IDHwt_NPC_UP)), color = "red", size = 6) + 
		labs(title = "IDHwt vs NPC", x = "log2 FC", y = "-log10 BH adj.p.value") + 
		theme_bw())
ggsave(IDHwt_NPC_volcano, file = "./limma/IDHwt_NPC_volcano.pdf", height = 5, width = 5)
(IDHwt_NPC_UP_DAVID <- enrich(name = "UP.IDHwt_NPC_limma", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 1e-5, p = "FDR", erminej = F, height = 7, width = 8))
(IDHwt_NPC_DN_DAVID <- enrich(name = "DN.IDHwt_NPC_limma", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 7, width = 8))
limma_venn_UP <- venn.diagram(list(IDHmut = IDHmut_NPC_UP$ENSG, IDHwt = IDHwt_NPC_UP$ENSG), filename = NULL, fill = c("red", "blue"), name = "UP compared to NPC")
grid.newpage()
grid.draw(limma_venn_UP)
limma_venn_DN <- venn.diagram(list(IDHmut = IDHmut_NPC_DN$ENSG, IDHwt = IDHwt_NPC_DN$ENSG), filename = NULL, fill = c("red", "blue"), name = "DN compared to NPC")
grid.newpage()
grid.draw(limma_venn_DN)
IDHmut_NPC_UP_specific <- IDHmut_NPC_UP %>% filter(!(ENSG %in% IDHwt_NPC_UP$ENSG))
write.table(IDHmut_NPC_UP_specific, file = "./limma/UP.IDHmut_NPC_limma_specific.txt", sep = "\t", quote = F, row.names = F, col.names = T)
IDHmut_NPC_DN_specific <- IDHmut_NPC_DN %>% filter(!(ENSG %in% IDHwt_NPC_DN$ENSG))
write.table(IDHmut_NPC_DN_specific, file = "./limma/DN.IDHmut_NPC_limma_specific.txt", sep = "\t", quote = F, row.names = F, col.names = T)
(IDHmut_NPC_UP_specific_DAVID <- enrich(name = "UP.IDHmut_NPC_limma_specific", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 7, width = 8))
(IDHmut_NPC_DN_specific_DAVID <- enrich(name = "DN.IDHmut_NPC_limma_specific", dirIn = "./limma/enrich/", dirOut = "./limma/enrich/", fdr = 0.01, p = "FDR", erminej = F, height = 7, width = 8))

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

save(list = c("spearman_RPKM_dend", ls(pattern = "DE"), ls(pattern = "UP"), ls(pattern = "DN"), ls(pattern = "limma"), ls(pattern = "volcano"), ls(pattern = "venn"), ls(pattern = "DAVID")), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/RNAseq.Rdata")
