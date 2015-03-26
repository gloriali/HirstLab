# Fetal Brain - Histone modifications FindER 
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(VennDiagram)
library(gridExtra)
source("/home/lli/HirstLab/Pipeline/R/enrich_GREAT.R")

setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/")
load("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FetalBrain_FindER.Rdata")
## ============ Sanity check ===============
### Summary 
FindER_summary <- read.delim("FindER.summary", head = F, as.is = T, col.names = c("Sample", "Mark", "N_peak", "N_base")) %>%
  mutate(Sample = str_split_fixed(Sample, "\\.", 4)[, 4], Ave_length = N_base/N_peak)
FindER_summary_tall <- melt(FindER_summary, id = c("Sample", "Mark"))
(FindER_summary_figure <- ggplot(FindER_summary_tall, aes(Sample, value)) + 
   geom_line(aes(group = 1, color = Mark)) + 
   facet_grid(variable ~ Mark, scales = "free") + 
   guides(color=FALSE) + 
   xlab("") + 
   ylab("") + 
   theme_bw() + 
   theme(axis.text.x = element_text(angle = 90)))
ggsave(FindER_summary_figure, file = "FindER_summary_figure.pdf", height = 10, width = 10)

## =========== Correlation with RPKM ================
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/")
colname <- c("chr", "start", "end", "id", "gene")
gene_FetalBrain <- read.delim("/home/lli/FetalBrain/RNAseq/rpkm_pc.txt", as.is = T, row.names = 1)
colnames(gene_FetalBrain) <- gsub("\\.HuFNSC", "", colnames(gene_FetalBrain))
HisMod_RPKM <- melt(gene_FetalBrain %>% select(Brain01, Brain02, Cortex01, Cortex02, GE01, GE02, GE04, id, name, coord), id = c("id", "name", "coord"), value.name = "RPKM", variable.name = "Sample") %>%
  mutate(K4_K27 = "neither", K36 = "no_H3K36me3")
# Brain01
H3K4me3_promoter_Brain01 <- read.delim("./H3K4me3/A03486.H3K4me3.Brain01.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain01"])
H3K27me3_promoter_Brain01 <- read.delim("./H3K27me3/A03488.H3K27me3.Brain01.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain01"])
H3K36me3_genebody_Brain01 <- read.delim("./H3K36me3/A03489.H3K36me3.Brain01.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain01"])
HisMod_RPKM[HisMod_RPKM$Sample == "Brain01" & HisMod_RPKM$id %in% H3K36me3_genebody_Brain01$gene, "K36"] <- "H3K36me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Brain01" & HisMod_RPKM$id %in% H3K4me3_promoter_Brain01$gene, "K4_K27"] <- "H3K4me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Brain01" & HisMod_RPKM$id %in% H3K27me3_promoter_Brain01$gene, "K4_K27"] <- "H3K27me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Brain01" & HisMod_RPKM$id %in% H3K4me3_promoter_Brain01$gene & HisMod_RPKM$id %in% H3K27me3_promoter_Brain01$gene, "K4_K27"] <- "bivalent"
# Brain02
H3K4me3_promoter_Brain02 <- read.delim("./H3K4me3/A03494.H3K4me3.Brain02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain02"])
H3K27me3_promoter_Brain02 <- read.delim("./H3K27me3/A03496.H3K27me3.Brain02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain02"])
H3K36me3_genebody_Brain02 <- read.delim("./H3K36me3/A03497.H3K36me3.Brain02.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain02"])
HisMod_RPKM[HisMod_RPKM$Sample == "Brain02" & HisMod_RPKM$id %in% H3K36me3_genebody_Brain02$gene, "K36"] <- "H3K36me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Brain02" & HisMod_RPKM$id %in% H3K4me3_promoter_Brain02$gene, "K4_K27"] <- "H3K4me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Brain02" & HisMod_RPKM$id %in% H3K27me3_promoter_Brain02$gene, "K4_K27"] <- "H3K27me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Brain02" & HisMod_RPKM$id %in% H3K4me3_promoter_Brain02$gene & HisMod_RPKM$id %in% H3K27me3_promoter_Brain02$gene, "K4_K27"] <- "bivalent"
# Cortex01
H3K27me3_promoter_Cortex01 <- read.delim("./H3K27me3/A03272.H3K27me3.Cortex01.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex01"])
H3K36me3_genebody_Cortex01 <- read.delim("./H3K36me3/A03273.H3K36me3.Cortex01.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex01"])
HisMod_RPKM[HisMod_RPKM$Sample == "Cortex01" & HisMod_RPKM$id %in% H3K36me3_genebody_Cortex01$gene, "K36"] <- "H3K36me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Cortex01" & HisMod_RPKM$id %in% H3K27me3_promoter_Cortex01$gene, "K4_K27"] <- "H3K27me3"
# Cortex02
H3K4me3_promoter_Cortex02 <- read.delim("./H3K4me3/A03282.H3K4me3.Cortex02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex02"])
H3K27me3_promoter_Cortex02 <- read.delim("./H3K27me3/A03284.H3K27me3.Cortex02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex02"])
H3K36me3_genebody_Cortex02 <- read.delim("./H3K36me3/A03285.H3K36me3.Cortex02.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex02"])
HisMod_RPKM[HisMod_RPKM$Sample == "Cortex02" & HisMod_RPKM$id %in% H3K36me3_genebody_Cortex02$gene, "K36"] <- "H3K36me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Cortex02" & HisMod_RPKM$id %in% H3K4me3_promoter_Cortex02$gene, "K4_K27"] <- "H3K4me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Cortex02" & HisMod_RPKM$id %in% H3K27me3_promoter_Cortex02$gene, "K4_K27"] <- "H3K27me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Cortex02" & HisMod_RPKM$id %in% H3K4me3_promoter_Cortex02$gene & HisMod_RPKM$id %in% H3K27me3_promoter_Cortex02$gene, "K4_K27"] <- "bivalent"
# GE01
H3K27me3_promoter_GE01 <- read.delim("./H3K27me3/A03278.H3K27me3.GE01.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE01"])
H3K36me3_genebody_GE01 <- read.delim("./H3K36me3/A03279.H3K36me3.GE01.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE01"])
HisMod_RPKM[HisMod_RPKM$Sample == "GE01" & HisMod_RPKM$id %in% H3K36me3_genebody_GE01$gene, "K36"] <- "H3K36me3"
HisMod_RPKM[HisMod_RPKM$Sample == "GE01" & HisMod_RPKM$id %in% H3K27me3_promoter_GE01$gene, "K4_K27"] <- "H3K27me3"
# GE02
H3K4me3_promoter_GE02 <- read.delim("./H3K4me3/A03478.H3K4me3.GE02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE02"])
H3K27me3_promoter_GE02 <- read.delim("./H3K27me3/A03480.H3K27me3.GE02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE02"])
H3K36me3_genebody_GE02 <- read.delim("./H3K36me3/A03481.H3K36me3.GE02.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE02"])
HisMod_RPKM[HisMod_RPKM$Sample == "GE02" & HisMod_RPKM$id %in% H3K36me3_genebody_GE02$gene, "K36"] <- "H3K36me3"
HisMod_RPKM[HisMod_RPKM$Sample == "GE02" & HisMod_RPKM$id %in% H3K4me3_promoter_GE02$gene, "K4_K27"] <- "H3K4me3"
HisMod_RPKM[HisMod_RPKM$Sample == "GE02" & HisMod_RPKM$id %in% H3K27me3_promoter_GE02$gene, "K4_K27"] <- "H3K27me3"
HisMod_RPKM[HisMod_RPKM$Sample == "GE02" & HisMod_RPKM$id %in% H3K4me3_promoter_GE02$gene & HisMod_RPKM$id %in% H3K27me3_promoter_GE02$gene, "K4_K27"] <- "bivalent"
# GE04
H3K4me3_promoter_GE04 <- read.delim("./H3K4me3/A19304.H3K4me3.GE04.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE04"])
H3K27me3_promoter_GE04 <- read.delim("./H3K27me3/A19306.H3K27me3.GE04.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE04"])
H3K36me3_genebody_GE04 <- read.delim("./H3K36me3/A19307.H3K36me3.GE04.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE04"])
HisMod_RPKM[HisMod_RPKM$Sample == "GE04" & HisMod_RPKM$id %in% H3K36me3_genebody_GE04$gene, "K36"] <- "H3K36me3"
HisMod_RPKM[HisMod_RPKM$Sample == "GE04" & HisMod_RPKM$id %in% H3K4me3_promoter_GE04$gene, "K4_K27"] <- "H3K4me3"
HisMod_RPKM[HisMod_RPKM$Sample == "GE04" & HisMod_RPKM$id %in% H3K27me3_promoter_GE04$gene, "K4_K27"] <- "H3K27me3"
HisMod_RPKM[HisMod_RPKM$Sample == "GE04" & HisMod_RPKM$id %in% H3K4me3_promoter_GE04$gene & HisMod_RPKM$id %in% H3K27me3_promoter_GE04$gene, "K4_K27"] <- "bivalent"
HisMod_RPKM[HisMod_RPKM$RPKM < 0.005, "RPKM"] <- 0.005
HisMod_RPKM_stat <- HisMod_RPKM %>% filter(!(Sample %in% c("Cortex01", "GE01"))) %>%
  mutate(RPKM = log10(RPKM), K4_K27 = factor(K4_K27, levels = c("H3K4me3", "bivalent", "H3K27me3", "neither")), K36 = factor(K36, levels = c("H3K36me3", "no_H3K36me3"))) %>%
  group_by(Sample, K4_K27, K36) %>% 
  summarise(count = n(), lower = quantile(RPKM, 0.25), middle = median(RPKM), upper = quantile(RPKM, 0.75), min = min(RPKM), max = max(RPKM)) %>% 
  mutate(Cell = factor(gsub("\\d+", "", Sample), levels = c("Brain", "Cortex", "GE")))
HisMod_RPKM_stat$ymin <- with(HisMod_RPKM_stat, apply(cbind(min, lower - 1.5*(upper-lower)), 1, max))
HisMod_RPKM_stat$ymax <- with(HisMod_RPKM_stat, apply(cbind(max, upper + 1.5*(upper-lower)), 1, min))
(HisMod_RPKM_figure <- ggplot(HisMod_RPKM_stat, aes(x = Sample, fill = Cell)) + 
   geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
   geom_hline(yintercept = 0) + 
   geom_text(aes(label = count, y = 4, color = Cell), size = 4, angle = 90) + 
   facet_grid(K36 ~ K4_K27, scales = "free_x", space = "free_x") + 
   coord_cartesian(ylim = c(-3, 5)) + 
   scale_fill_manual(values = c("green", "red", "blue"), name = "") + 
   guides(colour=FALSE) +  
   xlab("") + 
   ylab("log10(RPKM)") + 
   theme_bw() + 
   theme(axis.text.x = element_text(angle = 90)))
ggsave(HisMod_RPKM_figure, file = "HisMod_RPKM_figure.pdf", height = 8, width = 10)

## =========== Differential marked genes ================
# promoter signal file generated from ~/HirstLab/FetalBrain/HisMod/TSSsignal.sh and sinal_combine.sh
# normalize promoter signal against the total No. of unique mapped reads in wig files, from bam2wig log files in /home/lli/FetalBrain/HisMod/wigs/
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/signal/")
H3K4me3_TSS2000 <- read.delim("hg19v65_genes_TSS_2000.H3K4me3", head = F, as.is = T, row.names = 1)
colnames(H3K4me3_TSS2000) <- c("Brain01", "Brain02", "Cortex02", "GE02", "GE04")
Nreads_H3K4me3 <- c(Brain01 = 24536698, Brain02 = 26080706, Cortex02 = 16267106, GE02 = 17466495, GE04 = 34305005)
H3K4me3_TSS2000 <- as.data.frame(t(t(H3K4me3_TSS2000) * (Nreads_H3K4me3[1]/Nreads_H3K4me3)))
H3K27me3_TSS2000 <- read.delim("hg19v65_genes_TSS_2000.H3K27me3", head = F, as.is = T, row.names = 1)
colnames(H3K27me3_TSS2000) <- c("Brain01", "Brain02", "Cortex01", "Cortex02", "GE01", "GE02", "GE04")
Nreads_H3K27me3 <- c(Brain01 = 20882658, Brain02 = 37157778, Cortex01 = 28040379, Cortex02 = 44904699, GE01 = 76859708, GE02 = 38306253, GE04 = 35991053)
H3K27me3_TSS2000 <- as.data.frame(t(t(H3K27me3_TSS2000) * (Nreads_H3K27me3[1]/Nreads_H3K27me3)))
e <- 1e-4
fold <- 1 # log2 fold change cutoff (2-fold difference in normailzed signal)
cut <- 0.3 # min sum of signal of two samples, from ecdf plot
# Between MZ twins
Brain01_Brain02_DM_K4me3 <- rbind(data.frame(gene = rownames(H3K4me3_TSS2000[log2((H3K4me3_TSS2000$Brain01 + e)/(H3K4me3_TSS2000$Brain02 + e)) > fold & H3K4me3_TSS2000$Brain01 + H3K4me3_TSS2000$Brain02 > cut, ]), DM = "hyper"), 
                                  data.frame(gene = rownames(H3K4me3_TSS2000[log2((H3K4me3_TSS2000$Brain01 + e)/(H3K4me3_TSS2000$Brain02 + e)) < -fold & H3K4me3_TSS2000$Brain01 + H3K4me3_TSS2000$Brain02 > cut, ]), DM = "hypo"))
rownames(Brain01_Brain02_DM_K4me3) <- Brain01_Brain02_DM_K4me3$gene
Brain01_Brain02_DM_K27me3 <- rbind(data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$Brain01 + e)/(H3K27me3_TSS2000$Brain02 + e)) > fold & H3K27me3_TSS2000$Brain01 + H3K27me3_TSS2000$Brain02 > cut, ]), DM = "hyper"), 
                                   data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$Brain01 + e)/(H3K27me3_TSS2000$Brain02 + e)) < -fold & H3K27me3_TSS2000$Brain01 + H3K27me3_TSS2000$Brain02 > cut, ]), DM = "hypo"))
rownames(Brain01_Brain02_DM_K27me3) <- Brain01_Brain02_DM_K27me3$gene
Cortex01_Cortex02_DM_K27me3 <- rbind(data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$Cortex01 + e)/(H3K27me3_TSS2000$Cortex02 + e)) > fold & H3K27me3_TSS2000$Cortex01 + H3K27me3_TSS2000$Cortex02 > cut, ]), DM = "hyper"), 
                                     data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$Cortex01 + e)/(H3K27me3_TSS2000$Cortex02 + e)) < -fold & H3K27me3_TSS2000$Cortex01 + H3K27me3_TSS2000$Cortex02 > cut, ]), DM = "hypo"))
rownames(Cortex01_Cortex02_DM_K27me3) <- Cortex01_Cortex02_DM_K27me3$gene
GE01_GE02_DM_K27me3 <- rbind(data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$GE01 + e)/(H3K27me3_TSS2000$GE02 + e)) > fold & H3K27me3_TSS2000$GE01 + H3K27me3_TSS2000$GE02 > cut, ]), DM = "hyper"), 
                             data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$GE01 + e)/(H3K27me3_TSS2000$GE02 + e)) < -fold & H3K27me3_TSS2000$GE01 + H3K27me3_TSS2000$GE02 > cut, ]), DM = "hypo"))
rownames(GE01_GE02_DM_K27me3) <- GE01_GE02_DM_K27me3$gene
His_DM_MZ <- rbind(Brain01_Brain02_DM_K4me3 %>% mutate(Sample = "Brain", Mark = "H3K4me3"), 
                   Brain01_Brain02_DM_K27me3 %>% mutate(Sample = "Brain", Mark = "H3K27me3"), 
                   Cortex01_Cortex02_DM_K27me3 %>% mutate(Sample = "Cortex", Mark = "H3K27me3"), 
                   GE01_GE02_DM_K27me3 %>% mutate(Sample = "GE", Mark = "H3K27me3")) %>% mutate(Comparison = "MZ")
# Between neurospheres
Cortex01_GE01_DM_K27me3 <- rbind(data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$Cortex01 + e)/(H3K27me3_TSS2000$GE01 + e)) > fold & H3K27me3_TSS2000$Cortex01 + H3K27me3_TSS2000$GE01 > cut, ]), DM = "hyper"), 
                                 data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$Cortex01 + e)/(H3K27me3_TSS2000$GE01 + e)) < -fold & H3K27me3_TSS2000$Cortex01 + H3K27me3_TSS2000$GE01 > cut, ]), DM = "hypo"))
rownames(Cortex01_GE01_DM_K27me3) <- Cortex01_GE01_DM_K27me3$gene
Cortex02_GE02_DM_K4me3 <- rbind(data.frame(gene = rownames(H3K4me3_TSS2000[log2((H3K4me3_TSS2000$Cortex02 + e)/(H3K4me3_TSS2000$GE02 + e)) > fold & H3K4me3_TSS2000$Cortex02 + H3K4me3_TSS2000$GE02 > cut, ]), DM = "hyper"), 
                                data.frame(gene = rownames(H3K4me3_TSS2000[log2((H3K4me3_TSS2000$Cortex02 + e)/(H3K4me3_TSS2000$GE02 + e)) < -fold & H3K4me3_TSS2000$Cortex02 + H3K4me3_TSS2000$GE02 > cut, ]), DM = "hypo"))
rownames(Cortex02_GE02_DM_K4me3) <- Cortex02_GE02_DM_K4me3$gene
Cortex02_GE02_DM_K27me3 <- rbind(data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$Cortex02 + e)/(H3K27me3_TSS2000$GE02 + e)) > fold & H3K27me3_TSS2000$Cortex02 + H3K27me3_TSS2000$GE02 > cut, ]), DM = "hyper"), 
                                 data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$Cortex02 + e)/(H3K27me3_TSS2000$GE02 + e)) < -fold & H3K27me3_TSS2000$Cortex02 + H3K27me3_TSS2000$GE02 > cut, ]), DM = "hypo"))
rownames(Cortex02_GE02_DM_K27me3) <- Cortex02_GE02_DM_K27me3$gene
His_DM_neurospheres <- rbind(Cortex01_GE01_DM_K27me3 %>% mutate(Sample = "HuFNSC01", Mark = "H3K27me3"), 
                             Cortex02_GE02_DM_K4me3 %>% mutate(Sample = "HuFNSC02", Mark = "H3K4me3"), 
                             Cortex02_GE02_DM_K27me3 %>% mutate(Sample = "HuFNSC02", Mark = "H3K27me3")) %>% mutate(Comparison = "neurospheres")
# Between GW
GE01_GE04_DM_K27me3 <- rbind(data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$GE01 + e)/(H3K27me3_TSS2000$GE04 + e)) > fold & H3K27me3_TSS2000$GE01 + H3K27me3_TSS2000$GE04 > cut, ]), DM = "hyper"), 
                             data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$GE01 + e)/(H3K27me3_TSS2000$GE04 + e)) < -fold & H3K27me3_TSS2000$GE01 + H3K27me3_TSS2000$GE04 > cut, ]), DM = "hypo"))
rownames(GE01_GE04_DM_K27me3) <- GE01_GE04_DM_K27me3$gene
GE02_GE04_DM_K4me3 <- rbind(data.frame(gene = rownames(H3K4me3_TSS2000[log2((H3K4me3_TSS2000$GE02 + e)/(H3K4me3_TSS2000$GE04 + e)) > fold & H3K4me3_TSS2000$GE02 + H3K4me3_TSS2000$GE04 > cut, ]), DM = "hyper"), 
                            data.frame(gene = rownames(H3K4me3_TSS2000[log2((H3K4me3_TSS2000$GE02 + e)/(H3K4me3_TSS2000$GE04 + e)) < -fold & H3K4me3_TSS2000$GE02 + H3K4me3_TSS2000$GE04 > cut, ]), DM = "hypo"))
rownames(GE02_GE04_DM_K4me3) <- GE02_GE04_DM_K4me3$gene
GE02_GE04_DM_K27me3 <- rbind(data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$GE02 + e)/(H3K27me3_TSS2000$GE04 + e)) > fold & H3K27me3_TSS2000$GE02 + H3K27me3_TSS2000$GE04 > cut, ]), DM = "hyper"), 
                             data.frame(gene = rownames(H3K27me3_TSS2000[log2((H3K27me3_TSS2000$GE02 + e)/(H3K27me3_TSS2000$GE04 + e)) < -fold & H3K27me3_TSS2000$GE02 + H3K27me3_TSS2000$GE04 > cut, ]), DM = "hypo"))
rownames(GE02_GE04_DM_K27me3) <- GE02_GE04_DM_K27me3$gene
His_DM_GW <- rbind(GE01_GE04_DM_K27me3 %>% mutate(Sample = "GE01_GE04", Mark = "H3K27me3"), 
                   GE02_GE04_DM_K4me3 %>% mutate(Sample = "GE02_GE04", Mark = "H3K4me3"), 
                   GE02_GE04_DM_K27me3 %>% mutate(Sample = "GE02_GE04", Mark = "H3K27me3")) %>% mutate(Comparison = "GW")
His_DM_promoter <- rbind(His_DM_MZ, His_DM_neurospheres, His_DM_GW)
(His_DM_promoter_figure <- ggplot(His_DM_promoter, aes(x = Sample, fill = DM)) + 
   geom_bar(position = "dodge") + 
   facet_grid(Mark ~ Comparison, scales = "free") + 
   ylab("No. of genes") + 
   theme_bw())
ggsave(His_DM_promoter_figure, file = "His_DM_promoter_figure.pdf")

## =========== Correlation with DE genes ================
load("/home/lli/FetalBrain/RNAseq/DEfine/gene/FetalBrain_DEgenes.Rdata")
load("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR_MZ.Rdata")
load("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR_neurospheres.Rdata")
load("/projects/epigenomics/users/lli/FetalBrain/GW/GW.Rdata")
# Brain01 vs Brain02
Brain01_Brain02DE_epi <- brain01_brain02DE %>% mutate(K4me3 = "not_DM", K27me3 = "not_DM", DMR = "not_DM")
rownames(Brain01_Brain02DE_epi) <- Brain01_Brain02DE_epi$V1
Brain01_Brain02DE_epi[intersect(Brain01_Brain02_DM_K4me3$gene, Brain01_Brain02DE_epi$V1), "K4me3"] <- Brain01_Brain02_DM_K4me3[intersect(Brain01_Brain02_DM_K4me3$gene, Brain01_Brain02DE_epi$V1), "DM"]
Brain01_Brain02DE_epi[intersect(Brain01_Brain02_DM_K27me3$gene, Brain01_Brain02DE_epi$V1), "K27me3"] <- Brain01_Brain02_DM_K27me3[intersect(Brain01_Brain02_DM_K27me3$gene, Brain01_Brain02DE_epi$V1), "DM"]
Brain01_Brain02DE_epi[DMR_DE_Brain01_Brain02_hyper$DMR_gene_DE$id, "DMR"] <- "hyper"
Brain01_Brain02DE_epi[DMR_DE_Brain01_Brain02_hypo$DMR_gene_DE$id, "DMR"] <- "hypo"
Brain01_Brain02DE_epi_list <- list(H3K4me3 = as.character(filter(Brain01_Brain02DE_epi, K4me3 != "not_DM")[, "V1"]), H3K27me3 = as.character(filter(Brain01_Brain02DE_epi, K27me3 != "not_DM")[, "V1"]), DMR = as.character(filter(Brain01_Brain02DE_epi, DMR != "not_DM")[, "V1"]))
venn_Brain01_Brain02DE_epi <- venn.diagram(Brain01_Brain02DE_epi_list, filename = NULL, fill = c("green", "red", "blue"), main = paste0("Brain01 vs Brain02: total ", nrow(Brain01_Brain02DE_epi)))
plot.new()
grid.draw(venn_Brain01_Brain02DE_epi)
# Cortex01 vs Cortex02
Cortex01_Cortex02DE_epi <- cortex01_cortex02DE %>% mutate(K27me3 = "not_DM", DMR = "not_DM")
rownames(Cortex01_Cortex02DE_epi) <- Cortex01_Cortex02DE_epi$V1
Cortex01_Cortex02DE_epi[intersect(Cortex01_Cortex02_DM_K27me3$gene, Cortex01_Cortex02DE_epi$V1), "K27me3"] <- Cortex01_Cortex02_DM_K27me3[intersect(Cortex01_Cortex02_DM_K27me3$gene, Cortex01_Cortex02DE_epi$V1), "DM"]
Cortex01_Cortex02DE_epi[DMR_DE_Cortex01_Cortex02_hyper$DMR_gene_DE$id, "DMR"] <- "hyper"
Cortex01_Cortex02DE_epi[DMR_DE_Cortex01_Cortex02_hypo$DMR_gene_DE$id, "DMR"] <- "hypo"
Cortex01_Cortex02DE_epi_list <- list(H3K27me3 = as.character(filter(Cortex01_Cortex02DE_epi, K27me3 != "not_DM")[, "V1"]), DMR = as.character(filter(Cortex01_Cortex02DE_epi, DMR != "not_DM")[, "V1"]))
venn_Cortex01_Cortex02DE_epi <- venn.diagram(Cortex01_Cortex02DE_epi_list, filename = NULL, fill = c("red", "blue"), main = paste0("Cortex01 vs Cortex02: total ", nrow(Cortex01_Cortex02DE_epi)))
plot.new()
grid.draw(venn_Cortex01_Cortex02DE_epi)
# GE01 vs GE02
GE01_GE02DE_epi <- GE01_GE02DE %>% mutate(K27me3 = "not_DM", DMR = "not_DM")
rownames(GE01_GE02DE_epi) <- GE01_GE02DE_epi$V1
GE01_GE02DE_epi[intersect(GE01_GE02_DM_K27me3$gene, GE01_GE02DE_epi$V1), "K27me3"] <- GE01_GE02_DM_K27me3[intersect(GE01_GE02_DM_K27me3$gene, GE01_GE02DE_epi$V1), "DM"]
GE01_GE02DE_epi[DMR_DE_GE01_GE02_hyper$DMR_gene_DE$id, "DMR"] <- "hyper"
GE01_GE02DE_epi[DMR_DE_GE01_GE02_hypo$DMR_gene_DE$id, "DMR"] <- "hypo"
GE01_GE02DE_epi_list <- list(H3K27me3 = as.character(filter(GE01_GE02DE_epi, K27me3 != "not_DM")[, "V1"]), DMR = as.character(filter(GE01_GE02DE_epi, DMR != "not_DM")[, "V1"]))
venn_GE01_GE02DE_epi <- venn.diagram(GE01_GE02DE_epi_list, filename = NULL, fill = c("red", "blue"), main = paste0("GE01 vs GE02: total ", nrow(GE01_GE02DE_epi)))
plot.new()
grid.draw(venn_GE01_GE02DE_epi)
# Cortex01 vs GE01
Cortex01_GE01DE_epi <- cortex01_GE01DE %>% mutate(K27me3 = "not_DM", DMR = "not_DM")
rownames(Cortex01_GE01DE_epi) <- Cortex01_GE01DE_epi$V1
Cortex01_GE01DE_epi[intersect(Cortex01_GE01_DM_K27me3$gene, Cortex01_GE01DE_epi$V1), "K27me3"] <- Cortex01_GE01_DM_K27me3[intersect(Cortex01_GE01_DM_K27me3$gene, Cortex01_GE01DE_epi$V1), "DM"]
Cortex01_GE01DE_epi[DMR_DE_Cortex01_GE01_hyper$DMR_gene_DE$id, "DMR"] <- "hyper"
Cortex01_GE01DE_epi[DMR_DE_Cortex01_GE01_hypo$DMR_gene_DE$id, "DMR"] <- "hypo"
Cortex01_GE01DE_epi_list <- list(H3K27me3 = as.character(filter(Cortex01_GE01DE_epi, K27me3 != "not_DM")[, "V1"]), DMR = as.character(filter(Cortex01_GE01DE_epi, DMR != "not_DM")[, "V1"]))
venn_Cortex01_GE01DE_epi <- venn.diagram(Cortex01_GE01DE_epi_list, filename = NULL, fill = c("red", "blue"), main = paste0("Cortex01 vs GE01: total ", nrow(Cortex01_GE01DE_epi)))
plot.new()
grid.draw(venn_Cortex01_GE01DE_epi)
# Cortex02 vs GE02
Cortex02_GE02DE_epi <- cortex02_GE02DE %>% mutate(K4me3 = "not_DM", K27me3 = "not_DM", DMR = "not_DM")
rownames(Cortex02_GE02DE_epi) <- Cortex02_GE02DE_epi$V1
Cortex02_GE02DE_epi[intersect(Cortex02_GE02_DM_K4me3$gene, Cortex02_GE02DE_epi$V1), "K4me3"] <- Cortex02_GE02_DM_K4me3[intersect(Cortex02_GE02_DM_K4me3$gene, Cortex02_GE02DE_epi$V1), "DM"]
Cortex02_GE02DE_epi[intersect(Cortex02_GE02_DM_K27me3$gene, Cortex02_GE02DE_epi$V1), "K27me3"] <- Cortex02_GE02_DM_K27me3[intersect(Cortex02_GE02_DM_K27me3$gene, Cortex02_GE02DE_epi$V1), "DM"]
Cortex02_GE02DE_epi[DMR_DE_Cortex02_GE02_hyper$DMR_gene_DE$id, "DMR"] <- "hyper"
Cortex02_GE02DE_epi[DMR_DE_Cortex02_GE02_hypo$DMR_gene_DE$id, "DMR"] <- "hypo"
Cortex02_GE02DE_epi_list <- list(H3K4me3 = as.character(filter(Cortex02_GE02DE_epi, K4me3 != "not_DM")[, "V1"]), H3K27me3 = as.character(filter(Cortex02_GE02DE_epi, K27me3 != "not_DM")[, "V1"]), DMR = as.character(filter(Cortex02_GE02DE_epi, DMR != "not_DM")[, "V1"]))
venn_Cortex02_GE02DE_epi <- venn.diagram(Cortex02_GE02DE_epi_list, filename = NULL, fill = c("green", "red", "blue"), main = paste0("Cortex02 vs GE02: total ", nrow(Cortex02_GE02DE_epi)))
plot.new()
grid.draw(venn_Cortex02_GE02DE_epi)
# GE02 vs GE04
GE02_GE04DE_epi <- GE02_GE04DE %>% mutate(K4me3 = "not_DM", K27me3 = "not_DM", DMR = "not_DM")
rownames(GE02_GE04DE_epi) <- GE02_GE04DE_epi$ID
GE02_GE04DE_epi[intersect(GE02_GE04_DM_K4me3$gene, GE02_GE04DE_epi$ID), "K4me3"] <- GE02_GE04_DM_K4me3[intersect(GE02_GE04_DM_K4me3$gene, GE02_GE04DE_epi$ID), "DM"]
GE02_GE04DE_epi[intersect(GE02_GE04_DM_K27me3$gene, GE02_GE04DE_epi$ID), "K27me3"] <- GE02_GE04_DM_K27me3[intersect(GE02_GE04_DM_K27me3$gene, GE02_GE04DE_epi$ID), "DM"]
GE02_GE04DE_epi[DMR_DE_GE02_GE04_hyper$DMR_gene_DE$id, "DMR"] <- "hyper"
GE02_GE04DE_epi[DMR_DE_GE02_GE04_hypo$DMR_gene_DE$id, "DMR"] <- "hypo"
GE02_GE04DE_epi_list <- list(H3K4me3 = as.character(filter(GE02_GE04DE_epi, K4me3 != "not_DM")[, "ID"]), H3K27me3 = as.character(filter(GE02_GE04DE_epi, K27me3 != "not_DM")[, "ID"]), DMR = as.character(filter(GE02_GE04DE_epi, DMR != "not_DM")[, "ID"]))
venn_GE02_GE04DE_epi <- venn.diagram(GE02_GE04DE_epi_list, filename = NULL, fill = c("green", "red", "blue"), main = paste0("GE02 vs GE04: total ", nrow(GE02_GE04DE_epi)))
plot.new()
grid.draw(venn_GE02_GE04DE_epi)
grid.arrange(gTree(children = venn_Brain01_Brain02DE_epi), gTree(children = venn_Cortex01_Cortex02DE_epi), gTree(children = venn_GE01_GE02DE_epi), 
             gTree(children = venn_Cortex01_GE01DE_epi), gTree(children = venn_Cortex02_GE02DE_epi), gTree(children = venn_GE02_GE04DE_epi), nrow = 3)

## ============= enhancers: H3K4me1 =============
### closest genes
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1")
colname <- c("chr", "start", "end", "id", "gene", "dis")
closest_gene_Cortex01 <- read.delim("A03269.H3K4me1.Cortex01.closest.gene.pc", head = F, as.is = T, col.names = colname)
closest_gene_io_Cortex01 <- read.delim("A03269.H3K4me1.Cortex01.closest.gene.pc.io", head = F, as.is = T, col.names = colname)
length(unique(closest_gene_Cortex01$gene))     # 14660
length(unique(closest_gene_io_Cortex01$gene))  # 14630
length(intersect(unique(closest_gene_Cortex01$gene), unique(closest_gene_io_Cortex01$gene)))   # 11978
# ignore overlap doesn't have much effect on closest genes, use including overlap version. 
# closest genes included more than half of all pc genes, no point in enrichment analysis.  

### overlap with WGBS UMRs
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1")
UMR_enhancer <- read.delim("WGBS_UMR_enhancers.summary", as.is = T) %>%
  mutate(percent = No.enhancerUMR/No.UMR)
UMR_enhancer_enrich <- read.delim("./CpG/WGBS_UMR_enhancers_enrich.summary", as.is = T)
UMR_enhancer_summary <- data.frame(Sample = rep(c("Cortex02_GE02", "Cortex04_GE04", "Cortex02_Cortex04", "GE02_GE04"), each = 2), 
                                   compare = rep(c("Neurospheres", "GW"), each = 4), 
                                   UMR = UMR_enhancer$UMR, percentage = UMR_enhancer$percent*100, log2_enrichment = log2(UMR_enhancer_enrich$enrich)) %>%
  mutate(UMR = as.character(UMR), compare = as.character(compare))
UMR_enhancer_summary[UMR_enhancer_summary$compare == "GW" & UMR_enhancer_summary$UMR == "hyper", "UMR"] <- "GW13/Cortex UMRs"
UMR_enhancer_summary[UMR_enhancer_summary$compare == "GW" & UMR_enhancer_summary$UMR == "hypo", "UMR"] <- "GW17/GE UMRs"
UMR_enhancer_summary[UMR_enhancer_summary$compare == "Neurospheres" & UMR_enhancer_summary$UMR == "hypo", "UMR"] <- "GW13/Cortex UMRs"
UMR_enhancer_summary[UMR_enhancer_summary$compare == "Neurospheres" & UMR_enhancer_summary$UMR == "hyper", "UMR"] <- "GW17/GE UMRs"
UMR_enhancer_summary_tall <- melt(UMR_enhancer_summary, id = c("Sample", "compare", "UMR")) %>%
  mutate(compare = factor(compare, levels = c("Neurospheres", "GW")), UMR = factor(UMR, levels = c("GW13/Cortex UMRs", "GW17/GE UMRs")))
(UMR_enhancer_summary_figure <- ggplot(UMR_enhancer_summary_tall, aes(x = Sample, y = value, fill = UMR)) + 
   geom_bar(stat = "identity", position = position_dodge()) + 
   facet_grid(variable ~ compare, scales = "free") + 
   scale_fill_manual(values = c("red", "blue"), name = "") + 
   xlab("") + 
   ylab("") + 
   theme_bw())
ggsave(UMR_enhancer_summary_figure, file = "UMR_enhancer_summary_figure.pdf")
(GREAT_Cortex_UMR_enhancer_HuFNSC02 <- enrich_GREAT(file = "Cortex_UMR_enhancer_HuFNSC02", name = "Cortex_UMR_enhancer_HuFNSC02", height = 15))
(GREAT_GE_UMR_enhancer_HuFNSC02 <- enrich_GREAT(file = "GE_UMR_enhancer_HuFNSC02", name = "GE_UMR_enhancer_HuFNSC02", height = 15))
(GREAT_Cortex_UMR_enhancer_HuFNSC04 <- enrich_GREAT(file = "Cortex_UMR_enhancer_HuFNSC04", name = "Cortex_UMR_enhancer_HuFNSC04", height = 15))
(GREAT_GW17_UMR_enhancer_Cortex <- enrich_GREAT(file = "GW17_UMR_enhancer_Cortex", name = "GW17_UMR_enhancer_Cortex", height = 15))
(GREAT_GW17_UMR_enhancer_GE <- enrich_GREAT(file = "GW17_UMR_enhancer_GE", name = "GW17_UMR_enhancer_GE", height = 15))
(GREAT_GW13_UMR_enhancer_GE <- enrich_GREAT(file = "GW13_UMR_enhancer_GE", name = "GW13_UMR_enhancer_GE", height = 15))

save(FindER_summary, FindER_summary_figure, HisMod_RPKM, HisMod_RPKM_figure, 
     UMR_enhancer, UMR_enhancer_enrich, UMR_enhancer_summary_figure, 
     GREAT_Cortex_UMR_enhancer_HuFNSC02, GREAT_GE_UMR_enhancer_HuFNSC02, GREAT_Cortex_UMR_enhancer_HuFNSC04, 
     GREAT_GW17_UMR_enhancer_Cortex, GREAT_GW17_UMR_enhancer_GE, GREAT_GW13_UMR_enhancer_GE, 
     His_DM_promoter, His_DM_promoter_figure, 
     Brain01_Brain02DE_epi, Cortex01_Cortex02DE_epi, GE01_GE02DE_epi, Cortex01_GE01DE_epi, Cortex02_GE02DE_epi, GE02_GE04DE_epi, 
     venn_Brain01_Brain02DE_epi, venn_Cortex01_Cortex02DE_epi, venn_GE01_GE02DE_epi, venn_Cortex01_GE01DE_epi, venn_Cortex02_GE02DE_epi, venn_GE02_GE04DE_epi, 
     file = "/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FetalBrain_FindER.Rdata")
