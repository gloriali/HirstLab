# Fetal Brain - Histone modifications FindER 
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)

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

### correlation with expression 
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

save(FindER_summary, FindER_summary_figure, HisMod_RPKM, HisMod_RPKM_figure, 
     UMR_enhancer, UMR_enhancer_enrich, UMR_enhancer_summary_figure, 
     file = "/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FetalBrain_FindER.Rdata")
