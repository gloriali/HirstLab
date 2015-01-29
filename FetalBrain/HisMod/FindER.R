# Fetal Brain - Histone modifications FindER 
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)

setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/")
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
colname <- c("chr", "start", "end", "id", "gene")
gene_FetalBrain <- read.delim("/home/lli/FetalBrain/RNAseq/rpkm_pc.txt", as.is = T, row.names = 1)
colnames(gene_FetalBrain) <- gsub("\\.HuFNSC", "", colnames(gene_FetalBrain))
# Brain01
H3K4me3_promoter_Brain01 <- read.delim("./H3K4me3/A03486.H3K4me3.Brain01.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain01"])
H3K4me3_T_RPKM_Brain01_stat <- summarise(na.omit(gene_FetalBrain[H3K4me3_promoter_Brain01$gene, ]), count = n(), lower = quantile(Brain01, 0.25), middle = median(Brain01), upper = quantile(Brain01, 0.75), min = min(Brain01), max = max(Brain01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain01", Mark = "H3K4me3", Marked = T)
H3K4me3_F_RPKM_Brain01_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K4me3_promoter_Brain01$gene), ], count = n(), lower = quantile(Brain01, 0.25), middle = median(Brain01), upper = quantile(Brain01, 0.75), min = min(Brain01), max = max(Brain01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain01", Mark = "H3K4me3", Marked = F)
H3K27me3_promoter_Brain01 <- read.delim("./H3K27me3/A03488.H3K27me3.Brain01.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain01"])
H3K27me3_T_RPKM_Brain01_stat <- summarise(na.omit(gene_FetalBrain[H3K27me3_promoter_Brain01$gene, ]), count = n(), lower = quantile(Brain01, 0.25), middle = median(Brain01), upper = quantile(Brain01, 0.75), min = min(Brain01), max = max(Brain01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain01", Mark = "H3K27me3", Marked = T)
H3K27me3_F_RPKM_Brain01_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K27me3_promoter_Brain01$gene), ], count = n(), lower = quantile(Brain01, 0.25), middle = median(Brain01), upper = quantile(Brain01, 0.75), min = min(Brain01), max = max(Brain01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain01", Mark = "H3K27me3", Marked = F)
H3K36me3_genebody_Brain01 <- read.delim("./H3K36me3/A03489.H3K36me3.Brain01.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain01"])
H3K36me3_T_RPKM_Brain01_stat <- summarise(na.omit(gene_FetalBrain[H3K36me3_genebody_Brain01$gene, ]), count = n(), lower = quantile(Brain01, 0.25), middle = median(Brain01), upper = quantile(Brain01, 0.75), min = min(Brain01), max = max(Brain01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain01", Mark = "H3K36me3", Marked = T)
H3K36me3_F_RPKM_Brain01_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K36me3_genebody_Brain01$gene), ], count = n(), lower = quantile(Brain01, 0.25), middle = median(Brain01), upper = quantile(Brain01, 0.75), min = min(Brain01), max = max(Brain01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain01", Mark = "H3K36me3", Marked = F)
H3K4me3_H3K27me3_RPKM_Brain01_stat <- summarise(na.omit(gene_FetalBrain[intersect(H3K4me3_promoter_Brain01$gene, H3K27me3_promoter_Brain01$gene), ]), count = n(), lower = quantile(Brain01, 0.25), middle = median(Brain01), upper = quantile(Brain01, 0.75), min = min(Brain01), max = max(Brain01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain01", Mark = "H3K4me3+H3K27me3", Marked = T)
HisMod_RPKM_Brain01_stat <- rbind(H3K4me3_T_RPKM_Brain01_stat, H3K4me3_F_RPKM_Brain01_stat, H3K27me3_T_RPKM_Brain01_stat, H3K27me3_F_RPKM_Brain01_stat, H3K36me3_T_RPKM_Brain01_stat, H3K36me3_F_RPKM_Brain01_stat, H3K4me3_H3K27me3_RPKM_Brain01_stat)
# Brain02
H3K4me3_promoter_Brain02 <- read.delim("./H3K4me3/A03494.H3K4me3.Brain02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain02"])
H3K4me3_T_RPKM_Brain02_stat <- summarise(na.omit(gene_FetalBrain[H3K4me3_promoter_Brain02$gene, ]), count = n(), lower = quantile(Brain02, 0.25), middle = median(Brain02), upper = quantile(Brain02, 0.75), min = min(Brain02), max = max(Brain02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain02", Mark = "H3K4me3", Marked = T)
H3K4me3_F_RPKM_Brain02_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K4me3_promoter_Brain02$gene), ], count = n(), lower = quantile(Brain02, 0.25), middle = median(Brain02), upper = quantile(Brain02, 0.75), min = min(Brain02), max = max(Brain02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain02", Mark = "H3K4me3", Marked = F)
H3K27me3_promoter_Brain02 <- read.delim("./H3K27me3/A03496.H3K27me3.Brain02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain02"])
H3K27me3_T_RPKM_Brain02_stat <- summarise(na.omit(gene_FetalBrain[H3K27me3_promoter_Brain02$gene, ]), count = n(), lower = quantile(Brain02, 0.25), middle = median(Brain02), upper = quantile(Brain02, 0.75), min = min(Brain02), max = max(Brain02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain02", Mark = "H3K27me3", Marked = T)
H3K27me3_F_RPKM_Brain02_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K27me3_promoter_Brain02$gene), ], count = n(), lower = quantile(Brain02, 0.25), middle = median(Brain02), upper = quantile(Brain02, 0.75), min = min(Brain02), max = max(Brain02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain02", Mark = "H3K27me3", Marked = F)
H3K36me3_genebody_Brain02 <- read.delim("./H3K36me3/A03497.H3K36me3.Brain02.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Brain02"])
H3K36me3_T_RPKM_Brain02_stat <- summarise(na.omit(gene_FetalBrain[H3K36me3_genebody_Brain02$gene, ]), count = n(), lower = quantile(Brain02, 0.25), middle = median(Brain02), upper = quantile(Brain02, 0.75), min = min(Brain02), max = max(Brain02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain02", Mark = "H3K36me3", Marked = T)
H3K36me3_F_RPKM_Brain02_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K36me3_genebody_Brain02$gene), ], count = n(), lower = quantile(Brain02, 0.25), middle = median(Brain02), upper = quantile(Brain02, 0.75), min = min(Brain02), max = max(Brain02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain02", Mark = "H3K36me3", Marked = F)
H3K4me3_H3K27me3_RPKM_Brain02_stat <- summarise(na.omit(gene_FetalBrain[intersect(H3K4me3_promoter_Brain02$gene, H3K27me3_promoter_Brain02$gene), ]), count = n(), lower = quantile(Brain02, 0.25), middle = median(Brain02), upper = quantile(Brain02, 0.75), min = min(Brain02), max = max(Brain02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Brain02", Mark = "H3K4me3+H3K27me3", Marked = T)
HisMod_RPKM_Brain02_stat <- rbind(H3K4me3_T_RPKM_Brain02_stat, H3K4me3_F_RPKM_Brain02_stat, H3K27me3_T_RPKM_Brain02_stat, H3K27me3_F_RPKM_Brain02_stat, H3K36me3_T_RPKM_Brain02_stat, H3K36me3_F_RPKM_Brain02_stat, H3K4me3_H3K27me3_RPKM_Brain02_stat)
# Cortex01
H3K27me3_promoter_Cortex01 <- read.delim("./H3K27me3/A03272.H3K27me3.Cortex01.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex01"])
H3K27me3_T_RPKM_Cortex01_stat <- summarise(na.omit(gene_FetalBrain[H3K27me3_promoter_Cortex01$gene, ]), count = n(), lower = quantile(Cortex01, 0.25), middle = median(Cortex01), upper = quantile(Cortex01, 0.75), min = min(Cortex01), max = max(Cortex01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex01", Mark = "H3K27me3", Marked = T)
H3K27me3_F_RPKM_Cortex01_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K27me3_promoter_Cortex01$gene), ], count = n(), lower = quantile(Cortex01, 0.25), middle = median(Cortex01), upper = quantile(Cortex01, 0.75), min = min(Cortex01), max = max(Cortex01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex01", Mark = "H3K27me3", Marked = F)
H3K36me3_genebody_Cortex01 <- read.delim("./H3K36me3/A03273.H3K36me3.Cortex01.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex01"])
H3K36me3_T_RPKM_Cortex01_stat <- summarise(na.omit(gene_FetalBrain[H3K36me3_genebody_Cortex01$gene, ]), count = n(), lower = quantile(Cortex01, 0.25), middle = median(Cortex01), upper = quantile(Cortex01, 0.75), min = min(Cortex01), max = max(Cortex01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex01", Mark = "H3K36me3", Marked = T)
H3K36me3_F_RPKM_Cortex01_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K36me3_genebody_Cortex01$gene), ], count = n(), lower = quantile(Cortex01, 0.25), middle = median(Cortex01), upper = quantile(Cortex01, 0.75), min = min(Cortex01), max = max(Cortex01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex01", Mark = "H3K36me3", Marked = F)
HisMod_RPKM_Cortex01_stat <- rbind(H3K27me3_T_RPKM_Cortex01_stat, H3K27me3_F_RPKM_Cortex01_stat, H3K36me3_T_RPKM_Cortex01_stat, H3K36me3_F_RPKM_Cortex01_stat)
# Cortex02
H3K4me3_promoter_Cortex02 <- read.delim("./H3K4me3/A03282.H3K4me3.Cortex02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex02"])
H3K4me3_T_RPKM_Cortex02_stat <- summarise(na.omit(gene_FetalBrain[H3K4me3_promoter_Cortex02$gene, ]), count = n(), lower = quantile(Cortex02, 0.25), middle = median(Cortex02), upper = quantile(Cortex02, 0.75), min = min(Cortex02), max = max(Cortex02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex02", Mark = "H3K4me3", Marked = T)
H3K4me3_F_RPKM_Cortex02_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K4me3_promoter_Cortex02$gene), ], count = n(), lower = quantile(Cortex02, 0.25), middle = median(Cortex02), upper = quantile(Cortex02, 0.75), min = min(Cortex02), max = max(Cortex02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex02", Mark = "H3K4me3", Marked = F)
H3K27me3_promoter_Cortex02 <- read.delim("./H3K27me3/A03284.H3K27me3.Cortex02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex02"])
H3K27me3_T_RPKM_Cortex02_stat <- summarise(na.omit(gene_FetalBrain[H3K27me3_promoter_Cortex02$gene, ]), count = n(), lower = quantile(Cortex02, 0.25), middle = median(Cortex02), upper = quantile(Cortex02, 0.75), min = min(Cortex02), max = max(Cortex02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex02", Mark = "H3K27me3", Marked = T)
H3K27me3_F_RPKM_Cortex02_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K27me3_promoter_Cortex02$gene), ], count = n(), lower = quantile(Cortex02, 0.25), middle = median(Cortex02), upper = quantile(Cortex02, 0.75), min = min(Cortex02), max = max(Cortex02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex02", Mark = "H3K27me3", Marked = F)
H3K36me3_genebody_Cortex02 <- read.delim("./H3K36me3/A03285.H3K36me3.Cortex02.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex02"])
H3K36me3_T_RPKM_Cortex02_stat <- summarise(na.omit(gene_FetalBrain[H3K36me3_genebody_Cortex02$gene, ]), count = n(), lower = quantile(Cortex02, 0.25), middle = median(Cortex02), upper = quantile(Cortex02, 0.75), min = min(Cortex02), max = max(Cortex02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex02", Mark = "H3K36me3", Marked = T)
H3K36me3_F_RPKM_Cortex02_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K36me3_genebody_Cortex02$gene), ], count = n(), lower = quantile(Cortex02, 0.25), middle = median(Cortex02), upper = quantile(Cortex02, 0.75), min = min(Cortex02), max = max(Cortex02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex02", Mark = "H3K36me3", Marked = F)
H3K4me3_H3K27me3_RPKM_Cortex02_stat <- summarise(na.omit(gene_FetalBrain[intersect(H3K4me3_promoter_Cortex02$gene, H3K27me3_promoter_Cortex02$gene), ]), count = n(), lower = quantile(Cortex02, 0.25), middle = median(Cortex02), upper = quantile(Cortex02, 0.75), min = min(Cortex02), max = max(Cortex02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "Cortex02", Mark = "H3K4me3+H3K27me3", Marked = T)
HisMod_RPKM_Cortex02_stat <- rbind(H3K4me3_T_RPKM_Cortex02_stat, H3K4me3_F_RPKM_Cortex02_stat, H3K27me3_T_RPKM_Cortex02_stat, H3K27me3_F_RPKM_Cortex02_stat, H3K36me3_T_RPKM_Cortex02_stat, H3K36me3_F_RPKM_Cortex02_stat, H3K4me3_H3K27me3_RPKM_Cortex02_stat)
# GE01
H3K27me3_promoter_GE01 <- read.delim("./H3K27me3/A03278.H3K27me3.GE01.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE01"])
H3K27me3_T_RPKM_GE01_stat <- summarise(na.omit(gene_FetalBrain[H3K27me3_promoter_GE01$gene, ]), count = n(), lower = quantile(GE01, 0.25), middle = median(GE01), upper = quantile(GE01, 0.75), min = min(GE01), max = max(GE01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE01", Mark = "H3K27me3", Marked = T)
H3K27me3_F_RPKM_GE01_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K27me3_promoter_GE01$gene), ], count = n(), lower = quantile(GE01, 0.25), middle = median(GE01), upper = quantile(GE01, 0.75), min = min(GE01), max = max(GE01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE01", Mark = "H3K27me3", Marked = F)
H3K36me3_genebody_GE01 <- read.delim("./H3K36me3/A03279.H3K36me3.GE01.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE01"])
H3K36me3_T_RPKM_GE01_stat <- summarise(na.omit(gene_FetalBrain[H3K36me3_genebody_GE01$gene, ]), count = n(), lower = quantile(GE01, 0.25), middle = median(GE01), upper = quantile(GE01, 0.75), min = min(GE01), max = max(GE01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE01", Mark = "H3K36me3", Marked = T)
H3K36me3_F_RPKM_GE01_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K36me3_genebody_GE01$gene), ], count = n(), lower = quantile(GE01, 0.25), middle = median(GE01), upper = quantile(GE01, 0.75), min = min(GE01), max = max(GE01)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE01", Mark = "H3K36me3", Marked = F)
HisMod_RPKM_GE01_stat <- rbind(H3K27me3_T_RPKM_GE01_stat, H3K27me3_F_RPKM_GE01_stat, H3K36me3_T_RPKM_GE01_stat, H3K36me3_F_RPKM_GE01_stat)
# GE02
H3K4me3_promoter_GE02 <- read.delim("./H3K4me3/A03478.H3K4me3.GE02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE02"])
H3K4me3_T_RPKM_GE02_stat <- summarise(na.omit(gene_FetalBrain[H3K4me3_promoter_GE02$gene, ]), count = n(), lower = quantile(GE02, 0.25), middle = median(GE02), upper = quantile(GE02, 0.75), min = min(GE02), max = max(GE02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE02", Mark = "H3K4me3", Marked = T)
H3K4me3_F_RPKM_GE02_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K4me3_promoter_GE02$gene), ], count = n(), lower = quantile(GE02, 0.25), middle = median(GE02), upper = quantile(GE02, 0.75), min = min(GE02), max = max(GE02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE02", Mark = "H3K4me3", Marked = F)
H3K27me3_promoter_GE02 <- read.delim("./H3K27me3/A03480.H3K27me3.GE02.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE02"])
H3K27me3_T_RPKM_GE02_stat <- summarise(na.omit(gene_FetalBrain[H3K27me3_promoter_GE02$gene, ]), count = n(), lower = quantile(GE02, 0.25), middle = median(GE02), upper = quantile(GE02, 0.75), min = min(GE02), max = max(GE02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE02", Mark = "H3K27me3", Marked = T)
H3K27me3_F_RPKM_GE02_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K27me3_promoter_GE02$gene), ], count = n(), lower = quantile(GE02, 0.25), middle = median(GE02), upper = quantile(GE02, 0.75), min = min(GE02), max = max(GE02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE02", Mark = "H3K27me3", Marked = F)
H3K36me3_genebody_GE02 <- read.delim("./H3K36me3/A03481.H3K36me3.GE02.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE02"])
H3K36me3_T_RPKM_GE02_stat <- summarise(na.omit(gene_FetalBrain[H3K36me3_genebody_GE02$gene, ]), count = n(), lower = quantile(GE02, 0.25), middle = median(GE02), upper = quantile(GE02, 0.75), min = min(GE02), max = max(GE02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE02", Mark = "H3K36me3", Marked = T)
H3K36me3_F_RPKM_GE02_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K36me3_genebody_GE02$gene), ], count = n(), lower = quantile(GE02, 0.25), middle = median(GE02), upper = quantile(GE02, 0.75), min = min(GE02), max = max(GE02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE02", Mark = "H3K36me3", Marked = F)
H3K4me3_H3K27me3_RPKM_GE02_stat <- summarise(na.omit(gene_FetalBrain[intersect(H3K4me3_promoter_GE02$gene, H3K27me3_promoter_GE02$gene), ]), count = n(), lower = quantile(GE02, 0.25), middle = median(GE02), upper = quantile(GE02, 0.75), min = min(GE02), max = max(GE02)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE02", Mark = "H3K4me3+H3K27me3", Marked = T)
HisMod_RPKM_GE02_stat <- rbind(H3K4me3_T_RPKM_GE02_stat, H3K4me3_F_RPKM_GE02_stat, H3K27me3_T_RPKM_GE02_stat, H3K27me3_F_RPKM_GE02_stat, H3K36me3_T_RPKM_GE02_stat, H3K36me3_F_RPKM_GE02_stat, H3K4me3_H3K27me3_RPKM_GE02_stat)
# GE04
H3K4me3_promoter_GE04 <- read.delim("./H3K4me3/A19304.H3K4me3.GE04.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE04"])
H3K4me3_T_RPKM_GE04_stat <- summarise(na.omit(gene_FetalBrain[H3K4me3_promoter_GE04$gene, ]), count = n(), lower = quantile(GE04, 0.25), middle = median(GE04), upper = quantile(GE04, 0.75), min = min(GE04), max = max(GE04)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE04", Mark = "H3K4me3", Marked = T)
H3K4me3_F_RPKM_GE04_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K4me3_promoter_GE04$gene), ], count = n(), lower = quantile(GE04, 0.25), middle = median(GE04), upper = quantile(GE04, 0.75), min = min(GE04), max = max(GE04)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE04", Mark = "H3K4me3", Marked = F)
H3K27me3_promoter_GE04 <- read.delim("./H3K27me3/A19306.H3K27me3.GE04.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE04"])
H3K27me3_T_RPKM_GE04_stat <- summarise(na.omit(gene_FetalBrain[H3K27me3_promoter_GE04$gene, ]), count = n(), lower = quantile(GE04, 0.25), middle = median(GE04), upper = quantile(GE04, 0.75), min = min(GE04), max = max(GE04)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE04", Mark = "H3K27me3", Marked = T)
H3K27me3_F_RPKM_GE04_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K27me3_promoter_GE04$gene), ], count = n(), lower = quantile(GE04, 0.25), middle = median(GE04), upper = quantile(GE04, 0.75), min = min(GE04), max = max(GE04)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE04", Mark = "H3K27me3", Marked = F)
H3K36me3_genebody_GE04 <- read.delim("./H3K36me3/A19307.H3K36me3.GE04.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE04"])
H3K36me3_T_RPKM_GE04_stat <- summarise(na.omit(gene_FetalBrain[H3K36me3_genebody_GE04$gene, ]), count = n(), lower = quantile(GE04, 0.25), middle = median(GE04), upper = quantile(GE04, 0.75), min = min(GE04), max = max(GE04)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE04", Mark = "H3K36me3", Marked = T)
H3K36me3_F_RPKM_GE04_stat <- summarise(gene_FetalBrain[!(rownames(gene_FetalBrain) %in% H3K36me3_genebody_GE04$gene), ], count = n(), lower = quantile(GE04, 0.25), middle = median(GE04), upper = quantile(GE04, 0.75), min = min(GE04), max = max(GE04)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE04", Mark = "H3K36me3", Marked = F)
H3K4me3_H3K27me3_RPKM_GE04_stat <- summarise(na.omit(gene_FetalBrain[intersect(H3K4me3_promoter_GE04$gene, H3K27me3_promoter_GE04$gene), ]), count = n(), lower = quantile(GE04, 0.25), middle = median(GE04), upper = quantile(GE04, 0.75), min = min(GE04), max = max(GE04)) %>% mutate(ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)), Sample = "GE04", Mark = "H3K4me3+H3K27me3", Marked = T)
HisMod_RPKM_GE04_stat <- rbind(H3K4me3_T_RPKM_GE04_stat, H3K4me3_F_RPKM_GE04_stat, H3K27me3_T_RPKM_GE04_stat, H3K27me3_F_RPKM_GE04_stat, H3K36me3_T_RPKM_GE04_stat, H3K36me3_F_RPKM_GE04_stat, H3K4me3_H3K27me3_RPKM_GE04_stat)
HisMod_RPKM_stat <- rbind(HisMod_RPKM_Brain01_stat, HisMod_RPKM_Brain02_stat, HisMod_RPKM_Cortex01_stat, HisMod_RPKM_Cortex02_stat, HisMod_RPKM_GE01_stat, HisMod_RPKM_GE02_stat, HisMod_RPKM_GE04_stat) %>% 
  mutate(Mark = factor(Mark, levels = c("H3K4me3", "H3K27me3", "H3K4me3+H3K27me3", "H3K36me3")), Marked = factor(Marked, levels = c(T, F)))
(HisMod_RPKM_figure <- ggplot(HisMod_RPKM_stat, aes(x = Sample, fill = Marked)) + 
   geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", position = position_dodge(), outlier.shape = NA, width = 0.8) + 
   geom_hline(yintercept = 1) + 
   geom_text(aes(label = count, y = 50, color = Marked), size = 4, angle = 90, position = position_dodge(width = 1)) + 
   xlab("") + 
   ylab("") + 
   facet_wrap(~ Mark, nrow = 1, scales = "free_x") + 
   theme_bw() + 
   theme(axis.text.x = element_text(angle = 90)))
ggsave(HisMod_RPKM_figure, file = "HisMod_RPKM_figure.pdf", height = 6, width = 10)

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
                                   UMR = UMR_enhancer$UMR, percentage = UMR_enhancer$percent*100, log2_enrichment = log2(UMR_enhancer_enrich$enrich))
UMR_enhancer_summary_tall <- melt(UMR_enhancer_summary, id = c("Sample", "compare", "UMR"))
(UMR_enhancer_summary_figure <- ggplot(UMR_enhancer_summary_tall, aes(x = Sample, y = value, fill = UMR)) + 
   geom_bar(stat = "identity", position = position_dodge()) + 
   facet_grid(variable ~ compare, scales = "free") + 
   xlab("") + 
   ylab("") + 
   theme_bw())
ggsave(UMR_enhancer_summary_figure, file = "UMR_enhancer_summary_figure.pdf")

save(FindER_summary, FindER_summary_figure, HisMod_RPKM_stat, HisMod_RPKM_figure, H3K4me3_promoter_GE04, H3K27me3_promoter_GE04, H3K36me3_genebody_GE04, 
     H3K4me3_promoter_Brain01, H3K27me3_promoter_Brain01, H3K36me3_genebody_Brain01, H3K4me3_promoter_Brain02, H3K27me3_promoter_Brain02, H3K36me3_genebody_Brain02, 
     H3K27me3_promoter_Cortex01, H3K36me3_genebody_Cortex01, H3K4me3_promoter_Cortex02, H3K27me3_promoter_Cortex02, H3K36me3_genebody_Cortex02, 
     H3K27me3_promoter_GE01, H3K36me3_genebody_GE01, H3K4me3_promoter_GE02, H3K27me3_promoter_GE02, H3K36me3_genebody_GE02, 
     UMR_enhancer, UMR_enhancer_enrich, UMR_enhancer_summary_figure, 
     file = "/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FetalBrain_FindER.Rdata")
