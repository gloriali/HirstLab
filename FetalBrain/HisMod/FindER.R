# Fetal Brain - Histone modifications FindER 
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(VennDiagram)
library(gridExtra)
library(preprocessCore)
source("/home/lli/HirstLab/Pipeline/R/enrich.R")
source("/home/lli/HirstLab/Pipeline/R/enrich_GREAT.R")
load("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FetalBrain_FindER.Rdata")
split_on_last <- function(string, pattern){
  library(stringr)
  split <- unlist(str_split(string, pattern))
  first <- paste(split[1:length(split)-1], collapse = pattern)
  second <- split[length(split)]
  return(c(first, second))
}
split_on_last <- Vectorize(split_on_last, vectorize.args = "string", USE.NAMES = F)

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

## =========== Correlation with RPKM ================
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/")
colname <- c("chr", "start", "end", "id", "gene")
gene_FetalBrain <- read.delim("/home/lli/FetalBrain/RNAseq/rpkm_pc.txt", as.is = T, row.names = 1)
colnames(gene_FetalBrain) <- gsub("\\.HuFNSC", "", colnames(gene_FetalBrain))
HisMod_RPKM <- melt(gene_FetalBrain %>% select(Brain01, Brain02, Cortex01, Cortex02, GE01, GE02, GE04, id, name, coord), id = c("id", "name", "coord"), value.name = "RPKM", variable.name = "Sample") %>%
  mutate(K4_K27 = "neither", K36 = "no_H3K36me3")
### Brain01
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
### Brain02
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
### Cortex01
H3K27me3_promoter_Cortex01 <- read.delim("./H3K27me3/A03272.H3K27me3.Cortex01.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex01"])
H3K36me3_genebody_Cortex01 <- read.delim("./H3K36me3/A03273.H3K36me3.Cortex01.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "Cortex01"])
HisMod_RPKM[HisMod_RPKM$Sample == "Cortex01" & HisMod_RPKM$id %in% H3K36me3_genebody_Cortex01$gene, "K36"] <- "H3K36me3"
HisMod_RPKM[HisMod_RPKM$Sample == "Cortex01" & HisMod_RPKM$id %in% H3K27me3_promoter_Cortex01$gene, "K4_K27"] <- "H3K27me3"
### Cortex02
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
### GE01
H3K27me3_promoter_GE01 <- read.delim("./H3K27me3/A03278.H3K27me3.GE01.promoter.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE01"])
H3K36me3_genebody_GE01 <- read.delim("./H3K36me3/A03279.H3K36me3.GE01.genebody.bed", head = F, as.is = T, col.names = colname) %>% 
  mutate(type = gsub("ENSG\\d+_", "", gene), gene = gsub("_.+", "", gene), RPKM = gene_FetalBrain[gene, "GE01"])
HisMod_RPKM[HisMod_RPKM$Sample == "GE01" & HisMod_RPKM$id %in% H3K36me3_genebody_GE01$gene, "K36"] <- "H3K36me3"
HisMod_RPKM[HisMod_RPKM$Sample == "GE01" & HisMod_RPKM$id %in% H3K27me3_promoter_GE01$gene, "K4_K27"] <- "H3K27me3"
### GE02
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
### GE04
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
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/")
### Between MZ twins
Brain01_Brain02_DM_K4me3 <- rbind(H3K4me3_promoter_Brain01[H3K4me3_promoter_Brain01$gene %in% setdiff(H3K4me3_promoter_Brain01$gene, H3K4me3_promoter_Brain02$gene), ] %>% mutate(Marked = "HuFNSC01"), 
                                  H3K4me3_promoter_Brain02[H3K4me3_promoter_Brain02$gene %in% setdiff(H3K4me3_promoter_Brain02$gene, H3K4me3_promoter_Brain01$gene), ] %>% mutate(Marked = "HuFNSC02"))
Brain01_Brain02_DM_K27me3 <- rbind(H3K27me3_promoter_Brain01[H3K27me3_promoter_Brain01$gene %in% setdiff(H3K27me3_promoter_Brain01$gene, H3K27me3_promoter_Brain02$gene), ] %>% mutate(Marked = "HuFNSC01"), 
                                   H3K27me3_promoter_Brain02[H3K27me3_promoter_Brain02$gene %in% setdiff(H3K27me3_promoter_Brain02$gene, H3K27me3_promoter_Brain01$gene), ] %>% mutate(Marked = "HuFNSC02"))
Cortex01_Cortex02_DM_K27me3 <- rbind(H3K27me3_promoter_Cortex01[H3K27me3_promoter_Cortex01$gene %in% setdiff(H3K27me3_promoter_Cortex01$gene, H3K27me3_promoter_Cortex02$gene), ] %>% mutate(Marked = "HuFNSC01"), 
                                     H3K27me3_promoter_Cortex02[H3K27me3_promoter_Cortex02$gene %in% setdiff(H3K27me3_promoter_Cortex02$gene, H3K27me3_promoter_Cortex01$gene), ] %>% mutate(Marked = "HuFNSC02"))
GE01_GE02_DM_K27me3 <- rbind(H3K27me3_promoter_GE01[H3K27me3_promoter_GE01$gene %in% setdiff(H3K27me3_promoter_GE01$gene, H3K27me3_promoter_GE02$gene), ] %>% mutate(Marked = "HuFNSC01"), 
                             H3K27me3_promoter_GE02[H3K27me3_promoter_GE02$gene %in% setdiff(H3K27me3_promoter_GE02$gene, H3K27me3_promoter_GE01$gene), ] %>% mutate(Marked = "HuFNSC02"))
His_DM_MZ <- rbind(Brain01_Brain02_DM_K4me3 %>% mutate(Sample = "Brain", Mark = "H3K4me3"), 
                   Brain01_Brain02_DM_K27me3 %>% mutate(Sample = "Brain", Mark = "H3K27me3"), 
                   Cortex01_Cortex02_DM_K27me3 %>% mutate(Sample = "Cortex", Mark = "H3K27me3"), 
                   GE01_GE02_DM_K27me3 %>% mutate(Sample = "GE", Mark = "H3K27me3")) %>% mutate(Comparison = "MZ")
### Between neurospheres
Cortex01_GE01_DM_K27me3 <- rbind(H3K27me3_promoter_Cortex01[H3K27me3_promoter_Cortex01$gene %in% setdiff(H3K27me3_promoter_Cortex01$gene, H3K27me3_promoter_GE01$gene), ] %>% mutate(Marked = "Cortex"), 
                                 H3K27me3_promoter_GE01[H3K27me3_promoter_GE01$gene %in% setdiff(H3K27me3_promoter_GE01$gene, H3K27me3_promoter_Cortex01$gene), ] %>% mutate(Marked = "GE"))
Cortex02_GE02_DM_K4me3 <- rbind(H3K4me3_promoter_Cortex02[H3K4me3_promoter_Cortex02$gene %in% setdiff(H3K4me3_promoter_Cortex02$gene, H3K4me3_promoter_GE02$gene), ] %>% mutate(Marked = "Cortex"), 
                                H3K4me3_promoter_GE02[H3K4me3_promoter_GE02$gene %in% setdiff(H3K4me3_promoter_GE02$gene, H3K4me3_promoter_Cortex02$gene), ] %>% mutate(Marked = "GE"))
Cortex02_GE02_DM_K27me3 <- rbind(H3K27me3_promoter_Cortex02[H3K27me3_promoter_Cortex02$gene %in% setdiff(H3K27me3_promoter_Cortex02$gene, H3K27me3_promoter_GE02$gene), ] %>% mutate(Marked = "Cortex"), 
                                 H3K27me3_promoter_GE02[H3K27me3_promoter_GE02$gene %in% setdiff(H3K27me3_promoter_GE02$gene, H3K27me3_promoter_Cortex02$gene), ] %>% mutate(Marked = "GE"))
His_DM_neurospheres <- rbind(Cortex01_GE01_DM_K27me3 %>% mutate(Sample = "HuFNSC01", Mark = "H3K27me3"), 
                             Cortex02_GE02_DM_K4me3 %>% mutate(Sample = "HuFNSC02", Mark = "H3K4me3"), 
                             Cortex02_GE02_DM_K27me3 %>% mutate(Sample = "HuFNSC02", Mark = "H3K27me3")) %>% mutate(Comparison = "neurospheres")
### Between GW
GE01_GE04_DM_K27me3 <- rbind(H3K27me3_promoter_GE01[H3K27me3_promoter_GE01$gene %in% setdiff(H3K27me3_promoter_GE01$gene, H3K27me3_promoter_GE04$gene), ] %>% mutate(Marked = "GW17"), 
                             H3K27me3_promoter_GE04[H3K27me3_promoter_GE04$gene %in% setdiff(H3K27me3_promoter_GE04$gene, H3K27me3_promoter_GE01$gene), ] %>% mutate(Marked = "GW13"))
GE02_GE04_DM_K4me3 <- rbind(H3K4me3_promoter_GE02[H3K4me3_promoter_GE02$gene %in% setdiff(H3K4me3_promoter_GE02$gene, H3K4me3_promoter_GE04$gene), ] %>% mutate(Marked = "GW17"), 
                            H3K4me3_promoter_GE04[H3K4me3_promoter_GE04$gene %in% setdiff(H3K4me3_promoter_GE04$gene, H3K4me3_promoter_GE02$gene), ] %>% mutate(Marked = "GW13"))
GE02_GE04_DM_K27me3 <- rbind(H3K27me3_promoter_GE02[H3K27me3_promoter_GE02$gene %in% setdiff(H3K27me3_promoter_GE02$gene, H3K27me3_promoter_GE04$gene), ] %>% mutate(Marked = "GW17"), 
                             H3K27me3_promoter_GE04[H3K27me3_promoter_GE04$gene %in% setdiff(H3K27me3_promoter_GE04$gene, H3K27me3_promoter_GE02$gene), ] %>% mutate(Marked = "GW13"))
His_DM_GW <- rbind(GE01_GE04_DM_K27me3 %>% mutate(Sample = "GE01_GE04", Mark = "H3K27me3"), 
                   GE02_GE04_DM_K4me3 %>% mutate(Sample = "GE02_GE04", Mark = "H3K4me3"), 
                   GE02_GE04_DM_K27me3 %>% mutate(Sample = "GE02_GE04", Mark = "H3K27me3")) %>% mutate(Comparison = "GW")
His_DM_promoter <- rbind(His_DM_MZ, His_DM_neurospheres, His_DM_GW)
(His_DM_promoter_figure <- ggplot(His_DM_promoter, aes(x = Sample, fill = Marked)) + 
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
### Brain01 vs Brain02
Brain01_Brain02DE_epi <- brain01_brain02DE %>% mutate(K4_1 = F, K4_2 = F, K4 = F, K27_1 = F, K27_2 = F, K27 = F, DMR = F)
rownames(Brain01_Brain02DE_epi) <- Brain01_Brain02DE_epi$V1
Brain01_Brain02DE_epi[intersect(brain01_brain02DE$V1, H3K4me3_promoter_Brain01$gene), "K4_1"] <- T
Brain01_Brain02DE_epi[intersect(brain01_brain02DE$V1, H3K4me3_promoter_Brain02$gene), "K4_2"] <- T
Brain01_Brain02DE_epi[intersect(brain01_brain02DE$V1, H3K27me3_promoter_Brain01$gene), "K27_1"] <- T
Brain01_Brain02DE_epi[intersect(brain01_brain02DE$V1, H3K27me3_promoter_Brain02$gene), "K27_2"] <- T
Brain01_Brain02DE_epi[(Brain01_Brain02DE_epi$DE == "UP" & !Brain01_Brain02DE_epi$K4_1 & Brain01_Brain02DE_epi$K4_2), "K4"] <- T
Brain01_Brain02DE_epi[(Brain01_Brain02DE_epi$DE == "DN" & Brain01_Brain02DE_epi$K4_1 & !Brain01_Brain02DE_epi$K4_2), "K4"] <- T
Brain01_Brain02DE_epi[(Brain01_Brain02DE_epi$DE == "UP" & Brain01_Brain02DE_epi$K27_1 & !Brain01_Brain02DE_epi$K27_2), "K27"] <- T
Brain01_Brain02DE_epi[(Brain01_Brain02DE_epi$DE == "DN" & !Brain01_Brain02DE_epi$K27_1 & Brain01_Brain02DE_epi$K27_2), "K27"] <- T
Brain01_Brain02DE_epi[c(DMR_DE_Brain01_Brain02_hyper$DMR_gene_DE$id, DMR_DE_Brain01_Brain02_hypo$DMR_gene_DE$id), "DMR"] <- T
Brain01_Brain02DE_epi_list <- list(H3K4me3 = as.character(filter(Brain01_Brain02DE_epi, K4 == T)[, "V1"]), H3K27me3 = as.character(filter(Brain01_Brain02DE_epi, K27 == T)[, "V1"]), DMR = as.character(filter(Brain01_Brain02DE_epi, DMR == T)[, "V1"]))
venn_Brain01_Brain02DE_epi <- venn.diagram(Brain01_Brain02DE_epi_list, filename = NULL, fill = c("green", "red", "blue"), main = "Brain01 vs Brain02")
### Cortex01 vs Cortex02
Cortex01_Cortex02DE_epi <- cortex01_cortex02DE %>% mutate(K27_1 = F, K27_2 = F, K27 = F, DMR = F)
rownames(Cortex01_Cortex02DE_epi) <- Cortex01_Cortex02DE_epi$V1
Cortex01_Cortex02DE_epi[intersect(cortex01_cortex02DE$V1, H3K27me3_promoter_Cortex01$gene), "K27_1"] <- T
Cortex01_Cortex02DE_epi[intersect(cortex01_cortex02DE$V1, H3K27me3_promoter_Cortex02$gene), "K27_2"] <- T
Cortex01_Cortex02DE_epi[(Cortex01_Cortex02DE_epi$DE == "UP" & Cortex01_Cortex02DE_epi$K27_1 & !Cortex01_Cortex02DE_epi$K27_2), "K27"] <- T
Cortex01_Cortex02DE_epi[(Cortex01_Cortex02DE_epi$DE == "DN" & !Cortex01_Cortex02DE_epi$K27_1 & Cortex01_Cortex02DE_epi$K27_2), "K27"] <- T
Cortex01_Cortex02DE_epi[c(DMR_DE_Cortex01_Cortex02_hyper$DMR_gene_DE$id, DMR_DE_Cortex01_Cortex02_hypo$DMR_gene_DE$id), "DMR"] <- T
Cortex01_Cortex02DE_epi_list <- list(H3K27me3 = as.character(filter(Cortex01_Cortex02DE_epi, K27 == T)[, "V1"]), DMR = as.character(filter(Cortex01_Cortex02DE_epi, DMR == T)[, "V1"]))
venn_Cortex01_Cortex02DE_epi <- venn.diagram(Cortex01_Cortex02DE_epi_list, filename = NULL, fill = c("red", "blue"), main = "Cortex01 vs Cortex02")
### GE01 vs GE02
GE01_GE02DE_epi <- GE01_GE02DE %>% mutate(K27_1 = F, K27_2 = F, K27 = F, DMR = F)
rownames(GE01_GE02DE_epi) <- GE01_GE02DE_epi$V1
GE01_GE02DE_epi[intersect(GE01_GE02DE$V1, H3K27me3_promoter_GE01$gene), "K27_1"] <- T
GE01_GE02DE_epi[intersect(GE01_GE02DE$V1, H3K27me3_promoter_GE02$gene), "K27_2"] <- T
GE01_GE02DE_epi[(GE01_GE02DE_epi$DE == "UP" & GE01_GE02DE_epi$K27_1 & !GE01_GE02DE_epi$K27_2), "K27"] <- T
GE01_GE02DE_epi[(GE01_GE02DE_epi$DE == "DN" & !GE01_GE02DE_epi$K27_1 & GE01_GE02DE_epi$K27_2), "K27"] <- T
GE01_GE02DE_epi[c(DMR_DE_GE01_GE02_hyper$DMR_gene_DE$id, DMR_DE_GE01_GE02_hypo$DMR_gene_DE$id), "DMR"] <- T
GE01_GE02DE_epi_list <- list(H3K27me3 = as.character(filter(GE01_GE02DE_epi, K27 == T)[, "V1"]), DMR = as.character(filter(GE01_GE02DE_epi, DMR == T)[, "V1"]))
venn_GE01_GE02DE_epi <- venn.diagram(GE01_GE02DE_epi_list, filename = NULL, fill = c("red", "blue"), main = "GE01 vs GE02")
### Cortex01 vs GE01
Cortex01_GE01DE_epi <- cortex01_GE01DE %>% mutate(K27_1 = F, K27_2 = F, K27 = F, DMR = F)
rownames(Cortex01_GE01DE_epi) <- Cortex01_GE01DE_epi$V1
Cortex01_GE01DE_epi[intersect(cortex01_GE01DE$V1, H3K27me3_promoter_Cortex01$gene), "K27_1"] <- T
Cortex01_GE01DE_epi[intersect(cortex01_GE01DE$V1, H3K27me3_promoter_GE01$gene), "K27_2"] <- T
Cortex01_GE01DE_epi[(Cortex01_GE01DE_epi$DE == "UP" & Cortex01_GE01DE_epi$K27_1 & !Cortex01_GE01DE_epi$K27_2), "K27"] <- T
Cortex01_GE01DE_epi[(Cortex01_GE01DE_epi$DE == "DN" & !Cortex01_GE01DE_epi$K27_1 & Cortex01_GE01DE_epi$K27_2), "K27"] <- T
Cortex01_GE01DE_epi[c(DMR_DE_Cortex01_GE01_hyper$DMR_gene_DE$id, DMR_DE_Cortex01_GE01_hypo$DMR_gene_DE$id), "DMR"] <- T
Cortex01_GE01DE_epi_list <- list(H3K27me3 = as.character(filter(Cortex01_GE01DE_epi, K27 == T)[, "V1"]), DMR = as.character(filter(Cortex01_GE01DE_epi, DMR == T)[, "V1"]))
venn_Cortex01_GE01DE_epi <- venn.diagram(Cortex01_GE01DE_epi_list, filename = NULL, fill = c("red", "blue"), main = "Cortex01 vs GE01")
### Cortex02 vs GE02
Cortex02_GE02DE_epi <- cortex02_GE02DE %>% mutate(K4_1 = F, K4_2 = F, K4 = F, K27_1 = F, K27_2 = F, K27 = F, DMR = F)
rownames(Cortex02_GE02DE_epi) <- Cortex02_GE02DE_epi$V1
Cortex02_GE02DE_epi[intersect(cortex02_GE02DE$V1, H3K4me3_promoter_Cortex02$gene), "K4_1"] <- T
Cortex02_GE02DE_epi[intersect(cortex02_GE02DE$V1, H3K4me3_promoter_GE02$gene), "K4_2"] <- T
Cortex02_GE02DE_epi[intersect(cortex02_GE02DE$V1, H3K27me3_promoter_Cortex02$gene), "K27_1"] <- T
Cortex02_GE02DE_epi[intersect(cortex02_GE02DE$V1, H3K27me3_promoter_GE02$gene), "K27_2"] <- T
Cortex02_GE02DE_epi[(Cortex02_GE02DE_epi$DE == "UP" & !Cortex02_GE02DE_epi$K4_1 & Cortex02_GE02DE_epi$K4_2), "K4"] <- T
Cortex02_GE02DE_epi[(Cortex02_GE02DE_epi$DE == "DN" & Cortex02_GE02DE_epi$K4_1 & !Cortex02_GE02DE_epi$K4_2), "K4"] <- T
Cortex02_GE02DE_epi[(Cortex02_GE02DE_epi$DE == "UP" & Cortex02_GE02DE_epi$K27_1 & !Cortex02_GE02DE_epi$K27_2), "K27"] <- T
Cortex02_GE02DE_epi[(Cortex02_GE02DE_epi$DE == "DN" & !Cortex02_GE02DE_epi$K27_1 & Cortex02_GE02DE_epi$K27_2), "K27"] <- T
Cortex02_GE02DE_epi[c(DMR_DE_Cortex02_GE02_hyper$DMR_gene_DE$id, DMR_DE_Cortex02_GE02_hypo$DMR_gene_DE$id), "DMR"] <- T
Cortex02_GE02DE_epi_list <- list(H3K4me3 = as.character(filter(Cortex02_GE02DE_epi, K4 == T)[, "V1"]), H3K27me3 = as.character(filter(Cortex02_GE02DE_epi, K27 == T)[, "V1"]), DMR = as.character(filter(Cortex02_GE02DE_epi, DMR == T)[, "V1"]))
venn_Cortex02_GE02DE_epi <- venn.diagram(Cortex02_GE02DE_epi_list, filename = NULL, fill = c("green", "red", "blue"), main = "Cortex02 vs GE02")
### GE02 vs GE04
GE02_GE04DE_epi <- GE02_GE04DE %>% mutate(K4_1 = F, K4_2 = F, K4 = F, K27_1 = F, K27_2 = F, K27 = F, DMR = F)
rownames(GE02_GE04DE_epi) <- GE02_GE04DE_epi$ID
GE02_GE04DE_epi[intersect(GE02_GE04DE$ID, H3K4me3_promoter_GE02$gene), "K4_1"] <- T
GE02_GE04DE_epi[intersect(GE02_GE04DE$ID, H3K4me3_promoter_GE04$gene), "K4_2"] <- T
GE02_GE04DE_epi[intersect(GE02_GE04DE$ID, H3K27me3_promoter_GE02$gene), "K27_1"] <- T
GE02_GE04DE_epi[intersect(GE02_GE04DE$ID, H3K27me3_promoter_GE04$gene), "K27_2"] <- T
GE02_GE04DE_epi[(GE02_GE04DE_epi$DE == "UP" & !GE02_GE04DE_epi$K4_1 & GE02_GE04DE_epi$K4_2), "K4"] <- T
GE02_GE04DE_epi[(GE02_GE04DE_epi$DE == "DN" & GE02_GE04DE_epi$K4_1 & !GE02_GE04DE_epi$K4_2), "K4"] <- T
GE02_GE04DE_epi[(GE02_GE04DE_epi$DE == "UP" & GE02_GE04DE_epi$K27_1 & !GE02_GE04DE_epi$K27_2), "K27"] <- T
GE02_GE04DE_epi[(GE02_GE04DE_epi$DE == "DN" & !GE02_GE04DE_epi$K27_1 & GE02_GE04DE_epi$K27_2), "K27"] <- T
GE02_GE04DE_epi[c(DMR_DE_GE02_GE04_hyper$DMR_gene_DE$id, DMR_DE_GE02_GE04_hypo$DMR_gene_DE$id), "DMR"] <- T
GE02_GE04DE_epi_list <- list(H3K4me3 = as.character(filter(GE02_GE04DE_epi, K4 == T)[, "ID"]), H3K27me3 = as.character(filter(GE02_GE04DE_epi, K27 == T)[, "ID"]), DMR = as.character(filter(GE02_GE04DE_epi, DMR == T)[, "ID"]))
venn_GE02_GE04DE_epi <- venn.diagram(GE02_GE04DE_epi_list, filename = NULL, fill = c("green", "red", "blue"), main = "GE02 vs GE04")
grid.arrange(gTree(children = venn_Brain01_Brain02DE_epi), gTree(children = venn_Cortex01_Cortex02DE_epi), gTree(children = venn_GE01_GE02DE_epi), 
             gTree(children = venn_Cortex01_GE01DE_epi), gTree(children = venn_Cortex02_GE02DE_epi), gTree(children = venn_GE02_GE04DE_epi), nrow = 2)
pdf("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/venn_DE_epi.pdf")
grid.arrange(gTree(children = venn_Brain01_Brain02DE_epi), gTree(children = venn_Cortex01_Cortex02DE_epi), gTree(children = venn_GE01_GE02DE_epi), 
             gTree(children = venn_Cortex01_GE01DE_epi), gTree(children = venn_Cortex02_GE02DE_epi), gTree(children = venn_GE02_GE04DE_epi), nrow = 2)
dev.off()

## ============= core enhancers =============
### closest genes
# setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/closest_gene/")
# colname <- c("chr", "start", "end", "id", "gene", "dis")
# closest_gene_Cortex01 <- read.delim("A03269.H3K4me1.Cortex01.closest.gene.pc", head = F, as.is = T, col.names = colname)
# closest_gene_io_Cortex01 <- read.delim("A03269.H3K4me1.Cortex01.closest.gene.pc.io", head = F, as.is = T, col.names = colname)
# length(unique(closest_gene_Cortex01$gene))     # 14660
# length(unique(closest_gene_io_Cortex01$gene))  # 14630
# length(intersect(unique(closest_gene_Cortex01$gene), unique(closest_gene_io_Cortex01$gene)))   # 11978
# ignore overlap doesn't have much effect on closest genes, use including overlap version. 
# closest genes included more than half of all pc genes, no point in enrichment analysis.  
### GWAS
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/core/")
GWAS_ref <- read.delim("/projects/epigenomics/resources/UCSC_hg19/GWAS/gwasCatalog_July2014.table", as.is = T)
colnames <- c("enhancerChr", "enhancerStart", "enhancerEnd", "enhancerID", "gwasChr", "gwasStart", "gwasEnd", "gwasID", "trait", "genes")
GWAS_core_enhancer <- read.delim("./core_enhancers.GWAS.txt", head = F, col.names = colnames) 
GWAS_core_enhancer_TFBS <- read.delim("./core_enhancers.GWAS.TFBS.txt", head = F)[, 1:10]
colnames(GWAS_core_enhancer_TFBS) <- c("gwasChr", "gwasStart", "gwasEnd", "gwasID", "trait", "enhancerID", "TFBSchr", "TFBSstart", "TFBSend", "TF")
GWAS_core_enhancer_TFBS_brain <- GWAS_core_enhancer_TFBS %>%
  filter(as.character(trait) %in% c("Schizophrenia", "Response to antipsychotic treatment in schizophrenia (reasoning)", "Parkinson's disease", "Attention deficit hyperactivity disorder and conduct disorder", "Response to antipsychotic treatment in schizophrenia (working memory)", "Alzheimer's disease", "Bipolar disorder and schizophrenia", "Depression (quantitative trait)", "Optic disc parameters", "Optic nerve measurement (disc area)", "Optic nerve measurement (cup area)", 
                                    "Bipolar disorder", "Attention deficit hyperactivity disorder", "Migraine", "Alzheimer's disease (cognitive decline)", "Glaucoma (primary open-angle)", "Response to antipsychotic treatment", "Brain structure", "Glioma", "Myopia (pathological)", "Autism spectrum disorder, attention deficit-hyperactivity disorder, bipolar disorder, major depressive disorder, and schizophrenia (combined)", "Migraine without aura", "Migraine - clinic-based", "Central corneal thickness", 
                                    "Cognitive performance", "Non-word repetition", "Word reading", "Corneal astigmatism", "Corneal curvature", "Corneal structure", "Glaucoma (exfoliation)", "Behavioural disinhibition (generation interaction)", "Epilepsy (generalized)", "White matter hyperintensity burden", "Alzheimer's disease biomarkers", "Bipolar disorder (mood-incongruent)", "Attention deficit hyperactivity disorder (hyperactivity-impulsivity symptoms)", "Attention deficit hyperactivity disorder (combined symptoms)", 
                                    "Neuroblastoma (high-risk)", "Brain imaging in schizophrenia (interaction)", "Hippocampal atrophy", "Intelligence", "Conduct disorder (interaction)", "Response to antidepressant treatment", "Conduct disorder", "Autism", "Normalized brain volume", "Neuroblastoma", "Brain lesion load", "Alzheimer's disease (late onset)", "Psychosis (atypical)", "Personality dimensions", "Schizophrenia, bipolar disorder and depression (combined)", "Hypersomnia (HLA-DQB1*06:02 negative)"))
GWAS_core_enhancer_trait <- GWAS_core_enhancer %>% group_by(trait) %>% 
  summarise(Nsample_trait = n()) %>%                                                         # No. of records of the specific trait in the sample
  mutate(Nsample = nrow(GWAS_core_enhancer),                                                 # No. of records in the sample
         Ntrait = sapply(trait, function(x){return(nrow(GWAS_ref[GWAS_ref$trait == x, ]))}), # No. of records of the specific trait in the reference
         Ntotal = nrow(GWAS_ref),                                                            # total No. of records in the reference
         phyper = 1 - phyper(Nsample_trait, Nsample, Ntotal-Nsample, Ntrait), 
         FDR = p.adjust(phyper, method = "fdr")) %>% arrange(FDR)
GWAS_core_enhancer_trait_sig <- droplevels(GWAS_core_enhancer_trait[GWAS_core_enhancer_trait$FDR < 0.01, ])
GWAS_core_enhancer_trait_sig_brain <- GWAS_core_enhancer %>% filter(trait %in% c("Cortical structure", "Glaucoma (exfoliation)", "Glioblastoma", "Neuranatomic and neurocognitive phenotypes", "Neuroblastoma (high-risk)", "Odorant perception (isobutyraldehyde)", "Schizophrenia (cytomegalovirus infection interaction)", "Alzheimer's disease biomarkers")) %>% 
  arrange(trait)
write.table(GWAS_core_enhancer_TFBS_brain, file = "./GWAS_core_enhancer_TFBS_brain.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(GWAS_core_enhancer_trait, file = "./GWAS_core_enhancer_trait.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(GWAS_core_enhancer_trait_sig, file = "./GWAS_core_enhancer_trait_sig.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(GWAS_core_enhancer_trait_sig_brain, file = "./GWAS_core_enhancer_trait_sig_brain.txt", sep = "\t", col.names = T, row.names = F, quote = F)
### homer
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/core/homer/")
homer_core_enhancer <- read.delim("./knownResults.txt", head = T, as.is = T) %>% 
  mutate(Motif.Name = gsub("/Homer", "", Motif.Name), Motif.Name = gsub(" ", "_", Motif.Name), motif.file = paste0("./knownResults/known", 1:n(), ".motif"), 
         TF = gsub("\\(.+\\)", "", split_on_last(Motif.Name, "/")[1,]), Assay = split_on_last(Motif.Name, "/")[2,], 
         Percent.of.Target.Sequences.with.Motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif)), 
         Percent.of.Background.Sequences.with.Motif = as.numeric(gsub("%", "", X..of.Background.Sequences.with.Motif)), 
         FC = Percent.of.Target.Sequences.with.Motif/Percent.of.Background.Sequences.with.Motif) %>% 
  filter(q.value..Benjamini. < 0.01) %>% arrange(- Percent.of.Target.Sequences.with.Motif) %>%
  mutate(Motif.Name = factor(Motif.Name, levels = rev(Motif.Name)))
homer_core_enhancer_top <- filter(homer_core_enhancer, Percent.of.Target.Sequences.with.Motif >= 20)
write.table(homer_core_enhancer_top, file = "homer_core_enhancer_top.txt", sep = "\t", quote = F, row.names = F, col.names = T)
(homer_core_enhancer_figure <- ggplot(data = homer_core_enhancer_top, aes(Motif.Name, Percent.of.Target.Sequences.with.Motif)) +
  geom_bar(stat = "identity", width = .5, fill = "blue") + 
  coord_flip() + 
  ylab("Percent of Enhancers with Motif") + 
  xlab("") + 
  theme_bw())
ggsave(homer_core_enhancer_figure, file = "homer_core_enhancer_figure.pdf")
### overlap with WGBS UMRs
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/core/UMR/")
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

## ============= unique enhancers =============
### Intersections
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique")
unique_enhancers_summary <- read.delim("unique_enhancer.summary", as.is = T)
venn_unique_enhancer_HuFNSC01 <- draw.pairwise.venn(as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.MZ.Cortex01.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.MZ.GE01.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.MZ.HuFNSC01.bed", intern = T))), category = c("Cortex", "GE"), fill = c("red", "blue"), main = "MZ unique enhancers in HuFNSC01")
venn_unique_enhancer_HuFNSC02 <- draw.pairwise.venn(as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.MZ.Cortex02.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.MZ.GE02.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.MZ.HuFNSC02.bed", intern = T))), category = c("Cortex", "GE"), fill = c("red", "blue"), main = "MZ unique enhancers in HuFNSC02")
venn_unique_enhancer_Cortex <- draw.pairwise.venn(as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.Neurospheres.Cortex01.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.Neurospheres.Cortex02.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.Neurospheres.Cortex.bed", intern = T))), category = c("HuFNSC01", "HuFNSC02"), fill = c("red", "blue"), main = "Neurospheres unique enhancers in Cortex")
venn_unique_enhancer_GE <- draw.pairwise.venn(as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.Neurospheres.GE01.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.Neurospheres.GE02.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.Neurospheres.GE.bed", intern = T))), category = c("HuFNSC01", "HuFNSC02"), fill = c("red", "blue"), main = "Neurospheres unique enhancers in GE")
venn_unique_enhancer_GW13 <- draw.pairwise.venn(as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.GW.GW13_01.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.GW.GW13_02.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.GW.GW13.bed", intern = T))), category = c("HuFNSC01_04", "HuFNSC02_04"), fill = c("red", "blue"), main = "GW unique enhancers in GW13")
venn_unique_enhancer_GW17 <- draw.pairwise.venn(as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.GW.GW17_01.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.GW.GW17_02.bed", intern = T))), as.numeric(gsub(" .*", "", system("wc -l unique_enhancer.GW.GW17.bed", intern = T))), category = c("HuFNSC01_04", "HuFNSC02_04"), fill = c("red", "blue"), main = "GW unique enhancers in GW17")
pdf("venn_unique_enhancer.pdf", height = 8, width = 8)
grid.arrange(gTree(children = venn_unique_enhancer_HuFNSC01), gTree(children = venn_unique_enhancer_HuFNSC02),  
             gTree(children = venn_unique_enhancer_Cortex), gTree(children = venn_unique_enhancer_GE), 
             gTree(children = venn_unique_enhancer_GW13), gTree(children = venn_unique_enhancer_GW17), nrow = 3)
grid.text("MZ unique enhancers\n in HuFNSC01", 0.25, 0.95)
grid.text("MZ unique enhancers\n in HuFNSC02", 0.75, 0.95)
grid.text("Neurospheres unique enhancers\n in Cortex", 0.25, 0.65)
grid.text("Neurospheres unique enhancers\n in GE", 0.75, 0.65)
grid.text("GW unique enhancers\n in GW13", 0.25, 0.35)
grid.text("GW unique enhancers\n in GW17", 0.75, 0.35)
dev.off()
### Functional analysis - GREAT
(GREAT_unique_enhancer_HuFNSC01 <- enrich_GREAT(file = "unique_enhancer.MZ.HuFNSC01", name = "unique_enhancer.MZ.HuFNSC01", height = 12, top = 10))
(GREAT_unique_enhancer_HuFNSC02 <- enrich_GREAT(file = "unique_enhancer.MZ.HuFNSC02", name = "unique_enhancer.MZ.HuFNSC02", height = 8))
(GREAT_unique_enhancer_Cortex <- enrich_GREAT(file = "unique_enhancer.Neurospheres.Cortex", name = "unique_enhancer.Neurospheres.Cortex"))
(GREAT_unique_enhancer_GE <- enrich_GREAT(file = "unique_enhancer.Neurospheres.GE", name = "unique_enhancer.Neurospheres.GE", height = 4))
(GREAT_unique_enhancer_GW13 <- enrich_GREAT(file = "unique_enhancer.GW.GW13", name = "unique_enhancer.GW.GW13", height = 4))
(GREAT_unique_enhancer_GW17 <- enrich_GREAT(file = "unique_enhancer.GW.GW17", name = "unique_enhancer.GW.GW17", height = 4))
### GWAS
GWAS_ref <- read.delim("/projects/epigenomics/resources/UCSC_hg19/GWAS/gwasCatalog_July2014.table", as.is = T)
colnames <- c("enhancerChr", "enhancerStart", "enhancerEnd", "enhancerID", "gwasChr", "gwasStart", "gwasEnd", "gwasID", "trait", "genes")
GWAS_unique_enhancer_Brain01 <- read.delim("./GWAS/unique_enhancer.MZ.Brain01.GWAS.bed", head = F, col.names = colnames) %>% mutate(comparison = "MZ", Sample = "Brain01")
GWAS_unique_enhancer_Brain02 <- read.delim("./GWAS/unique_enhancer.MZ.Brain02.GWAS.bed", head = F, col.names = colnames) %>% mutate(comparison = "MZ", Sample = "Brain02")
GWAS_unique_enhancer_Cortex01 <- read.delim("./GWAS/unique_enhancer.MZ.Cortex01.GWAS.bed", head = F, col.names = colnames) %>% mutate(comparison = "MZ", Sample = "Cortex01")
GWAS_unique_enhancer_Cortex02 <- read.delim("./GWAS/unique_enhancer.MZ.Cortex02.GWAS.bed", head = F, col.names = colnames) %>% mutate(comparison = "MZ", Sample = "Cortex02")
GWAS_unique_enhancer_GE01 <- read.delim("./GWAS/unique_enhancer.MZ.GE01.GWAS.bed", head = F, col.names = colnames) %>% mutate(comparison = "MZ", Sample = "GE01")
GWAS_unique_enhancer_GE02 <- read.delim("./GWAS/unique_enhancer.MZ.GE02.GWAS.bed", head = F, col.names = colnames) %>% mutate(comparison = "MZ", Sample = "GE02")
GWAS_unique_enhancer_Cortex <- read.delim("./GWAS/unique_enhancer.Neurospheres.Cortex.GWAS.bed", head = F, col.names = colnames) %>% mutate(comparison = "Neurospheres", Sample = "Cortex")
GWAS_unique_enhancer_GE <- read.delim("./GWAS/unique_enhancer.Neurospheres.GE.GWAS.bed", head = F, col.names = colnames) %>% mutate(comparison = "Neurospheres", Sample = "GE")
GWAS_unique_enhancer_GW13 <- read.delim("./GWAS/unique_enhancer.GW.GW13.GWAS.bed", head = F, col.names = colnames) %>% mutate(comparison = "GW", Sample = "GW13")
GWAS_unique_enhancer_GW17 <- read.delim("./GWAS/unique_enhancer.GW.GW17.GWAS.bed", head = F, col.names = colnames) %>% mutate(comparison = "GW", Sample = "GW17")
GWAS_unique_enhancer_all <- rbind(GWAS_unique_enhancer_Brain01, GWAS_unique_enhancer_Brain02, GWAS_unique_enhancer_Cortex01, GWAS_unique_enhancer_Cortex02, GWAS_unique_enhancer_GE01, GWAS_unique_enhancer_GE02, 
                                  GWAS_unique_enhancer_Cortex, GWAS_unique_enhancer_GE, GWAS_unique_enhancer_GW13, GWAS_unique_enhancer_GW17)
GWAS_unique_enhancer_trait <- GWAS_unique_enhancer_all %>% group_by(comparison, Sample, trait) %>% summarise(Nsample_trait = n()) %>% 
  mutate(Nsample = sapply(Sample, function(x){return(nrow(GWAS_unique_enhancer_all[GWAS_unique_enhancer_all$Sample == x, ]))}), 
         Ntrait = sapply(trait, function(x){return(nrow(GWAS_ref[GWAS_ref$trait == x, ]))}), 
         Ntotal = nrow(GWAS_ref), 
         phyper = 1 - phyper(Nsample_trait, Nsample, Ntotal-Nsample, Ntrait), 
         FDR = p.adjust(phyper, method = "fdr")) %>% arrange(comparison, Sample, FDR)
GWAS_unique_enhancer_trait_sig <- droplevels(GWAS_unique_enhancer_trait[GWAS_unique_enhancer_trait$FDR < 0.01, ])
GWAS_unique_enhancer_trait_sig_brain <- GWAS_unique_enhancer_all %>% 
  filter(trait %in% c("Gliomas", "Alzheimer's disease (cognitive decline)", "Response to cholinesterase inhibitors in Alzheimer's disease", "Autism", "Schizophrenia or bipolar disorder", "Attention deficit hyperactivity disorder", 
                      "Cognitive test performance", "Information processing speed", "Corneal astigmatism", "Hippocampal atrophy", "Central corneal thickness", "Neuroblastoma (high-risk)", "Psychosis (methamphetamine induced)", 
                      "Hearing impairment", "Hippocampal volume", "Optic nerve measurement (cup-to-disc ratio)", "Cerebrospinal AB1-42 levels", "Cerebrospinal P-tau181p levels", "Optic nerve measurement (cup area)", "Myopia (pathological)", 
                      "Cognitive performance", "Migraine - clinic-based", "Intelligence (childhood)", "Inattentive symptoms", "Behavioural disinhibition (generation interaction)", "Normalized brain volume", "Corneal curvature", "Intelligence"))
write.table(GWAS_unique_enhancer_all, file = "./GWAS/GWAS_unique_enhancer_all.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(GWAS_unique_enhancer_all[, c(5:9, 4, 11:12)], file = "./GWAS/GWAS_unique_enhancer_all.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(GWAS_unique_enhancer_trait, file = "./GWAS/GWAS_unique_enhancer_trait.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(GWAS_unique_enhancer_trait_sig, file = "./GWAS/GWAS_unique_enhancer_trait_sig.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(GWAS_unique_enhancer_trait_sig_brain, file = "./GWAS/GWAS_unique_enhancer_trait_sig_brain.txt", sep = "\t", col.names = T, row.names = F, quote = F)
GWAS_unique_enhancer_all_TFBS <- read.delim("./GWAS/GWAS_unique_enhancer_all_TFBS.txt", head = F)[, 1:12]
colnames(GWAS_unique_enhancer_all_TFBS) <- c("gwasChr", "gwasStart", "gwasEnd", "gwasID", "trait", "enhancerID", "comparison", "sample", "TFBSchr", "TFBSstart", "TFBSend", "TF")
GWAS_unique_enhancer_brain_TFBS <- GWAS_unique_enhancer_all_TFBS %>% 
  filter(as.character(trait) %in% c("Amyotrophic lateral sclerosis (age of onset)", "Response to antipsychotic treatment in schizophrenia (reasoning)", "Bipolar disorder", "Schizophrenia", "Intraocular pressure", "Attention deficit hyperactivity disorder (combined symptoms)", "Neuroblastoma (high-risk)", "Attention deficit hyperactivity disorder", "Brain lesion load", "Personality dimensions", "Cognitive performance", "Bipolar disorder and schizophrenia", 
                                    "White matter hyperintensity burden", "Major depressive disorder", "Autism spectrum disorder, attention deficit-hyperactivity disorder, bipolar disorder, major depressive disorder, and schizophrenia (combined)", "Alzheimer's disease biomarkers", "Intelligence", "Attention deficit hyperactivity disorder symptoms (interaction)", "Intelligence (childhood)", "Word reading", "Parkinson's disease (age of onset)", "Glaucoma", 
                                    "Alzheimer's disease (late onset)", "Parkinson's disease", "Alzheimer's disease (cognitive decline)", "Brain structure", "Behavioural disinhibition (generation interaction)", "Panic disorder", "Alzheimer's disease", "Corneal structure", "Central corneal thickness"))
write.table(GWAS_unique_enhancer_brain_TFBS, file = "./GWAS/GWAS_unique_enhancer_brain_TFBS.txt", sep = "\t", col.names = T, row.names = F, quote = F)
### homer
load("/home/lli/FetalBrain/RNAseq/DEfine/gene/FetalBrain_DEgenes.Rdata")
load("/projects/epigenomics/users/lli/FetalBrain/GW/GW.Rdata")
load("/home/lli/hg19/hg19v65_genes.Rdata")
#### GW
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/homer/")
DE_GW13_GW17_cortex <- data.frame(ID = intersect(cortex01_cortex04DE$ID, cortex02_cortex04DE$ID)) %>%
  mutate(cortex01_cortex04 = cortex01_cortex04DE[as.character(ID), "DE"], cortex02_cortex04 = cortex02_cortex04DE[as.character(ID), "DE"], 
         cortex01 = cortex01_cortex04DE[as.character(ID), "rpkm1"], cortex02 = cortex02_cortex04DE[as.character(ID), "rpkm1"], cortex04 = cortex01_cortex04DE[as.character(ID), "rpkm2"], 
         GW13 = cortex04, GW17 = (cortex01 + cortex02)/2, log2FC = log2(GW13/GW17 + 1e-6))
DE_GW13_GW17_GE <- data.frame(ID = intersect(GE01_GE04DE$ID, GE02_GE04DE$ID)) %>%
  mutate(GE01_GE04 = GE01_GE04DE[as.character(ID), "DE"], GE02_GE04 = GE02_GE04DE[as.character(ID), "DE"], 
         GE01 = GE01_GE04DE[as.character(ID), "rpkm1"], GE02 = GE02_GE04DE[as.character(ID), "rpkm1"], GE04 = GE01_GE04DE[as.character(ID), "rpkm2"], 
         GW13 = GE04, GW17 = (GE01 + GE02)/2, log2FC = log2(GW13/GW17 + 1e-6))
homer_unique_enhancer_GW_GW13 <- read.delim("./GW.GW13/knownResults.txt", head = T, as.is = T) %>% 
  mutate(Motif.Name = gsub("/Homer", "", Motif.Name), Motif.Name = gsub(" ", "_", Motif.Name), motif.file = paste0("./knownResults/known", 1:n(), ".motif"), 
         TF = gsub("\\(.+\\)", "", split_on_last(Motif.Name, "/")[1,]), Assay = split_on_last(Motif.Name, "/")[2,], 
         Percent.of.Target.Sequences.with.Motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif)), 
         Percent.of.Background.Sequences.with.Motif = as.numeric(gsub("%", "", X..of.Background.Sequences.with.Motif)), 
         FC = Percent.of.Target.Sequences.with.Motif/Percent.of.Background.Sequences.with.Motif) %>% 
  filter(q.value..Benjamini. < 0.01) %>% 
  mutate(Motif.Name = factor(Motif.Name, levels = rev(Motif.Name)))
homer_unique_enhancer_GW_GW17 <- read.delim("./GW.GW17/knownResults.txt", head = T, as.is = T) %>% 
  mutate(Motif.Name = gsub("/Homer", "", Motif.Name), Motif.Name = gsub(" ", "_", Motif.Name), motif.file = paste0("./knownResults/known", 1:n(), ".motif"), 
         TF = gsub("\\(.+\\)", "", split_on_last(Motif.Name, "/")[1,]), Assay = split_on_last(Motif.Name, "/")[2,], 
         Percent.of.Target.Sequences.with.Motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif)), 
         Percent.of.Background.Sequences.with.Motif = as.numeric(gsub("%", "", X..of.Background.Sequences.with.Motif)), 
         FC = Percent.of.Target.Sequences.with.Motif/Percent.of.Background.Sequences.with.Motif) %>% 
  filter(q.value..Benjamini. < 0.01) %>% 
  mutate(Motif.Name = factor(Motif.Name, levels = rev(Motif.Name)))
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/homer/GW/")
##### Common TFs
homer_unique_enhancer_GW_common <- merge(homer_unique_enhancer_GW_GW13, homer_unique_enhancer_GW_GW17, by = "Motif.Name") %>% 
  mutate(Consensus = Consensus.x, TF = TF.x, Assay = Assay.x, 
         GW13_logP = Log.P.value.x, GW13_q = q.value..Benjamini..x, GW13_motif_file = motif.file.x, GW13_Percent_Target = Percent.of.Target.Sequences.with.Motif.x, GW13_FC = FC.x, 
         GW17_logP = Log.P.value.y, GW17_q = q.value..Benjamini..y, GW17_motif_file = motif.file.y, GW17_Percent_Target = Percent.of.Target.Sequences.with.Motif.y, GW17_FC = FC.y, delta_Percent_Target = abs(GW13_Percent_Target - GW17_Percent_Target)) %>%
  select(Motif.Name, Consensus, TF, Assay, GW13_logP, GW13_q, GW13_motif_file, GW13_Percent_Target, GW13_FC, GW17_logP, GW17_q, GW17_motif_file, GW17_Percent_Target, GW17_FC, delta_Percent_Target) %>% 
  arrange(- delta_Percent_Target)
write.table(homer_unique_enhancer_GW_common, file = "homer_unique_enhancer_GW_common.txt", sep = "\t", col.names = T, row.names = F, quote = F)
homer_unique_enhancer_GW_common_top <- homer_unique_enhancer_GW_common %>% filter(GW13_Percent_Target >= 20 | GW17_Percent_Target >= 20) %>% arrange(GW13_Percent_Target + GW17_Percent_Target) %>% mutate(Motif.Name = factor(Motif.Name, levels = as.character(Motif.Name)))
homer_unique_enhancer_GW_common_top <- data.frame(Motif.Name = rep(homer_unique_enhancer_GW_common_top$Motif.Name, 2), GW = rep(c("GW13", "GW17"), each = nrow(homer_unique_enhancer_GW_common_top)), Percent_Target = c(homer_unique_enhancer_GW_common_top$GW13_Percent_Target, homer_unique_enhancer_GW_common_top$GW17_Percent_Target))
(homer_unique_enhancer_GW_common_figure <- ggplot(data = homer_unique_enhancer_GW_common_top, aes(Motif.Name, Percent_Target)) +
   geom_bar(aes(fill = GW), stat = "identity", width = .5, position = position_dodge()) + 
   coord_flip() + 
   scale_fill_manual(values = c("red", "blue"), name = "") + 
   ylab("Percent of Enhancers with Motif") + 
   xlab("") + 
   theme_bw())
ggsave(homer_unique_enhancer_GW_common_figure, file = "homer_unique_enhancer_GW_common_figure.pdf")
homer_unique_enhancer_GW_GW13_common_Sox3 <- read.delim("Common_GW13_Sox3.annotate", head = T, as.is = T) %>% filter(Gene.Type == "protein-coding", Nearest.Ensembl != "") %>% distinct(Nearest.Ensembl)
row.names(homer_unique_enhancer_GW_GW13_common_Sox3) <- homer_unique_enhancer_GW_GW13_common_Sox3$Nearest.Ensembl
homer_unique_enhancer_GW_GW17_common_Sox3 <- read.delim("Common_GW17_Sox3.annotate", head = T, as.is = T) %>% filter(Gene.Type == "protein-coding", Nearest.Ensembl != "") %>% distinct(Nearest.Ensembl)
row.names(homer_unique_enhancer_GW_GW17_common_Sox3) <- homer_unique_enhancer_GW_GW17_common_Sox3$Nearest.Ensembl
homer_unique_enhancer_GW_common_Sox3_list <- list(GW13 = homer_unique_enhancer_GW_GW13_common_Sox3$Nearest.Ensembl, GW17 = homer_unique_enhancer_GW_GW17_common_Sox3$Nearest.Ensembl, DE = c(as.character(DE_GW13_GW17_cortex$ID), as.character(DE_GW13_GW17_GE$ID)))
Venn_homer_unique_enhancer_GW_common_Sox3 <- venn.diagram(homer_unique_enhancer_GW_common_Sox3_list, filename = NULL, fill = c("red", "blue", "green"), main = "GW common Sox3", force.unique = T)
plot.new()
grid.draw(Venn_homer_unique_enhancer_GW_common_Sox3)
homer_unique_enhancer_GW_common_Sox3 <- data.frame(ID = unique(c(homer_unique_enhancer_GW_GW13_common_Sox3$Nearest.Ensembl, homer_unique_enhancer_GW_GW17_common_Sox3$Nearest.Ensembl))) %>%
  mutate(GW13 = homer_unique_enhancer_GW_GW13_common_Sox3[as.character(ID), 1], GW17 = homer_unique_enhancer_GW_GW17_common_Sox3[as.character(ID), 1])
homer_unique_enhancer_GW_common_Sox3_DE_cortex <- merge(homer_unique_enhancer_GW_common_Sox3, DE_GW13_GW17_cortex, by = "ID") %>% mutate(name = ensembl[as.character(ID), "name"], description = ensembl[as.character(ID), "description"]) %>% filter(GW13.y + GW13.y >= 0.2) %>% arrange(-abs(log2FC))
homer_unique_enhancer_GW_common_Sox3_DE_GE <- merge(homer_unique_enhancer_GW_common_Sox3, DE_GW13_GW17_GE, by = "ID") %>% mutate(name = ensembl[as.character(ID), "name"], description = ensembl[as.character(ID), "description"]) %>% filter(GW13.y + GW13.y >= 0.2) %>% arrange(-abs(log2FC))
write.table(homer_unique_enhancer_GW_common_Sox3_DE_cortex, file = "homer_unique_enhancer_GW_common_Sox3_DE_cortex.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(homer_unique_enhancer_GW_common_Sox3_DE_GE, file = "homer_unique_enhancer_GW_common_Sox3_DE_GE.txt", sep = "\t", col.names = T, row.names = F, quote = F)
##### GW13 TFs
homer_unique_enhancer_GW_GW13only <- homer_unique_enhancer_GW_GW13[!(homer_unique_enhancer_GW_GW13$Motif.Name %in% homer_unique_enhancer_GW_GW17$Motif.Name), ] %>%
  select(Motif.Name, TF, Assay, motif.file, Consensus, Log.P.value, q.value..Benjamini., Percent.of.Target.Sequences.with.Motif, Percent.of.Background.Sequences.with.Motif, FC) %>% 
  arrange(- Percent.of.Target.Sequences.with.Motif) %>% mutate(Motif.Name = factor(Motif.Name, levels = rev(as.character(Motif.Name))))
write.table(homer_unique_enhancer_GW_GW13only, file = "homer_unique_enhancer_GW_GW13only.txt", sep = "\t", col.names = T, row.names = F, quote = F)
(homer_unique_enhancer_GW_GW13only_figure <- ggplot(data = homer_unique_enhancer_GW_GW13only, aes(Motif.Name, Percent.of.Target.Sequences.with.Motif)) +
   geom_bar(fill = "blue", stat = "identity", width = .5) + 
   coord_flip() + 
   ylab("Percent of Enhancers with Motif") + 
   xlab("") + 
   theme_bw())
ggsave(homer_unique_enhancer_GW_GW13only_figure, file = "homer_unique_enhancer_GW_GW13only_figure.pdf")
##### GW17 TFs
homer_unique_enhancer_GW_GW17only <- homer_unique_enhancer_GW_GW17[!(homer_unique_enhancer_GW_GW17$Motif.Name %in% homer_unique_enhancer_GW_GW13$Motif.Name), ] %>%
  select(Motif.Name, TF, Assay, motif.file, Consensus, Log.P.value, q.value..Benjamini., Percent.of.Target.Sequences.with.Motif, Percent.of.Background.Sequences.with.Motif, FC) %>% 
  arrange(- Percent.of.Target.Sequences.with.Motif) %>% mutate(Motif.Name = factor(Motif.Name, levels = rev(as.character(Motif.Name))))
write.table(homer_unique_enhancer_GW_GW17only, file = "homer_unique_enhancer_GW_GW17only.txt", sep = "\t", col.names = T, row.names = F, quote = F)
(homer_unique_enhancer_GW_GW17only_figure <- ggplot(data = filter(homer_unique_enhancer_GW_GW17only, Percent.of.Target.Sequences.with.Motif >= 20), aes(Motif.Name, Percent.of.Target.Sequences.with.Motif)) +
   geom_bar(fill = "blue", stat = "identity", width = .5) + 
   coord_flip() + 
   ylab("Percent of Enhancers with Motif") + 
   xlab("") + 
   theme_bw())
ggsave(homer_unique_enhancer_GW_GW17only_figure, file = "homer_unique_enhancer_GW_GW17only_figure.pdf", height = 4)
homer_unique_enhancer_GW_GW17only_Olig2 <- read.delim("GW17only_Olig2.annotate", head = T, as.is = T) %>% filter(Gene.Type == "protein-coding", Nearest.Ensembl != "") %>% distinct(Nearest.Ensembl) %>%
  mutate(ID = Nearest.Ensembl, Enhancer = paste0(Chr, ":", Start, "-", End)) %>% select(ID, Enhancer, Annotation, Distance.to.TSS)
row.names(homer_unique_enhancer_GW_GW17only_Olig2) <- homer_unique_enhancer_GW_GW17only_Olig2$ID
homer_unique_enhancer_GW_GW17only_Olig2_list <- list(GW17 = homer_unique_enhancer_GW_GW17only_Olig2$ID, DE = c(as.character(DE_GW13_GW17_cortex$ID), as.character(DE_GW13_GW17_GE$ID)))
Venn_homer_unique_enhancer_GW_GW17only_Olig2 <- venn.diagram(homer_unique_enhancer_GW_GW17only_Olig2_list, filename = NULL, fill = c("red", "blue"), main = "GW GW17only Olig2", force.unique = T)
plot.new()
grid.draw(Venn_homer_unique_enhancer_GW_GW17only_Olig2)
homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex <- merge(homer_unique_enhancer_GW_GW17only_Olig2, DE_GW13_GW17_cortex, by = "ID") %>% 
  mutate(name = ensembl[as.character(ID), "name"], description = ensembl[as.character(ID), "description"]) %>% filter(GW13 + GW13 >= 0.2) %>% arrange(-abs(log2FC))
write.table(homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex, file = "homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex.txt", sep = "\t", col.names = T, row.names = F, quote = F)
homer_unique_enhancer_GW_GW17only_Olig2_DE_GE <- merge(homer_unique_enhancer_GW_GW17only_Olig2, DE_GW13_GW17_GE, by = "ID") %>% 
  mutate(name = ensembl[as.character(ID), "name"], description = ensembl[as.character(ID), "description"]) %>% filter(GW13 + GW13 >= 0.2) %>% arrange(-abs(log2FC))
write.table(homer_unique_enhancer_GW_GW17only_Olig2_DE_GE, file = "homer_unique_enhancer_GW_GW17only_Olig2_DE_GE.txt", sep = "\t", col.names = T, row.names = F, quote = F)
(homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_DAVID <- enrich("homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex", erminej = F, fdr = 0.05, height = 2))
(homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_DAVID <- enrich("homer_unique_enhancer_GW_GW17only_Olig2_DE_GE", erminej = F, fdr = 0.05, height = 8))

homer_unique_enhancer_GW_GW17only_Foxo1 <- read.delim("GW17only_Foxo1.annotate", head = T, as.is = T) %>% filter(Gene.Type == "protein-coding", Nearest.Ensembl != "") %>% distinct(Nearest.Ensembl) %>%
  mutate(ID = Nearest.Ensembl, Enhancer = paste0(Chr, ":", Start, "-", End)) %>% select(ID, Enhancer, Annotation, Distance.to.TSS)
row.names(homer_unique_enhancer_GW_GW17only_Foxo1) <- homer_unique_enhancer_GW_GW17only_Foxo1$ID
homer_unique_enhancer_GW_GW17only_Foxo1_list <- list(GW17 = homer_unique_enhancer_GW_GW17only_Foxo1$ID, DE = c(as.character(DE_GW13_GW17_cortex$ID), as.character(DE_GW13_GW17_GE$ID)))
Venn_homer_unique_enhancer_GW_GW17only_Foxo1 <- venn.diagram(homer_unique_enhancer_GW_GW17only_Foxo1_list, filename = NULL, fill = c("red", "blue"), main = "GW GW17only Foxo1", force.unique = T)
plot.new()
grid.draw(Venn_homer_unique_enhancer_GW_GW17only_Foxo1)
homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex <- merge(homer_unique_enhancer_GW_GW17only_Foxo1, DE_GW13_GW17_cortex, by = "ID") %>% 
  mutate(name = ensembl[as.character(ID), "name"], description = ensembl[as.character(ID), "description"]) %>% filter(GW13 + GW13 >= 0.2) %>% arrange(-abs(log2FC))
write.table(homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex, file = "homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex.txt", sep = "\t", col.names = T, row.names = F, quote = F)
homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE <- merge(homer_unique_enhancer_GW_GW17only_Foxo1, DE_GW13_GW17_GE, by = "ID") %>% 
  mutate(name = ensembl[as.character(ID), "name"], description = ensembl[as.character(ID), "description"]) %>% filter(GW13 + GW13 >= 0.2) %>% arrange(-abs(log2FC))
write.table(homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE, file = "homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE.txt", sep = "\t", col.names = T, row.names = F, quote = F)
#### Neurospheres
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/homer/")
cortex_GE_UP_duplicated <- c(cortex01_GE01UP$V1, cortex02_GE02UP$V1, cortex03_GE03UP$V1, cortex04_GE04UP$V1)
cortex_GE_UP_duplicated <- ensembl[unique(cortex_GE_UP_duplicated[duplicated(cortex_GE_UP_duplicated)]), ]
cortex_GE_DN_duplicated <- c(cortex01_GE01DN$V1, cortex02_GE02DN$V1, cortex03_GE03DN$V1, cortex04_GE04DN$V1)
cortex_GE_DN_duplicated <- ensembl[unique(cortex_GE_DN_duplicated[duplicated(cortex_GE_DN_duplicated)]), ]
cortex_GE_DE_duplicated <- rbind(data.frame(cortex_GE_UP_duplicated, DE = "UP"), data.frame(cortex_GE_DN_duplicated, DE = "DN"))
homer_unique_enhancer_Neurospheres_Cortex <- read.delim("./Neurospheres.Cortex/knownResults.txt", head = T, as.is = T) %>% 
  mutate(Motif.Name = gsub("/Homer", "", Motif.Name), Motif.Name = gsub(" ", "_", Motif.Name), motif.file = paste0("./knownResults/known", 1:n(), ".motif"), 
         TF = gsub("\\(.+\\)", "", split_on_last(Motif.Name, "/")[1,]), Assay = split_on_last(Motif.Name, "/")[2,], 
         Percent.of.Target.Sequences.with.Motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif)), 
         Percent.of.Background.Sequences.with.Motif = as.numeric(gsub("%", "", X..of.Background.Sequences.with.Motif)), 
         FC = Percent.of.Target.Sequences.with.Motif/Percent.of.Background.Sequences.with.Motif) %>% 
  filter(q.value..Benjamini. < 0.01) %>% 
  mutate(Motif.Name = factor(Motif.Name, levels = rev(Motif.Name)))
homer_unique_enhancer_Neurospheres_GE <- read.delim("./Neurospheres.GE/knownResults.txt", head = T, as.is = T) %>% 
  mutate(Motif.Name = gsub("/Homer", "", Motif.Name), Motif.Name = gsub(" ", "_", Motif.Name), motif.file = paste0("./knownResults/known", 1:n(), ".motif"), 
         TF = gsub("\\(.+\\)", "", split_on_last(Motif.Name, "/")[1,]), Assay = split_on_last(Motif.Name, "/")[2,], 
         Percent.of.Target.Sequences.with.Motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif)), 
         Percent.of.Background.Sequences.with.Motif = as.numeric(gsub("%", "", X..of.Background.Sequences.with.Motif)), 
         FC = Percent.of.Target.Sequences.with.Motif/Percent.of.Background.Sequences.with.Motif) %>% 
  filter(q.value..Benjamini. < 0.01) %>% 
  mutate(Motif.Name = factor(Motif.Name, levels = rev(Motif.Name)))
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/homer/Neurospheres/")
##### Common TFs
homer_unique_enhancer_Neurospheres_common <- merge(homer_unique_enhancer_Neurospheres_Cortex, homer_unique_enhancer_Neurospheres_GE, by = "Motif.Name") %>% 
  mutate(Consensus = Consensus.x, TF = TF.x, Assay = Assay.x, 
         Cortex_logP = Log.P.value.x, Cortex_q = q.value..Benjamini..x, Cortex_motif_file = motif.file.x, Cortex_Percent_Target = Percent.of.Target.Sequences.with.Motif.x, Cortex_FC = FC.x, 
         GE_logP = Log.P.value.y, GE_q = q.value..Benjamini..y, GE_motif_file = motif.file.y, GE_Percent_Target = Percent.of.Target.Sequences.with.Motif.y, GE_FC = FC.y, delta_Percent_Target = abs(Cortex_Percent_Target - GE_Percent_Target)) %>%
  select(Motif.Name, Consensus, TF, Assay, Cortex_logP, Cortex_q, Cortex_motif_file, Cortex_Percent_Target, Cortex_FC, GE_logP, GE_q, GE_motif_file, GE_Percent_Target, GE_FC, delta_Percent_Target) %>% 
  arrange(- delta_Percent_Target)
write.table(homer_unique_enhancer_Neurospheres_common, file = "homer_unique_enhancer_Neurospheres_common.txt", sep = "\t", col.names = T, row.names = F, quote = F)
homer_unique_enhancer_Neurospheres_common_top <- homer_unique_enhancer_Neurospheres_common %>% filter(Cortex_Percent_Target >= 20 | GE_Percent_Target >= 20) %>% arrange(Cortex_Percent_Target + GE_Percent_Target) %>% mutate(Motif.Name = factor(Motif.Name, levels = as.character(Motif.Name)))
homer_unique_enhancer_Neurospheres_common_top <- data.frame(Motif.Name = rep(homer_unique_enhancer_Neurospheres_common_top$Motif.Name, 2), Neurospheres = rep(c("Cortex", "GE"), each = nrow(homer_unique_enhancer_Neurospheres_common_top)), Percent_Target = c(homer_unique_enhancer_Neurospheres_common_top$Cortex_Percent_Target, homer_unique_enhancer_Neurospheres_common_top$GE_Percent_Target))
(homer_unique_enhancer_Neurospheres_common_figure <- ggplot(data = homer_unique_enhancer_Neurospheres_common_top, aes(Motif.Name, Percent_Target)) +
   geom_bar(aes(fill = Neurospheres), stat = "identity", width = .5, position = position_dodge()) + 
   coord_flip() + 
   scale_fill_manual(values = c("red", "blue"), name = "") + 
   ylab("Percent of Enhancers with Motif") + 
   xlab("") + 
   theme_bw())
ggsave(homer_unique_enhancer_Neurospheres_common_figure, file = "homer_unique_enhancer_Neurospheres_common_figure.pdf")
homer_unique_enhancer_Neurospheres_Cortex_common_Lhx3 <- read.delim("Common_Cortex_Lhx3.annotate", head = T, as.is = T) %>% filter(Gene.Type == "protein-coding", Nearest.Ensembl != "") %>% distinct(Nearest.Ensembl)
row.names(homer_unique_enhancer_Neurospheres_Cortex_common_Lhx3) <- homer_unique_enhancer_Neurospheres_Cortex_common_Lhx3$Nearest.Ensembl
homer_unique_enhancer_Neurospheres_GE_common_Lhx3 <- read.delim("Common_GE_Lhx3.annotate", head = T, as.is = T) %>% filter(Gene.Type == "protein-coding", Nearest.Ensembl != "") %>% distinct(Nearest.Ensembl)
row.names(homer_unique_enhancer_Neurospheres_GE_common_Lhx3) <- homer_unique_enhancer_Neurospheres_GE_common_Lhx3$Nearest.Ensembl
homer_unique_enhancer_Neurospheres_common_Lhx3_list <- list(Cortex = homer_unique_enhancer_Neurospheres_Cortex_common_Lhx3$Nearest.Ensembl, GE = homer_unique_enhancer_Neurospheres_GE_common_Lhx3$Nearest.Ensembl, DE = cortex_GE_DE_duplicated$id)
Venn_homer_unique_enhancer_Neurospheres_common_Lhx3 <- venn.diagram(homer_unique_enhancer_Neurospheres_common_Lhx3_list, filename = NULL, fill = c("red", "blue", "green"), main = "Neurospheres common Lhx3", force.unique = T)
plot.new()
grid.draw(Venn_homer_unique_enhancer_Neurospheres_common_Lhx3)
homer_unique_enhancer_Neurospheres_common_Lhx3 <- data.frame(id = unique(c(homer_unique_enhancer_Neurospheres_Cortex_common_Lhx3$Nearest.Ensembl, homer_unique_enhancer_Neurospheres_GE_common_Lhx3$Nearest.Ensembl))) %>%
  mutate(Cortex = homer_unique_enhancer_Neurospheres_Cortex_common_Lhx3[as.character(id), 1], GE = homer_unique_enhancer_Neurospheres_GE_common_Lhx3[as.character(id), 1])
homer_unique_enhancer_Neurospheres_common_Lhx3_DE <- merge(homer_unique_enhancer_Neurospheres_common_Lhx3, cortex_GE_DE_duplicated, by = "id") 
write.table(homer_unique_enhancer_Neurospheres_common_Lhx3_DE, file = "homer_unique_enhancer_Neurospheres_common_Lhx3_DE.txt", sep = "\t", col.names = T, row.names = F, quote = F)
##### Cortex TFs
homer_unique_enhancer_Neurospheres_CortexOnly <- homer_unique_enhancer_Neurospheres_Cortex[!(homer_unique_enhancer_Neurospheres_Cortex$Motif.Name %in% homer_unique_enhancer_Neurospheres_GE$Motif.Name), ] %>%
  select(Motif.Name, TF, Assay, motif.file, Consensus, Log.P.value, q.value..Benjamini., Percent.of.Target.Sequences.with.Motif, Percent.of.Background.Sequences.with.Motif, FC) %>% 
  arrange(- Percent.of.Target.Sequences.with.Motif) %>% mutate(Motif.Name = factor(Motif.Name, levels = rev(as.character(Motif.Name))))
write.table(homer_unique_enhancer_Neurospheres_CortexOnly, file = "homer_unique_enhancer_Neurospheres_CortexOnly.txt", sep = "\t", col.names = T, row.names = F, quote = F)
(homer_unique_enhancer_Neurospheres_CortexOnly_figure <- ggplot(data = homer_unique_enhancer_Neurospheres_CortexOnly, aes(Motif.Name, Percent.of.Target.Sequences.with.Motif)) +
   geom_bar(fill = "blue", stat = "identity", width = .5) + 
   coord_flip() + 
   ylab("Percent of Enhancers with Motif") + 
   xlab("") + 
   theme_bw())
ggsave(homer_unique_enhancer_Neurospheres_CortexOnly_figure, file = "homer_unique_enhancer_Neurospheres_CortexOnly_figure.pdf")
##### GE TFs
homer_unique_enhancer_Neurospheres_GEonly <- homer_unique_enhancer_Neurospheres_GE[!(homer_unique_enhancer_Neurospheres_GE$Motif.Name %in% homer_unique_enhancer_Neurospheres_Cortex$Motif.Name), ] %>%
  select(Motif.Name, TF, Assay, motif.file, Consensus, Log.P.value, q.value..Benjamini., Percent.of.Target.Sequences.with.Motif, Percent.of.Background.Sequences.with.Motif, FC) %>% 
  arrange(- Percent.of.Target.Sequences.with.Motif) %>% mutate(Motif.Name = factor(Motif.Name, levels = rev(as.character(Motif.Name))))
write.table(homer_unique_enhancer_Neurospheres_GEonly, file = "homer_unique_enhancer_Neurospheres_GEonly.txt", sep = "\t", col.names = T, row.names = F, quote = F)
(homer_unique_enhancer_Neurospheres_GEonly_figure <- ggplot(data = filter(homer_unique_enhancer_Neurospheres_GEonly, Percent.of.Target.Sequences.with.Motif >= 20), aes(Motif.Name, Percent.of.Target.Sequences.with.Motif)) +
   geom_bar(fill = "blue", stat = "identity", width = .5) + 
   coord_flip() + 
   ylab("Percent of Enhancers with Motif") + 
   xlab("") + 
   theme_bw())
ggsave(homer_unique_enhancer_Neurospheres_GEonly_figure, file = "homer_unique_enhancer_Neurospheres_GEonly_figure.pdf", height = 4)
homer_unique_enhancer_Neurospheres_GEonly_Olig2 <- read.delim("GEonly_Olig2.annotate", head = T, as.is = T) %>% filter(Gene.Type == "protein-coding", Nearest.Ensembl != "") %>% distinct(Nearest.Ensembl) %>%
  mutate(id = Nearest.Ensembl, Enhancer = paste0(Chr, ":", Start, "-", End)) %>% select(id, Enhancer, Annotation, Distance.to.TSS)
row.names(homer_unique_enhancer_Neurospheres_GEonly_Olig2) <- homer_unique_enhancer_Neurospheres_GEonly_Olig2$id
homer_unique_enhancer_Neurospheres_GEonly_Olig2_list <- list(GE = homer_unique_enhancer_Neurospheres_GEonly_Olig2$id, DE = cortex_GE_DE_duplicated$id)
Venn_homer_unique_enhancer_Neurospheres_GEonly_Olig2 <- venn.diagram(homer_unique_enhancer_Neurospheres_GEonly_Olig2_list, filename = NULL, fill = c("red", "blue"), main = "Neurospheres GE_only Olig2", force.unique = T)
plot.new()
grid.draw(Venn_homer_unique_enhancer_Neurospheres_GEonly_Olig2)
homer_unique_enhancer_Neurospheres_GEonly_Olig2_DE <- merge(homer_unique_enhancer_Neurospheres_GEonly_Olig2, cortex_GE_DE_duplicated, by = "id") 
write.table(homer_unique_enhancer_Neurospheres_GEonly_Olig2_DE, file = "homer_unique_enhancer_Neurospheres_GEonly_Olig2_DE.txt", sep = "\t", col.names = T, row.names = F, quote = F)
(homer_unique_enhancer_Neurospheres_GEonly_Olig2_DE_DAVID <- enrich("homer_unique_enhancer_Neurospheres_GEonly_Olig2_DE", erminej = F, fdr = 0.05, height = 8))
#### cytoscape
setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/homer/GW/")
homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_node <- homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex %>% 
  filter(abs(Distance.to.TSS) < 10000) %>% select(name, log2FC, cortex01_cortex04) %>% mutate(log2FC = abs(log2FC))
homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_network <- homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_node %>% mutate(Source = "OLIG2", Interaction = "activate")
homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_network[homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_network$cortex01_cortex04 == "DN", "Interaction"] <- "repress"
homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_node <- rbind(data.frame(name = "OLIG2", log2FC = abs(DE_GW13_GW17_cortex[DE_GW13_GW17_cortex$ID == "ENSG00000205927", "log2FC"]), cortex01_cortex04 = DE_GW13_GW17_cortex[DE_GW13_GW17_cortex$ID == "ENSG00000205927", "cortex01_cortex04"]), homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_node)
write.table(homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_node, file = "homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_node.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_network, file = "homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_network.txt", sep = "\t", quote = F, col.names = T, row.names = F)
homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_node <- homer_unique_enhancer_GW_GW17only_Olig2_DE_GE %>% 
  filter(abs(Distance.to.TSS) < 10000) %>% select(name, log2FC, GE01_GE04) %>% mutate(log2FC = abs(log2FC))
homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_network <- homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_node %>% mutate(Source = "OLIG2", Interaction = "activate")
homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_network[homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_network$GE01_GE04 == "DN", "Interaction"] <- "repress"
homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_node <- rbind(data.frame(name = "OLIG2", log2FC = abs(DE_GW13_GW17_GE[DE_GW13_GW17_GE$ID == "ENSG00000205927", "log2FC"]), GE01_GE04 = DE_GW13_GW17_GE[DE_GW13_GW17_GE$ID == "ENSG00000205927", "GE01_GE04"]), homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_node)
write.table(homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_node, file = "homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_node.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_network, file = "homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_network.txt", sep = "\t", quote = F, col.names = T, row.names = F)
homer_unique_enhancer_GW_GW17only_Olig2_DE_node <- rbind(homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_node %>% mutate(DE = GE01_GE04) %>% select(-GE01_GE04), homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_node %>% mutate(DE = cortex01_cortex04) %>% select(-cortex01_cortex04)) %>% 
  mutate(name = as.character(name)) %>% distinct(name)
homer_unique_enhancer_GW_GW17only_Olig2_DE_network <- homer_unique_enhancer_GW_GW17only_Olig2_DE_node[2:nrow(homer_unique_enhancer_GW_GW17only_Olig2_DE_node),] %>% mutate(Source = "OLIG2", Interaction = "activate")
write.table(homer_unique_enhancer_GW_GW17only_Olig2_DE_node, file = "homer_unique_enhancer_GW_GW17only_Olig2_DE_node.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(homer_unique_enhancer_GW_GW17only_Olig2_DE_network, file = "homer_unique_enhancer_GW_GW17only_Olig2_DE_network.txt", sep = "\t", quote = F, col.names = T, row.names = F)

homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_node <- homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex %>% 
  filter(abs(Distance.to.TSS) < 10000) %>% select(name, log2FC, cortex01_cortex04) %>% mutate(log2FC = abs(log2FC))
homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_network <- homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_node %>% mutate(Source = "Foxo1", Interaction = "activate")
homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_network[homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_network$cortex01_cortex04 == "DN", "Interaction"] <- "repress"
homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_node <- rbind(data.frame(name = "Foxo1", log2FC = abs(DE_GW13_GW17_cortex[DE_GW13_GW17_cortex$ID == "ENSG00000205927", "log2FC"]), cortex01_cortex04 = DE_GW13_GW17_cortex[DE_GW13_GW17_cortex$ID == "ENSG00000205927", "cortex01_cortex04"]), homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_node)
write.table(homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_node, file = "homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_node.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_network, file = "homer_unique_enhancer_GW_GW17only_Foxo1_DE_cortex_network.txt", sep = "\t", quote = F, col.names = T, row.names = F)
homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_node <- homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE %>% 
  filter(abs(Distance.to.TSS) < 10000) %>% select(name, log2FC, GE01_GE04) %>% mutate(log2FC = abs(log2FC))
homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_network <- homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_node %>% mutate(Source = "Foxo1", Interaction = "activate")
homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_network[homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_network$GE01_GE04 == "DN", "Interaction"] <- "repress"
homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_node <- rbind(data.frame(name = "Foxo1", log2FC = abs(DE_GW13_GW17_GE[DE_GW13_GW17_GE$ID == "ENSG00000205927", "log2FC"]), GE01_GE04 = DE_GW13_GW17_GE[DE_GW13_GW17_GE$ID == "ENSG00000205927", "GE01_GE04"]), homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_node)
write.table(homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_node, file = "homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_node.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_network, file = "homer_unique_enhancer_GW_GW17only_Foxo1_DE_GE_network.txt", sep = "\t", quote = F, col.names = T, row.names = F)

save(FindER_summary, FindER_summary_figure, HisMod_RPKM, HisMod_RPKM_figure, 
     H3K4me3_TSS1500, H3K27me3_TSS1500, His_DM_promoter, His_DM_promoter_figure, 
     Brain01_Brain02DE_epi, Cortex01_Cortex02DE_epi, GE01_GE02DE_epi, Cortex01_GE01DE_epi, Cortex02_GE02DE_epi, GE02_GE04DE_epi, 
     venn_Brain01_Brain02DE_epi, venn_Cortex01_Cortex02DE_epi, venn_GE01_GE02DE_epi, venn_Cortex01_GE01DE_epi, venn_Cortex02_GE02DE_epi, venn_GE02_GE04DE_epi, 
     GWAS_core_enhancer, GWAS_core_enhancer_TFBS, GWAS_core_enhancer_TFBS_brain, GWAS_core_enhancer_trait, GWAS_core_enhancer_trait_sig, GWAS_core_enhancer_trait_sig_brain, 
     homer_core_enhancer, homer_core_enhancer_top, homer_core_enhancer_figure, homer_unique_enhancer_GW_common_figure, 
     UMR_enhancer, UMR_enhancer_enrich, UMR_enhancer_summary_figure, 
     GREAT_Cortex_UMR_enhancer_HuFNSC02, GREAT_GE_UMR_enhancer_HuFNSC02, GREAT_Cortex_UMR_enhancer_HuFNSC04, 
     GREAT_GW17_UMR_enhancer_Cortex, GREAT_GW17_UMR_enhancer_GE, GREAT_GW13_UMR_enhancer_GE, 
     unique_enhancers_summary, GWAS_unique_enhancer_all, GWAS_unique_enhancer_all_TFBS, GWAS_unique_enhancer_brain_TFBS, GWAS_unique_enhancer_trait, GWAS_unique_enhancer_trait_sig, GWAS_unique_enhancer_trait_sig_brain, 
     venn_unique_enhancer_HuFNSC01, venn_unique_enhancer_HuFNSC02, venn_unique_enhancer_Cortex, venn_unique_enhancer_GE, venn_unique_enhancer_GW13, venn_unique_enhancer_GW17, 
     GREAT_unique_enhancer_HuFNSC01, GREAT_unique_enhancer_HuFNSC02, GREAT_unique_enhancer_Cortex, GREAT_unique_enhancer_GE, GREAT_unique_enhancer_GW13, GREAT_unique_enhancer_GW17, 
     homer_unique_enhancer_GW_GW13, homer_unique_enhancer_GW_GW17, homer_unique_enhancer_GW_common, Venn_homer_unique_enhancer_GW_common_Sox3, 
     DE_GW13_GW17_cortex, DE_GW13_GW17_GE, homer_unique_enhancer_GW_common_Sox3_DE_cortex, homer_unique_enhancer_GW_common_Sox3_DE_GE, 
     homer_unique_enhancer_GW_GW13only, homer_unique_enhancer_GW_GW13only_figure, homer_unique_enhancer_GW_GW17only, homer_unique_enhancer_GW_GW17only_Olig2, Venn_homer_unique_enhancer_GW_GW17only_Olig2, 
     homer_unique_enhancer_GW_GW17only_figure, homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex, homer_unique_enhancer_GW_GW17only_Olig2_DE_GE, homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_DAVID, homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_DAVID, 
     homer_unique_enhancer_Neurospheres_Cortex, homer_unique_enhancer_Neurospheres_GE, homer_unique_enhancer_Neurospheres_common, homer_unique_enhancer_Neurospheres_common_figure, 
     cortex_GE_DE_duplicated, Venn_homer_unique_enhancer_Neurospheres_common_Lhx3, homer_unique_enhancer_Neurospheres_common_Lhx3_DE, 
     homer_unique_enhancer_Neurospheres_CortexOnly, homer_unique_enhancer_Neurospheres_CortexOnly_figure, homer_unique_enhancer_Neurospheres_GEonly, homer_unique_enhancer_Neurospheres_GEonly_figure, 
     homer_unique_enhancer_Neurospheres_GEonly_Olig2, Venn_homer_unique_enhancer_Neurospheres_GEonly_Olig2, homer_unique_enhancer_Neurospheres_GEonly_Olig2_DE, homer_unique_enhancer_Neurospheres_GEonly_Olig2_DE_DAVID, 
     homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_network, homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_node, homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_network, homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_node, 
     file = "/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FetalBrain_FindER.Rdata")
