# Fetal Brain - Histone modifications FindER 
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)

setwd("/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER")
## ============ Sanity check ===============
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


save(FindER_summary, FindER_summary_figure, 
     file = "FetalBrain_FindER.Rdata")
