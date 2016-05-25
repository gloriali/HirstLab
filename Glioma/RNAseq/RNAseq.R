# Glioma RNA-seq expression analysis
library(ggplot2)
library(dplyr)
library(reshape2)

## ----- IDH1/2 expression levels ------
setwd("/projects/epigenomics2/users/lli/glioma/RNAseq/RPKM")
IDH_glioma <- read.delim("IDH.RPKM", head = F, col.names = c("ID", "Gene", "Sample", "RPKM"))
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
ggsave(IDH_RPKM_figure, file = "IDH_RPKM_figure.pdf", height = 5, width = 6)
