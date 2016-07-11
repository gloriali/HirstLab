# Glioma RNA-seq expression analysis
library(ggplot2)
library(dplyr)
library(reshape2)
source('~/HirstLab/Pipeline/R/enrich.R')
load("/projects/epigenomics2/users/lli/glioma/RNAseq/RNAseq.Rdata")
setwd("/projects/epigenomics2/users/lli/glioma/RNAseq/")
libs <- c("CEMT_19", "CEMT_21", "CEMT_22", "CEMT_23", "CEMT_47")

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

save(list = c("IDH_RPKM_figure", "DE_pc_summary_figure", "libs", ls(pattern = "DAVID")), 
		 file = "/projects/epigenomics2/users/lli/glioma/RNAseq/RNAseq.Rdata")
