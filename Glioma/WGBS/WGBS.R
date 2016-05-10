# Glioma - WGBS analysis

library(ggplot2)
library(plyr)
library(VennDiagram)
library(gridExtra)
library(gplots)
library(dendextend)
library(reshape)
library(wq)
library(dplyr)
library(RCircos)
source('~/HirstLab/Pipeline/R/DMR.figures.R')
source('~/HirstLab/Pipeline/R/enrich.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')

## ------- DMR glioma vs NPC -------
setwd("/projects/epigenomics2/users/lli/glioma/WGBS/DMR")
libs <- c("CEMT_19", "CEMT_21", "CEMT_22", "CEMT_23")
### -------- summary ---------
DMR_summary <- read.delim("DMR.summary.stats", head = F, as.is = T, col.names = c("Sample", "size", "cut", "N_hyper", "N_hypo", "len_hyper", "len_hypo"))
DMR_summary_tall <- data.frame(Sample = rep(gsub("_NPCs", "", DMR_summary$Sample), 2), DM = rep(c("hyper", "hypo"), each = nrow(DMR_summary)), length = c(DMR_summary$len_hyper, -DMR_summary$len_hypo)/10^6)
DMR_summary_figure <- ggplot(DMR_summary_tall, aes(x = Sample, y = length, fill = DM)) + 
	geom_bar(stat="identity", width=.5, position = "identity") + 
	xlab("") + 
	ylab("Total DMR length (Mb)") + 
	scale_fill_manual(name = "", values = c("red", "blue")) + 
	theme_bw()
ggsave(DMR_summary_figure, file = "DMR_summary_figure.pdf", height = 5, width = 6)

### ------- visualization -----
colname <- c("chr", "start", "end", "ID", "DM", "len")
for(lib in libs){
	print(lib)
	assign(paste0("DMR_", lib), read.delim(paste0("DMR.", lib, "_NPC"), head = F, as.is = T, col.names = colname) %>% mutate(chr = paste0("chr", chr)))
	assign(paste0("DMR_", lib, "_figure"), DMR_figures(get(paste0("DMR_", lib)), lib, "NPCs", figures = c("length", "frequency", "circos"), colname = c("chr", "start", "end", "ID", "DM", "length"), hist_width = 10))
}

### -------- GREAT --------
(GREAT_DMR_CEMT_19_hyper <- enrich_GREAT("CEMT_19_hyper", "CEMT_19_hyper", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))
(GREAT_DMR_CEMT_19_hypo <- enrich_GREAT("CEMT_19_hypo", "CEMT_19_hypo", categories = c("GOBP", "MSigMotif", "MSigPerturbation"), height = 12, width = 9))
(GREAT_DMR_CEMT_21_hyper <- enrich_GREAT("CEMT_21_hyper", "CEMT_21_hyper", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))
(GREAT_DMR_CEMT_21_hypo <- enrich_GREAT("CEMT_21_hypo", "CEMT_21_hypo", categories = c("GOBP", "MSigPerturbation"), height = 3, width = 6))
(GREAT_DMR_CEMT_22_hyper <- enrich_GREAT("CEMT_22_hyper", "CEMT_22_hyper", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))
(GREAT_DMR_CEMT_22_hypo <- enrich_GREAT("CEMT_22_hypo", "CEMT_22_hypo", categories = c("MSigMotif", "MSigPerturbation"), height = 12, width = 9))
(GREAT_DMR_CEMT_23_hyper <- enrich_GREAT("CEMT_23_hyper", "CEMT_23_hyper", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 9, width = 7))
(GREAT_DMR_CEMT_23_hypo <- enrich_GREAT("CEMT_23_hypo", "CEMT_23_hypo", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))



save(list = c("DMR_summary_figure", ls(pattern = "DMR_CEMT_.*_figure"), ls(pattern = "GREAT_DMR_*")), 
		 file = "/projects/epigenomics2/users/lli/glioma/WGBS/WGBS.Rdata")
