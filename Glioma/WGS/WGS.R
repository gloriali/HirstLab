# Glioma - WGS analysis

library(ggplot2)
library(dplyr)
library(ggbio)
library(GenomicRanges)
library(RCircos)
load("/projects/epigenomics2/users/lli/glioma/WGS/WGS.Rdata")
setwd("/projects/epigenomics2/users/lli/glioma/WGS/")
libs <- c("CEMT_19", "CEMT_21", "CEMT_22", "CEMT_23", "CEMT_47")

## Common genetic variations in glioma
glioma_mutations <- read.delim("/home/lli/HirstLab/Glioma/WGS/CEMT_gliomaMutations.tsv", as.is = T)
IDH <- read.delim("/home/lli/HirstLab/Glioma/WGS/IDH.tsv", as.is = T)

## CNVs
colname <- c("chr", "start", "end", "ID", "CN", "type")
for(lib in libs){
	print(lib)
	assign(paste0("CNV_", lib), read.delim(paste0("./CNV/", lib, ".CNV.bed"), head = F, as.is = T, col.names = colname) %>% 
				 	mutate(chr = paste0("chr", chr), xstart = start, xend = end, CN = CN * 5 + 11))
	CNV_gr <- keepSeqlevels(as(get(paste0("CNV_", lib)), "GRanges"), paste0("chr", c(1:22, "X")))
	data(hg19IdeogramCyto, package = "biovizBase")
	hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X")))
	assign(paste0("CNV_", lib, "_figure"), 
				 ggplot(hg19) + 
				 	layout_karyogram(cytoband = TRUE) + 
				 	layout_karyogram(CNV_gr, geom = "errorbarh", aes(x = xstart, xmin = xstart, xmax = xend, y = CN, color = type, height = 0)) + 
				 	xlab("") + 
				 	ylab("Copy Number") +
				 	ggtitle(paste("CNV -", lib)) + 
				 	coord_cartesian(ylim = c(0, 35)) + 
				 	scale_y_continuous(breaks = c(11, 21, 31), labels = c(0, 2, 10)) + 
				 	scale_color_manual(values = c("red", "blue"), name = "") + 
				 	guides(fill = "none") + 
				 	theme_bw() + 
				 	theme(axis.text.y = element_text(size = 5), strip.text = element_text(size = 8)))
	ggsave(get(paste0("CNV_", lib, "_figure"))@ggplot, file = paste0("./CNV/CNV_", lib, "_figure.pdf"), height = 8, width = 6)
}


save(list = c("glioma_mutations", "IDH", ls(pattern = "CNV_.*_figure")), 
		 file = "/projects/epigenomics2/users/lli/glioma/WGS/WGS.Rdata")


