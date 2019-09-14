# Glioma - WGS analysis

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggbio)
library(GenomicRanges)
library(RCircos)
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/glioma/")

## ------- IDH mut rate --------------
mut_rate <- read.delim("mut_rate.txt", as.is = T) %>% mutate(rate = rate * 100, category = gsub("\\..*", "", sample), sample = gsub(".*\\.", "", sample))
(mut_rate_figure <- ggplot(mut_rate, aes(sample, rate, fill = category)) + 
		geom_bar(stat = "identity") + 
		geom_hline(yintercept = 30) + 
		scale_fill_manual(values = c(hcl(h = seq(15, 375, length = 5 + 1)[1], l = 65, c = 100), hcl(h = seq(15, 375, length = 5 + 1)[2], l = 65, c = 100),hcl(h = seq(15, 375, length = 5 + 1)[4], l = 65, c = 100))) + 
		xlab("") + 
		ylab("Ratio of mutant IDH read count") + 
		coord_flip() + 
		theme_bw())
ggsave(mut_rate_figure, file = "mut_rate.pdf", height = 4, width = 5)

## -------- sample genetic mutations ------
mutation <- read.delim("mutation.txt", as.is = T) %>% melt(id = "Gene") %>% mutate(Gene = factor(Gene, levels = c("TP53", "CIC", "1p19q LOH", "IDH2", "IDH1")))
(mutation_figure <- ggplot(mutation, aes(variable, Gene, fill = value)) +
		geom_tile() +
		scale_fill_manual(values = c("IDH1 R132H" = "red", "IDH2 R172K" = "dark red", "IDH2 R172M" = "dark red", "CIC R215W" = "purple", "CIC R202W" = "purple", "TP53 R273L" = "blue", "TP53 R273C" = "blue", "1p19q LOH" = "dark green", "N" = "white"), guide = "none") +
		xlab("") +
		ylab("") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, size = 15, hjust = 1), axis.text.y = element_text(size = 15)))
ggsave(mutation_figure, file = "mutation.pdf", width = 8, height = 2.5)

##################################################################

setwd("/projects/epigenomics2/users/lli/glioma/WGS/")
load("/projects/epigenomics2/users/lli/glioma/WGS/WGS.Rdata")
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


