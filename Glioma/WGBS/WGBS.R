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
load("/projects/epigenomics2/users/lli/glioma/WGBS/WGBS.Rdata")
setwd("/projects/epigenomics2/users/lli/glioma/WGBS/DMR")
libs <- c("CEMT_19", "CEMT_21", "CEMT_22", "CEMT_23", "CEMT_47")

## ------- 5mC distribution------
quantile_5mC <- read.delim("../qc.5mC.quantile", as.is = T) %>% mutate(type = gsub(".*\\d+_", "", sample), sample = gsub("_[gC].*", "", sample), category = gsub("[_\\.].*", "", sample))
quantile_5mC_figure <- ggplot(quantile_5mC, aes(x = sample, fill = category)) + 
	geom_boxplot(aes(lower = lower, middle = median, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.5, color = "grey") + 
	facet_grid(type ~ .) + 
	xlab("") + 
	ylab("Fractional methylation") + 
	scale_fill_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw()
ggsave(quantile_5mC_figure, file = "../quantile_5mC_figure.pdf", height = 5, width = 6)

## ------- DMR glioma vs NPC -------
### -------- summary ---------
DMR_summary <- read.delim("./intermediate/DMR.summary.stats", head = T, as.is = T)
DMR_summary_tall <- data.frame(glioma = rep(gsub("_NPC.*", "", DMR_summary$sample), 2), NPC = rep(gsub(".*_NPC\\.", "", DMR_summary$sample), 2), DM = rep(c("hyper", "hypo"), each = nrow(DMR_summary)), length = c(DMR_summary$hyper, -DMR_summary$hypo)/10^6)
DMR_summary_figure <- ggplot(DMR_summary_tall, aes(x = glioma, y = length, color = DM, shape = NPC)) + 
	geom_point(position = position_jitter(width = 0.1), size = 3) + 
	geom_hline(yintercept = 0) + 
	xlab("") + 
	ylab("Total DMR length (Mb)") + 
	scale_color_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw()
ggsave(DMR_summary_figure, file = "DMR_summary_figure.pdf", height = 5, width = 6)

### ------- visualization -----
colname <- c("chr", "start", "end", "ID", "DM", "length")
for(lib in libs){
	print(lib)
	assign(paste0("DMR_", lib), read.delim(paste0("DMR.", lib, "_NPC"), head = F, as.is = T, col.names = colname) %>% mutate(chr = paste0("chr", chr)))
	assign(paste0("DMR_", lib, "_figure"), DMR_figures(get(paste0("DMR_", lib)), lib, "NPCs", figures = c("length", "frequency", "circos"), colname = colname, hist_width = 10))
}

### -------- enrichment in genomic regions ------
genomic_breakdown <- read.delim("./CpG/genomic.breakdown.summary", as.is = T) %>% 
	mutate(sample = gsub("_NPC.*", "", Name), DM = gsub(".*NPC\\.", "", Name), NCpG = NULL, Name = NULL)
genomic_breakdown_tall <- melt(genomic_breakdown, id = c("sample", "DM")) %>% 
	mutate(value = ifelse(DM == "hyper", value, -value))
genomic_breakdown_figure <- ggplot(genomic_breakdown_tall, aes(variable, value, fill = DM)) + 
	geom_bar(position = "identity", stat = "identity", width = 0.5) + 
	facet_wrap(~sample) + 
	xlab("") + 
	ylab("Fold enrichment") + 
	scale_fill_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw()
ggsave(genomic_breakdown_figure, file = "genomic_breakdown_figure.pdf", height = 5, width = 6)

### -------- distance to closest CGI --------
colname <- c("chr", "start", "end", "ID", "CGI_chr", "CGI_start", "CGI_end", "CGI_ID", "distance", "norm_dis")
for(lib in libs){
	print(lib)
	assign(paste0("DMR_", lib, "_CGI_dis"), rbind((read.delim(paste0("./CGI_dis/DMR.", lib, "_NPC.hyper.CGI.dis"), head = F, as.is = T, col.names = colname) %>% mutate(DM = "hyper")), 
																								(read.delim(paste0("./CGI_dis/DMR.", lib, "_NPC.hypo.CGI.dis"), head = F, as.is = T, col.names = colname) %>% mutate(DM = "hypo"))))
	assign(paste0("DMR_", lib, "_CGI_dis_figure"), ggplot(get(paste0("DMR_", lib, "_CGI_dis")) %>% filter(abs(norm_dis) <= 3), aes(norm_dis, color = DM)) + 
				 	geom_density() + 
				 	geom_vline(xintercept = c(-1, 1)) + 
				 	scale_color_manual(name = "", values = c("red", "blue")) + 
				 	xlab("Normalized distance to CGI") + 
				 	ylab("density") + 
				 	ggtitle(lib) + 
				 	theme_bw())
	ggsave(get(paste0("DMR_", lib, "_CGI_dis_figure")), file = paste0("DMR_", lib, "_CGI_dis_figure.pdf"), height = 5, width = 6)
}

### -------- GREAT --------
(GREAT_DMR_CEMT_19_hyper <- enrich_GREAT("CEMT_19_hyper", "CEMT_19_hyper", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))
(GREAT_DMR_CEMT_19_hypo <- enrich_GREAT("CEMT_19_hypo", "CEMT_19_hypo", categories = c("GOBP"), height = 3, width = 7))
(GREAT_DMR_CEMT_21_hyper <- enrich_GREAT("CEMT_21_hyper", "CEMT_21_hyper", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))
(GREAT_DMR_CEMT_21_hypo <- enrich_GREAT("CEMT_21_hypo", "CEMT_21_hypo", categories = c("MSigPerturbation"), height = 4, width = 7))
(GREAT_DMR_CEMT_22_hyper <- enrich_GREAT("CEMT_22_hyper", "CEMT_22_hyper", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))
# GREAT_DMR_CEMT_22_hypo : no enrichment
(GREAT_DMR_CEMT_23_hyper <- enrich_GREAT("CEMT_23_hyper", "CEMT_23_hyper", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 9, width = 7))
(GREAT_DMR_CEMT_23_hypo <- enrich_GREAT("CEMT_23_hypo", "CEMT_23_hypo", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))
(GREAT_DMR_CEMT_47_hyper <- enrich_GREAT("CEMT_47_hyper", "CEMT_47_hyper", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 9, width = 7))
(GREAT_DMR_CEMT_47_hypo <- enrich_GREAT("CEMT_47_hypo", "CEMT_47_hypo", categories = c("MSigPerturbation"), height = 2, width = 7))

### -------- Venn Diagram --------
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
for(i1 in 1:(length(libs)-1)){
	for(i2 in (i1+1):length(libs)){
		print(c(i1, i2, libs[i1], libs[i2]))
		assign(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hyper"), 
					 draw.pairwise.venn(as.numeric(system(paste0("less DMR.", libs[i1], "_NPC.hyper.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T)), 
					 									 as.numeric(system(paste0("less DMR.", libs[i2], "_NPC.hyper.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T)), 
					 									 as.numeric(system(paste0(BEDTOOLS, "/intersectBed -a  DMR.", libs[i1], "_NPC.hyper.bed -b DMR.", libs[i2], "_NPC.hyper.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T)), 
					 									 category = c(libs[i1], libs[i2]), cat.pos = 0, ext.pos = 180))
		assign(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hypo"), 
					 draw.pairwise.venn(as.numeric(system(paste0("less DMR.", libs[i1], "_NPC.hypo.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T)), 
					 									 as.numeric(system(paste0("less DMR.", libs[i2], "_NPC.hypo.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T)), 
					 									 as.numeric(system(paste0(BEDTOOLS, "/intersectBed -a  DMR.", libs[i1], "_NPC.hypo.bed -b DMR.", libs[i2], "_NPC.hypo.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T)), 
					 									 category = c(libs[i1], libs[i2]), cat.pos = 0, ext.pos = 180))
	}
}
pdf("Venn_DMR.pdf")
for(i1 in 1:(length(libs)-1)){
	for(i2 in (i1+1):length(libs)){
		grid.arrange(gTree(children = get(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hyper"))), gTree(children = get(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hypo"))), nrow = 2)
		grid.text("hyper", x = unit(0.1, "npc"), y = unit(0.98, "npc"))
		grid.text("hypo", x = unit(0.1, "npc"), y = unit(0.48, "npc"))
	}
}
dev.off()


save(list = c("quantile_5mC_figure", "DMR_summary_figure", "genomic_breakdown_figure", 
							ls(pattern = "DMR_CEMT_\\d+_figure"), ls(pattern = "DMR_CEMT_\\d+_CGI_dis_figure"), 
							ls(pattern = "GREAT_DMR_*"), ls(pattern = "Venn_DMR_*")),
		 file = "/projects/epigenomics2/users/lli/glioma/WGBS/WGBS.Rdata")

