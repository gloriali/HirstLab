# Glioma - WGBS analysis

library(ggplot2)
library(plyr)
library(VennDiagram)
library(gridExtra)
library(gplots)
library(dendextend)
library(reshape2)
library(wq)
library(dplyr)
library(RCircos)
library(stringr)
source('~/HirstLab/Pipeline/R/DMR.figures.R')
source('~/HirstLab/Pipeline/R/enrich.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
load("/projects/epigenomics2/users/lli/glioma/WGBS/WGBS.Rdata")
setwd("/projects/epigenomics2/users/lli/glioma/WGBS/DMR")
libs <- c("CEMT_19", "CEMT_21", "CEMT_22", "CEMT_23", "CEMT_47")

## ------- 5mC distribution------
quantile_5mC <- read.delim("../qc.5mC.quantile", as.is = T) %>% mutate(type = gsub(".*\\d+_", "", sample), sample = gsub("_[gC].*", "", sample), category = gsub("[_\\.].*", "", sample))
(quantile_5mC_figure <- ggplot(quantile_5mC, aes(x = sample, fill = category)) + 
	geom_boxplot(aes(lower = lower, middle = median, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.5, color = "grey") + 
	facet_grid(type ~ .) + 
	xlab("") + 
	ylab("Fractional methylation") + 
	scale_fill_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw())
ggsave(quantile_5mC_figure, file = "../quantile_5mC_figure.pdf", height = 5, width = 6)

## ------- changes at CGI edges
CGI_edge <- read.delim("../CGI_edge/CGI.edge.all", head = F, col.names = c("chr", "start", "end", "ID", "mC", "CGI", "edge", "dis", "sample")) %>% 
	mutate(edge = revalue(edge, c("L" = "5-prime", "R" = "3-prime")))
CGI_edge_CpG_ID <- CGI_edge %>% group_by(ID) %>% summarise(n = n()) %>% filter(n == 9)
CGI_edge <- CGI_edge %>% filter(ID %in% CGI_edge_CpG_ID$ID)
(CGI_edge_figure <- ggplot(CGI_edge, aes(-dis, mC, color = sample)) + 
	geom_smooth() + 
	facet_wrap(~ edge, scales = "free_x") + 
	xlab("Distance to CGI edge (bp)") + 
	ylab("Fractional methylation") + 
	theme_bw())
ggsave(CGI_edge_figure, file = "../CGI_edge/CGI_edge_figure.pdf", width = 8, height = 5)
CGI_edge_delta <- read.delim("../CGI_edge/CGI.edge.delta.all", head = F, col.names = c("ID", "CGI", "edge", "dis", "delta", "sample", "NPC")) %>% 
	mutate(NPC = gsub("CEMT_[0-9]+-", "", NPC), edge = revalue(edge, c("L" = "5-prime", "R" = "3-prime")))
(CGI_edge_delta_figure <- ggplot(CGI_edge_delta, aes(-dis, delta, color = NPC)) + 
		geom_smooth() + 
		facet_grid(sample ~ edge, scales = "free_x") + 
		xlab("Distance to CGI edge (bp)") + 
		ylab("Difference in fractional methylation\nglioma - NPC") + 
		theme_bw())
ggsave(CGI_edge_delta_figure, file = "../CGI_edge/CGI_edge_delta_figure.pdf", width = 8, height = 8)

## ------- 5mC at CTCF loss regions -------
CTCF_loss_5mC <- read.delim("../CTCF/CTCF.loss.5mC", as.is = T) %>% mutate(type = ifelse(grepl("NPC", sample), "NPC", ifelse(sample == "CEMT_23", "IDHwt", "IDHmut")), sample = gsub("NPC.", "", sample))
(CTCF_loss_5mC_figure <- ggplot(CTCF_loss_5mC, aes(sample, fractional, fill = type)) + 
		geom_boxplot() + 
		xlab("") + 
		ylab("Fractional methylation") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CTCF_loss_5mC_figure, file = "../CTCF/CTCF_loss_5mC_figure.pdf", height = 4, width = 5)

## ------- DMR glioma vs NPC -------
### -------- summary ---------
DMR_summary <- read.delim("./intermediate/DMR.summary.stats", head = T, as.is = T)
DMR_summary_tall <- data.frame(glioma = rep(gsub("_NPC.*", "", DMR_summary$sample), 2), NPC = rep(gsub(".*_NPC\\.", "", DMR_summary$sample), 2), DM = rep(c("hyper", "hypo"), each = nrow(DMR_summary)), length = c(DMR_summary$hyper, -DMR_summary$hypo)/10^6)
#DMR_summary_tall <- DMR_summary_tall %>% filter(glioma != "CEMT_21") %>% mutate(glioma = ifelse(glioma == "CEMT_23", paste0("IDHwt\n", glioma), paste0("IDHmut\n", glioma)))
DMR_summary_figure <- ggplot(DMR_summary_tall, aes(x = glioma, y = length, color = DM, shape = NPC)) + 
	geom_point(position = position_jitter(width = 0.1), size = 3) + 
	geom_hline(yintercept = 0) + 
	xlab("") + 
	ylab("Total DMR length (Mb)") + 
	scale_color_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw() + 
	theme(axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15))
ggsave(DMR_summary_figure, file = "DMR_summary_figure.pdf", height = 5, width = 6)

### ------- visualization ----- 
colname <- c("chr", "start", "end", "ID", "DM", "length")
for(lib in libs){
	print(lib)
	assign(paste0("DMR_", lib), read.delim(paste0("DMR.", lib, "_NPC"), head = F, as.is = T, col.names = colname))
	assign(paste0("DMR_", lib, "_figure"), DMR_figures(get(paste0("DMR_", lib)), lib, "NPCs", figures = c("length", "frequency", "circos"), colname = colname, hist_width = 3))
}

### -------- enrichment in genomic regions ------
genomic_breakdown <- read.delim("./CpG/genomic.breakdown.summary", as.is = T) %>% 
	mutate(sample = gsub("_NPC.*", "", Name), DM = gsub(".*NPC\\.", "", Name), NCpG = NULL, Name = NULL)
genomic_breakdown_tall <- melt(genomic_breakdown, id = c("sample", "DM")) %>% 
	mutate(value = ifelse(DM == "hyper", value, -value))
genomic_breakdown_figure <- ggplot(genomic_breakdown_tall, aes(variable, value, fill = DM)) + 
	geom_bar(position = "identity", stat = "identity", width = 0.5) + 
	geom_hline(yintercept = c(-2, 2)) + 
	facet_wrap(~sample) + 
	xlab("") + 
	ylab("Fold enrichment") + 
	scale_fill_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw()
ggsave(genomic_breakdown_figure, file = "genomic_breakdown_figure.pdf", height = 5, width = 6)

### -------- hyper CGI with K36me3 ---------
CGI_DMR_hyper_summary <- read.delim("./CGI/CGI.DMR.hyper.summary", as.is = T)
CGI_hyper_summary <- read.delim("./CGI/CGI.hyper.H3K36me3.summary", as.is = T)

### -------- % of hyper CpGs in hyper CGIs -----------
CGI_DMR_hyper_DM_all <- read.delim("./CGI/CGI.DMR.hyper.DM.all", as.is = T)
(CGI_DMR_hyper_DM_figure <- ggplot(CGI_DMR_hyper_DM_all, aes(glioma, percent, fill = NPC)) + 
		geom_boxplot(position = position_dodge()) + 
		xlab("")+
		ylab("Percent of hyper CpGs in each hyper CGIs") + 
		facet_grid(Type ~.) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(CGI_DMR_hyper_DM_figure, file = "./CGI/CGI_DMR_hyper_DM_figure.pdf")

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
(GREAT_DMR_CEMT_19_hyper <- enrich_GREAT("CEMT_19_hyper", "CEMT_19_hyper", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 8, width = 7))
(GREAT_DMR_CEMT_19_hypo <- enrich_GREAT("CEMT_19_hypo", "CEMT_19_hypo", categories = c("GOBP"), height = 3, width = 7))
(GREAT_DMR_CEMT_21_hyper <- enrich_GREAT("CEMT_21_hyper", "CEMT_21_hyper", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 8, width = 7))
(GREAT_DMR_CEMT_21_hypo <- enrich_GREAT("CEMT_21_hypo", "CEMT_21_hypo", categories = c("MSigPerturbation"), height = 4, width = 7))
(GREAT_DMR_CEMT_22_hyper <- enrich_GREAT("CEMT_22_hyper", "CEMT_22_hyper", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 8, width = 7))
# GREAT_DMR_CEMT_22_hypo : no enrichment
(GREAT_DMR_CEMT_23_hyper <- enrich_GREAT("CEMT_23_hyper", "CEMT_23_hyper", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 9, width = 7))
(GREAT_DMR_CEMT_23_hypo <- enrich_GREAT("CEMT_23_hypo", "CEMT_23_hypo", categories = c("GOBP", "DiseaseOntology", "InterPro"), height = 8, width = 7))
(GREAT_DMR_CEMT_47_hyper <- enrich_GREAT("CEMT_47_hyper", "CEMT_47_hyper", categories = c("GOBP", "DiseaseOntology", "GOCC"), height = 9, width = 7))
(GREAT_DMR_CEMT_47_hypo <- enrich_GREAT("CEMT_47_hypo", "CEMT_47_hypo", categories = c("MSigPerturbation"), height = 2, width = 7))

### -------- Intersect --------
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
DMR_intersect <- data.frame(sample1 = NA, sample2 = NA, hyper1 = NA, hyper2 = NA, hyper_intersect = NA, hyper_percent1 = NA, hyper_percent2 = NA, hyper_jaccard = NA, hypo1 = NA, hypo2 = NA, hypo_intersect = NA, hypo_percent1 = NA, hypo_percent2 = NA, hypo_jaccard = NA)
for(i1 in 1:(length(libs)-1)){
	for(i2 in (i1+1):length(libs)){
		print(c(i1, i2, libs[i1], libs[i2]))
		a1_hyper <- as.numeric(system(paste0("less DMR.", libs[i1], "_NPC.hyper.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		a2_hyper <- as.numeric(system(paste0("less DMR.", libs[i2], "_NPC.hyper.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		cross_hyper <- as.numeric(system(paste0(BEDTOOLS, "/intersectBed -a  DMR.", libs[i1], "_NPC.hyper.bed -b DMR.", libs[i2], "_NPC.hyper.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		jaccard_hyper <- as.numeric(str_split_fixed(system(paste0(BEDTOOLS, "bedtools jaccard -a DMR.", libs[i1], "_NPC.hyper.bed -b DMR.", libs[i2], "_NPC.hyper.bed"), intern = T)[2], pattern = "\\t", n = 4)[1,3])
		assign(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hyper"), 
					 draw.pairwise.venn(a1_hyper, a2_hyper, cross_hyper, category = c(libs[i1], libs[i2]), cat.pos = 0, ext.pos = 180))
		a1_hypo <- as.numeric(system(paste0("less DMR.", libs[i1], "_NPC.hypo.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		a2_hypo <- as.numeric(system(paste0("less DMR.", libs[i2], "_NPC.hypo.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		cross_hypo <- as.numeric(system(paste0(BEDTOOLS, "/intersectBed -a  DMR.", libs[i1], "_NPC.hypo.bed -b DMR.", libs[i2], "_NPC.hypo.bed | awk '{s=s+$3-$2}END{print s}'"), intern = T))
		jaccard_hypo <- as.numeric(str_split_fixed(system(paste0(BEDTOOLS, "bedtools jaccard -a DMR.", libs[i1], "_NPC.hypo.bed -b DMR.", libs[i2], "_NPC.hypo.bed"), intern = T)[2], pattern = "\\t", n = 4)[1,3])
		assign(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hypo"), 
					 draw.pairwise.venn(a1_hypo, a2_hypo, cross_hypo, category = c(libs[i1], libs[i2]), cat.pos = 0, ext.pos = 180))
		DMR_intersect <- rbind(DMR_intersect, data.frame(sample1 = libs[i1], sample2 = libs[i2], hyper1 = a1_hyper, hyper2 = a2_hyper, hyper_intersect = cross_hyper, hyper_percent1 = cross_hyper/a1_hyper, hyper_percent2 = cross_hyper/a2_hyper, hyper_jaccard = jaccard_hyper, hypo1 = a1_hypo, hypo2 = a2_hypo, hypo_intersect = cross_hypo, hypo_percent1 = cross_hypo/a1_hypo, hypo_percent2 = cross_hypo/a2_hypo, hypo_jaccard = jaccard_hypo))
	}
}
DMR_intersect <- na.omit(DMR_intersect)
DMR_jaccard_hyper <- data.frame(sample1 = c(DMR_intersect$sample1, DMR_intersect$sample2, libs), sample2 = c(DMR_intersect$sample2, DMR_intersect$sample1, libs), jaccard = c(rep(DMR_intersect$hyper_jaccard, 2), rep(1, length(libs)))) 
DMR_jaccard_hyper_figure <- ggplot(DMR_jaccard_hyper, aes(sample1, sample2, fill = jaccard)) + 
	geom_tile() + 
	xlab("") + 
	ylab("") + 
	ggtitle("Jaccard similarity - hyper DMRs") + 
	theme_bw()
ggsave(DMR_jaccard_hyper_figure, file = "DMR_jaccard_hyper_figure.pdf", height = 5, width = 6)
DMR_jaccard_hypo <- data.frame(sample1 = c(DMR_intersect$sample1, DMR_intersect$sample2, libs), sample2 = c(DMR_intersect$sample2, DMR_intersect$sample1, libs), jaccard = c(rep(DMR_intersect$hypo_jaccard, 2), rep(1, length(libs)))) 
DMR_jaccard_hypo_figure <- ggplot(DMR_jaccard_hypo, aes(sample1, sample2, fill = jaccard)) + 
	geom_tile() + 
	xlab("") + 
	ylab("") + 
	ggtitle("Jaccard similarity - hypo DMRs") + 
	theme_bw()
ggsave(DMR_jaccard_hypo_figure, file = "DMR_jaccard_hypo_figure.pdf", height = 5, width = 6)
pdf("Venn_DMR.pdf")
for(i1 in 1:(length(libs)-1)){
	for(i2 in (i1+1):length(libs)){
		grid.arrange(gTree(children = get(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hyper"))), gTree(children = get(paste0("Venn_DMR_", libs[i1], ".", libs[i2], "_hypo"))), nrow = 2)
		grid.text("hyper", x = unit(0.1, "npc"), y = unit(0.98, "npc"))
		grid.text("hypo", x = unit(0.1, "npc"), y = unit(0.48, "npc"))
	}
}
dev.off()

### -------- enrichment in chromatin states --------
DMR_ChromHMM_summary <- read.delim("./CpG/DMR.chromHMM.enrich.summary", as.is = T) 
(DMR_ChromHMM_summary_figure <- ggplot(DMR_ChromHMM_summary, aes(Name, Enrichment, fill = Sample)) + 
	geom_bar(stat = "identity", position = position_dodge()) + 
	facet_grid(. ~ DM) + 
	coord_flip() + 
	xlab("") + 
	ylab("Fold enrichment") +
	theme_bw())
ggsave(DMR_ChromHMM_summary_figure, file = "DMR_ChromHMM_summary_figure.pdf")

### -------- enrichment in differentially marked histone mods -------
DMR_DHM_enrich <- read.delim("./DMR.uniqueHM.enrichment.summary", as.is = T) %>% 
	mutate(HM = gsub("CEMT.*.unique", "Histone modification gain", HM), HM = gsub("NPC.*.unique", "Histone modification loss", HM), sig = ifelse(p.value <= 0.01, "p-value <= 0.01", "p-value > 0.01"))
(DMR_DHM_enrich_figure <- ggplot(DMR_DHM_enrich, aes(Sample, log2(Fold), fill = Mark, color = sig)) + 
	geom_bar(stat = "identity", position = position_dodge(), width = 0.8) + 
	scale_color_manual(name = "", values = c("red", "transparent")) + 
	facet_grid(HM ~ DMR) + 
	coord_flip() + 
	xlab("") + 
	ylab("log2 Fold Enrichemnt") + 
	theme_bw())
ggsave(DMR_DHM_enrich_figure, file = "DMR_DHM_enrich_figure.pdf")
DMR_enhancer_enrich <- DMR_DHM_enrich %>% filter(Sample != "CEMT_21", Mark == "H3K27ac")
pdf("DMR_enhancer_enrich_venn.pdf")
DMR_enhancer_enrich_venn <- draw.pairwise.venn(area1 = 9128, area2 = 10940, cross.area = 826, category = c("Hypermethylation", "Loss of H3K27ac"), fill = c("orange", "blue"), cat.pos = c(200, 160), cat.cex = 1.5)
dev.off()

### -------- associated with DE genes ----------
DMR_DE <- read.delim("./DE/DMR.DE.summary", as.is = T) %>% mutate(Significant = p_Fisher < 0.01)
(DMR_DE_figure <- ggplot(DMR_DE, aes(DM, Percent_intersect, fill = DE, color = Significant)) + 
	geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
	scale_fill_manual(name = "", values = c("blue", "red")) + 
	scale_color_manual(values = c("transparent", "green")) + 
	facet_wrap(~ Sample) + 
	xlab("") + 
	ylab("Fraction of DE genes") + 
	theme_bw())
ggsave(DMR_DE_figure, file = "DMR_DE_figure.pdf")

save(list = c("DMR_intersect",  
							ls(pattern = "summary"), ls(pattern = "figure"), 
							ls(pattern = "GREAT_DMR_*"), ls(pattern = "Venn_DMR_*")),
		 file = "/projects/epigenomics2/users/lli/glioma/WGBS/WGBS.Rdata")

