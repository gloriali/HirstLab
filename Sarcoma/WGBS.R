# Sarcoma WGBS analysis

library(ggplot2)
library(plyr)
library(VennDiagram)
library(gridExtra)
library(gplots)
library(reshape2)
library(wq)
library(dplyr)
library(RCircos)
library(stringr)
source('~/HirstLab/Pipeline/R/DMR.figures.R')
source('~/HirstLab/Pipeline/R/enrich.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
load("/projects/epigenomics3/epigenomics3_results/alorzadeh/Sarcoma/WGBS/WGBS.Rdata")
setwd("/projects/epigenomics3/epigenomics3_results/alorzadeh/Sarcoma/WGBS/")

## ========== global profile ==========
quantile_5mC <- read.delim("qc.5mC.quantile", as.is = T) %>% mutate(category = gsub("\\..*", "", sample), category = ifelse(grepl("DG", category), "Sarcoma", category))
(quantile_5mC_figure <- ggplot(quantile_5mC, aes(x = sample, fill = category, lower = lower, middle = median, upper = upper, ymin = ymin, ymax = ymax)) + 
		geom_boxplot(stat = "identity") + 
		facet_wrap(~type) + 
		scale_fill_manual(values = c("MSC" = "salmon", "H1" = "gray", "NPC" = "sky blue", "Sarcoma" = "orange"), name = "") + 
		xlab("") + 
		ylab("Fractional methylation") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(quantile_5mC_figure, file = "quantile_5mC_figure.pdf", height = 4, width = 5)
CGI_promoter_5mC <- read.delim("CGI.promoter.5mC", as.is = T, col.names = c("ID", "fractional", "sample")) %>% dcast(ID ~ sample, value.var = "fractional") %>% na.omit
write.table(CGI_promoter_5mC, file = "CGI_promoter_5mC.txt", sep = "\t", quote = F, row.names = F, col.names = T)

## ============ DMRs ================
### -------- summary ---------
DMR_summary <- read.delim("./DMR/DMR.summary.stats", head = T, as.is = T) %>% 
	mutate(sarcoma = gsub("_.*", "", sample), normal = gsub(".*_", "", sample), hypo = -hypo) %>%
	select(sarcoma, normal, hyper, hypo) %>% melt(id = c("sarcoma", "normal"), value.name = "length", variable.name = "DM")
(DMR_summary_figure <- ggplot(DMR_summary, aes(x = sarcoma, y = length/1e6, color = DM, shape = normal)) + 
	geom_point(position = position_jitter(width = 0.1), size = 3) + 
	geom_hline(yintercept = 0) + 
	xlab("") + 
	ylab("Total DMR length (Mb)") + 
	scale_color_manual(name = "", values = c("red", "blue")) + 
	coord_flip(ylim = c(-180, 180)) + 
	theme_bw())
ggsave(DMR_summary_figure, file = "./DMR/DMR_summary_figure.pdf", height = 5, width = 6)

### ------- visualization ----- 
colname <- c("chr", "start", "end", "ID", "DM", "count", "length")
pairs <- data.frame(s1 = rep(c("DG1343a", "DG1344", "DG1346", "DG1348", "DG1349"), each = 3), s2 = rep(c("H1", "MSC", "NPC"), 5)) %>% 
	mutate(lib = paste0(s1, "_", s2))
for(i in 1:nrow(pairs)){
	s1 <- as.character(pairs[i, "s1"])
	s2 <- as.character(pairs[i, "s2"])
	lib <- as.character(pairs[i, "lib"])
	print(lib)
	assign(paste0("DMR_", lib), read.delim(paste0("./DMR/DMR.", lib, ".s500.c3"), head = F, as.is = T, col.names = colname) %>% mutate(chr = paste0("chr", chr)))
	assign(paste0("DMR_", lib, "_figure"), DMR_figures(get(paste0("DMR_", lib)), s1, s2, dirOut = "./DMR/", figures = c("length", "frequency", "circos"), colname = colname, hist_width = 3))
}

### -------- enrichment in genomic regions ------
genomic_breakdown <- read.delim("./DMR/CpG/genomic.breakdown.summary", as.is = T) %>% 
	mutate(sarcoma = gsub("_.*", "", Name), normal = gsub(".*_", "", gsub("\\..*", "", Name)), DM = gsub(".*c3\\.", "", Name), NCpG = NULL, Name = NULL)
genomic_breakdown_tall <- melt(genomic_breakdown, id = c("sarcoma", "normal", "DM")) %>% 
	mutate(value = ifelse(DM == "hyper", value, -value))
(genomic_breakdown_figure <- ggplot(genomic_breakdown_tall, aes(variable, value, fill = DM)) + 
	geom_bar(position = "identity", stat = "identity", width = 0.5) + 
	geom_hline(yintercept = c(-2, 2)) + 
	facet_grid(sarcoma ~ normal) + 
	xlab("") + 
	ylab("Fold enrichment") + 
	scale_fill_manual(name = "", values = c("red", "blue")) + 
	coord_flip() + 
	theme_bw())
ggsave(genomic_breakdown_figure, file = "./DMR/genomic_breakdown_figure.pdf", height = 8, width = 6)



save(list = c(ls(pattern = "figure"), "CGI_promoter_5mC"), 
		 file = "/projects/epigenomics3/epigenomics3_results/alorzadeh/Sarcoma/WGBS/WGBS.Rdata")

