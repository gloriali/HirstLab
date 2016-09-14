# MeDIP and MRE data from Kongkham's Lab
library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)
setwd("/projects/epigenomics2/users/lli/glioma/Kongkham/")
load("/projects/epigenomics2/users/lli/glioma/Kongkham/Kongkham.Rdata")

## ====== QC ======
sample_info <- read.delim("Sample_info.txt", as.is = T)
rownames(sample_info) <- sample_info$ID
qc_rep <- read.delim("QC_repetitive", as.is = T) 
rownames(qc_rep) <- qc_rep$Library
qc <- read.delim("QC_summary.txt", as.is = T) %>% 
	mutate(Number_Multiple_Alignment = qc_rep[Library, "Number_Repetitive_Mapping"], Library = gsub("_[0-9]+_R1_trimmed", "", Library), Type = gsub("hmedip", "hMC", gsub(".*_", "", Library)), 
				 Library = gsub("_(MC|hMC|Input)", "", Library), SampleType = sample_info[Library, "SampleType"], IDH1mut = sample_info[Library, "IDH1mutation"], Percent_Multiple_Alignment = Number_Multiple_Alignment/Total_Number_Of_Reads*100)
qc[is.na(qc$SampleType), "SampleType"] <- "Cell line - neural stem cell"
qc[is.na(qc$IDH1mut), "IDH1mut"] <- "negative"
qc_long <- qc %>% select(-Estim_X_coverage) %>% melt(id = c("Library", "Type", "SampleType", "IDH1mut")) %>% 
	mutate(variable = gsub("Uniquely_Aligned_Reads_without_Dups_and_Q_.._10", "Reads_After_Filter", variable), variable = factor(variable, levels = c("Total_Number_Of_Reads", "Number_Reads_Aligned", "Mapping_Efficiency", "Number_Duplicates", "Percent_of_dups", "Number_Multiple_Alignment", "Percent_Multiple_Alignment", "Number_Reads_After_Filter", "Percent_Reads_After_Filter")))
(qc_figure <- ggplot(qc_long, aes(Library, value, fill = SampleType)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	facet_grid(variable ~ Type, scale = "free") + 
	xlab("") + 
	ylab("") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, size = 4)))
ggsave(qc_figure, file = "qc_figure.pdf", height = 20, width = 12)
(qc_Nmap_figure <- ggplot(qc_long %>% filter(variable %in% c("Total_Number_Of_Reads", "Number_Reads_Aligned")), aes(Library, value, fill = SampleType)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	facet_grid(variable ~ Type, scale = "free") + 
	scale_y_continuous(breaks = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8), labels = c(0, 20, 40, 60, 80, 100, 120)) + 
	xlab("") + 
	ylab("# million reads") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, size = 4)))
ggsave(qc_Nmap_figure, file = "qc_Nmap_figure.pdf", height = 6, width = 12)
(qc_Percentmap_figure <- ggplot(qc_long %>% filter(variable %in% c("Mapping_Efficiency")), aes(Library, value, fill = SampleType)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	facet_grid(variable ~ Type, scale = "free") + 
	xlab("") + 
	ylab("Mapped%") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, size = 4)))
ggsave(qc_Percentmap_figure, file = "qc_Percentmap_figure.pdf", height = 4, width = 12)
(qc_Ndup_figure <- ggplot(qc_long %>% filter(variable %in% c("Number_Duplicates")), aes(Library, value, fill = SampleType)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	facet_grid(variable ~ Type, scale = "free") + 
	scale_y_continuous(breaks = c(0, 1e7, 2e7, 3e7, 4e7), labels = c(0, 10, 20, 30, 40)) + 
	xlab("") + 
	ylab("# million reads") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, size = 4)))
ggsave(qc_Ndup_figure, file = "qc_Ndup_figure.pdf", height = 4, width = 12)
(qc_Percentdup_figure <- ggplot(qc_long %>% filter(variable %in% c("Percent_of_dups")), aes(Library, value, fill = SampleType)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	facet_grid(variable ~ Type, scale = "free") + 
	xlab("") + 
	ylab("Duplicate%") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, size = 4)))
ggsave(qc_Percentdup_figure, file = "qc_Percentdup_figure.pdf", height = 4, width = 12)
(qc_Nuniq_figure <- ggplot(qc_long %>% filter(variable %in% c("Number_Reads_After_Filter")), aes(Library, value, fill = SampleType)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	facet_grid(variable ~ Type, scale = "free") + 
	scale_y_continuous(breaks = c(0, 2e7, 4e7, 6e7), labels = c(0, 20, 40, 60)) + 
	xlab("") + 
	ylab("# million reads") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, size = 4)))
ggsave(qc_Nuniq_figure, file = "qc_Nuniq_figure.pdf", height = 4, width = 12)
(qc_Percentuniq_figure <- ggplot(qc_long %>% filter(variable %in% c("Percent_Reads_After_Filter")), aes(Library, value, fill = SampleType)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	facet_grid(variable ~ Type, scale = "free") + 
	xlab("") + 
	ylab("Uniquely aligned%") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, size = 4)))
ggsave(qc_Percentuniq_figure, file = "qc_Percentuniq_figure.pdf", height = 4, width = 12)

## ====== CpG coverage ======
CG_coverage_MC <- read.delim("./bam/CG.coverage.MC", as.is = T) %>% melt(id = "Coverage") %>% filter(Coverage != ">50") %>%
	mutate(Coverage = as.numeric(Coverage), Library = gsub("X", "", variable), SampleType = sample_info[Library, "SampleType"], IDH1mut = sample_info[Library, "IDH1mutation"])
(CG_coverage_MC_figure <- ggplot(CG_coverage_MC, aes(Coverage, value/10^6, color = SampleType)) + 
	geom_point() + 
	geom_line(aes(group = Library)) + 
	coord_cartesian(xlim = c(0, 30)) + 
	ylab("# million CpGs") + 
	ggtitle("MC CpG Coverage Distribution") + 
	theme_bw())
ggsave(CG_coverage_MC_figure, file = "CG_coverage_MC_figure.pdf", height = 6, width = 8)
CG_coverage_MC_filter_summary <- CG_coverage_MC %>% filter(Coverage >= 3) %>% group_by(Library) %>% 
	summarize(SampleType = SampleType[1], N = sum(value))
(CG_coverage_MC_filter_summary_figure <- ggplot(CG_coverage_MC_filter_summary, aes(Library, N/10^6, fill = SampleType)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	scale_y_continuous(breaks = seq(0, 16, by = 2)) + 
	xlab("") + 
	ylab("# million CpGs") + 
	ggtitle("MC CpG Coverage >= 3") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, size = 4)))
ggsave(CG_coverage_MC_filter_summary_figure, file = "CG_coverage_MC_filter_summary_figure.pdf", height = 6, width = 10)

CG_coverage_hMC <- read.delim("./bam/CG.coverage.hMC", as.is = T) %>% melt(id = "Coverage") %>% filter(Coverage != ">50") %>%
	mutate(Coverage = as.numeric(Coverage), Library = gsub("X", "", variable), SampleType = sample_info[Library, "SampleType"], IDH1mut = sample_info[Library, "IDH1mutation"])
(CG_coverage_hMC_figure <- ggplot(CG_coverage_hMC, aes(Coverage, value/10^6, color = SampleType)) + 
	geom_point() + 
	geom_line(aes(group = Library)) + 
	coord_cartesian(xlim = c(0, 30)) + 
	ylab("# million CpGs") + 
	ggtitle("hMC CpG Coverage Distribution") + 
	theme_bw())
ggsave(CG_coverage_hMC_figure, file = "CG_coverage_hMC_figure.pdf", height = 6, width = 8)
CG_coverage_hMC_filter_summary <- CG_coverage_hMC %>% filter(Coverage >= 3) %>% group_by(Library) %>% 
	summarize(SampleType = SampleType[1], N = sum(value))
(CG_coverage_hMC_filter_summary_figure <- ggplot(CG_coverage_hMC_filter_summary, aes(Library, N/10^6, fill = SampleType)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	scale_y_continuous(breaks = seq(0, 16, by = 2)) + 
	xlab("") + 
	ylab("# million CpGs") + 
	ggtitle("hMC CpG Coverage >= 3") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, size = 4)))
ggsave(CG_coverage_hMC_filter_summary_figure, file = "CG_coverage_hMC_filter_summary_figure.pdf", height = 6, width = 10)

save(list = c("qc", ls(pattern = "_figure")),
		 file = "/projects/epigenomics2/users/lli/glioma/Kongkham/Kongkham.Rdata")

