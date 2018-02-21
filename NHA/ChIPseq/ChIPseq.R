# NHA VitC - ChIP-seq analysis

library(ggplot2)
library(plyr)
library(VennDiagram)
library(grid)
library(gridExtra)
library(gplots)
library(dendextend)
library(reshape2)
library(dplyr)
library(RCircos)
library(stringr)
library(scales)
source('~/HirstLab/Pipeline/R/enrich.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
load("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/ChIPseq.Rdata")
setwd("/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/")
marks <- c("H3K27me3", "H3K9me3", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3")
RPKM <- read.delim("../RNAseq/RPKM/vitc.RPKM", as.is = T)
rownames(RPKM) <- RPKM$ENSG

## ============= 5hmC stain ===========
stain_5hmC <- read.delim("../5hmC.stain", as.is = T) %>% melt(id = "Sample") %>% mutate(Sample = gsub("_", "\n", Sample))
(stain_5hmC_figure <- ggplot(stain_5hmC, aes(Sample, value, fill = variable)) + 
		geom_bar(stat = "identity", position = position_dodge()) + 
		scale_fill_manual(values = c("blue", "red"), name = "") + 
		xlab("") + 
		ylab("Mean intensity") + 
		ggtitle("5hmC positive staining") + 
		theme_bw() + 
		theme(axis.text.x = element_text(size = 12)))
ggsave(stain_5hmC_figure, file = "../stain_5hmC_figure.pdf", height = 3, width = 4)

## ============= QC =============
QC_summary <- read.delim("./bam/QC_summary.txt", as.is = T) %>% mutate(Mark = gsub("_.*", "", Library), Sample = gsub(".*_NHA", "NHA", Library), Sample = gsub(".*_MGG.*_", "MGG119_", Sample))
QC_summary_reads <- QC_summary %>% select(Mark, Sample, Total_Number_Of_Reads, Number_Reads_Aligned, Number_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = c("Mark", "Sample"))
(QC_summary_reads_figure <- ggplot(QC_summary_reads, aes(Sample, value/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ Mark, scales = "free_x") + 
		xlab("") + 
		ylab("No. of million reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_reads_figure, file = "./bam/QC_summary_reads_figure.pdf", height = 8, width = 8)
QC_summary_percent <- QC_summary %>% select(Mark, Sample, Mapping_Efficiency, Percent_of_dups, Percent_of_Paired_Alignments, Percent_Uniquely_Aligned_Reads_After_Filter) %>% melt(id.var = c("Mark", "Sample"))
(QC_summary_percent_figure <- ggplot(QC_summary_percent, aes(Sample, value, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(variable ~ Mark, scales = "free") + 
		xlab("") + 
		ylab("Percent of total number of reads") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_percent_figure, file = "./bam/QC_summary_percent_figure.pdf", height = 10, width = 8)
QC_target_region <- read.delim("./bam/QC.target_regions.txt", as.is = T)
QC_summary_domain <- merge(QC_summary, QC_target_region) %>% mutate(Percent_Domain_Reads = Number_Target_Regions*100/Number_Reads_Aligned)  %>% filter(!grepl("input", Mark))
(QC_summary_domain_figure <- ggplot(QC_summary_domain, aes(Sample, Percent_Domain_Reads, fill = Sample)) + 
		geom_bar(stat = "identity", width = 0.5) + 
		facet_grid(Mark ~ ., scales = "free_y") + 
		xlab("") + 
		ylab("Domain reads (% of mapped)") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(QC_summary_domain_figure, file = "./bam/QC_summary_domain_figure.pdf", height = 7, width = 4)
### wrong IP for NHAR_H3K9me3
NHAR_H3K9me3_H3K4me3 <- read.delim("./bam/NHAR.H3K9me3.H3K4me3.bed", as.is = T, head = F, col.names = c("chr", "start", "end", "ID", "N", "Sample", "RPKM")) %>% select(ID, Sample, RPKM) %>%
	reshape(direction = "wide", idvar = "ID", timevar = "Sample")
(NHAR_H3K9me3_H3K4me3_figure <- ggplot(NHAR_H3K9me3_H3K4me3, aes(log10(RPKM.H3K9me3_NHAR_vitc), log10(RPKM.H3K4me3_NHAR_vitc))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.1, 0.2, 0.3, 0.04, 0.5, 0.8, 1), guide = "none") +	
		coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) + 
		ggtitle(cor(NHAR_H3K9me3_H3K4me3$RPKM.H3K9me3_NHAR_vitc, NHAR_H3K9me3_H3K4me3$RPKM.H3K4me3_NHAR_vitc)) + 
		theme_bw())
ggsave(NHAR_H3K9me3_H3K4me3_figure, file = "./bam/NHAR_H3K9me3_H3K4me3_figure.pdf", width = 5, height = 5)
(H3K9me3_NHAR_vitC_NTC_figure <- ggplot(NHAR_H3K9me3_H3K4me3, aes(log10(RPKM.H3K9me3_NHAR_vitc), log10(RPKM.H3K9me3_NHAR_control))) + 
		stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) +
		scale_fill_gradientn(colors = c("black", "darkblue", "lightblue", "green", "yellow", "orange", "red", "darkred"), values = c(0, 0.1, 0.2, 0.3, 0.04, 0.5, 0.8, 1), guide = "none") +	
		coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) + 
		ggtitle(cor(NHAR_H3K9me3_H3K4me3$RPKM.H3K9me3_NHAR_vitc, NHAR_H3K9me3_H3K4me3$RPKM.H3K9me3_NHAR_control)) + 
		theme_bw())
ggsave(H3K9me3_NHAR_vitC_NTC_figure, file = "./bam/H3K9me3_NHAR_vitC_NTC_figure.pdf", width = 5, height = 5)
NHAR_H3K9me3_H3K4me3_venn <- draw.pairwise.venn(area1 = 44496864, area2 = 54317631, cross.area = 41973315, category = c("H3K9me3_NHAR_vitc", "H3K4me3_NHAR_vitc"), fill = c("red", "blue"))
H3K9me3_NHAR_vitC_NTC_venn <- draw.pairwise.venn(area1 = 44496864, area2 = 121812813, cross.area = 41856862, category = c("H3K9me3_NHAR_vitc", "H3K9me3_NHAR_control"), fill = c("red", "blue"))

## ============= ER =============
ER_summary_MACS2 <- read.delim("./MACS2/ER_summary.txt", as.is = T)
(ER_summary_MACS2_figure <- ggplot(ER_summary_MACS2, aes(Sample, Total_length/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(Mark ~., scales = "free_y") + 
		xlab("") + 
		ylab("Total length (MB)") + 
		ggtitle("MACS2 ER") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(ER_summary_MACS2_figure, file = "./MACS2/ER_summary_MACS2_figure.pdf", height = 6, width = 5)
ER_summary_FindER <- read.delim("./FindER/ER_summary.txt", as.is = T) %>% mutate(Sample = gsub("MGG_", "MGG119_", Sample))
(ER_summary_FindER_figure <- ggplot(ER_summary_FindER %>% filter(Mark != "H3K9me3"), aes(Sample, Total_length/10^6, fill = Sample)) + 
		geom_bar(stat = "identity") + 
		facet_grid(Mark ~., scales = "free_y") + 
		xlab("") + 
		ylab("Total length (MB)") + 
		ggtitle("FindER ER") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(ER_summary_FindER_figure, file = "./FindER/ER_summary_FindER_figure.pdf", height = 6, width = 5)

## ============= unique ER =============
ER_unique_summary <- read.delim("./unique2/ER.unique.summary", as.is = T) %>% 
	melt(id = c("Sample1", "Sample2", "Mark")) %>% 
	mutate(Sample = gsub("MGG_", "MGG119_", paste0(Sample1, ".", Sample2)), category = ifelse(grepl("len", variable), "Total No. of bases", "Total No. of regions"), unique = ifelse(grepl("unique", variable), T, F), value = ifelse(grepl("1", variable), value, -value))
(ER_unique_summary_figure <- ggplot(ER_unique_summary %>% filter(unique == T, Sample %in% c("NHAR_vitc.NHAR_control", "NHAR_control.NHA_control", "MGG119_vitc.MGG119_control")), aes(Mark, value, fill = Sample)) + 
		geom_bar(stat = "identity", position = position_dodge()) +
		geom_hline(yintercept = 0) + 
		facet_grid(category ~., scales = "free_y") + 
		xlab("") + 
		ylab("") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90)))
ggsave(ER_unique_summary_figure, file = "./unique2/ER_unique_summary_figure.pdf", width = 7, height = 5)

### unique K27ac homer
TF_ENSG <- read.delim("./unique2/TF.ENSG", as.is = T) %>% merge(RPKM)
homer_unique_enhancer_NHARcontrol_NHARvitc <- read.delim("./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/NHAR_control/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) %>% merge(TF_ENSG)
(homer_unique_enhancer_NHARcontrol_NHARvitc_figure <- ggplot(homer_unique_enhancer_NHARcontrol_NHARvitc, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("NHARvitc_NHARcontrol NHARcontrol unique") + 
		theme_bw())
ggsave(homer_unique_enhancer_NHARcontrol_NHARvitc_figure, file = "./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_NHARcontrol_NHARvitc_figure.pdf", height = 3, width = 5)
homer_unique_enhancer_NHARcontrol_NHARvitc_RPKM <- homer_unique_enhancer_NHARcontrol_NHARvitc %>% select(TF, MGG119_control:NHA_vitc) %>% melt(id = "TF") %>% mutate(TF = factor(TF, levels = levels(homer_unique_enhancer_NHARcontrol_NHARvitc$TF)))
(homer_unique_enhancer_NHARcontrol_NHARvitc_RPKM_figure <- ggplot(homer_unique_enhancer_NHARcontrol_NHARvitc_RPKM, aes(TF, log10(value), color = variable)) + 
		geom_point(size = 2, position = position_jitter(width = 0.1)) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("NHARvitc_NHARcontrol NHARcontrol unique") + 
		guides(color = guide_legend(title = "")) + 
		theme_bw() + 
		theme(axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(homer_unique_enhancer_NHARcontrol_NHARvitc_RPKM_figure, file = "./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_NHARcontrol_NHARvitc_RPKM_figure.pdf", height = 3, width = 5)
pdf("./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_NHARcontrol_NHARvitc_figure.pdf", height = 3, width = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(homer_unique_enhancer_NHARcontrol_NHARvitc_figure + theme(plot.margin = unit(c(1,0,1,1), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(homer_unique_enhancer_NHARcontrol_NHARvitc_RPKM_figure + ggtitle("") + theme(plot.margin = unit(c(1,1,1,-0.5), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
homer_unique_enhancer_NHARvitc_NHARcontrol <- read.delim("./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/NHAR_vitc/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) %>% merge(TF_ENSG)
(homer_unique_enhancer_NHARvitc_NHARcontrol_figure <- ggplot(homer_unique_enhancer_NHARvitc_NHARcontrol, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("NHARvitc_NHARcontrol NHARvitc unique") + 
		theme_bw())
ggsave(homer_unique_enhancer_NHARvitc_NHARcontrol_figure, file = "./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_NHARvitc_NHARcontrol_figure.pdf", height = 2, width = 5)
homer_unique_enhancer_NHARvitc_NHARcontrol_RPKM <- homer_unique_enhancer_NHARvitc_NHARcontrol %>% select(TF, MGG119_control:NHA_vitc) %>% melt(id = "TF") %>% mutate(TF = factor(TF, levels = levels(homer_unique_enhancer_NHARvitc_NHARcontrol$TF)))
(homer_unique_enhancer_NHARvitc_NHARcontrol_RPKM_figure <- ggplot(homer_unique_enhancer_NHARvitc_NHARcontrol_RPKM, aes(TF, log10(value), color = variable)) + 
		geom_point(size = 2, position = position_jitter(width = 0.1)) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("NHARvitc_NHARcontrol NHARcontrol unique") + 
		guides(color = guide_legend(title = "")) + 
		theme_bw() + 
		theme(axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(homer_unique_enhancer_NHARvitc_NHARcontrol_RPKM_figure, file = "./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_NHARvitc_NHARcontrol_RPKM_figure.pdf", height = 2, width = 5)
pdf("./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_NHARvitc_NHARcontrol_figure.pdf", height = 2, width = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(homer_unique_enhancer_NHARvitc_NHARcontrol_figure + theme(plot.margin = unit(c(1,0,1,1), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(homer_unique_enhancer_NHARvitc_NHARcontrol_RPKM_figure + ggtitle("") + theme(plot.margin = unit(c(1,1,1,-0.5), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
homer_unique_enhancer_NHARcontrol_NHAcontrol <- read.delim("./unique2/NHAR_control.NHA_control/H3K27ac/Homer2/NHAR_control/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) %>% merge(TF_ENSG)
(homer_unique_enhancer_NHARcontrol_NHAcontrol_figure <- ggplot(homer_unique_enhancer_NHARcontrol_NHAcontrol, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("NHARcontrol_NHAcontrol NHARcontrol unique") + 
		theme_bw())
ggsave(homer_unique_enhancer_NHARcontrol_NHAcontrol_figure, file = "./unique2/NHAR_control.NHA_control/H3K27ac/Homer2/homer_unique_enhancer_NHARcontrol_NHAcontrol_figure.pdf", height = 3, width = 5)
homer_unique_enhancer_NHARcontrol_NHAcontrol_RPKM <- homer_unique_enhancer_NHARcontrol_NHAcontrol %>% select(TF, MGG119_control:NHA_vitc) %>% melt(id = "TF") %>% mutate(TF = factor(TF, levels = levels(homer_unique_enhancer_NHARcontrol_NHAcontrol$TF)))
(homer_unique_enhancer_NHARcontrol_NHAcontrol_RPKM_figure <- ggplot(homer_unique_enhancer_NHARcontrol_NHAcontrol_RPKM, aes(TF, log10(value), color = variable)) + 
		geom_point(size = 2, position = position_jitter(width = 0.1)) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("NHARvitc_NHARcontrol NHARcontrol unique") + 
		guides(color = guide_legend(title = "")) + 
		theme_bw() + 
		theme(axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(homer_unique_enhancer_NHARcontrol_NHAcontrol_RPKM_figure, file = "./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_NHARcontrol_NHAcontrol_RPKM_figure.pdf", height = 3, width = 5)
pdf("./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_NHARcontrol_NHAcontrol_figure.pdf", height = 3, width = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(homer_unique_enhancer_NHARcontrol_NHAcontrol_figure + theme(plot.margin = unit(c(1,0,1,1), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(homer_unique_enhancer_NHARcontrol_NHAcontrol_RPKM_figure + ggtitle("") + theme(plot.margin = unit(c(1,1,1,-0.5), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
homer_unique_enhancer_NHAcontrol_NHARcontrol <- read.delim("./unique2/NHAR_control.NHA_control/H3K27ac/Homer2/NHA_control/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) %>% merge(TF_ENSG)
(homer_unique_enhancer_NHAcontrol_NHARcontrol_figure <- ggplot(homer_unique_enhancer_NHAcontrol_NHARcontrol, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("NHARcontrol_NHAcontrol NHAcontrol unique") + 
		theme_bw())
ggsave(homer_unique_enhancer_NHAcontrol_NHARcontrol_figure, file = "./unique2/NHAR_control.NHA_control/H3K27ac/Homer2/homer_unique_enhancer_NHAcontrol_NHARcontrol_figure.pdf", height = 2, width = 5)
homer_unique_enhancer_NHAcontrol_NHARcontrol_RPKM <- homer_unique_enhancer_NHAcontrol_NHARcontrol %>% select(TF, MGG119_control:NHA_vitc) %>% melt(id = "TF") %>% mutate(TF = factor(TF, levels = levels(homer_unique_enhancer_NHAcontrol_NHARcontrol$TF)))
(homer_unique_enhancer_NHAcontrol_NHARcontrol_RPKM_figure <- ggplot(homer_unique_enhancer_NHAcontrol_NHARcontrol_RPKM, aes(TF, log10(value), color = variable)) + 
		geom_point(size = 2, position = position_jitter(width = 0.1)) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("NHARvitc_NHARcontrol NHARcontrol unique") + 
		guides(color = guide_legend(title = "")) + 
		theme_bw() + 
		theme(axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(homer_unique_enhancer_NHAcontrol_NHARcontrol_RPKM_figure, file = "./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_NHAcontrol_NHARcontrol_RPKM_figure.pdf", height = 2, width = 5)
pdf("./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_NHAcontrol_NHARcontrol_figure.pdf", height = 2, width = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(homer_unique_enhancer_NHAcontrol_NHARcontrol_figure + theme(plot.margin = unit(c(1,0,1,1), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(homer_unique_enhancer_NHAcontrol_NHARcontrol_RPKM_figure + ggtitle("") + theme(plot.margin = unit(c(1,1,1,-0.5), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
homer_unique_enhancer_MGG119control_MGG119vitc <- read.delim("./unique2/MGG_vitc.MGG119_control/H3K27ac/Homer2/MGG119_control/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) %>% merge(TF_ENSG)
(homer_unique_enhancer_MGG119control_MGG119vitc_figure <- ggplot(homer_unique_enhancer_MGG119control_MGG119vitc, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("MGG119vitc_MGG119control MGG119control unique") + 
		theme_bw())
ggsave(homer_unique_enhancer_MGG119control_MGG119vitc_figure, file = "./unique2/MGG_vitc.MGG119_control/H3K27ac/Homer2/homer_unique_enhancer_MGG119control_MGG119vitc_figure.pdf", height = 2, width = 5)
homer_unique_enhancer_MGG119control_MGG119vitc_RPKM <- homer_unique_enhancer_MGG119control_MGG119vitc %>% select(TF, MGG119_control:NHA_vitc) %>% melt(id = "TF") %>% mutate(TF = factor(TF, levels = levels(homer_unique_enhancer_MGG119control_MGG119vitc$TF)))
(homer_unique_enhancer_MGG119control_MGG119vitc_RPKM_figure <- ggplot(homer_unique_enhancer_MGG119control_MGG119vitc_RPKM, aes(TF, log10(value), color = variable)) + 
		geom_point(size = 2, position = position_jitter(width = 0.1)) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("NHARvitc_NHARcontrol NHARcontrol unique") + 
		guides(color = guide_legend(title = "")) + 
		theme_bw() + 
		theme(axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(homer_unique_enhancer_MGG119control_MGG119vitc_RPKM_figure, file = "./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_MGG119control_MGG119vitc_RPKM_figure.pdf", height = 2, width = 5)
pdf("./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_MGG119control_MGG119vitc_RPKM_figure.pdf", height = 2, width = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(homer_unique_enhancer_MGG119control_MGG119vitc_figure + theme(plot.margin = unit(c(1,0,1,1), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(homer_unique_enhancer_MGG119control_MGG119vitc_RPKM_figure + ggtitle("") + theme(plot.margin = unit(c(1,1,1,-0.5), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
homer_unique_enhancer_MGG119vitc_MGG119control <- read.delim("./unique2/MGG_vitc.MGG119_control/H3K27ac/Homer2/MGG_vitc/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) %>% merge(TF_ENSG)
(homer_unique_enhancer_MGG119vitc_MGG119control_figure <- ggplot(homer_unique_enhancer_MGG119vitc_MGG119control, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("MGG119vitc_MGG119control MGG119vitc unique") + 
		theme_bw())
ggsave(homer_unique_enhancer_MGG119vitc_MGG119control_figure, file = "./unique2/MGG_vitc.MGG119_control/H3K27ac/Homer2/homer_unique_enhancer_MGG119vitc_MGG119control_figure.pdf", height = 5, width = 5)
homer_unique_enhancer_MGG119vitc_MGG119control_RPKM <- homer_unique_enhancer_MGG119vitc_MGG119control %>% select(TF, MGG119_control:NHA_vitc) %>% melt(id = "TF") %>% mutate(TF = factor(TF, levels = levels(homer_unique_enhancer_MGG119vitc_MGG119control$TF)))
(homer_unique_enhancer_MGG119vitc_MGG119control_RPKM_figure <- ggplot(homer_unique_enhancer_MGG119vitc_MGG119control_RPKM, aes(TF, log10(value), color = variable)) + 
		geom_point(size = 2, position = position_jitter(width = 0.1)) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("NHARvitc_NHARcontrol NHARcontrol unique") + 
		guides(color = guide_legend(title = "")) + 
		theme_bw() + 
		theme(axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(homer_unique_enhancer_MGG119vitc_MGG119control_RPKM_figure, file = "./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_MGG119vitc_MGG119control_RPKM_figure.pdf", height = 5, width = 5)
pdf("./unique2/NHAR_vitc.NHAR_control/H3K27ac/Homer2/homer_unique_enhancer_MGG119vitc_MGG119control_figure.pdf", height = 5, width = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(homer_unique_enhancer_MGG119vitc_MGG119control_figure + theme(plot.margin = unit(c(1,0,1,1), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(homer_unique_enhancer_MGG119vitc_MGG119control_RPKM_figure + ggtitle("") + theme(plot.margin = unit(c(1,1,1,-0.5), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

### vitc reversed histone mods
#### enrichment in genomic features
vitc_reversed_genomic_breakdown <- read.delim("./unique2/intersect/genomic.breakdown.summary", as.is = T) %>% mutate(Mark = gsub("\\..*", "", Name), Category = gsub(".*mut", "mut", Name), X = NULL, NCpG = NULL, Name = NULL)
vitc_reversed_genomic_breakdown_tall <- melt(vitc_reversed_genomic_breakdown, id = c("Mark", "Category")) %>% 
	mutate(value = ifelse(Category == "mut_loss.vitc_gain", value, -value))
(vitc_reversed_genomic_breakdown_figure <- ggplot(vitc_reversed_genomic_breakdown_tall, aes(variable, value, fill = Category)) + 
		geom_bar(position = "identity", stat = "identity", width = 0.5) + 
		geom_hline(yintercept = c(-2, 2)) + 
		facet_wrap(~Mark, scales = "free_x") + 
		xlab("") + 
		ylab("Fold enrichment") + 
		scale_fill_manual(name = "", values = c("red", "blue")) + 
		coord_flip() + 
		theme_bw())
ggsave(vitc_reversed_genomic_breakdown_figure, file = "./unique2/vitc_reversed_genomic_breakdown_figure.pdf", height = 6, width = 7)
#### H3K27ac
K27ac_MGG119_vitc_unique <- as.numeric(system("less ./unique2/MGG_vitc.MGG119_control/H3K27ac/H3K27ac.MGG_vitc.MGG119_control.MGG_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_MGG119_control_unique <- as.numeric(system("less ./unique2/MGG_vitc.MGG119_control/H3K27ac/H3K27ac.MGG_vitc.MGG119_control.MGG119_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_NHAR_vitc_unique <- as.numeric(system("less ./unique2/NHAR_vitc.NHAR_control/H3K27ac/H3K27ac.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_NHAR_control_unique <- as.numeric(system("less ./unique2/NHAR_vitc.NHAR_control/H3K27ac/H3K27ac.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_NHA_mut_unique <- as.numeric(system("less ./unique2/NHAR_control.NHA_control/H3K27ac/H3K27ac.NHAR_control.NHA_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_NHA_wt_unique <- as.numeric(system("less ./unique2/NHAR_control.NHA_control/H3K27ac/H3K27ac.NHAR_control.NHA_control.NHA_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_MGG119_NHAR_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/MGG_vitc.MGG119_control/H3K27ac/H3K27ac.MGG_vitc.MGG119_control.MGG_vitc.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K27ac/H3K27ac.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_MGG119_NHAR_vitc_venn <- draw.pairwise.venn(area1 = K27ac_MGG119_vitc_unique, area2 = K27ac_NHAR_vitc_unique, cross.area = K27ac_MGG119_NHAR_vitc, category = c("MGG119", "NHAR"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique2/K27ac_MGG119_NHAR_vitc_venn.pdf", height = 3, width = 4)
grid.draw(K27ac_MGG119_NHAR_vitc_venn)
dev.off()
K27ac_MGG119_NHAR_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/MGG_vitc.MGG119_control/H3K27ac/H3K27ac.MGG_vitc.MGG119_control.MGG119_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K27ac/H3K27ac.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_MGG119_NHAR_control_venn <- draw.pairwise.venn(area1 = K27ac_MGG119_control_unique, area2 = K27ac_NHAR_control_unique, cross.area = K27ac_MGG119_NHAR_control, category = c("MGG119", "NHAR"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique2/K27ac_MGG119_NHAR_control_venn.pdf", height = 3, width = 4)
grid.draw(K27ac_MGG119_NHAR_control_venn)
dev.off()
K27ac_NHA_mut_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K27ac/H3K27ac.NHAR_control.NHA_control.NHAR_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K27ac/H3K27ac.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_NHA_mut_control_venn <- draw.pairwise.venn(area1 = K27ac_NHA_mut_unique, area2 = K27ac_NHAR_control_unique, cross.area = K27ac_NHA_mut_control, category = c("NHAR_NHA\n mut unique", "NHAR vitc_control\n control unique"), fill = c("red", "blue"), cex = 1, cat.pos = c(-90, 90), cat.dist = 0.05, cat.cex = 1, margin = 0.1)
pdf("./unique2/K27ac_NHA_mut_control_venn.pdf", height = 2, width = 7)
grid.draw(K27ac_NHA_mut_control_venn)
dev.off()
K27ac_NHA_wt_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K27ac/H3K27ac.NHAR_control.NHA_control.NHA_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K27ac/H3K27ac.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_NHA_wt_vitc_venn <- draw.pairwise.venn(area1 = K27ac_NHA_wt_unique, area2 = K27ac_NHAR_vitc_unique, cross.area = K27ac_NHA_wt_vitc, category = c("NHAR_NHA\n wt unique", "NHAR vitc_control\n vitc unique"), fill = c("red", "blue"), cex = 1, cat.pos = c(-90, 90), cat.dist = 0.05, cat.cex = 1, margin = 0.1)
pdf("./unique2/K27ac_NHA_wt_vitc_venn.pdf", height = 2, width = 7)
grid.draw(K27ac_NHA_wt_vitc_venn)
dev.off()
K27ac_NHA_mut_MGG119_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K27ac/H3K27ac.NHAR_control.NHA_control.NHAR_control.unique -b ./unique2/MGG_vitc.MGG119_control/H3K27ac/H3K27ac.MGG_vitc.MGG119_control.MGG119_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_NHA_mut_MGG119_control_venn <- draw.pairwise.venn(area1 = K27ac_NHA_mut_unique, area2 = K27ac_MGG119_control_unique, cross.area = K27ac_NHA_mut_MGG119_control, category = c("NHAR_NHA\n mut unique", "MGG119 vitc_control\n control unique"), fill = c("red", "blue"), cex = 1.3, cat.pos = c(-10, 170), cat.dist = 0.05, cat.cex = 1)
pdf("./unique2/K27ac_NHA_mut_MGG119_control_venn.pdf", height = 3, width = 4)
grid.draw(K27ac_NHA_mut_MGG119_control_venn)
dev.off()
K27ac_NHA_wt_MGG119_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K27ac/H3K27ac.NHAR_control.NHA_control.NHA_control.unique -b ./unique2/MGG_vitc.MGG119_control/H3K27ac/H3K27ac.MGG_vitc.MGG119_control.MGG_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27ac_NHA_wt_MGG119_vitc_venn <- draw.pairwise.venn(area1 = K27ac_NHA_wt_unique, area2 = K27ac_MGG119_vitc_unique, cross.area = K27ac_NHA_wt_MGG119_vitc, category = c("NHAR_NHA\n wt unique", "MGG119 vitc_control\n vitc unique"), fill = c("red", "blue"), cex = 1.3, cat.pos = c(-10, 170), cat.dist = 0.05, cat.cex = 1)
pdf("./unique2/K27ac_NHA_wt_MGG119_vitc_venn.pdf", height = 3, width = 4)
grid.draw(K27ac_NHA_wt_MGG119_vitc_venn)
dev.off()
(GREAT_K27ac_mutGain_vitcLoss_figure <- enrich_GREAT("H3K27ac.mut_gain.vitc_loss", "H3K27ac.mut_gain.vitc_loss", dirIn = "./unique2/enrich/", dirOut = "./unique2/enrich/", height = 8, width = 7, categories = c("GOBP", "MousePhenotype", "MSigPathway")))
(GREAT_K27ac_mutLoss_vitcGain_figure <- enrich_GREAT("H3K27ac.mut_loss.vitc_gain", "H3K27ac.mut_loss.vitc_gain", dirIn = "./unique2/enrich/", dirOut = "./unique2/enrich/", height = 8, width = 7))
TF_ENSG <- read.delim("./unique2/TF.ENSG", as.is = T) %>% merge(RPKM)
homer_K27ac_mutGain_vitcLoss <- read.delim("./unique2/homer/H3K27ac.mut_gain.vitc_loss/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) %>% merge(TF_ENSG)
(homer_K27ac_mutGain_vitcLoss_figure <- ggplot(homer_K27ac_mutGain_vitcLoss, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("Mutant gain & VitC loss H3K27ac") + 
		theme_bw())
ggsave(homer_K27ac_mutGain_vitcLoss_figure, file = "./unique2/homer/homer_K27ac_mutGain_vitcLoss_figure.pdf", height = 3, width = 5)
homer_K27ac_mutGain_vitcLoss_RPKM <- homer_K27ac_mutGain_vitcLoss %>% select(TF, MGG119_control:NHA_vitc) %>% melt(id = "TF") %>% mutate(TF = factor(TF, levels = levels(homer_K27ac_mutGain_vitcLoss$TF)))
(homer_K27ac_mutGain_vitcLoss_RPKM_figure <- ggplot(homer_K27ac_mutGain_vitcLoss_RPKM, aes(TF, log10(value), color = variable)) + 
		geom_point(size = 2, position = position_jitter(width = 0.1)) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("Mutant gain & VitC loss H3K27ac") + 
		guides(color = guide_legend(title = "")) + 
		theme_bw() + 
		theme(axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(homer_K27ac_mutGain_vitcLoss_RPKM_figure, file = "./unique2/homer/homer_K27ac_mutGain_vitcLoss_RPKM_figure.pdf", height = 3, width = 5)
pdf("./unique2/homer/homer_K27ac_mutGain_vitcLoss_figure.pdf", height = 3, width = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(homer_K27ac_mutGain_vitcLoss_figure + theme(plot.margin = unit(c(1,0,1,1), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(homer_K27ac_mutGain_vitcLoss_RPKM_figure + ggtitle("") + theme(plot.margin = unit(c(1,1,1,-0.5), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
homer_K27ac_mutLoss_vitcGain <- read.delim("./unique2/homer/H3K27ac.mut_loss.vitc_gain/knownResults.txt", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_enhancer_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>%
	filter(Percent_enhancer_with_motif >= 15, P.value<= 0.05) %>% arrange(Percent_enhancer_with_motif) %>% mutate(TF = factor(TF, levels = TF)) %>% merge(TF_ENSG)
(homer_K27ac_mutLoss_vitcGain_figure <- ggplot(homer_K27ac_mutLoss_vitcGain, aes(TF, Percent_enhancer_with_motif)) + 
		geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
		coord_flip() + 
		xlab("") + 
		ylab("Percent of enhancers with motif") + 
		ggtitle("Mutant loss & VitC gain H3K27ac") + 
		theme_bw())
ggsave(homer_K27ac_mutLoss_vitcGain_figure, file = "./unique2/homer/homer_K27ac_mutLoss_vitcGain_figure.pdf", height = 2, width = 5)
homer_K27ac_mutLoss_vitcGain_RPKM <- homer_K27ac_mutLoss_vitcGain %>% select(TF, MGG119_control:NHA_vitc) %>% melt(id = "TF") %>% mutate(TF = factor(TF, levels = levels(homer_K27ac_mutLoss_vitcGain$TF)))
(homer_K27ac_mutLoss_vitcGain_RPKM_figure <- ggplot(homer_K27ac_mutLoss_vitcGain_RPKM, aes(TF, log10(value), color = variable)) + 
		geom_point(size = 2, position = position_jitter(width = 0.1)) + 
		coord_flip() + 
		xlab("") + 
		ylab("log10 RPKM") + 
		ggtitle("Mutant loss & VitC gain H3K27ac") + 
		guides(color = guide_legend(title = "")) + 
		theme_bw() + 
		theme(axis.text.y = element_text(size = 0), plot.background = element_rect(fill = "transparent")))
ggsave(homer_K27ac_mutLoss_vitcGain_RPKM_figure, file = "./unique2/homer/homer_K27ac_mutLoss_vitcGain_RPKM_figure.pdf", height = 2, width = 5)
pdf("./unique2/homer/homer_K27ac_mutLoss_vitcGain_figure.pdf", height = 2, width = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(homer_K27ac_mutLoss_vitcGain_figure + theme(plot.margin = unit(c(1,0,1,1), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(homer_K27ac_mutLoss_vitcGain_RPKM_figure + ggtitle("") + theme(plot.margin = unit(c(1,1,1,-0.5), "cm")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
(enrich_K27ac_mutGain_vitcLoss_GATA3 <- enrich(name = "H3K27ac.mut_gain.vitc_loss.GATA3", dirIn = "./unique2/enrich/", dirOut = "./unique2/enrich/", fdr = 0.1, p = "FDR", erminej = F, height = 2, width = 8))
(enrich_K27ac_mutGain_vitcLoss_Foxo1 <- enrich(name = "H3K27ac.mut_gain.vitc_loss.Foxo1", dirIn = "./unique2/enrich/", dirOut = "./unique2/enrich/", fdr = 0.1, p = "FDR", erminej = F, height = 2, width = 8))
(enrich_K27ac_mutGain_vitcLoss_BMYB <- enrich(name = "H3K27ac.mut_gain.vitc_loss.BMYB", dirIn = "./unique2/enrich/", dirOut = "./unique2/enrich/", fdr = 0.1, p = "FDR", erminej = F, height = 2, width = 8))
(enrich_K27ac_mutLoss_vitcGain_Bcl6 <- enrich(name = "H3K27ac.mut_loss.vitc_gain.Bcl6", dirIn = "./unique2/enrich/", dirOut = "./unique2/enrich/", fdr = 0.1, p = "FDR", erminej = F, height = 2, width = 8))
e <- 1e-5
K27ac_mutGain_vitcLoss_GATA3_RPKM <- read.delim("./unique2/H3K27ac.mut_gain.vitc_loss.GATA3.closest.gene.RPKM", as.is = T) %>% select(ENSG, Name:NHA_vitc) 
heatmap.2(as.matrix(K27ac_mutGain_vitcLoss_GATA3_RPKM %>% select(-ENSG, -Name)), scale = "row", dendrogram = "column", trace = "none", margins = c(10,3), labRow = NA, keysize = 1, density.info = "none", key.title = "", key.xlab = "RPKM Z-score", main = "H3K27ac.mut_gain.vitc_loss.GATA3")
K27ac_mutGain_vitcLoss_Foxo1_RPKM <- read.delim("./unique2/H3K27ac.mut_gain.vitc_loss.Foxo1.closest.gene.RPKM", as.is = T) %>% select(ENSG, Name:NHA_vitc) 
heatmap.2(as.matrix(K27ac_mutGain_vitcLoss_Foxo1_RPKM %>% select(-ENSG, -Name)), scale = "row", dendrogram = "column", trace = "none", margins = c(10,3), labRow = NA, keysize = 1, density.info = "none", key.title = "", key.xlab = "RPKM Z-score", main = "H3K27ac.mut_gain.vitc_loss.Foxo1")
UP_NHAR_control_NHA_control <- read.delim("../RNAseq/DEfine/UP.NHAR_control_NHA_control.FDR_0.01.rmin_0.005.Nmin_25", as.is = T, col.names = c("ENSG", "RPKM1", "RPKM2", "FDR")) %>% select(ENSG) %>% mutate(DE = "UP", Sample = "NHAR_control_NHA_control")
DN_NHAR_control_NHA_control <- read.delim("../RNAseq/DEfine/DN.NHAR_control_NHA_control.FDR_0.01.rmin_0.005.Nmin_25", as.is = T, col.names = c("ENSG", "RPKM1", "RPKM2", "FDR")) %>% select(ENSG) %>% mutate(DE = "DN", Sample = "NHAR_control_NHA_control")
UP_NHAR_vitc_NHAR_control <- read.delim("../RNAseq/DEfine/UP.NHAR_vitc_NHAR_control.FDR_0.01.rmin_0.005.Nmin_25", as.is = T, col.names = c("ENSG", "RPKM1", "RPKM2", "FDR")) %>% select(ENSG) %>% mutate(DE = "UP", Sample = "NHAR_vitc_NHAR_control")
DN_NHAR_vitc_NHAR_control <- read.delim("../RNAseq/DEfine/DN.NHAR_vitc_NHAR_control.FDR_0.01.rmin_0.005.Nmin_25", as.is = T, col.names = c("ENSG", "RPKM1", "RPKM2", "FDR")) %>% select(ENSG) %>% mutate(DE = "DN", Sample = "NHAR_vitc_NHAR_control")
DE <- rbind(UP_NHAR_control_NHA_control, DN_NHAR_control_NHA_control, UP_NHAR_vitc_NHAR_control, DN_NHAR_vitc_NHAR_control)
K27ac_mutGain_vitcLoss_GATA3_DE <- merge(K27ac_mutGain_vitcLoss_GATA3_RPKM, DE) %>% unique()
K27ac_mutGain_vitcLoss_Foxo1_DE <- merge(K27ac_mutGain_vitcLoss_Foxo1_RPKM, DE) %>% unique()
#### H3K4me1
K4me1_MGG119_vitc_unique <- as.numeric(system("less ./unique2/MGG_vitc.MGG119_control/H3K4me1/H3K4me1.MGG_vitc.MGG119_control.MGG_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_MGG119_control_unique <- as.numeric(system("less ./unique2/MGG_vitc.MGG119_control/H3K4me1/H3K4me1.MGG_vitc.MGG119_control.MGG119_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_NHAR_vitc_unique <- as.numeric(system("less ./unique2/NHAR_vitc.NHAR_control/H3K4me1/H3K4me1.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_NHAR_control_unique <- as.numeric(system("less ./unique2/NHAR_vitc.NHAR_control/H3K4me1/H3K4me1.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_NHA_mut_unique <- as.numeric(system("less ./unique2/NHAR_control.NHA_control/H3K4me1/H3K4me1.NHAR_control.NHA_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_NHA_wt_unique <- as.numeric(system("less ./unique2/NHAR_control.NHA_control/H3K4me1/H3K4me1.NHAR_control.NHA_control.NHA_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_MGG119_NHAR_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/MGG_vitc.MGG119_control/H3K4me1/H3K4me1.MGG_vitc.MGG119_control.MGG_vitc.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K4me1/H3K4me1.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_MGG119_NHAR_vitc_venn <- draw.pairwise.venn(area1 = K4me1_MGG119_vitc_unique, area2 = K4me1_NHAR_vitc_unique, cross.area = K4me1_MGG119_NHAR_vitc, category = c("MGG119", "NHAR"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique2/K4me1_MGG119_NHAR_vitc_venn.pdf", height = 3, width = 4)
grid.draw(K4me1_MGG119_NHAR_vitc_venn)
dev.off()
K4me1_MGG119_NHAR_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/MGG_vitc.MGG119_control/H3K4me1/H3K4me1.MGG_vitc.MGG119_control.MGG119_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K4me1/H3K4me1.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_MGG119_NHAR_control_venn <- draw.pairwise.venn(area1 = K4me1_MGG119_control_unique, area2 = K4me1_NHAR_control_unique, cross.area = K4me1_MGG119_NHAR_control, category = c("MGG119", "NHAR"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique2/K4me1_MGG119_NHAR_control_venn.pdf", height = 3, width = 4)
grid.draw(K4me1_MGG119_NHAR_control_venn)
dev.off()
K4me1_NHA_mut_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K4me1/H3K4me1.NHAR_control.NHA_control.NHAR_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K4me1/H3K4me1.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_NHA_mut_control_venn <- draw.pairwise.venn(area1 = K4me1_NHA_mut_unique, area2 = K4me1_NHAR_control_unique, cross.area = K4me1_NHA_mut_control, category = c("NHAR_NHA\n mut unique", "NHAR vitc_control\n control unique"), fill = c("red", "blue"), cex = 1, cat.pos = c(-90, 90), cat.dist = 0.05, cat.cex = 1, margin = 0.1)
pdf("./unique2/K4me1_NHA_mut_control_venn.pdf", height = 2, width = 7)
grid.draw(K4me1_NHA_mut_control_venn)
dev.off()
K4me1_NHA_wt_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K4me1/H3K4me1.NHAR_control.NHA_control.NHA_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K4me1/H3K4me1.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_NHA_wt_vitc_venn <- draw.pairwise.venn(area1 = K4me1_NHA_wt_unique, area2 = K4me1_NHAR_vitc_unique, cross.area = K4me1_NHA_wt_vitc, category = c("NHAR_NHA\n wt unique", "NHAR vitc_control\n vitc unique"), fill = c("red", "blue"), cex = 1, cat.pos = c(-90, 90), cat.dist = 0.05, cat.cex = 1, margin = 0.1)
pdf("./unique2/K4me1_NHA_wt_vitc_venn.pdf", height = 2, width = 7)
grid.draw(K4me1_NHA_wt_vitc_venn)
dev.off()
K4me1_NHA_mut_MGG119_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K4me1/H3K4me1.NHAR_control.NHA_control.NHAR_control.unique -b ./unique2/MGG_vitc.MGG119_control/H3K4me1/H3K4me1.MGG_vitc.MGG119_control.MGG119_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_NHA_mut_MGG119_control_venn <- draw.pairwise.venn(area1 = K4me1_NHA_mut_unique, area2 = K4me1_MGG119_control_unique, cross.area = K4me1_NHA_mut_MGG119_control, category = c("NHAR_NHA\n mut unique", "MGG119 vitc_control\n control unique"), fill = c("red", "blue"), cex = 1.3, cat.pos = c(-10, 170), cat.dist = 0.05, cat.cex = 1)
pdf("./unique2/K4me1_NHA_mut_MGG119_control_venn.pdf", height = 3, width = 4)
grid.draw(K4me1_NHA_mut_MGG119_control_venn)
dev.off()
K4me1_NHA_wt_MGG119_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K4me1/H3K4me1.NHAR_control.NHA_control.NHA_control.unique -b ./unique2/MGG_vitc.MGG119_control/H3K4me1/H3K4me1.MGG_vitc.MGG119_control.MGG_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me1_NHA_wt_MGG119_vitc_venn <- draw.pairwise.venn(area1 = K4me1_NHA_wt_unique, area2 = K4me1_MGG119_vitc_unique, cross.area = K4me1_NHA_wt_MGG119_vitc, category = c("NHAR_NHA\n wt unique", "MGG119 vitc_control\n vitc unique"), fill = c("red", "blue"), cex = 1.3, cat.pos = c(-10, 170), cat.dist = 0.05, cat.cex = 1)
pdf("./unique2/K4me1_NHA_wt_MGG119_vitc_venn.pdf", height = 3, width = 4)
grid.draw(K4me1_NHA_wt_MGG119_vitc_venn)
dev.off()
(GREAT_K4me1_mutGain_vitcLoss_figure <- enrich_GREAT("H3K4me1.mut_gain.vitc_loss", "H3K4me1.mut_gain.vitc_loss", dirIn = "./unique2/enrich/", dirOut = "./unique2/enrich/", height = 8, width = 7, categories = c("GOBP", "GOMF", "MSigPathway")))
(GREAT_K4me1_mutLoss_vitcGain_figure <- enrich_GREAT("H3K4me1.mut_loss.vitc_gain", "H3K4me1.mut_loss.vitc_gain", dirIn = "./unique2/enrich/", dirOut = "./unique2/enrich/", height = 8, width = 7, categories = c("GOBP", "HumanPhenotype", "MSigPathway")))
#### H3K4me3
K4me3_MGG119_vitc_unique <- as.numeric(system("less ./unique2/MGG_vitc.MGG119_control/H3K4me3/H3K4me3.MGG_vitc.MGG119_control.MGG_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_MGG119_control_unique <- as.numeric(system("less ./unique2/MGG_vitc.MGG119_control/H3K4me3/H3K4me3.MGG_vitc.MGG119_control.MGG119_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_NHAR_vitc_unique <- as.numeric(system("less ./unique2/NHAR_vitc.NHAR_control/H3K4me3/H3K4me3.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_NHAR_control_unique <- as.numeric(system("less ./unique2/NHAR_vitc.NHAR_control/H3K4me3/H3K4me3.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_NHA_mut_unique <- as.numeric(system("less ./unique2/NHAR_control.NHA_control/H3K4me3/H3K4me3.NHAR_control.NHA_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_NHA_wt_unique <- as.numeric(system("less ./unique2/NHAR_control.NHA_control/H3K4me3/H3K4me3.NHAR_control.NHA_control.NHA_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_MGG119_NHAR_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/MGG_vitc.MGG119_control/H3K4me3/H3K4me3.MGG_vitc.MGG119_control.MGG_vitc.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K4me3/H3K4me3.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_MGG119_NHAR_vitc_venn <- draw.pairwise.venn(area1 = K4me3_MGG119_vitc_unique, area2 = K4me3_NHAR_vitc_unique, cross.area = K4me3_MGG119_NHAR_vitc, category = c("MGG119", "NHAR"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique2/K4me3_MGG119_NHAR_vitc_venn.pdf", height = 3, width = 4)
grid.draw(K4me3_MGG119_NHAR_vitc_venn)
dev.off()
K4me3_MGG119_NHAR_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/MGG_vitc.MGG119_control/H3K4me3/H3K4me3.MGG_vitc.MGG119_control.MGG119_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K4me3/H3K4me3.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_MGG119_NHAR_control_venn <- draw.pairwise.venn(area1 = K4me3_MGG119_control_unique, area2 = K4me3_NHAR_control_unique, cross.area = K4me3_MGG119_NHAR_control, category = c("MGG119", "NHAR"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique2/K4me3_MGG119_NHAR_control_venn.pdf", height = 3, width = 4)
grid.draw(K4me3_MGG119_NHAR_control_venn)
dev.off()
K4me3_NHA_mut_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K4me3/H3K4me3.NHAR_control.NHA_control.NHAR_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K4me3/H3K4me3.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_NHA_mut_control_venn <- draw.pairwise.venn(area1 = K4me3_NHA_mut_unique, area2 = K4me3_NHAR_control_unique, cross.area = K4me3_NHA_mut_control, category = c("NHAR_NHA\n mut unique", "NHAR vitc_control\n control unique"), fill = c("red", "blue"), cex = 1, cat.pos = c(-90, 90), cat.dist = 0.05, cat.cex = 1, margin = 0.1)
pdf("./unique2/K4me3_NHA_mut_control_venn.pdf", height = 2, width = 7)
grid.draw(K4me3_NHA_mut_control_venn)
dev.off()
K4me3_NHA_wt_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K4me3/H3K4me3.NHAR_control.NHA_control.NHA_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K4me3/H3K4me3.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_NHA_wt_vitc_venn <- draw.pairwise.venn(area1 = K4me3_NHA_wt_unique, area2 = K4me3_NHAR_vitc_unique, cross.area = K4me3_NHA_wt_vitc, category = c("NHAR_NHA\n wt unique", "NHAR vitc_control\n vitc unique"), fill = c("red", "blue"), cex = 1, cat.pos = c(-90, 90), cat.dist = 0.05, cat.cex = 1, margin = 0.1)
pdf("./unique2/K4me3_NHA_wt_vitc_venn.pdf", height = 2, width = 7)
grid.draw(K4me3_NHA_wt_vitc_venn)
dev.off()
K4me3_NHA_mut_MGG119_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K4me3/H3K4me3.NHAR_control.NHA_control.NHAR_control.unique -b ./unique2/MGG_vitc.MGG119_control/H3K4me3/H3K4me3.MGG_vitc.MGG119_control.MGG119_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_NHA_mut_MGG119_control_venn <- draw.pairwise.venn(area1 = K4me3_NHA_mut_unique, area2 = K4me3_MGG119_control_unique, cross.area = K4me3_NHA_mut_MGG119_control, category = c("NHAR_NHA\n mut unique", "MGG119 vitc_control\n control unique"), fill = c("red", "blue"), cex = 1.3, cat.pos = c(-10, 170), cat.dist = 0.05, cat.cex = 1)
pdf("./unique2/K4me3_NHA_mut_MGG119_control_venn.pdf", height = 3, width = 4)
grid.draw(K4me3_NHA_mut_MGG119_control_venn)
dev.off()
K4me3_NHA_wt_MGG119_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K4me3/H3K4me3.NHAR_control.NHA_control.NHA_control.unique -b ./unique2/MGG_vitc.MGG119_control/H3K4me3/H3K4me3.MGG_vitc.MGG119_control.MGG_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K4me3_NHA_wt_MGG119_vitc_venn <- draw.pairwise.venn(area1 = K4me3_NHA_wt_unique, area2 = K4me3_MGG119_vitc_unique, cross.area = K4me3_NHA_wt_MGG119_vitc, category = c("NHAR_NHA\n wt unique", "MGG119 vitc_control\n vitc unique"), fill = c("red", "blue"), cex = 1.3, cat.pos = c(-10, 170), cat.dist = 0.05, cat.cex = 1)
pdf("./unique2/K4me3_NHA_wt_MGG119_vitc_venn.pdf", height = 3, width = 4)
grid.draw(K4me3_NHA_wt_MGG119_vitc_venn)
dev.off()
(GREAT_K4me3_mutGain_vitcLoss_figure <- enrich_GREAT("H3K4me3.mut_gain.vitc_loss", "H3K4me3.mut_gain.vitc_loss", dirIn = "./unique2/enrich/", dirOut = "./unique2/enrich/", height = 8, width = 7, categories = c("GOBP", "GOMF", "MSigPathway")))
(GREAT_K4me3_mutLoss_vitcGain_figure <- enrich_GREAT("H3K4me3.mut_loss.vitc_gain", "H3K4me3.mut_loss.vitc_gain", dirIn = "./unique2/enrich/", dirOut = "./unique2/enrich/", height = 8, width = 7))
#### H3K27me3
K27me3_MGG119_vitc_unique <- as.numeric(system("less ./unique2/MGG_vitc.MGG119_control/H3K27me3/H3K27me3.MGG_vitc.MGG119_control.MGG_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_MGG119_control_unique <- as.numeric(system("less ./unique2/MGG_vitc.MGG119_control/H3K27me3/H3K27me3.MGG_vitc.MGG119_control.MGG119_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_NHAR_vitc_unique <- as.numeric(system("less ./unique2/NHAR_vitc.NHAR_control/H3K27me3/H3K27me3.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_NHAR_control_unique <- as.numeric(system("less ./unique2/NHAR_vitc.NHAR_control/H3K27me3/H3K27me3.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_NHA_mut_unique <- as.numeric(system("less ./unique2/NHAR_control.NHA_control/H3K27me3/H3K27me3.NHAR_control.NHA_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_NHA_wt_unique <- as.numeric(system("less ./unique2/NHAR_control.NHA_control/H3K27me3/H3K27me3.NHAR_control.NHA_control.NHA_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_MGG119_NHAR_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/MGG_vitc.MGG119_control/H3K27me3/H3K27me3.MGG_vitc.MGG119_control.MGG_vitc.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K27me3/H3K27me3.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_MGG119_NHAR_vitc_venn <- draw.pairwise.venn(area1 = K27me3_MGG119_vitc_unique, area2 = K27me3_NHAR_vitc_unique, cross.area = K27me3_MGG119_NHAR_vitc, category = c("MGG119", "NHAR"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique2/K27me3_MGG119_NHAR_vitc_venn.pdf", height = 3, width = 4)
grid.draw(K27me3_MGG119_NHAR_vitc_venn)
dev.off()
K27me3_MGG119_NHAR_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/MGG_vitc.MGG119_control/H3K27me3/H3K27me3.MGG_vitc.MGG119_control.MGG119_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K27me3/H3K27me3.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_MGG119_NHAR_control_venn <- draw.pairwise.venn(area1 = K27me3_MGG119_control_unique, area2 = K27me3_NHAR_control_unique, cross.area = K27me3_MGG119_NHAR_control, category = c("MGG119", "NHAR"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique2/K27me3_MGG119_NHAR_control_venn.pdf", height = 3, width = 4)
grid.draw(K27me3_MGG119_NHAR_control_venn)
dev.off()
K27me3_NHA_mut_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K27me3/H3K27me3.NHAR_control.NHA_control.NHAR_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K27me3/H3K27me3.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_NHA_mut_control_venn <- draw.pairwise.venn(area1 = K27me3_NHA_mut_unique, area2 = K27me3_NHAR_control_unique, cross.area = K27me3_NHA_mut_control, category = c("NHAR_NHA\n mut unique", "NHAR vitc_control\n control unique"), fill = c("red", "blue"), cex = 1, cat.pos = c(-90, 90), cat.dist = 0.05, cat.cex = 1, margin = 0.1)
pdf("./unique2/K27me3_NHA_mut_control_venn.pdf", height = 2, width = 7)
grid.draw(K27me3_NHA_mut_control_venn)
dev.off()
K27me3_NHA_wt_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K27me3/H3K27me3.NHAR_control.NHA_control.NHA_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K27me3/H3K27me3.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_NHA_wt_vitc_venn <- draw.pairwise.venn(area1 = K27me3_NHA_wt_unique, area2 = K27me3_NHAR_vitc_unique, cross.area = K27me3_NHA_wt_vitc, category = c("NHAR_NHA\n wt unique", "NHAR vitc_control\n vitc unique"), fill = c("red", "blue"), cex = 1, cat.pos = c(-90, 90), cat.dist = 0.05, cat.cex = 1, margin = 0.1)
pdf("./unique2/K27me3_NHA_wt_vitc_venn.pdf", height = 2, width = 7)
grid.draw(K27me3_NHA_wt_vitc_venn)
dev.off()
K27me3_NHA_mut_MGG119_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K27me3/H3K27me3.NHAR_control.NHA_control.NHAR_control.unique -b ./unique2/MGG_vitc.MGG119_control/H3K27me3/H3K27me3.MGG_vitc.MGG119_control.MGG119_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_NHA_mut_MGG119_control_venn <- draw.pairwise.venn(area1 = K27me3_NHA_mut_unique, area2 = K27me3_MGG119_control_unique, cross.area = K27me3_NHA_mut_MGG119_control, category = c("NHAR_NHA\n mut unique", "MGG119 vitc_control\n control unique"), fill = c("red", "blue"), cex = 1.3, cat.pos = c(-10, 170), cat.dist = 0.05, cat.cex = 1)
pdf("./unique2/K27me3_NHA_mut_MGG119_control_venn.pdf", height = 3, width = 4)
grid.draw(K27me3_NHA_mut_MGG119_control_venn)
dev.off()
K27me3_NHA_wt_MGG119_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K27me3/H3K27me3.NHAR_control.NHA_control.NHA_control.unique -b ./unique2/MGG_vitc.MGG119_control/H3K27me3/H3K27me3.MGG_vitc.MGG119_control.MGG_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K27me3_NHA_wt_MGG119_vitc_venn <- draw.pairwise.venn(area1 = K27me3_NHA_wt_unique, area2 = K27me3_MGG119_vitc_unique, cross.area = K27me3_NHA_wt_MGG119_vitc, category = c("NHAR_NHA\n wt unique", "MGG119 vitc_control\n vitc unique"), fill = c("red", "blue"), cex = 1.3, cat.pos = c(-10, 170), cat.dist = 0.05, cat.cex = 1)
pdf("./unique2/K27me3_NHA_wt_MGG119_vitc_venn.pdf", height = 3, width = 4)
grid.draw(K27me3_NHA_wt_MGG119_vitc_venn)
dev.off()
#### H3K36me3
K36me3_MGG119_vitc_unique <- as.numeric(system("less ./unique2/MGG_vitc.MGG119_control/H3K36me3/H3K36me3.MGG_vitc.MGG119_control.MGG_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_MGG119_control_unique <- as.numeric(system("less ./unique2/MGG_vitc.MGG119_control/H3K36me3/H3K36me3.MGG_vitc.MGG119_control.MGG119_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_NHAR_vitc_unique <- as.numeric(system("less ./unique2/NHAR_vitc.NHAR_control/H3K36me3/H3K36me3.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_NHAR_control_unique <- as.numeric(system("less ./unique2/NHAR_vitc.NHAR_control/H3K36me3/H3K36me3.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_NHA_mut_unique <- as.numeric(system("less ./unique2/NHAR_control.NHA_control/H3K36me3/H3K36me3.NHAR_control.NHA_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_NHA_wt_unique <- as.numeric(system("less ./unique2/NHAR_control.NHA_control/H3K36me3/H3K36me3.NHAR_control.NHA_control.NHA_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_MGG119_NHAR_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/MGG_vitc.MGG119_control/H3K36me3/H3K36me3.MGG_vitc.MGG119_control.MGG_vitc.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K36me3/H3K36me3.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_MGG119_NHAR_vitc_venn <- draw.pairwise.venn(area1 = K36me3_MGG119_vitc_unique, area2 = K36me3_NHAR_vitc_unique, cross.area = K36me3_MGG119_NHAR_vitc, category = c("MGG119", "NHAR"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique2/K36me3_MGG119_NHAR_vitc_venn.pdf", height = 3, width = 4)
grid.draw(K36me3_MGG119_NHAR_vitc_venn)
dev.off()
K36me3_MGG119_NHAR_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/MGG_vitc.MGG119_control/H3K36me3/H3K36me3.MGG_vitc.MGG119_control.MGG119_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K36me3/H3K36me3.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_MGG119_NHAR_control_venn <- draw.pairwise.venn(area1 = K36me3_MGG119_control_unique, area2 = K36me3_NHAR_control_unique, cross.area = K36me3_MGG119_NHAR_control, category = c("MGG119", "NHAR"), fill = c("red", "blue"), cex = 1.5, cat.pos = c(-20, 160), cat.dist = 0.05, cat.cex = 1.5)
pdf("./unique2/K36me3_MGG119_NHAR_control_venn.pdf", height = 3, width = 4)
grid.draw(K36me3_MGG119_NHAR_control_venn)
dev.off()
K36me3_NHA_mut_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K36me3/H3K36me3.NHAR_control.NHA_control.NHAR_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K36me3/H3K36me3.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_NHA_mut_control_venn <- draw.pairwise.venn(area1 = K36me3_NHA_mut_unique, area2 = K36me3_NHAR_control_unique, cross.area = K36me3_NHA_mut_control, category = c("NHAR_NHA\n mut unique", "NHAR vitc_control\n control unique"), fill = c("red", "blue"), cex = 1, cat.pos = c(-90, 90), cat.dist = 0.05, cat.cex = 1, margin = 0.1)
pdf("./unique2/K36me3_NHA_mut_control_venn.pdf", height = 2, width = 7)
grid.draw(K36me3_NHA_mut_control_venn)
dev.off()
K36me3_NHA_wt_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K36me3/H3K36me3.NHAR_control.NHA_control.NHA_control.unique -b ./unique2/NHAR_vitc.NHAR_control/H3K36me3/H3K36me3.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_NHA_wt_vitc_venn <- draw.pairwise.venn(area1 = K36me3_NHA_wt_unique, area2 = K36me3_NHAR_vitc_unique, cross.area = K36me3_NHA_wt_vitc, category = c("NHAR_NHA\n wt unique", "NHAR vitc_control\n vitc unique"), fill = c("red", "blue"), cex = 1, cat.pos = c(-90, 90), cat.dist = 0.05, cat.cex = 1, margin = 0.1)
pdf("./unique2/K36me3_NHA_wt_vitc_venn.pdf", height = 2, width = 7)
grid.draw(K36me3_NHA_wt_vitc_venn)
dev.off()
K36me3_NHA_mut_MGG119_control <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K36me3/H3K36me3.NHAR_control.NHA_control.NHAR_control.unique -b ./unique2/MGG_vitc.MGG119_control/H3K36me3/H3K36me3.MGG_vitc.MGG119_control.MGG119_control.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_NHA_mut_MGG119_control_venn <- draw.pairwise.venn(area1 = K36me3_NHA_mut_unique, area2 = K36me3_MGG119_control_unique, cross.area = K36me3_NHA_mut_MGG119_control, category = c("NHAR_NHA\n mut unique", "MGG119 vitc_control\n control unique"), fill = c("red", "blue"), cex = 1.3, cat.pos = c(-10, 170), cat.dist = 0.05, cat.cex = 1)
pdf("./unique2/K36me3_NHA_mut_MGG119_control_venn.pdf", height = 3, width = 4)
grid.draw(K36me3_NHA_mut_MGG119_control_venn)
dev.off()
K36me3_NHA_wt_MGG119_vitc <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/intersectBed -a ./unique2/NHAR_control.NHA_control/H3K36me3/H3K36me3.NHAR_control.NHA_control.NHA_control.unique -b ./unique2/MGG_vitc.MGG119_control/H3K36me3/H3K36me3.MGG_vitc.MGG119_control.MGG_vitc.unique | awk '{s=s+$3-$2}END{print s}'", intern = T))
K36me3_NHA_wt_MGG119_vitc_venn <- draw.pairwise.venn(area1 = K36me3_NHA_wt_unique, area2 = K36me3_MGG119_vitc_unique, cross.area = K36me3_NHA_wt_MGG119_vitc, category = c("NHAR_NHA\n wt unique", "MGG119 vitc_control\n vitc unique"), fill = c("red", "blue"), cex = 1.3, cat.pos = c(-10, 170), cat.dist = 0.05, cat.cex = 1)
pdf("./unique2/K36me3_NHA_wt_MGG119_vitc_venn.pdf", height = 3, width = 4)
grid.draw(K36me3_NHA_wt_MGG119_vitc_venn)
dev.off()


## ============= save ===========
save(list = c(ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "enrich"), ls(pattern = "DAVID"), ls(pattern = "GREAT"), ls(pattern = "venn"), 
							"K27ac_mutGain_vitcLoss_GATA3_DE", "K27ac_mutGain_vitcLoss_Foxo1_DE"), 
		 file = "/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/ChIPseq.Rdata")

