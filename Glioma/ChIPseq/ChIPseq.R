# Glioma - ChIP-seq analysis

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
library(scales)
source('~/HirstLab/Pipeline/R/enrich.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
load("/projects/epigenomics2/users/lli/glioma/ChIPseq/ChIPseq.Rdata")
setwd("/projects/epigenomics2/users/lli/glioma/ChIPseq/")
libs <- c("CEMT_19", "CEMT_21", "CEMT_22", "CEMT_23", "CEMT_47", "GE04")

## -------- ER summary -----------
### adjust for differences in sequencing depth: subsampling deep sequenced bam file to the average sequence depth of all other samples
glioma_qc <- read.delim("~/HirstLab/Glioma/GliomaLibrariesDetails.tsv", as.is = T) %>% 
	mutate(Assay = gsub("ChIP-Seq\\/", "", Assay))
NPC_qc <- read.delim("~/HirstLab/FetalBrain/FetalBrainLibrariesDetail.tsv", as.is = T) %>% 
	mutate(Sample = gsub("HuFNSC", "", paste0(Cell.Type, Donor)), Sample = gsub("Primary Cell Culture Neurospheres, ", "", Sample), Sample = gsub("Ganglionic Eminence Derived", "GE", Sample))
sequencing_depth <- data.frame(ID = c(paste0(NPC_qc$Sample, "_", NPC_qc$Library.Strategy), paste0(glioma_qc$Donor, "_", glioma_qc$Assay)), depth = c(NPC_qc$Total_Number_Of_Reads_After_Filtering, glioma_qc$Number_Uniquely_Aligned_Reads_without_Dups_and_Q_.._10))
rownames(sequencing_depth) <- sequencing_depth$ID
ER_adjust <- read.delim("./FindER/subsample/ER.summary.adjust", as.is = T)
ER_adjust_summary <- ER_adjust %>% group_by(mark, sample) %>% 
	summarise(frac = frac[1], n_original = n_original[1], len_original = len_original[1], 
						n_adjust_ave = mean(n_adjust), n_adjust_sd = sd(n_adjust), len_adjust_ave = mean(len_adjust), len_adjust_sd = sd(len_adjust), 
						n_intersect_ave = mean(n_intersect), n_intersect_sd = sd(n_intersect), len_intersect_ave = mean(len_intersect), len_intersect_sd = sd(len_intersect))
rownames(ER_adjust_summary) <- paste0(ER_adjust_summary$sample, "_", ER_adjust_summary$mark)
ER_summary <- read.delim("./FindER/ER.summary.stats", as.is = T) %>% 
	mutate(Sample = gsub("NPC_", "", Sample), Seq_depth = sequencing_depth[paste0(Sample, "_", Mark), "depth"], N_region_adjust = N_region, Total_length_adjust = Total_length)
rownames(ER_summary) <- paste0(ER_summary$Sample, "_", ER_summary$Mark)
ER_summary[rownames(ER_adjust_summary), "N_region_adjust"] <- ER_adjust_summary$n_adjust_ave
ER_summary[rownames(ER_adjust_summary), "Total_length_adjust"] <- ER_adjust_summary$len_adjust_ave
ER_summary_tall <- ER_summary %>% melt(id = c("Mark", "Sample")) %>% mutate(sd = NA, variable = factor(variable, levels = c("N_region", "N_region_adjust", "Total_length", "Total_length_adjust", "Seq_depth")))
rownames(ER_summary_tall) <- paste0(ER_summary_tall$Sample, "_", ER_summary_tall$Mark, "_", ER_summary_tall$variable)
ER_summary_tall[paste0(rownames(ER_adjust_summary), "_N_region_adjust"), "sd"] <- ER_adjust_summary$n_adjust_sd
ER_summary_tall[paste0(rownames(ER_adjust_summary), "_Total_length_adjust"), "sd"] <- ER_adjust_summary$len_adjust_sd
(ER_summary_figure <- ggplot(ER_summary_tall, aes(Mark, value, fill = Sample)) + 
	geom_bar(position = position_dodge(width=0.9), stat = "identity") + 
	geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width=0.9), width=0.5) + 
	facet_grid(variable ~ ., scales = "free_y") +
	scale_y_continuous(label=scientific_format()) + 
	xlab("") + 
	ylab("") + 
	theme_bw())
ggsave(ER_summary_figure, file = "./FindER/ER_summary_figure.pdf", width = 8, height = 8)

for(i in 1:nrow(ER_adjust_summary)){
	assign(paste0("Venn_adjust_N_region", rownames(ER_adjust_summary)[i]), draw.pairwise.venn(as.integer(ER_adjust_summary[i, "n_adjust_ave"]), as.integer(ER_adjust_summary[i, "n_original"]), as.integer(ER_adjust_summary[i, "n_intersect_ave"]), category = c("adjusted", "original")))
	assign(paste0("Venn_adjust_Total_length", rownames(ER_adjust_summary)[i]), draw.pairwise.venn(as.integer(ER_adjust_summary[i, "len_adjust_ave"]), as.integer(ER_adjust_summary[i, "len_original"]), as.integer(ER_adjust_summary[i, "len_intersect_ave"]), category = c("adjusted", "original")))
}
pdf("./FindER/subsample/subsample_summary.pdf", height = 8, width = 4)
for(i in 1:nrow(ER_adjust_summary)){
	grid.arrange(gTree(children = get(paste0("Venn_adjust_N_region", rownames(ER_adjust_summary)[i]))), gTree(children = get(paste0("Venn_adjust_Total_length", rownames(ER_adjust_summary)[i]))), nrow = 2)
	grid.text(paste0(rownames(ER_adjust_summary)[i], "\nfraction_", ER_adjust_summary[i, "frac"]), x = unit(0.8, "npc"), y = unit(0.95, "npc"))
	grid.text("No. of regions", x = unit(0.2, "npc"), y = unit(0.98, "npc"))
	grid.text("Total length", x = unit(0.2, "npc"), y = unit(0.48, "npc"))
}
dev.off()

### genic/intergenic distribution
HM_distrbution <- read.delim("./bam/HM.distribution.summary", as.is = T) %>%
	mutate(Type = rep(c("NPC", rep("glioma", 5)), 2))
(HM_distrbution_figure <- ggplot(HM_distrbution, aes(Sample, genic_intergenic, fill = Type)) + 
	geom_bar(stat = "identity", width = 0.5) + 
	facet_wrap(~ Mark) + 
	xlab("") + 
	ylab("# gene-associated reads / # intergenic reads") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90)))
ggsave(HM_distrbution_figure, file = "HM_distrbution_figure.pdf", width = 8, height = 6)

## -------- Differentially marked regions --------
### test: how to set the background coverage cutoff? CEMT_19 vs GE04 H3K27ac
signal <- read.delim("./unique/test/63.CEMT_19.vs.NPC_GE04.CEMT_19.signal.coverage", as.is = T, head = F, col.names = c("chr", "start", "end", "ID", "cov", "max")) 
pdf("./unique/test/test.pdf")
plot(c(0, 20), c(0, 1), type = "n", main = "ER.unique cutoff", xlab = "average coverage", ylab = "ecdf")
lines(ecdf(signal$cov), col = "red")
for(i in 1:10){
	background <- read.delim(paste0("./unique/test/CEMT_19.vs.NPC_GE04.CEMT_19.background.", i, ".bed.CEMT_19.vs.NPC_GE04.CEMT_19.background.", i, ".coverage"), as.is = T, head = F, col.names = c("chr", "start", "end", "i", "ID", "cov", "max"))
	lines(ecdf(background$cov), col = i + 2)
	abline(v = quantile(background$cov, 0.9), col = i + 2)
}
dev.off()

### summary
ER_unique_summary <- read.delim("./unique/ER.unique.summary", as.is = T) %>% 
	melt(id = c("Sample", "Mark")) %>% 
	mutate(Sample = paste0("CEMT_", Sample), category = ifelse(grepl("len", variable), "Total No. of bases", "Total No. of regions"), type = ifelse(grepl("glioma", variable), "glioma", "NPC"), unique = ifelse(grepl("unique", variable), T, F), value = ifelse(type == "glioma", value, -value))
(ER_unique_summary_figure <- ggplot(ER_unique_summary %>% filter(unique == T), aes(Mark, value, fill = Sample)) + 
	geom_bar(stat = "identity", position = position_dodge()) +
	geom_hline(yintercept = 0) + 
	facet_grid(category ~., scales = "free_y") + 
	xlab("") + 
	ylab("") + 
	theme_bw())
ggsave(ER_unique_summary_figure, file = "./unique/ER_unique_summary_figure.pdf", width = 8, height = 6)

### associated with DE genes 
DHM_DE <- read.delim("./unique/DHM.DE.summary", as.is = T) %>% mutate(Significant = p_Fisher < 0.01, DHM = ifelse(Marked == "NPC_GE04", "Loss", "Gain"))
(DHM_DE_figure <- ggplot(DHM_DE, aes(DHM, Percent_intersect, fill = DE, color = Significant)) + 
	geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
	scale_fill_manual(name = "", values = c("blue", "red")) + 
	scale_color_manual(values = c("transparent", "green")) + 
	facet_grid(Mark ~ Sample) + 
	xlab("") + 
	ylab("Fraction of DE genes") + 
	theme_bw())
ggsave(DHM_DE_figure, file = "./unique/DHM_DE_figure.pdf")

## -------- Chromatin states ----------
chromHMM_summary <- read.delim("./ChromHMM/chromatin.states.summary", as.is = T) 
chromHMM_summary_tall <- melt(chromHMM_summary, id.vars = c("State", "Name"), variable.name = "Sample")
(chromHMM_summary_figure <- ggplot(chromHMM_summary_tall %>% filter(State != "E18"), aes(Name, value/10^6, fill = Sample)) + 
	geom_bar(position = position_dodge(), stat = "identity") + 
	xlab("") + 
	ylab("Total length (Mb)") + 
	coord_flip() + 
	theme_bw())
ggsave(chromHMM_summary_figure, file = "./ChromHMM/chromHMM_summary_figure.pdf")

## H3K36me3 methyltransferases expression
H3K36me3_methyltransferases_RPKM <- read.delim("/projects/epigenomics2/users/lli/glioma/RNAseq/H3K36me3.methyltransferases.RPKM", as.is = T) %>%
	select(-Brain01, -Brain02, -CEMT_21) %>% melt(id = c("ENSG", "Name")) %>% 
	mutate(Type = ifelse(variable == "CEMT_23", "IDHwt", ifelse(grepl("CEMT", variable), "IDHmut", "NPCs"))) %>% 
	group_by(ENSG, Type) %>% summarize(Name = Name[1], mean = mean(value), sd = sd(value)) %>% mutate(sd = ifelse(is.na(sd), 0, sd))
(H3K36me3_methyltransferases_RPKM_figure <- ggplot(H3K36me3_methyltransferases_RPKM, aes(Name, mean, ymax = mean + sd, ymin = mean - sd, color = Type)) + 
	geom_point(size = 3) + 
	geom_errorbar(width = 0.2) + 
	coord_flip() + 
	xlab("") + 
	ylab("RPKM") + 
	theme_bw())
ggsave(H3K36me3_methyltransferases_RPKM_figure, file = "H3K36me3_methyltransferases_RPKM_figure.pdf", height = 5, width = 6)


save(list = c("ER_summary", "ER_adjust_summary", "chromHMM_summary", 
							ls(pattern = "figure"), ls(pattern = "DAVID")), 
		 file = "/projects/epigenomics2/users/lli/glioma/ChIPseq/ChIPseq.Rdata")


