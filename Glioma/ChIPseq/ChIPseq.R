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
Ensembl <- read.delim("/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.EnsID_sorted.HUGO", as.is = T, head = F, row.names = 1)

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

### unique enhancers homer
#### H3K27ac
homer_unique_K27ac_NPC_tf <- read.delim("./unique/H3K27ac/homer/IDHmut_NPC.NPC.tf", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	filter(Percent_with_motif >= 20) %>% arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_unique_K27ac_NPC_tf_figure <- ggplot(homer_unique_K27ac_NPC_tf, aes(TF, Percent_with_motif)) + 
	geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
	coord_flip() + 
	xlab("") + 
	ylab("Percent of enhancers with motif") + 
	theme_bw())
ggsave(homer_unique_K27ac_NPC_tf_figure, file = "./unique/H3K27ac/homer/homer_unique_K27ac_NPC_tf_figure.pdf", height = 4, width = 3)
homer_unique_K27ac_NPC_tf_RPKM <- read.delim("./unique/H3K27ac/homer/IDHmut_NPC.NPC.tf.RPKM", as.is = T) %>%
	select(-Brain01, -Brain02, -CEMT_21) %>% melt(id = c("ENSG", "Motif")) %>% 
	mutate(Type = ifelse(variable == "CEMT_23", "IDHwt", ifelse(grepl("CEMT", variable), "IDHmut", "NPCs"))) %>% 
	group_by(ENSG, Type) %>% summarize(Motif = Motif[1], mean = mean(value), sd = sd(value)) %>% mutate(sd = ifelse(is.na(sd), 0, sd), Motif = factor(Motif, levels = gsub("-halfsite", "", levels(homer_unique_K27ac_NPC_tf$TF))))
(homer_unique_K27ac_NPC_tf_RPKM_figure <- ggplot(homer_unique_K27ac_NPC_tf_RPKM, aes(Motif, mean, ymax = mean + sd, ymin = mean - sd, color = Type)) + 
	geom_point(size = 3) + 
	geom_errorbar(width = 0.4) + 
	coord_flip() + 
	xlab("") + 
	ylab("RPKM") + 
	theme_bw())
ggsave(homer_unique_K27ac_NPC_tf_RPKM_figure, file = "./unique/H3K27ac/homer/homer_unique_K27ac_NPC_tf_RPKM_figure.pdf", height = 4, width = 4)
homer_unique_K27ac_NPC_Sox3_5mC <- read.delim("./unique/H3K27ac/homer/IDHmut_NPC.NPC.Sox3.annotate.5mC", as.is = T) %>% 
	mutate(IDHmut = (CEMT_19 + CEMT_22 + CEMT_47)/3, NPC = (NPC.Cortex02 + NPC.Cortex04 + NPC.GE02 + NPC.GE04)/4, delta = IDHmut - NPC) %>% 
	filter(delta >= 0.2) %>% arrange(delta) %>% mutate(ID = factor(ID, levels = ID))
homer_unique_K27ac_NPC_Sox3_5mC_tall <- homer_unique_K27ac_NPC_Sox3_5mC %>%	select(-(chr:end), -CEMT_21, -IDHmut, -NPC, -delta) %>% melt(id.var = "ID") %>% 
	mutate(variable = ifelse(variable == "CEMT_23", gsub("CEMT", "IDHwt_CEMT", variable), ifelse(grepl("CEMT", variable), gsub("CEMT", "IDHmut_CEMT", variable), as.character(variable))))
(homer_unique_K27ac_NPC_Sox3_5mC_heatmap <- ggplot(homer_unique_K27ac_NPC_Sox3_5mC_tall, aes(x = variable, y = ID, fill = value)) + 
	geom_tile() + 
	scale_fill_gradient(name = " Fractional\nmethylation", low = "lightblue", high = "black") + 
	xlab("") + 
	ylab("") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size = 0), axis.ticks.y = element_line(size = 0)))	
ggsave(homer_unique_K27ac_NPC_Sox3_5mC_heatmap, file = "./unique/H3K27ac/homer/homer_unique_K27ac_NPC_Sox3_5mC_heatmap.pdf", height = 6, width = 6)
homer_unique_K27ac_NPC_Sox3_DN <- read.delim("./unique/H3K27ac/homer/IDHmut_NPC.NPC.Sox3.annotate.closest.gene.RPKM.DN", as.is = T) %>%
	filter(ID %in% homer_unique_K27ac_NPC_Sox3_5mC$ID, abs(dis) <= 2000) %>% distinct(ENSG) %>% mutate(Name = Ensembl[ENSG, "V2"]) 
write.table(homer_unique_K27ac_NPC_Sox3_DN, file = "./unique/H3K27ac/homer/homer_unique_K27ac_NPC_Sox3_DN.txt", row.names = F, quote = F, sep = "\t")
(homer_unique_K27ac_NPC_Sox3_DN_DAVID <- enrich("IDHmut_NPC.NPC.Sox3.annotate.closest.gene.RPKM.DN", dirIn = "./unique/H3K27ac/homer/enrich/", dirOut = "./unique/H3K27ac/homer/enrich/", fdr = 0.05, p = "Benjamini", erminej = F, height = 8) + ggtitle("NPC-unique H3K27ac Sox3 DN"))

#### H3K4me1
homer_unique_K4me1_IDHmut_tf <- read.delim("./unique/H3K4me1/homer/IDHmut_NPC.IDHmut.tf", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_unique_K4me1_IDHmut_tf_figure <- ggplot(homer_unique_K4me1_IDHmut_tf, aes(TF, Percent_with_motif)) + 
	geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
	coord_flip() + 
	xlab("") + 
	ylab("Percent of enhancers with motif") + 
	theme_bw())
ggsave(homer_unique_K4me1_IDHmut_tf_figure, file = "./unique/H3K4me1/homer/homer_unique_K4me1_IDHmut_tf_figure.pdf", height = 6, width = 6)
homer_unique_K4me1_IDHmut_tf_RPKM <- read.delim("./unique/H3K4me1/homer/IDHmut_NPC.IDHmut.tf.RPKM", as.is = T) %>%
	select(-Brain01, -Brain02, -CEMT_21) %>% melt(id = c("ENSG", "Motif")) %>% 
	mutate(Type = ifelse(variable == "CEMT_23", "IDHwt", ifelse(grepl("CEMT", variable), "IDHmut", "NPCs"))) %>% 
	group_by(ENSG, Type) %>% summarize(Motif = Motif[1], mean = mean(value), sd = sd(value)) %>% mutate(sd = ifelse(is.na(sd), 0, sd), Motif = factor(Motif, levels = gsub("-halfsite", "", levels(homer_unique_K4me1_IDHmut_tf$TF))))
(homer_unique_K4me1_IDHmut_tf_RPKM_figure <- ggplot(homer_unique_K4me1_IDHmut_tf_RPKM, aes(Motif, mean, ymax = mean + sd, ymin = mean - sd, color = Type)) + 
	geom_point(size = 3) + 
	geom_errorbar(width = 0.4) + 
	coord_flip() + 
	xlab("") + 
	ylab("RPKM") + 
	theme_bw())
ggsave(homer_unique_K4me1_IDHmut_tf_RPKM_figure, file = "./unique/H3K4me1/homer/homer_unique_K4me1_IDHmut_tf_RPKM_figure.pdf", height = 6, width = 6)
(homer_unique_K4me1_IDHmut_Ascl1_UP_DAVID <- enrich("IDHmut_NPC.IDHmut.Ascl1.annotate.closest.gene.UP", dirIn = "./unique/H3K4me1/homer/enrich/", dirOut = "./unique/H3K4me1/homer/enrich/", fdr = 0.05, p = "Benjamini", erminej = F, height = 4) + ggtitle("glioma-unique H3K4me1 Ascl1 UP"))
(homer_unique_K4me1_IDHmut_Olig2_UP_DAVID <- enrich("IDHmut_NPC.IDHmut.Olig2.annotate.closest.gene.UP", dirIn = "./unique/H3K4me1/homer/enrich/", dirOut = "./unique/H3K4me1/homer/enrich/", fdr = 0.05, p = "Benjamini", erminej = F, height = 4) + ggtitle("glioma-unique H3K4me1 Olig2 UP"))
(homer_unique_K4me1_IDHmut_HEB_UP_DAVID <- enrich("IDHmut_NPC.IDHmut.HEB.annotate.closest.gene.UP", dirIn = "./unique/H3K4me1/homer/enrich/", dirOut = "./unique/H3K4me1/homer/enrich/", fdr = 0.05, p = "Benjamini", erminej = F) + ggtitle("glioma-unique H3K4me1 HEB UP"))
homer_unique_K4me1_NPC_tf <- read.delim("./unique/H3K4me1/homer/IDHmut_NPC.NPC.tf", as.is = T) %>%
	mutate(TF = gsub("\\(.*", "", Motif.Name), Percent_with_motif = as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))) %>% 
	arrange(Percent_with_motif) %>% mutate(TF = factor(TF, levels = TF))
(homer_unique_K4me1_NPC_tf_figure <- ggplot(homer_unique_K4me1_NPC_tf, aes(TF, Percent_with_motif)) + 
	geom_bar(stat = "identity", width = 0.5, fill = "blue") + 
	coord_flip() + 
	xlab("") + 
	ylab("Percent of enhancers with motif") + 
	theme_bw())
ggsave(homer_unique_K4me1_NPC_tf_figure, file = "./unique/H3K4me1/homer/homer_unique_K4me1_NPC_tf_figure.pdf", height = 6, width = 6)
homer_unique_K4me1_NPC_tf_RPKM <- read.delim("./unique/H3K4me1/homer/IDHmut_NPC.NPC.tf.RPKM", as.is = T) %>%
	select(-Brain01, -Brain02, -CEMT_21) %>% melt(id = c("ENSG", "Motif")) %>% 
	mutate(Type = ifelse(variable == "CEMT_23", "IDHwt", ifelse(grepl("CEMT", variable), "IDHmut", "NPCs"))) %>% 
	group_by(ENSG, Type) %>% summarize(Motif = Motif[1], mean = mean(value), sd = sd(value)) %>% mutate(sd = ifelse(is.na(sd), 0, sd), Motif = factor(Motif, levels = gsub("-halfsite", "", levels(homer_unique_K4me1_NPC_tf$TF))))
(homer_unique_K4me1_NPC_tf_RPKM_figure <- ggplot(homer_unique_K4me1_NPC_tf_RPKM, aes(Motif, mean, ymax = mean + sd, ymin = mean - sd, color = Type)) + 
	geom_point(size = 3) + 
	geom_errorbar(width = 0.4) + 
	coord_flip() + 
	xlab("") + 
	ylab("RPKM") + 
	theme_bw())
ggsave(homer_unique_K4me1_NPC_tf_RPKM_figure, file = "./unique/H3K4me1/homer/homer_unique_K4me1_NPC_tf_RPKM_figure.pdf", height = 6, width = 6)
(homer_unique_K4me1_NPC_Sox3_DN_DAVID <- enrich("IDHmut_NPC.NPC.Sox3.annotate.closest.gene.RPKM.DN", dirIn = "./unique/H3K4me1/homer/enrich/", dirOut = "./unique/H3K4me1/homer/enrich/", fdr = 0.05, p = "Benjamini", erminej = F, height = 2) + ggtitle("NPC-unique H3K4me1 Sox3 DN"))
(homer_unique_K4me1_NPC_Sox6_DN_DAVID <- enrich("IDHmut_NPC.NPC.Sox6.annotate.closest.gene.RPKM.DN", dirIn = "./unique/H3K4me1/homer/enrich/", dirOut = "./unique/H3K4me1/homer/enrich/", fdr = 0.05, p = "Benjamini", erminej = F, height = 2) + ggtitle("NPC-unique H3K4me1 Sox6 DN"))
(homer_unique_K4me1_NPC_Lhx2_DN_DAVID <- enrich("IDHmut_NPC.NPC.Lhx2.annotate.closest.gene.RPKM.DN", dirIn = "./unique/H3K4me1/homer/enrich/", dirOut = "./unique/H3K4me1/homer/enrich/", fdr = 0.05, p = "Benjamini", erminej = F, height = 2) + ggtitle("NPC-unique H3K4me1 Lhx2 DN"))

### loss of H3K36me3 
H3K36me3_unique_NPC_5mC <- read.delim("./unique/H3K36me3/IDH_NPC_NPC.unique.5mC", as.is = T) %>% 
	mutate(IDHmut = (CEMT_19 + CEMT_22 + CEMT_47)/3, NPC = (NPC.Cortex02 + NPC.Cortex04 + NPC.GE02 + NPC.GE04)/4, delta = IDHmut - NPC) %>% 
	arrange(delta) %>% mutate(ID = factor(ID, levels = ID))
H3K36me3_unique_NPC_5mC_tall <- H3K36me3_unique_NPC_5mC %>%	select(-(chr:end), -CEMT_21, -IDHmut, -NPC, -delta) %>% melt(id.var = "ID") %>% 
	mutate(variable = ifelse(variable == "CEMT_23", gsub("CEMT", "IDHwt_CEMT", variable), ifelse(grepl("CEMT", variable), gsub("CEMT", "IDHmut_CEMT", variable), as.character(variable))))
(H3K36me3_unique_NPC_5mC_heatmap <- ggplot(H3K36me3_unique_NPC_5mC_tall, aes(x = variable, y = ID, fill = value)) + 
	geom_tile() + 
	scale_fill_gradient(name = " Fractional\nmethylation", low = "lightblue", high = "black") + 
	xlab("") + 
	ylab("") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size = 0), axis.ticks.y = element_line(size = 0)))	
ggsave(H3K36me3_unique_NPC_5mC_heatmap, file = "./unique/H3K36me3/H3K36me3_unique_NPC_5mC_heatmap.pdf", height = 6, width = 6)

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
Histone_modifiers_RPKM <- read.delim("/projects/epigenomics2/users/lli/glioma/RNAseq/Histone.modifiers.RPKM", as.is = T) %>%
	select(-Brain01, -Brain02, -CEMT_21) %>% melt(id = c("ENSG", "Name")) %>% 
	mutate(Type = ifelse(variable == "CEMT_23", "IDHwt", ifelse(grepl("CEMT", variable), "IDHmut", "NPCs"))) %>% 
	group_by(ENSG, Type) %>% summarize(Name = Name[1], mean = mean(value), sd = sd(value)) %>% mutate(sd = ifelse(is.na(sd), 0, sd))
(Histone_modifiers_RPKM_figure <- ggplot(Histone_modifiers_RPKM, aes(Name, mean, ymax = mean + sd, ymin = mean - sd, color = Type)) + 
	geom_point(size = 3) + 
	geom_errorbar(width = 0.2) + 
	coord_flip() + 
	xlab("") + 
	ylab("RPKM") + 
	theme_bw())
ggsave(Histone_modifiers_RPKM_figure, file = "Histone_modifiers_RPKM_figure.pdf", height = 5, width = 6)


save(list = c("ER_summary", "ER_adjust_summary", "chromHMM_summary", "homer_unique_K27ac_NPC_Sox3_5mC_heatmap", "homer_unique_K27ac_NPC_Sox3_DN", 
							ls(pattern = "figure"), ls(pattern = "DAVID")), 
		 file = "/projects/epigenomics2/users/lli/glioma/ChIPseq/ChIPseq.Rdata")


