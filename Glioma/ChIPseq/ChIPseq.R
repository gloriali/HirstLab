# Glioma - ChIP-seq analysis

library(ggplot2)
library(plyr)
library(VennDiagram)
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
(ER_summary_figure <- ggplot(ER_summary_tall %>% filter(variable %in% c("N_region", "Total_length", "Seq_depth")), aes(Mark, value, fill = Sample)) + 
	geom_bar(position = position_dodge(width=0.9), stat = "identity") + 
#	geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width=0.9), width=0.5) + 
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
(ER_unique_summary_figure <- ggplot(ER_unique_summary %>% filter(unique == T, Sample %in% libs), aes(Mark, value, fill = Sample)) + 
	geom_bar(stat = "identity", position = position_dodge()) +
	geom_hline(yintercept = 0) + 
	facet_grid(category ~., scales = "free_y") + 
	xlab("") + 
	ylab("") + 
	theme_bw())
ggsave(ER_unique_summary_figure, file = "./unique/ER_unique_summary_figure.pdf", width = 8, height = 6)
H3K27me3_unique_summary <- ER_unique_summary %>% filter(Mark == "H3K27me3", variable %in% c("len_glioma_unique", "len_NPC_unique")) %>%
  mutate(Sample = ifelse(Sample == "CEMT_23", gsub("CEMT", "IDHwt_CEMT", Sample), gsub("CEMT", "IDHmut_CEMT", Sample)), value = abs(value)/10^6, type = paste0(type, "-specific"))
(H3K27me3_unique_summary_figure <- ggplot(H3K27me3_unique_summary, aes(Sample, value, fill = type)) + 
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) + 
  scale_fill_manual(values = c("red", "blue"), name = "", labels = c("Gain", "Loss")) + 
  xlab("") + 
  ylab("Total length (MB)") + 
  ggtitle("H3K27me3") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)))
ggsave(H3K27me3_unique_summary_figure, file = "./unique/H3K27me3_unique_summary_figure.pdf", height = 5, width = 5)

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

### ======== unique enhancers =========
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

### ============ loss of H3K36me3 ===========
H3K36me3_unique_summary <- ER_unique_summary %>% filter(Mark == "H3K36me3", variable %in% c("len_glioma_unique", "len_NPC_unique")) %>%
	mutate(Sample = ifelse(Sample == "CEMT_23", gsub("CEMT", "IDHwt_CEMT", Sample), ifelse(Sample %in% libs, gsub("CEMT", "IDHmut_CEMT", Sample), Sample)), value = abs(value)/10^6, type = paste0(type, "-specific"))
(H3K36me3_unique_summary_figure <- ggplot(H3K36me3_unique_summary, aes(Sample, value, fill = type)) + 
	geom_bar(stat = "identity", width = 0.5, position = position_dodge()) + 
	scale_fill_manual(values = c("red", "blue"), name = "", labels = c("Gain", "Loss")) + 
	xlab("") + 
	ylab("Total length (MB)") + 
	ggtitle("H3K36me3 unique regions") +  
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90)))
ggsave(H3K36me3_unique_summary_figure, file = "./unique/H3K36me3_unique_summary_figure.pdf", height = 5, width = 5)
(H3K36me3_unique_GREAT <- enrich_GREAT(file = "IDH_NPC_NPC.unique", name = "H3K36me3_loss_regions", dirIn = "./unique/H3K36me3/enrich/", dirOut = "./unique/H3K36me3/enrich/", categories = c("GOBP", "MousePhenotype", "MSigPathway")))
H3K36me3_unique_enrich_summary <- read.delim("./unique/DHM.enrich.summary", as.is = T) %>% filter(Mark == "H3K36me3") %>% mutate(type = ifelse(Marked == "NPC_GE04", "loss", "gain")) %>% 
	select(Sample, type, FC_genebody, FC_intergenic) %>% melt(id = c("Sample", "type")) %>% mutate(region = gsub("FC_", "", variable))
(H3K36me3_unique_enrich_summary_figure <- ggplot(H3K36me3_unique_enrich_summary, aes(Sample, log2(value), fill = region)) + 
	geom_bar(stat = "identity", position = position_dodge()) + 
	facet_wrap(~ type) + 
	coord_flip() + 
	xlab("") + 
	ylab("log2 Fold change") + 
	ggtitle("H3K36me3 unique regions enrichment") + 
	theme_bw())
ggsave(H3K36me3_unique_enrich_summary_figure, file = "./unique/H3K36me3/H3K36me3_unique_enrich_summary_figure.pdf")
H3K36me3_loss_intersect_summary <- read.delim("./unique/H3K36me3/H3K36me3.loss.intersect.summary", as.is = T) %>% mutate(samples = paste0(Sample1, "&", Sample2))
(H3K36me3_loss_intersect_figure <- ggplot(H3K36me3_loss_intersect_summary, aes(samples, log2(FC))) + 
	geom_bar(stat = "identity", fill = "blue") + 
	coord_flip() + 
	xlab("") + 
	ylab("log2 Fold enrichment") + 
	ggtitle("H3K36me3 loss regions \noverlap between samples") + 
	theme_bw())
ggsave(H3K36me3_loss_intersect_figure, file = "./unique/H3K36me3/H3K36me3_loss_intersect_figure.pdf")
H3K36me3_loss_intersect_venn <- draw.pairwise.venn(area1 = H3K36me3_loss_intersect_summary[1, "N1"], area2 = H3K36me3_loss_intersect_summary[1, "N2"], cross.area = H3K36me3_loss_intersect_summary[1, "Intersect"], category = c(H3K36me3_loss_intersect_summary[1, "Sample1"], H3K36me3_loss_intersect_summary[1, "Sample2"]), fill = c("red", "blue"))
#### RPKM
e <- 1e-6
colname <- c("ENSG", "chr", "start", "end", "ID", "geneChr", "geneStart", "geneEnd", "NPC.Cortex01", "NPC.GE01", "NPC.Cortex02", "NPC.GE02", "Brain01", "NPC.Cortex03", "Brain02", "NPC.GE03", "NPC.Cortex04", "NPC.GE04", "CEMT_19", "CEMT_21", "CEMT_22", "CEMT_23", "CEMT_47")
H3K36me3_unique_NPC_gene_RPKM <- read.delim("./unique/H3K36me3/IDH_NPC_NPC.unique.gene.RPKM", as.is = T, head = F, col.names = colname) %>% 
  mutate(IDHmut = (CEMT_19 + CEMT_22 + CEMT_47)/3, NPC = (NPC.Cortex02 + NPC.Cortex04 + NPC.GE02 + NPC.GE04)/4, logFC = log2((IDHmut + e)/(NPC + e))) %>% 
  arrange(logFC) %>% mutate(ENSG = factor(ENSG, levels = ENSG))
(H3K36me3_unique_NPC_gene_RPKM_figure <- ggplot(H3K36me3_unique_NPC_gene_RPKM, aes(logFC)) + 
	geom_density() + 
	coord_cartesian(xlim = c(-10, 10)) + 
	xlab("RPKM log2 IDHmut/NPC") + 
	theme_bw())
ggsave(H3K36me3_unique_NPC_gene_RPKM_figure, file = "./unique/H3K36me3_unique_NPC_gene_RPKM_figure.pdf", height = 6, width = 6)
H3K36me3_unique_NPC_gene_RPKM_tall <- H3K36me3_unique_NPC_gene_RPKM %>%	distinct(ENSG) %>% select(ENSG, NPC.Cortex02, NPC.Cortex04, NPC.GE02, NPC.GE04, CEMT_19, CEMT_22, CEMT_23, CEMT_47) %>% melt(id.var = "ENSG") %>% 
  mutate(variable = ifelse(variable == "CEMT_23", gsub("CEMT", "IDHwt_CEMT", variable), ifelse(grepl("CEMT", variable), gsub("CEMT", "IDHmut_CEMT", variable), as.character(variable))))
(H3K36me3_unique_NPC_gene_RPKM_heatmap <- ggplot(H3K36me3_unique_NPC_gene_RPKM_tall, aes(x = variable, y = ENSG, fill = log10(value))) + 
  geom_tile() + 
  scale_fill_gradient(name = "log10 RPKM", low = "lightgreen", high = "black") + 
  xlab("") + 
  ylab("") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), axis.text.y = element_text(size = 0), axis.ticks.y = element_line(size = 0)))	
ggsave(H3K36me3_unique_NPC_gene_RPKM_heatmap, file = "./unique/H3K36me3/H3K36me3_unique_NPC_gene_RPKM_heatmap.pdf", height = 6, width = 6)
H3K36me3_loss_DE_summary <- read.delim("./unique/H3K36me3/H3K36me3_loss_DE.summary", as.is = T) %>% 
	filter(grepl("CEMT", File)) %>% mutate(File = gsub("\\.vs.*", "", File))
(H3K36me3_loss_gene_DAVID <- enrich(name = "IDH_NPC_NPC.unique_gene", dirIn = "./unique/H3K36me3/enrich/", dirOut = "./unique/H3K36me3/enrich/", erminej = F, category = "GOBP"))
(H3K36me3_loss_gene_DN_DAVID <- enrich(name = "IDH_NPC_NPC.unique_gene_DN", dirIn = "./unique/H3K36me3/enrich/", dirOut = "./unique/H3K36me3/enrich/", erminej = F))
#### 5mC
H3K36me3_unique_NPC_5mC <- read.delim("./unique/H3K36me3/IDH_NPC_NPC.unique.5mC", as.is = T) %>% 
  mutate(IDHmut = (CEMT_19 + CEMT_22 + CEMT_47)/3, NPC = (NPC.Cortex02 + NPC.Cortex04 + NPC.GE02 + NPC.GE04)/4, delta = IDHmut - NPC) %>% 
  filter(ID %in% H3K36me3_unique_NPC_gene_RPKM$ID) %>% arrange(delta) %>% mutate(ID = factor(ID, levels = ID))
(H3K36me3_unique_NPC_5mC_figure <- ggplot(H3K36me3_unique_NPC_5mC, aes(delta)) + 
	geom_density() + 
	coord_cartesian(xlim = c(-0.5, 0.5)) + 
	xlab("delta fractional methylation IDHmut - NPC") + 
	theme_bw())
ggsave(H3K36me3_unique_NPC_5mC_figure, file = "./unique/H3K36me3_unique_NPC_5mC_figure.pdf", height = 6, width = 6)
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
#### H3K27me3
for(lib in libs[1:5]){
  print(lib)
  assign(paste0("H3K36me3_unique_NPC_K27me3_K36me3_", lib), read.delim(paste0("./unique/H3K36me3/", lib, ".vs.NPC_GE04.NPC_GE04.unique.H3K27me3.H3K36me3"), as.is = T))
  assign(paste0("H3K36me3_unique_NPC_K27me3_K36me3_", lib, "_density2d"), ggplot(get(paste0("H3K36me3_unique_NPC_K27me3_K36me3_", lib)), aes(diff_H3K27me3, diff_H3K36me3)) + 
           geom_point(alpha = 0.2, size = 0.5) + 
           stat_density2d(aes(fill = ..level..), geom="polygon") + 
           coord_cartesian(xlim = c(-25, 25), ylim = c(-30, 0)) + 
           ggtitle(lib) + 
           xlab("H3K27me3 normalized signal\nGlioma - NPC") + 
           ylab("H3K36me3 normalized signal\nGlioma - NPC") + 
           theme_bw())
  print(get(paste0("H3K36me3_unique_NPC_K27me3_K36me3_", lib, "_density2d")))
  ggsave(get(paste0("H3K36me3_unique_NPC_K27me3_K36me3_", lib, "_density2d")), file = paste0("./unique/H3K36me3/H3K36me3_unique_NPC_K27me3_K36me3_", lib, "_density2d.pdf"), height = 5, width = 6)
  assign(paste0("H3K36me3_unique_NPC_K27me3_K36me3_", lib, "_hist"), ggplot(get(paste0("H3K36me3_unique_NPC_K27me3_K36me3_", lib)), aes(diff_H3K27me3)) + 
           geom_histogram(fill = "blue", binwidth = 0.1) + 
           geom_vline(xintercept = 0) + 
           coord_cartesian(xlim = c(-10, 10), ylim = c(0, 4000)) + 
           ggtitle(lib) + 
           xlab("H3K27me3 normalized signal\nGlioma - NPC") + 
           ylab("No. of regions") + 
           theme_bw())
  print(get(paste0("H3K36me3_unique_NPC_K27me3_K36me3_", lib, "_hist")))
  ggsave(get(paste0("H3K36me3_unique_NPC_K27me3_K36me3_", lib, "_hist")), file = paste0("./unique/H3K36me3/H3K36me3_unique_NPC_K27me3_K36me3_", lib, "_hist.pdf"), height = 5, width = 6)
}
H3K36me3_loss_H3K27me3_summary <- read.delim("./unique/H3K36me3/H3K36me3.loss.H3K27me3.summary", as.is = T)
(H3K36me3_loss_H3K27me3_summary_figure <- ggplot(H3K36me3_loss_H3K27me3_summary, aes(Sample, FC_gene, fill = H3K27me3)) + 
	geom_bar(stat = "identity", position = position_dodge()) + 
	coord_flip() + 
	xlab("") + 
	ylab("Fold enrichment in genebody") + 
	theme_bw())
ggsave(H3K36me3_loss_H3K27me3_summary_figure, file = "./unique/H3K36me3/H3K36me3_loss_H3K27me3_summary_figure.pdf")
H3K36me3_loss_H3K27me3_intersect_summary <- read.delim("./unique/H3K36me3/H3K36me3.loss.H3K27me3.intersect.summary", as.is = T)
(H3K36me3_loss_H3K27me3gain_GREAT <- enrich_GREAT(file = "IDH_NPC_NPC.unique.H3K27me3.H3K36me3.K27gain", name = "H3K36me3_loss_H3K27me3gain", dirIn = "./unique/H3K36me3/enrich/", dirOut = "./unique/H3K36me3/enrich/", categories = c("GOBP", "MousePhenotype", "HumanPhenotype")))
(H3K36me3_loss_H3K27me3loss_GREAT <- enrich_GREAT(file = "IDH_NPC_NPC.unique.H3K27me3.H3K36me3.K27loss", name = "H3K36me3_loss_H3K27me3loss", dirIn = "./unique/H3K36me3/enrich/", dirOut = "./unique/H3K36me3/enrich/", categories = c("GOBP", "MousePhenotype", "HumanPhenotype")))
#### genic/intergenic distribution
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
Histone_modifiers_RPKM <- read.delim("../RNAseq/Histone.modifiers.RPKM", as.is = T) %>%
	select(-Brain01, -Brain02, -CEMT_21) %>% melt(id = c("ENSG", "Name")) %>% filter(Name %in% c("NSD1", "NSD2", "SETD2")) %>%
	mutate(Type = ifelse(variable == "CEMT_23", "IDHwt", ifelse(grepl("CEMT", variable), "IDHmut", "NPCs"))) %>% 
	group_by(ENSG, Type) %>% summarize(Name = Name[1], mean = mean(value), min = min(value), max = max(value)) 
(Histone_modifiers_RPKM_figure <- ggplot(Histone_modifiers_RPKM, aes(Name, mean, ymax = max, ymin = min, color = Type)) + 
	geom_point(size = 3) + 
	geom_errorbar(width = 0.2) + 
	coord_flip() + 
	xlab("") + 
	ylab("RPKM") + 
	theme_bw())
ggsave(Histone_modifiers_RPKM_figure, file = "Histone_modifiers_RPKM_figure.pdf", height = 5, width = 6)

## HOXA/B expression
HOXA_HOAXB_RPKM <- read.delim("../RNAseq/HOXA.HOXB.RPKM", as.is = T) %>%
	select(-Brain01, -Brain02, -CEMT_21) %>% melt(id = c("ENSG", "Name")) %>%
	mutate(Type = ifelse(variable == "CEMT_23", "IDHwt", ifelse(grepl("CEMT", variable), "IDHmut", "NPCs"))) %>% 
	group_by(ENSG, Type) %>% summarize(Name = Name[1], mean = mean(value), min = min(value), max = max(value)) 
(HOXA_HOAXB_RPKM_figure <- ggplot(HOXA_HOAXB_RPKM, aes(Name, mean, ymax = max, ymin = min, color = Type)) + 
	geom_point(size = 3) + 
	geom_errorbar(width = 0.2) + 
	coord_flip() + 
	scale_y_log10() + 
	xlab("") + 
	ylab("log10 RPKM") + 
	theme_bw())
ggsave(HOXA_HOAXB_RPKM_figure, file = "HOXA_HOAXB_RPKM_figure.pdf", height = 5, width = 6)

## Global comparisons
for(mark in c("H3K27me3", "H3K36me3", "H3K9me3", "H3K27ac", "H3K4me1", "H3K4me3", "Input")){
	assign(paste0("global_", mark), read.delim(paste0("./global/", mark, ".glioma_NPC.quantile"), as.is = T))
	assign(paste0("global_", mark, "_figure"), ggplot(get(paste0("global_", mark)), aes(S2, ymin = log2(ymin), lower = log2(bottom), middle = log2(middle), upper = log2(top), ymax = log2(ymax), fill = Type)) + 
				 	geom_boxplot(stat = "identity", width = 0.8, posistion = position_dodge()) + 
				 	facet_wrap(~ S1, scales = "free_y", ncol = 1) + 
				 	scale_fill_manual(name = "", values = c("red", "blue")) + 
				 	xlab("") + 
				 	ylab("log2 ratio glioma/NPC") + 
				 	ggtitle(mark) + 
				 	theme_bw() + 
				 	theme(axis.text.x = element_text(angle = 90, size = 8)))
	(get(paste0("global_", mark, "_figure")))
	ggsave(get(paste0("global_", mark, "_figure")), file = paste0("./global/global_", mark, "_figure.pdf"), height = 10, width = 10)
}
for(mark in c("H3K27me3", "H3K36me3", "H3K9me3", "H3K27ac", "H3K4me1", "H3K4me3", "Input")){
	e <- 1e-5
	assign(paste0("CEMT_19.NPC_GE04_", mark), read.delim(paste0("./global/", mark, "/CEMT_19.NPC_GE04.coverage"), as.is = T, head = F))
	assign(paste0("density2d_CEMT_19.NPC_GE04_", mark), ggplot(get(paste0("CEMT_19.NPC_GE04_", mark)), aes(log10(V2 + e), log10(V3 + e))) + 
				 	geom_point(data = get(paste0("CEMT_19.NPC_GE04_", mark)) %>% filter(V2 > 10, V3 > 10), aes(log10(V2 + e), log10(V3 + e)), alpha = 0.5, size = 0.5) + 
				 	stat_density2d(aes(fill = ..level..), geom="polygon") + 
				 	ggtitle(mark) + 
				 	xlab("log10 CEMT_19") + 
				 	ylab("log10 NPC_GE04") + 
				 	theme_bw())
	print(get(paste0("density2d_CEMT_19.NPC_GE04_", mark)))
	ggsave(get(paste0("density2d_CEMT_19.NPC_GE04_", mark)), file = paste0("./global/density2d_CEMT_19.NPC_GE04_", mark, ".pdf"), height = 5, width = 6)
}	

## ========== save ===========
save(list = c("homer_unique_K27ac_NPC_Sox3_DN", 
							ls(pattern = "summary"), ls(pattern = "figure"), ls(pattern = "heatmap"), ls(pattern = "density2d"), ls(pattern = "hist"), ls(pattern = "DAVID"), ls(pattern = "GREAT"), ls(pattern = "venn")), 
		 file = "/projects/epigenomics2/users/lli/glioma/ChIPseq/ChIPseq.Rdata")


