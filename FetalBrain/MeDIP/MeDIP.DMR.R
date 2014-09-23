# DMR analysis from MeDIP fractional calls  

# DMR identification:  
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.signal.sh
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.fractional.m
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.DM.sh

setwd("~/快盘/FetalBrain/MeDIP/DMR/")
library(ggplot2)

# testing DMR script parameters
DM_test_summary <- read.delim("DM.summary.stats.test", head = F, as.is = T, col.names = c("Sample", "delta", "Total.DM.CpGs", "Hyper.DM.CpGs", "Hypo.DM.CpGs"))
DMR_test_summary <- read.delim("DMR.summary.stats.test", head = F, as.is = T, col.names = c("Sample", "delta", "size", "CpG.cut", "Median.length", "Median.CpG", "Total.DMR", "Hyper.DMR", "Hypo.DMR"))
pdf("DMR_test_delta.pdf", width = 9)
(ggplot(DM_test_summary, aes(x = delta, y = Total.DM.CpGs)) + geom_point() + geom_line() + theme_bw() + xlab("delta cutoff") + ylab("Total No. of DM CpGs"))
(ggplot(DM_test_summary, aes(x = delta, y = log2(Hyper.DM.CpGs/Hypo.DM.CpGs))) + geom_point() + geom_line() + geom_hline(aes(yintercept = 0)) + theme_bw() + xlab("delta cutoff") + ylab("log2(Hyper/Hypo)"))
(ggplot(DMR_test_summary, aes(x = delta, y = Median.length, color = factor(size))) + geom_point() + geom_line() + geom_hline(aes(yintercept = 250)) + theme_bw() + xlab("Max distance between adjacent DM CpGs") + ylab("Median DMR length"))
(ggplot(DMR_test_summary, aes(x = delta, y = Total.DMR, , color = factor(size))) + geom_point() + geom_line() + theme_bw() + xlab("delta cutoff") + ylab("Total No. of DMRs"))
(ggplot(DMR_test_summary, aes(x = delta, y = log2(Hyper.DMR/Hypo.DMR), color = factor(size))) + geom_point() + geom_line() + geom_hline(aes(yintercept = 0)) + theme_bw() + xlab("delta cutoff") + ylab("log2(Hyper/Hypo)"))
(ggplot(DMR_test_summary, aes(x = size, y = Median.length, color = factor(delta))) + geom_point() + geom_line() + geom_hline(aes(yintercept = 250)) + theme_bw() + xlab("Max distance between adjacent DM CpGs") + ylab("Median DMR length"))
dev.off()

source("~/HirstLab/Pipeline/DMR.figures.R")
col <- c("chr", "start", "end", "ID", "DM", "CpG_count", "length") # format of DMR files
# Between MZ twins 
setwd("~/快盘/FetalBrain/MeDIP/DMR/MZ/")
Brain01_Brain02_DMR <- read.delim("DMR.Brain-HuFNSC01_Brain-HuFNSC02.d0.4.s200.c3.bed", as.is = T, head = F)
Cortex01_Cortex02_DMR <- read.delim("DMR.Cortex-HuFNSC01_Cortex-HuFNSC02.d0.4.s200.c3.bed", as.is = T, head = F)
GE01_GE02_DMR <- read.delim("DMR.GE-HuFNSC01_GE-HuFNSC02.d0.4.s200.c3.bed", as.is = T, head = F)
DMR_MZ_summary <- data.frame(Total = c(nrow(Brain01_Brain02_DMR), nrow(Cortex01_Cortex02_DMR), nrow(GE01_GE02_DMR)), 
                          Hyper = c(nrow(Brain01_Brain02_DMR[Brain01_Brain02_DMR$V5 == 1, ]), nrow(Cortex01_Cortex02_DMR[Cortex01_Cortex02_DMR$V5 == 1, ]), nrow(GE01_GE02_DMR[GE01_GE02_DMR$V5 == 1, ])), 
                          Hypo = c(nrow(Brain01_Brain02_DMR[Brain01_Brain02_DMR$V5 == -1, ]), nrow(Cortex01_Cortex02_DMR[Cortex01_Cortex02_DMR$V5 == -1, ]), nrow(GE01_GE02_DMR[GE01_GE02_DMR$V5 == -1, ])))
rownames(DMR_MZ_summary) <- c("Brain", "Cortex", "GE")

# sanity check and visualization 
Brain01_Brain02_DMR_figures <- DMR_figures(Brain01_Brain02_DMR, sample1 = "Brain-HuFNSC01", sample2 = "Brain-HuFNSC02")
(Brain01_Brain02_DMR_figures$length)
(Brain01_Brain02_DMR_figures$count)
(Brain01_Brain02_DMR_figures$dis)
(Brain01_Brain02_DMR_figures$freq)
(Brain01_Brain02_DMR_figures$pos)
Cortex01_Cortex02_DMR_figures <- DMR_figures(Cortex01_Cortex02_DMR, sample1 = "Cortex-HuFNSC01", sample2 = "Cortex-HuFNSC02")
(Cortex01_Cortex02_DMR_figures$length)
(Cortex01_Cortex02_DMR_figures$count)
(Cortex01_Cortex02_DMR_figures$dis)
(Cortex01_Cortex02_DMR_figures$freq)
(Cortex01_Cortex02_DMR_figures$pos)
GE01_GE02_DMR_figures <- DMR_figures(GE01_GE02_DMR, sample1 = "GE-HuFNSC01", sample2 = "GE-HuFNSC02")
(GE01_GE02_DMR_figures$length)
(GE01_GE02_DMR_figures$count)
(GE01_GE02_DMR_figures$dis)
(GE01_GE02_DMR_figures$freq)
(GE01_GE02_DMR_figures$pos)

DMR_length_MZ <- rbind(cbind(Brain01_Brain02_DMR_figures$length$data, cell = "Brain"), cbind(Cortex01_Cortex02_DMR_figures$length$data, cell = "Cortex"), cbind(GE01_GE02_DMR_figures$length$data, cell = "GE"))
DMR_length_MZ_figure <- ggplot(DMR_length_MZ, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("DMR length (bp)") + 
  facet_wrap(~ cell) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_length_MZ_figure, file = "DMRlength_MZ.pdf")
(DMR_length_MZ_figure + ggtitle("DMR length between monozygotic twins"))

DMR_count_MZ <- rbind(cbind(Brain01_Brain02_DMR_figures$count$data, cell = "Brain"), cbind(Cortex01_Cortex02_DMR_figures$count$data, cell = "Cortex"), cbind(GE01_GE02_DMR_figures$count$data, cell = "GE"))
DMR_count_MZ_figure <- ggplot(DMR_count_MZ, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("No. of CpGs per DMR") + 
  facet_wrap(~ cell) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_count_MZ_figure, file = "DMRcount_MZ.pdf")
(DMR_count_MZ_figure + ggtitle("No. of CpGs per DMR between monozygotic twins"))

DMR_dis_MZ <- rbind(cbind(Brain01_Brain02_DMR_figures$dis$data, cell = "Brain"), cbind(Cortex01_Cortex02_DMR_figures$dis$data, cell = "Cortex"), cbind(GE01_GE02_DMR_figures$dis$data, cell = "GE"))
DMR_dis_MZ_figure <- ggplot(DMR_dis_MZ, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("Distance between adjacent DMRs (bp)") + 
  facet_wrap(~ cell) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_dis_MZ_figure, file = "DMRdis_MZ.pdf")
(DMR_dis_MZ_figure + ggtitle("Distance between adjacent DMRs between monozygotic twins"))

DMR_freq_MZ <- rbind(cbind(Brain01_Brain02_DMR_figures$freq$data, cell = "Brain"), cbind(Cortex01_Cortex02_DMR_figures$freq$data, cell = "Cortex"), cbind(GE01_GE02_DMR_figures$freq$data, cell = "GE"))
DMR_freq_MZ_figure <- ggplot(DMR_freq_MZ, aes(x = chr, y = freq, fill = DM)) + 
  geom_bar(position = "identity", stat = "identity", width = 0.8) + 
  facet_wrap(~ cell) + 
  xlab("") + 
  ylab("DMR frequency (bp/MB)") + 
  coord_flip() + 
  scale_fill_manual(values = c("red", "blue"), labels = c("hyper", "hypo"), name = "") + 
  theme_bw()
ggsave(DMR_freq_MZ_figure, file = "DMRfreq_MZ.pdf")
(DMR_freq_MZ_figure + ggtitle("DMR frequency between monozygotic twins"))

chrs = c(paste0("chr", as.character(1:22)), "chrX")
chrlength <- read.csv("~/快盘/hg19/chrlen_hg19.csv", as.is = T, row.names = 1)
chrlength <- chrlength[chrlength$chr %in% chrs, ]
chrlength$chr <- factor(chrlength$chr, levels = chrs[1:length(chrs)])
DMR_pos_MZ <- rbind(cbind(Brain01_Brain02_DMR_figures$pos$data, cell = "Brain"), cbind(Cortex01_Cortex02_DMR_figures$pos$data, cell = "Cortex"), cbind(GE01_GE02_DMR_figures$pos$data, cell = "GE"))
DMR_pos_MZ_figure <- ggplot(DMR_pos_MZ) + 
  geom_linerange(aes(x = factor(chr, levels = chr[length(chr):1]), ymin = 0, ymax = length), data = chrlength, alpha = 0.5) + 
  geom_point(aes(x = (as.numeric(chr) + 0.25*DM), y = pos, color = factor(DM, levels = c("1", "-1"))), position = position_jitter(width = 0.05), size = 0.5, alpha = 0.2) +  
  xlab("") + 
  ylab("Position of DMRs on the chromosome") +
  coord_flip() + 
  facet_wrap(~ cell) + 
  scale_color_manual(values = c("red", "blue"), labels = c("hyper", "hypo"), name = "") + 
  theme_bw()
ggsave(DMR_pos_MZ_figure, file = "DMRpos_MZ.pdf")
(DMR_pos_MZ_figure + ggtitle("DMR positions on the chromosomes between monozygotic twins"))

# intersect DMRs with genomic regions 
# /home/lli/bin/shell/DMR.intersect.sh -d /projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR
genomicBreak_MZ <- read.delim("genomic.breakdown.summary", head = F, as.is = T, row.names = 1, col.names = c("Name", "Total", "Intergenic", "Intron", "Exon", "Gene", "Promoter"))
genomicBreak_MZ_tall <- genomicBreak_MZ[, -1]
genomicBreak_MZ_tall <- data.frame(Sample = rep(row.names(genomicBreak_MZ_tall), ncol(genomicBreak_MZ_tall)), Region = factor(rep(colnames(genomicBreak_MZ_tall), each = nrow(genomicBreak_MZ_tall)), levels = colnames(genomicBreak_MZ_tall)), CpG = as.vector(as.matrix(genomicBreak_MZ_tall)))
genomicBreak_MZ_figure <- ggplot(genomicBreak_MZ_tall, aes(x = Region, y = CpG, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5) + 
  xlab("") + 
  ylab("Fraction of CpG") + 
  coord_flip() + 
  scale_fill_manual(values = c("green", "red", "blue"), labels = c("Brain", "Cortex", "GE")) + 
  theme_bw()
ggsave(genomicBreak_MZ_figure, file = "genomicBreak_MZ.pdf")
(genomicBreak_MZ_figure + ggtitle("DM CpG breakdown between monozygotic twins"))

# proximal DMRs (promoter: TSS+/-1.5Kb) and associated genes
load("~/hg19/hg19v65_genes.Rdata")
load("~/快盘/FetalBrain/RNAseq/DEfine/gene/FetalBrain_DEgenes.Rdata")
source('~/HirstLab/Pipeline/enrich.R')
library(VennDiagram)
DMR_gene_MZ_summary <- data.frame(matrix(NA, ncol = 7, nrow = 3, dimnames = list(c("Brain", "Cortex", "GE"), c("pc.Genes", "unique.Genes", "pc.Promoters", "unique.Promoters", "proximal.DE.Genes", "same.direction", "unique.DE.Genes"))))
# Brain
Brain01_Brain02_DMR_pcGene <- read.delim("DMR.Brain-HuFNSC01_Brain-HuFNSC02.CpG_gene_pc.bed", head = F, as.is = T)
Brain01_Brain02_DMR_pcPromoter <- read.delim("DMR.Brain-HuFNSC01_Brain-HuFNSC02.CpG_promoter_pc.bed", head = F, as.is = T)
Brain01_Brain02_DMR_pcPromoter_DE <- Brain01_Brain02_DMR_pcPromoter[Brain01_Brain02_DMR_pcPromoter$V5 %in% brain01_brain02DE$V1, ]
Brain01_Brain02_DMR_pcPromoter_DE <- Brain01_Brain02_DMR_pcPromoter_DE[!duplicated(Brain01_Brain02_DMR_pcPromoter_DE$V4), ]
Brain01_Brain02_DMR_pcPromoter_DE$DM <- as.integer(gsub("chr[0-9X:-]+>", "", Brain01_Brain02_DMR_pcPromoter_DE$V4))
Brain01_Brain02_DMR_pcPromoter_DE$DE <- brain01_brain02DE[Brain01_Brain02_DMR_pcPromoter_DE$V5, "DE"]
DMR_gene_MZ_summary["Brain", "pc.Genes"] <- length(unique(Brain01_Brain02_DMR_pcGene$V4))
DMR_gene_MZ_summary["Brain", "unique.Genes"] <- length(unique(Brain01_Brain02_DMR_pcGene$V5))
DMR_gene_MZ_summary["Brain", "pc.Promoters"] <- length(unique(Brain01_Brain02_DMR_pcPromoter$V4))
DMR_gene_MZ_summary["Brain", "unique.Promoters"] <- length(unique(Brain01_Brain02_DMR_pcPromoter$V5))
DMR_gene_MZ_summary["Brain", "proximal.DE.Genes"] <- nrow(Brain01_Brain02_DMR_pcPromoter_DE)
DMR_gene_MZ_summary["Brain", "same.direction"] <- sum((Brain01_Brain02_DMR_pcPromoter_DE$DM == 1 & Brain01_Brain02_DMR_pcPromoter_DE$DE == "DN")|(Brain01_Brain02_DMR_pcPromoter_DE$DM == -1 & Brain01_Brain02_DMR_pcPromoter_DE$DE == "UP"))
DMR_gene_MZ_summary["Brain", "unique.DE.Genes"] <- length(unique(Brain01_Brain02_DMR_pcPromoter_DE$V5))
Brain01_Brain02_DMR_pcPromoter <- Brain01_Brain02_DMR_pcPromoter[!duplicated(Brain01_Brain02_DMR_pcPromoter$V5), ]
write.table(Brain01_Brain02_DMR_pcPromoter, file = "DMR_pcPromoter.Brain-HuFNSC01_Brain-HuFNSC02.txt", sep = "\t", quote = F, col.names = F, row.names = F)
(Brain01_Brain02_DMR_pcPromoter_DAVID <- enrich(name = "DMR_pcPromoter.Brain01_Brain02", erminej = F, fdr = 1e-6))
Brain01_Brain02_DMR_pcPromoter_DE <- Brain01_Brain02_DMR_pcPromoter_DE[!duplicated(Brain01_Brain02_DMR_pcPromoter_DE$V5), ]
write.table(Brain01_Brain02_DMR_pcPromoter_DE, file = "DMR_pcPromoter_DE.Brain-HuFNSC01_Brain-HuFNSC02.txt", sep = "\t", quote = F, col.names = F, row.names = F)
(Brain01_Brain02_DMR_pcPromoter_DE_DAVID <- enrich(name = "DMR_pcPromoter_DE.Brain01_Brain02", erminej = F))
# Cortex
Cortex01_Cortex02_DMR_pcGene <- read.delim("DMR.Cortex-HuFNSC01_Cortex-HuFNSC02.CpG_gene_pc.bed", head = F, as.is = T)
Cortex01_Cortex02_DMR_pcPromoter <- read.delim("DMR.Cortex-HuFNSC01_Cortex-HuFNSC02.CpG_promoter_pc.bed", head = F, as.is = T)
Cortex01_Cortex02_DMR_pcPromoter_DE <- Cortex01_Cortex02_DMR_pcPromoter[Cortex01_Cortex02_DMR_pcPromoter$V5 %in% cortex01_cortex02DE$V1, ]
Cortex01_Cortex02_DMR_pcPromoter_DE <- Cortex01_Cortex02_DMR_pcPromoter_DE[!duplicated(Cortex01_Cortex02_DMR_pcPromoter_DE$V4), ]
Cortex01_Cortex02_DMR_pcPromoter_DE$DM <- as.integer(gsub("chr[0-9X:-]+>", "", Cortex01_Cortex02_DMR_pcPromoter_DE$V4))
Cortex01_Cortex02_DMR_pcPromoter_DE$DE <- cortex01_cortex02DE[Cortex01_Cortex02_DMR_pcPromoter_DE$V5, "DE"]
DMR_gene_MZ_summary["Cortex", "pc.Genes"] <- length(unique(Cortex01_Cortex02_DMR_pcGene$V4))
DMR_gene_MZ_summary["Cortex", "unique.Genes"] <- length(unique(Cortex01_Cortex02_DMR_pcGene$V5))
DMR_gene_MZ_summary["Cortex", "pc.Promoters"] <- length(unique(Cortex01_Cortex02_DMR_pcPromoter$V4))
DMR_gene_MZ_summary["Cortex", "unique.Promoters"] <- length(unique(Cortex01_Cortex02_DMR_pcPromoter$V5))
DMR_gene_MZ_summary["Cortex", "proximal.DE.Genes"] <- nrow(Cortex01_Cortex02_DMR_pcPromoter_DE)
DMR_gene_MZ_summary["Cortex", "same.direction"] <- sum((Cortex01_Cortex02_DMR_pcPromoter_DE$DM == 1 & Cortex01_Cortex02_DMR_pcPromoter_DE$DE == "DN")|(Cortex01_Cortex02_DMR_pcPromoter_DE$DM == -1 & Cortex01_Cortex02_DMR_pcPromoter_DE$DE == "UP"))
DMR_gene_MZ_summary["Cortex", "unique.DE.Genes"] <- length(unique(Cortex01_Cortex02_DMR_pcPromoter_DE$V5))
Cortex01_Cortex02_DMR_pcPromoter <- Cortex01_Cortex02_DMR_pcPromoter[!duplicated(Cortex01_Cortex02_DMR_pcPromoter$V5), ]
write.table(Cortex01_Cortex02_DMR_pcPromoter, file = "DMR_pcPromoter.Cortex-HuFNSC01_Cortex-HuFNSC02.txt", sep = "\t", quote = F, col.names = F, row.names = F)
(Cortex01_Cortex02_DMR_pcPromoter_DAVID <- enrich(name = "DMR_pcPromoter.Cortex01_Cortex02", erminej = F, fdr = 1e-6))
Cortex01_Cortex02_DMR_pcPromoter_DE <- Cortex01_Cortex02_DMR_pcPromoter_DE[!duplicated(Cortex01_Cortex02_DMR_pcPromoter_DE$V5), ]
write.table(Cortex01_Cortex02_DMR_pcPromoter_DE, file = "DMR_pcPromoter_DE.Cortex-HuFNSC01_Cortex-HuFNSC02.txt", sep = "\t", quote = F, col.names = F, row.names = F)
(Cortex01_Cortex02_DMR_pcPromoter_DE_DAVID <- enrich(name = "DMR_pcPromoter_DE.Cortex01_Cortex02", erminej = F))
# GE
GE01_GE02_DMR_pcGene <- read.delim("DMR.GE-HuFNSC01_GE-HuFNSC02.CpG_gene_pc.bed", head = F, as.is = T)
GE01_GE02_DMR_pcPromoter <- read.delim("DMR.GE-HuFNSC01_GE-HuFNSC02.CpG_promoter_pc.bed", head = F, as.is = T)
GE01_GE02_DMR_pcPromoter_DE <- GE01_GE02_DMR_pcPromoter[GE01_GE02_DMR_pcPromoter$V5 %in% GE01_GE02DE$V1, ]
GE01_GE02_DMR_pcPromoter_DE <- GE01_GE02_DMR_pcPromoter_DE[!duplicated(GE01_GE02_DMR_pcPromoter_DE$V4), ]
GE01_GE02_DMR_pcPromoter_DE$DM <- as.integer(gsub("chr[0-9X:-]+>", "", GE01_GE02_DMR_pcPromoter_DE$V4))
GE01_GE02_DMR_pcPromoter_DE$DE <- GE01_GE02DE[GE01_GE02_DMR_pcPromoter_DE$V5, "DE"]
DMR_gene_MZ_summary["GE", "pc.Genes"] <- length(unique(GE01_GE02_DMR_pcGene$V4))
DMR_gene_MZ_summary["GE", "unique.Genes"] <- length(unique(GE01_GE02_DMR_pcGene$V5))
DMR_gene_MZ_summary["GE", "pc.Promoters"] <- length(unique(GE01_GE02_DMR_pcPromoter$V4))
DMR_gene_MZ_summary["GE", "unique.Promoters"] <- length(unique(GE01_GE02_DMR_pcPromoter$V5))
DMR_gene_MZ_summary["GE", "proximal.DE.Genes"] <- nrow(GE01_GE02_DMR_pcPromoter_DE)
DMR_gene_MZ_summary["GE", "same.direction"] <- sum((GE01_GE02_DMR_pcPromoter_DE$DM == 1 & GE01_GE02_DMR_pcPromoter_DE$DE == "DN")|(GE01_GE02_DMR_pcPromoter_DE$DM == -1 & GE01_GE02_DMR_pcPromoter_DE$DE == "UP"))
DMR_gene_MZ_summary["GE", "unique.DE.Genes"] <- length(unique(GE01_GE02_DMR_pcPromoter_DE$V5))
GE01_GE02_DMR_pcPromoter <- GE01_GE02_DMR_pcPromoter[!duplicated(GE01_GE02_DMR_pcPromoter$V5), ]
write.table(GE01_GE02_DMR_pcPromoter, file = "DMR_pcPromoter.GE-HuFNSC01_GE-HuFNSC02.txt", sep = "\t", quote = F, col.names = F, row.names = F)
(GE01_GE02_DMR_pcPromoter_DAVID <- enrich(name = "DMR_pcPromoter.GE01_GE02", erminej = F, fdr = 1e-6))
GE01_GE02_DMR_pcPromoter_DE <- GE01_GE02_DMR_pcPromoter_DE[!duplicated(GE01_GE02_DMR_pcPromoter_DE$V5), ]
write.table(GE01_GE02_DMR_pcPromoter_DE, file = "DMR_pcPromoter_DE.GE-HuFNSC01_GE-HuFNSC02.txt", sep = "\t", quote = F, col.names = F, row.names = F)
(GE01_GE02_DMR_pcPromoter_DE_DAVID <- enrich(name = "DMR_pcPromoter_DE.GE01_GE02", erminej = F))

DMR_pcPromoter_MZ <- list(Brain = Brain01_Brain02_DMR_pcPromoter$V5, Cortex = Cortex01_Cortex02_DMR_pcPromoter$V5, GE = GE01_GE02_DMR_pcPromoter$V5)
venn_DMR_pcPromoter_MZ <- venn.diagram(DMR_pcPromoter_MZ, filename = NULL, fill = c("green", "red", "blue"), main = "Venn diagram of DE genes with promoter DMR between MZ twins")
pdf("venn_DMR_pcPromoter_MZ.pdf")
plot.new()
grid.draw(venn_DMR_pcPromoter_MZ)
dev.off()
DMR_pcPromoter_DE_MZ <- list(Brain = Brain01_Brain02_DMR_pcPromoter_DE$V5, Cortex = Cortex01_Cortex02_DMR_pcPromoter_DE$V5, GE = GE01_GE02_DMR_pcPromoter_DE$V5)
venn_DMR_pcPromoter_DE_MZ <- venn.diagram(DMR_pcPromoter_DE_MZ, filename = NULL, fill = c("green", "red", "blue"), main = "Venn diagram of DE genes with promoter DMR between MZ twins")
pdf("venn_DMR_pcPromoter_DE_MZ.pdf")
plot.new()
grid.draw(venn_DMR_pcPromoter_DE_MZ)
dev.off()

save(DMR_length_MZ_figure, DMR_count_MZ_figure, DMR_dis_MZ_figure, DMR_freq_MZ_figure, DMR_pos_MZ_figure, 
     DMR_MZ_summary, genomicBreak_MZ, genomicBreak_MZ_figure, DMR_gene_MZ_summary, venn_DMR_pcPromoter_MZ, venn_DMR_pcPromoter_DE_MZ, 
     Brain01_Brain02_DMR_pcPromoter_DAVID, Cortex01_Cortex02_DMR_pcPromoter_DAVID, GE01_GE02_DMR_pcPromoter_DAVID, 
     Brain01_Brain02_DMR_pcPromoter_DE_DAVID, Cortex01_Cortex02_DMR_pcPromoter_DE_DAVID, GE01_GE02_DMR_pcPromoter_DE_DAVID, 
     file = "DMR_MZ.Rdata")

DMR_summary <- read.delim("DMR.summary.stats", head = F, as.is = T, row.names = 1, col.names = c("Sample", "delta", "size", "CpG.cut", "Median.length", "Median.CpG", "Total.DMR", "Hyper.DMR", "Hypo.DMR"))
# Between Cortex and GE
rm(list = ls())
setwd("~/快盘/FetalBrain/MeDIP/DMR/")
source("~/HirstLab/Pipeline/DMR.figures.R")
col <- c("chr", "start", "end", "ID", "DM", "CpG_count", "length") # format of DMR files
Cortex01_GE01_DMR <- read.delim("DMR.Cortex-HuFNSC01_GE-HuFNSC01.d0.6.s500.c3", head = F, as.is = T, col.names = col)
Cortex02_GE02_DMR <- read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.d0.6.s500.c3", head = F, as.is = T, col.names = col)
DMR_Cortex_GE_summary <- DMR_summary[grep("Cortex.*_GE", rownames(DMR_summary)), grep("DMR", colnames(DMR_summary))]
Cortex01_GE01_DMR_figures <- DMR_figures(Cortex01_GE01_DMR, sample1 = "Cortex-HuFNSC01", sample2 = "GE-HuFNSC01")
(Cortex01_GE01_DMR_figures$length)
(Cortex01_GE01_DMR_figures$count)
(Cortex01_GE01_DMR_figures$dis)
(Cortex01_GE01_DMR_figures$freq)
(Cortex01_GE01_DMR_figures$pos)
Cortex02_GE02_DMR_figures <- DMR_figures(Cortex02_GE02_DMR, sample1 = "Cortex-HuFNSC02", sample2 = "GE-HuFNSC02")
(Cortex02_GE02_DMR_figures$length)
(Cortex02_GE02_DMR_figures$count)
(Cortex02_GE02_DMR_figures$dis)
(Cortex02_GE02_DMR_figures$freq)
(Cortex02_GE02_DMR_figures$pos)



