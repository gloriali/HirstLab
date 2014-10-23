# DMR analysis from MeDIP fractional calls  

# DMR identification:  
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.signal.sh
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.fractional.m
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.DM.sh

setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/")
library(ggplot2)
library(labeling)
library(VennDiagram)
source('~/HirstLab/Pipeline/R/enrich.R')
source("~/HirstLab/Pipeline/R/DMR.figures.R")
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
source('~/HirstLab/Pipeline/R/DMR_DE.R')
load("~/FetalBrain/RNAseq/DEfine/gene/FetalBrain_DEgenes.Rdata")

#' testing DMR script parameters 
DM_test_summary <- read.delim("DM.summary.stats", head = F, as.is = T, col.names = c("Sample", "m", "delta", "Total.DM.CpGs", "Hyper.DM.CpGs", "Hypo.DM.CpGs"))
DMR_test_summary <- read.delim("DMR.summary.stats", head = F, as.is = T, col.names = c("Sample", "size", "CpG.cut", "Median.length", "Median.CpG", "Total.DMR", "Hyper.DMR", "Hypo.DMR"))
DMR_test_summary$m <- gsub(".*m", "", DMR_test_summary$Sample)
DMR_test_summary$m <- as.numeric(gsub(".d.*", "", DMR_test_summary$m))
DMR_test_summary$delta <- as.numeric(gsub(".*d", "", DMR_test_summary$Sample))
DMR_test_summary$Sample <- gsub(".m.*.d.*", "", DMR_test_summary$Sample)
pdf("DMR_parameters.pdf", width = 9)
(ggplot(DM_test_summary, aes(x = delta, y = Total.DM.CpGs, color = Sample)) + geom_point() + geom_line() + facet_wrap(~ m) + theme_bw() + xlab("delta cutoff") + ylab("Total No. of DM CpGs"))
(ggplot(DM_test_summary, aes(x = delta, y = log2(Hyper.DM.CpGs/Hypo.DM.CpGs), color = Sample)) + geom_point() + geom_line() + geom_hline(aes(yintercept = 0)) + facet_wrap(~ m) + theme_bw() + xlab("delta cutoff") + ylab("log2(Hyper/Hypo)"))
(ggplot(DMR_test_summary, aes(x = size, y = Median.length, color = Sample)) + geom_point() + geom_line() + geom_hline(aes(yintercept = 250)) + facet_wrap(~ delta) + theme_bw() + xlab("max distance between adjacent DM CpGs") + ylab("Median DMR length"))
(ggplot(DMR_test_summary, aes(x = size, y = Total.DMR, , color = Sample)) + geom_point() + geom_line() + facet_wrap(~ delta) + theme_bw() + xlab("max distance between adjacent DM CpGs") + ylab("Total No. of DMRs"))
(ggplot(DMR_test_summary, aes(x = CpG.cut, y = Total.DMR, , color = Sample)) + geom_point() + geom_line() + theme_bw() + xlab("min No. of CpGs per DMR") + ylab("Total No. of DMRs"))
dev.off()
# set parameters to m = 0.75, delta = 0.6, size = 300, CpG cut = 4

#' Between MZ twins 
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/")
DMR_summary <- read.delim("DMR.summary.stats", head = F, as.is = T, col.names = c("Sample", "size", "CpG.cut", "Median.length", "Median.CpG", "Total.DMR", "Hyper.DMR", "Hypo.DMR"))
DMR_MZ_summary <- DMR_summary[grepl("HuFNSC01.*HuFNSC02", DMR_summary$Sample) & DMR_summary$CpG.cut == 4, grepl("DMR", colnames(DMR_summary)) | colnames(DMR_summary) == "Sample"]
DMR_MZ_summary$Sample <- gsub(".m0.75.d0.6", "", DMR_MZ_summary$Sample)
col <- c("chr", "start", "end", "ID", "DM", "CpG_count", "length") # format of DMR files
Brain01_Brain02_DMR <- read.delim("DMR.Brain-HuFNSC01_Brain-HuFNSC02.m0.75.d0.6.s300.c4", as.is = T, head = F, col.names = col)
Cortex01_Cortex02_DMR <- read.delim("DMR.Cortex-HuFNSC01_Cortex-HuFNSC02.m0.75.d0.6.s300.c4", as.is = T, head = F, col.names = col)
GE01_GE02_DMR <- read.delim("DMR.GE-HuFNSC01_GE-HuFNSC02.m0.75.d0.6.s300.c4", as.is = T, head = F, col.names = col)

#' sanity check and visualization 
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
  facet_wrap(DM ~ cell) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_length_MZ_figure, file = "DMRlength_MZ.pdf")
(DMR_length_MZ_figure + ggtitle("DMR length between monozygotic twins"))

DMR_count_MZ <- rbind(cbind(Brain01_Brain02_DMR_figures$count$data, cell = "Brain"), cbind(Cortex01_Cortex02_DMR_figures$count$data, cell = "Cortex"), cbind(GE01_GE02_DMR_figures$count$data, cell = "GE"))
DMR_count_MZ_figure <- ggplot(DMR_count_MZ, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("No. of CpGs per DMR") + 
  facet_wrap(DM ~ cell) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_count_MZ_figure, file = "DMRcount_MZ.pdf")
(DMR_count_MZ_figure + ggtitle("No. of CpGs per DMR between monozygotic twins"))

DMR_dis_MZ <- rbind(cbind(Brain01_Brain02_DMR_figures$dis$data, cell = "Brain"), cbind(Cortex01_Cortex02_DMR_figures$dis$data, cell = "Cortex"), cbind(GE01_GE02_DMR_figures$dis$data, cell = "GE"))
DMR_dis_MZ_figure <- ggplot(DMR_dis_MZ, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("Distance between adjacent DMRs (bp)") + 
  facet_wrap(DM ~ cell) + 
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
  scale_fill_manual(values = c("red", "blue"), labels = c("HuFNSC02 UMRs", "HuFNSC01UMRs"), name = "") + 
  theme_bw()
ggsave(DMR_freq_MZ_figure, file = "DMRfreq_MZ.pdf")
(DMR_freq_MZ_figure + ggtitle("DMR frequency asymmetry between monozygotic twins"))

chrs = c(paste0("chr", as.character(1:22)), "chrX")
chrlength <- read.csv("~/hg19/chrlen_hg19.csv", as.is = T, row.names = 1)
chrlength <- chrlength[chrlength$chr %in% chrs, ]
chrlength$chr <- factor(chrlength$chr, levels = chrs[1:length(chrs)])
DMR_pos_MZ <- rbind(cbind(Brain01_Brain02_DMR_figures$pos$data, cell = "Brain"), cbind(Cortex01_Cortex02_DMR_figures$pos$data, cell = "Cortex"), cbind(GE01_GE02_DMR_figures$pos$data, cell = "GE"))
DMR_pos_MZ_figure <- ggplot(DMR_pos_MZ) + 
  geom_linerange(aes(x = factor(chr, levels = chr[length(chr):1]), ymin = 0, ymax = length), data = chrlength, alpha = 0.5) + 
  geom_point(aes(x = (as.numeric(chr) + 0.25*DM), y = pos, color = factor(DM, levels = c("1", "-1"))), position = position_jitter(width = 0.05), size = 1, alpha = 0.2) +  
  xlab("") + 
  ylab("Position of DMRs on the chromosome") +
  coord_flip() + 
  facet_wrap(~ cell) + 
  scale_color_manual(values = c("red", "blue"), labels = c("HuFNSC02 UMRs", "HuFNSC01UMRs"), name = "") + 
  theme_bw()
ggsave(DMR_pos_MZ_figure, file = "DMRpos_MZ.pdf")
(DMR_pos_MZ_figure + ggtitle("DMR positions on the chromosomes between monozygotic twins"))

#' GREAT enrichment 
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/")
(GREAT_HuFNSC01.UMR.Brain <- enrich_GREAT(file = "DMR.Brain-HuFNSC01_Brain-HuFNSC02.hypo", name = "Brain-HuFNSC01.UMRs", height = 4))
(GREAT_HuFNSC02.UMR.Brain <- enrich_GREAT(file = "DMR.Brain-HuFNSC01_Brain-HuFNSC02.hyper", name = "Brain-HuFNSC02.UMRs", height = 6))
(GREAT_HuFNSC01.UMR.Cortex <- enrich_GREAT(file = "DMR.Cortex-HuFNSC01_Cortex-HuFNSC02.hypo", name = "Cortex-HuFNSC01.UMRs", height = 6))
(GREAT_HuFNSC02.UMR.Cortex <- enrich_GREAT(file = "DMR.Cortex-HuFNSC01_Cortex-HuFNSC02.hyper", name = "Cortex-HuFNSC02.UMRs", height = 2))
(GREAT_HuFNSC01.UMR.GE <- enrich_GREAT(file = "DMR.GE-HuFNSC01_GE-HuFNSC02.hypo", name = "GE-HuFNSC01.UMRs", height = 3))
(GREAT_HuFNSC02.UMR.GE <- enrich_GREAT(file = "DMR.GE-HuFNSC01_GE-HuFNSC02.hyper", name = "GE-HuFNSC02.UMRs", height = 3))

#' intersect DMRs with genomic regions 
# /home/lli/HirstLab/Pipeline/shell/DMR.intersect.sh -d /projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/CpG/")
genomicBreak_MZ <- read.delim("genomic.breakdown.summary", head = F, as.is = T, row.names = 1, col.names = c("Name", "Total", "Intergenic", "Intron", "Exon", "Gene", "Promoter", "CGI"))
genomicBreak_MZ <- genomicBreak_MZ[grep("HuFNSC01.*HuFNSC02", rownames(genomicBreak_MZ)), -1]
genomicBreak_MZ_tall <- data.frame(Sample = rep(row.names(genomicBreak_MZ), ncol(genomicBreak_MZ)), Region = factor(rep(colnames(genomicBreak_MZ), each = nrow(genomicBreak_MZ)), levels = colnames(genomicBreak_MZ)), CpG = as.vector(as.matrix(genomicBreak_MZ)))
genomicBreak_MZ_tall$DM <- gsub(".*c4.", "", genomicBreak_MZ_tall$Sample)
genomicBreak_MZ_tall$Sample <- gsub(".m0.75.*", "", genomicBreak_MZ_tall$Sample)
genomicBreak_MZ_figure <- ggplot(genomicBreak_MZ_tall, aes(x = Region, y = CpG, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5) + 
  xlab("") + 
  ylab("Fraction of CpG") + 
  facet_wrap(~ DM) + 
  coord_flip() + 
  scale_fill_manual(values = c("green", "red", "blue"), labels = c("Brain", "Cortex", "GE")) + 
  theme_bw()
ggsave(genomicBreak_MZ_figure, file = "genomicBreak_MZ.pdf")
(genomicBreak_MZ_figure + ggtitle("DMR breakdown between monozygotic twins"))

#' proximal DMRs (promoter: TSS+/-1.5Kb) and associated genes
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/CpG/")
DMR_proximal_Brain01_Brain02_hyper <- read.delim("DMR.Brain-HuFNSC01_Brain-HuFNSC02.m0.75.d0.6.s300.c4.hyper.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_Brain01_Brain02_hyper <- DMR_DE(DMR_gene = DMR_proximal_Brain01_Brain02_hyper, DE = brain01_brain02DE, DM = "hyper", name = "Brain-HuFNSC01_Brain-HuFNSC02.hyper.proximal")
DMR_proximal_Brain01_Brain02_hypo <- read.delim("DMR.Brain-HuFNSC01_Brain-HuFNSC02.m0.75.d0.6.s300.c4.hypo.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_Brain01_Brain02_hypo <- DMR_DE(DMR_gene = DMR_proximal_Brain01_Brain02_hypo, DE = brain01_brain02DE, DM = "hypo", name = "Brain-HuFNSC01_Brain-HuFNSC02.hypo.proximal")
DMR_proximal_Cortex01_Cortex02_hyper <- read.delim("DMR.Cortex-HuFNSC01_Cortex-HuFNSC02.m0.75.d0.6.s300.c4.hyper.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_Cortex01_Cortex02_hyper <- DMR_DE(DMR_gene = DMR_proximal_Cortex01_Cortex02_hyper, DE = cortex01_cortex02DE, DM = "hyper", name = "Cortex-HuFNSC01_Cortex-HuFNSC02.hyper.proximal")
DMR_proximal_Cortex01_Cortex02_hypo <- read.delim("DMR.Cortex-HuFNSC01_Cortex-HuFNSC02.m0.75.d0.6.s300.c4.hypo.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_Cortex01_Cortex02_hypo <- DMR_DE(DMR_gene = DMR_proximal_Cortex01_Cortex02_hypo, DE = cortex01_cortex02DE, DM = "hypo", name = "Cortex-HuFNSC01_Cortex-HuFNSC02.hypo.proximal")
DMR_proximal_GE01_GE02_hyper <- read.delim("DMR.GE-HuFNSC01_GE-HuFNSC02.m0.75.d0.6.s300.c4.hyper.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_GE01_GE02_hyper <- DMR_DE(DMR_gene = DMR_proximal_GE01_GE02_hyper, DE = GE01_GE02DE, DM = "hyper", name = "GE-HuFNSC01_GE-HuFNSC02.hyper.proximal")
DMR_proximal_GE01_GE02_hypo <- read.delim("DMR.GE-HuFNSC01_GE-HuFNSC02.m0.75.d0.6.s300.c4.hypo.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_GE01_GE02_hypo <- DMR_DE(DMR_gene = DMR_proximal_GE01_GE02_hypo, DE = GE01_GE02DE, DM = "hypo", name = "GE-HuFNSC01_GE-HuFNSC02.hypo.proximal")
DMR_proximal_MZ_summary <- rbind(DMR_DE_Brain01_Brain02_hyper$summary, DMR_DE_Brain01_Brain02_hypo$summary, DMR_DE_Cortex01_Cortex02_hyper$summary, DMR_DE_Cortex01_Cortex02_hypo$summary, DMR_DE_GE01_GE02_hyper$summary, DMR_DE_GE01_GE02_hypo$summary)
rownames(DMR_proximal_MZ_summary) <- c("Brain01_Brain02_hyper", "Brain01_Brain02_hypo", "Cortex01_Cortex02_hyper", "Cortex01_Cortex02_hypo", "GE01_GE02_hyper", "GE01_GE02_hypo")

DMR_proximal_MZ_hyper <- list(Brain = DMR_DE_Brain01_Brain02_hyper$DMR_gene$id, Cortex = DMR_DE_Cortex01_Cortex02_hyper$DMR_gene$id, GE = DMR_DE_GE01_GE02_hyper$DMR_gene$id)
venn_DMR_proximal_MZ_hyper <- venn.diagram(DMR_proximal_MZ_hyper, filename = NULL, fill = c("green", "red", "blue"), main = "Venn diagram of proximal HuFNSC02 UMRs between MZ twins", force.unique = T)
pdf("venn_DMR_proximal_MZ_hyper.pdf")
plot.new()
grid.draw(venn_DMR_proximal_MZ_hyper)
dev.off()
DMR_proximal_MZ_hypo <- list(Brain = DMR_DE_Brain01_Brain02_hypo$DMR_gene$id, Cortex = DMR_DE_Cortex01_Cortex02_hypo$DMR_gene$id, GE = DMR_DE_GE01_GE02_hypo$DMR_gene$id)
venn_DMR_proximal_MZ_hypo <- venn.diagram(DMR_proximal_MZ_hypo, filename = NULL, fill = c("green", "red", "blue"), main = "Venn diagram of proximal HuFNSC01 UMRs between MZ twins", force.unique = T)
pdf("venn_DMR_proximal_MZ_hypo.pdf")
plot.new()
grid.draw(venn_DMR_proximal_MZ_hypo)
dev.off()
pcgene <- 19819 # wc -l /projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID
phyper(length(intersect(DMR_DE_Brain01_Brain02_hyper$DMR_gene$id, DMR_DE_Cortex01_Cortex02_hyper$DMR_gene$id)), DMR_proximal_MZ_summary["Brain01_Brain02_hyper", "unique.genes"], pcgene - DMR_proximal_MZ_summary["Brain01_Brain02_hyper", "unique.genes"], DMR_proximal_MZ_summary["Cortex01_Cortex02_hyper", "unique.genes"], lower.tail = F, log = T) 
phyper(length(intersect(DMR_DE_Brain01_Brain02_hypo$DMR_gene$id, DMR_DE_Cortex01_Cortex02_hypo$DMR_gene$id)), DMR_proximal_MZ_summary["Brain01_Brain02_hypo", "unique.genes"], pcgene - DMR_proximal_MZ_summary["Brain01_Brain02_hypo", "unique.genes"], DMR_proximal_MZ_summary["Cortex01_Cortex02_hypo", "unique.genes"], lower.tail = F, log = T) 
phyper(length(intersect(DMR_DE_Brain01_Brain02_hyper$DMR_gene$id, DMR_DE_GE01_GE02_hyper$DMR_gene$id)), DMR_proximal_MZ_summary["Brain01_Brain02_hyper", "unique.genes"], pcgene - DMR_proximal_MZ_summary["Brain01_Brain02_hyper", "unique.genes"], DMR_proximal_MZ_summary["GE01_GE02_hyper", "unique.genes"], lower.tail = F, log = T) 
phyper(length(intersect(DMR_DE_Brain01_Brain02_hypo$DMR_gene$id, DMR_DE_GE01_GE02_hypo$DMR_gene$id)), DMR_proximal_MZ_summary["Brain01_Brain02_hypo", "unique.genes"], pcgene - DMR_proximal_MZ_summary["Brain01_Brain02_hypo", "unique.genes"], DMR_proximal_MZ_summary["GE01_GE02_hypo", "unique.genes"], lower.tail = F, log = T) 
phyper(length(intersect(DMR_DE_Cortex01_Cortex02_hyper$DMR_gene$id, DMR_DE_GE01_GE02_hyper$DMR_gene$id)), DMR_proximal_MZ_summary["Cortex01_Cortex02_hyper", "unique.genes"], pcgene - DMR_proximal_MZ_summary["Cortex01_Cortex02_hyper", "unique.genes"], DMR_proximal_MZ_summary["GE01_GE02_hyper", "unique.genes"], lower.tail = F, log = T) 
phyper(length(intersect(DMR_DE_Cortex01_Cortex02_hypo$DMR_gene$id, DMR_DE_GE01_GE02_hypo$DMR_gene$id)), DMR_proximal_MZ_summary["Cortex01_Cortex02_hypo", "unique.genes"], pcgene - DMR_proximal_MZ_summary["Cortex01_Cortex02_hypo", "unique.genes"], DMR_proximal_MZ_summary["GE01_GE02_hypo", "unique.genes"], lower.tail = F, log = T) 
#' DAVID functional enrichment showed neuron development, but none reached statistical significance.  
 
#' overlap with TFBSs
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/TF/")
DMR_MZ_TF <- read.delim("DMR.MZ.m0.75.d0.6.s300.c4.TF.summary", head = F, as.is = T, col.names = c("TF", "Brain.hypo", "Brain.hyper", "Ratio.Brain", "Cortex.hypo", "Cortex.hyper", "Ratio.Cortex", "GE.hypo", "GE.hyper", "Ratio.GE"))
DMR_MZ_TF_tall <- data.frame(TF = rep(DMR_MZ_TF$TF, 3), Cell = rep(c("Brain", "Cortex", "GE"), each = nrow(DMR_MZ_TF)), hyper = c(DMR_MZ_TF$Brain.hyper, DMR_MZ_TF$Cortex.hyper, DMR_MZ_TF$GE.hyper), hypo = c(DMR_MZ_TF$Brain.hypo, DMR_MZ_TF$Cortex.hypo, DMR_MZ_TF$GE.hypo), Ratio = c(DMR_MZ_TF$Ratio.Brain, DMR_MZ_TF$Ratio.Cortex, DMR_MZ_TF$Ratio.GE))
DMR_MZ_TF_tall$Asymmetry <- NA
DMR_MZ_TF_tall[DMR_MZ_TF_tall$Ratio > 1,]$Asymmetry <- "HuFNSC01 enriched"
DMR_MZ_TF_tall[DMR_MZ_TF_tall$Ratio < 1,]$Asymmetry <- "HuFNSC02 enriched"
DMR_MZ_TF_tall[DMR_MZ_TF_tall$Ratio == 1,]$Asymmetry <- "Equal"
DMR_MZ_TF_tall$Asymmetry <- factor(DMR_MZ_TF_tall$Asymmetry, levels = c("HuFNSC01 enriched", "Equal", "HuFNSC02 enriched"))
(DMR_MZ_TF_figure <- ggplot(DMR_MZ_TF_tall, aes(x = hyper, y = hypo, label = TF, color = Asymmetry)) + 
   geom_point() + 
   geom_abline(intercept=0, slope=1) + 
   geom_text(data = subset(DMR_MZ_TF_tall, (Ratio >= 2 | Ratio <= 0.3)), angle = 45, size = 4, hjust = 0.2, vjust = 0.2) + 
   facet_wrap(~ Cell) + 
   scale_x_log10() + 
   scale_y_log10() + 
   scale_color_hue(l = 50) + 
   theme_bw())
ggsave(DMR_MZ_TF_figure, file = "DMR_MZ_TF.pdf", height = 6)

save(Brain01_Brain02_DMR, Cortex01_Cortex02_DMR, GE01_GE02_DMR, 
     DMR_length_MZ_figure, DMR_count_MZ_figure, DMR_dis_MZ_figure, DMR_freq_MZ_figure, DMR_pos_MZ_figure, DMR_MZ_summary, 
     GREAT_HuFNSC01.UMR.Brain, GREAT_HuFNSC02.UMR.Brain, GREAT_HuFNSC01.UMR.Cortex, GREAT_HuFNSC02.UMR.Cortex, GREAT_HuFNSC01.UMR.GE, GREAT_HuFNSC02.UMR.GE, 
     genomicBreak_MZ, genomicBreak_MZ_figure, DMR_proximal_MZ_summary, DMR_proximal_MZ_hyper, venn_DMR_proximal_MZ_hyper, DMR_proximal_MZ_hypo, venn_DMR_proximal_MZ_hypo, 
     DMR_DE_Brain01_Brain02_hyper, DMR_DE_Brain01_Brain02_hypo, DMR_DE_Cortex01_Cortex02_hyper, DMR_DE_Cortex01_Cortex02_hypo, DMR_DE_GE01_GE02_hyper, DMR_DE_GE01_GE02_hypo, 
     DMR_MZ_TF, DMR_MZ_TF_figure, 
     file = "/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR_MZ.Rdata")

#' Between Cortex and GE
rm(list = ls())
source('~/HirstLab/Pipeline/R/enrich.R')
source("~/HirstLab/Pipeline/R/DMR.figures.R")
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
source('~/HirstLab/Pipeline/R/DMR_DE.R')
load("~/FetalBrain/RNAseq/DEfine/gene/FetalBrain_DEgenes.Rdata")
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/")
DMR_summary <- read.delim("DMR.summary.stats", head = F, as.is = T, col.names = c("Sample", "size", "CpG.cut", "Median.length", "Median.CpG", "Total.DMR", "Hyper.DMR", "Hypo.DMR"))
DMR_neurospheres_summary <- DMR_summary[grepl("Cortex.*_GE", DMR_summary$Sample) & DMR_summary$CpG.cut == 4, grepl("DMR", colnames(DMR_summary)) | colnames(DMR_summary) == "Sample"]
DMR_neurospheres_summary$Sample <- gsub(".m0.75.d0.6", "", DMR_neurospheres_summary$Sample)
col <- c("chr", "start", "end", "ID", "DM", "CpG_count", "length") # format of DMR files
Cortex01_GE01_DMR <- read.delim("DMR.Cortex-HuFNSC01_GE-HuFNSC01.m0.75.d0.6.s300.c4", head = F, as.is = T, col.names = col)
Cortex02_GE02_DMR <- read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.s300.c4", head = F, as.is = T, col.names = col)
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

DMR_length_neurospheres <- rbind(cbind(Cortex01_GE01_DMR_figures$length$data, donor = "HuFNSC01"), cbind(Cortex02_GE02_DMR_figures$length$data, donor = "HuFNSC02"))
DMR_length_neurospheres_figure <- ggplot(DMR_length_neurospheres, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("DMR length (bp)") + 
  facet_wrap(DM ~ donor) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_length_neurospheres_figure, file = "DMRlength_neurospheres.pdf")
(DMR_length_neurospheres_figure + ggtitle("DMR length between Cortex and GE"))

DMR_count_neurospheres <- rbind(cbind(Cortex01_GE01_DMR_figures$count$data, donor = "HuFNSC01"), cbind(Cortex02_GE02_DMR_figures$count$data, donor = "HuFNSC02"))
DMR_count_neurospheres_figure <- ggplot(DMR_count_neurospheres, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("No. of CpGs per DMR") + 
  facet_wrap(DM ~ donor) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_count_neurospheres_figure, file = "DMRcount_neurospheres.pdf")
(DMR_count_neurospheres_figure + ggtitle("No. of CpGs per DMR between Cortex and GE"))

DMR_dis_neurospheres <- rbind(cbind(Cortex01_GE01_DMR_figures$dis$data, donor = "HuFNSC01"), cbind(Cortex02_GE02_DMR_figures$dis$data, donor = "HuFNSC02"))
DMR_dis_neurospheres_figure <- ggplot(DMR_dis_neurospheres, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("Distance between adjacent DMRs (bp)") + 
  facet_wrap(DM ~ donor) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_dis_neurospheres_figure, file = "DMRdis_neurospheres.pdf")
(DMR_dis_neurospheres_figure + ggtitle("Distance between adjacent DMRs between Cortex and GE"))

DMR_freq_neurospheres <- rbind(cbind(Cortex01_GE01_DMR_figures$freq$data, donor = "HuFNSC01"), cbind(Cortex02_GE02_DMR_figures$freq$data, donor = "HuFNSC02"))
DMR_freq_neurospheres_figure <- ggplot(DMR_freq_neurospheres, aes(x = chr, y = freq, fill = DM)) + 
  geom_bar(position = "identity", stat = "identity", width = 0.8) + 
  facet_wrap(~ donor) + 
  xlab("") + 
  ylab("DMR frequency (bp/MB)") + 
  coord_flip() + 
  scale_fill_manual(values = c("red", "blue"), labels = c("hyper", "hypo"), name = "") + 
  theme_bw()
ggsave(DMR_freq_neurospheres_figure, file = "DMRfreq_neurospheres.pdf")
(DMR_freq_neurospheres_figure + ggtitle("DMR frequency asymmetry between Cortex and GE"))

chrs = c(paste0("chr", as.character(1:22)), "chrX")
chrlength <- read.csv("~/hg19/chrlen_hg19.csv", as.is = T, row.names = 1)
chrlength <- chrlength[chrlength$chr %in% chrs, ]
chrlength$chr <- factor(chrlength$chr, levels = chrs[1:length(chrs)])
DMR_pos_neurospheres <- rbind(cbind(Cortex01_GE01_DMR_figures$pos$data, donor = "HuFNSC01"), cbind(Cortex02_GE02_DMR_figures$pos$data, donor = "HuFNSC02"))
DMR_pos_neurospheres_figure <- ggplot(DMR_pos_neurospheres) + 
  geom_linerange(aes(x = factor(chr, levels = chr[length(chr):1]), ymin = 0, ymax = length), data = chrlength, alpha = 0.5) + 
  geom_point(aes(x = (as.numeric(chr) + 0.25*DM), y = pos, color = factor(DM, levels = c("1", "-1"))), position = position_jitter(width = 0.05), size = 1, alpha = 0.2) +  
  xlab("") + 
  ylab("Position of DMRs on the chromosome") +
  coord_flip() + 
  facet_wrap(~ donor) + 
  scale_color_manual(values = c("red", "blue"), labels = c("hyper", "hypo"), name = "") + 
  theme_bw()
ggsave(DMR_pos_neurospheres_figure, file = "DMRpos_neurospheres.pdf")
(DMR_pos_neurospheres_figure + ggtitle("DMR positions on the chromosomes between Cortex and GE"))

#' GREAT enrichment 
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/")
(GREAT_Cortex01.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC01_GE-HuFNSC01.hypo", name = "HuFNSC01-Cortex.UMRs", height = 12))
(GREAT_GE01.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC01_GE-HuFNSC01.hyper", name = "HuFNSC01-GE.UMRs", height = 4))
(GREAT_Cortex02.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC02_GE-HuFNSC02.hypo", name = "HuFNSC02-Cortex.UMRs", height = 12))
(GREAT_GE02.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC02_GE-HuFNSC02.hyper", name = "HuFNSC02-GE.UMRs", height = 8))

#' intersect DMRs with genomic regions 
# /home/lli/HirstLab/Pipeline/shell/DMR.intersect.sh -d /projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/CpG/")
genomicBreak_neurospheres <- read.delim("genomic.breakdown.summary", head = F, as.is = T, row.names = 1, col.names = c("Name", "Total", "Intergenic", "Intron", "Exon", "Gene", "Promoter", "CGI"))
genomicBreak_neurospheres <- genomicBreak_neurospheres[grep("Cortex.*GE", rownames(genomicBreak_neurospheres)), -1]
genomicBreak_neurospheres_tall <- data.frame(Sample = rep(row.names(genomicBreak_neurospheres), ncol(genomicBreak_neurospheres)), Region = factor(rep(colnames(genomicBreak_neurospheres), each = nrow(genomicBreak_neurospheres)), levels = colnames(genomicBreak_neurospheres)), CpG = as.vector(as.matrix(genomicBreak_neurospheres)))
genomicBreak_neurospheres_tall$DM <- gsub(".*c4.", "", genomicBreak_neurospheres_tall$Sample)
genomicBreak_neurospheres_tall$Sample <- gsub(".m0.75.*", "", genomicBreak_neurospheres_tall$Sample)
genomicBreak_neurospheres_figure <- ggplot(genomicBreak_neurospheres_tall, aes(x = Region, y = CpG, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5) + 
  xlab("") + 
  ylab("Fraction of CpG") + 
  facet_wrap(~ DM) + 
  coord_flip() + 
  scale_fill_hue(l = 50, labels = c("HuFNSC01", "HuFNSC02")) + 
  theme_bw()
ggsave(genomicBreak_neurospheres_figure, file = "genomicBreak_neurospheres.pdf")
(genomicBreak_neurospheres_figure + ggtitle("DMR breakdown between Cortex and GE"))

#' proximal DMRs (promoter: TSS+/-1.5Kb) and associated genes
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/CpG/")
DMR_proximal_Cortex01_GE01_hyper <- read.delim("DMR.Cortex-HuFNSC01_GE-HuFNSC01.m0.75.d0.6.s300.c4.hyper.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_Cortex01_GE01_hyper <- DMR_DE(DMR_gene = DMR_proximal_Cortex01_GE01_hyper, DE = cortex01_GE01DE, DM = "hyper", name = "Cortex-HuFNSC01_GE-HuFNSC01.hyper.proximal")
DMR_proximal_Cortex01_GE01_hypo <- read.delim("DMR.Cortex-HuFNSC01_GE-HuFNSC01.m0.75.d0.6.s300.c4.hypo.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_Cortex01_GE01_hypo <- DMR_DE(DMR_gene = DMR_proximal_Cortex01_GE01_hypo, DE = cortex01_GE01DE, DM = "hypo", name = "Cortex-HuFNSC01_GE-HuFNSC01.hypo.proximal")
DMR_proximal_Cortex02_GE02_hyper <- read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.s300.c4.hyper.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_Cortex02_GE02_hyper <- DMR_DE(DMR_gene = DMR_proximal_Cortex02_GE02_hyper, DE = cortex02_GE02DE, DM = "hyper", name = "Cortex-HuFNSC02_GE-HuFNSC02.hyper.proximal")
DMR_proximal_Cortex02_GE02_hypo <- read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.s300.c4.hypo.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_Cortex02_GE02_hypo <- DMR_DE(DMR_gene = DMR_proximal_Cortex02_GE02_hypo, DE = cortex02_GE02DE, DM = "hypo", name = "Cortex-HuFNSC02_GE-HuFNSC02.hypo.proximal")
DMR_proximal_neurospheres_summary <- rbind(DMR_DE_Cortex01_GE01_hyper$summary, DMR_DE_Cortex01_GE01_hypo$summary, DMR_DE_Cortex02_GE02_hyper$summary, DMR_DE_Cortex02_GE02_hypo$summary)
rownames(DMR_proximal_neurospheres_summary) <- c("GE01.UMRs", "Cortex01.UMRs", "GE02.UMRs", "Cortex02.UMRs")

DMR_proximal_neurospheres_hyper <- list(HuFNSC01 = DMR_DE_Cortex01_GE01_hyper$DMR_gene$id, HuFNSC02 = DMR_DE_Cortex02_GE02_hyper$DMR_gene$id)
venn_DMR_proximal_neurospheres_hyper <- venn.diagram(DMR_proximal_neurospheres_hyper, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of proximal GE UMRs between neurospheres twins", force.unique = T)
pdf("venn_DMR_proximal_neurospheres_hyper.pdf")
plot.new()
grid.draw(venn_DMR_proximal_neurospheres_hyper)
dev.off()
DMR_proximal_neurospheres_hypo <- list(HuFNSC01 = DMR_DE_Cortex01_GE01_hypo$DMR_gene$id, HuFNSC02 = DMR_DE_Cortex02_GE02_hypo$DMR_gene$id)
venn_DMR_proximal_neurospheres_hypo <- venn.diagram(DMR_proximal_neurospheres_hypo, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of proximal Cortex UMRs between neurospheres twins", force.unique = T, na = "remove")
pdf("venn_DMR_proximal_neurospheres_hypo.pdf")
plot.new()
grid.draw(venn_DMR_proximal_neurospheres_hypo)
dev.off()
pcgene <- 19819 # wc -l /projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID
phyper(length(intersect(DMR_DE_Cortex01_GE01_hyper$DMR_gene$id, DMR_DE_Cortex02_GE02_hyper$DMR_gene$id)), DMR_proximal_neurospheres_summary["GE01.UMRs", "unique.genes"], pcgene - DMR_proximal_neurospheres_summary["GE01.UMRs", "unique.genes"], DMR_proximal_neurospheres_summary["GE02.UMRs", "unique.genes"], lower.tail = F, log = T) 
phyper(length(intersect(DMR_DE_Cortex01_GE01_hypo$DMR_gene$id, DMR_DE_Cortex02_GE02_hypo$DMR_gene$id)), DMR_proximal_neurospheres_summary["Cortex01.UMRs", "unique.genes"], pcgene - DMR_proximal_neurospheres_summary["Cortex01.UMRs", "unique.genes"], DMR_proximal_neurospheres_summary["Cortex02.UMRs", "unique.genes"], lower.tail = F, log = T) 

#' DAVID functional enrichment for HuFNSC01 have no significant terms.  
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/")
(Cortex02_UMR_proximal_DAVID <- enrich(name = "Cortex-HuFNSC02_GE-HuFNSC02.hypo.proximal", erminej = F, height = 2))
(GE02_UMR_proximal_DAVID <- enrich(name = "Cortex-HuFNSC02_GE-HuFNSC02.hyper.proximal", erminej = F, height = 5))

#' overlap with TFBSs
setwd("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/TF/")
DMR_neurospheres_TF <- read.delim("DMR.Cortex_GE.m0.75.d0.6.s300.c4.TF.summary", head = F, as.is = T, col.names = c("TF", "Cortex01UMR", "GE01UMR", "Ratio01", "Cortex02UMR", "GE02UMR", "Ratio02"))
DMR_neurospheres_TF_tall <- data.frame(TF = rep(DMR_neurospheres_TF$TF, 2), Donor = rep(c("HuFNSC01", "HuFNSC02"), each = nrow(DMR_neurospheres_TF)), Cortex_UMR = c(DMR_neurospheres_TF$Cortex01UMR, DMR_neurospheres_TF$Cortex02UMR), GE_UMR = c(DMR_neurospheres_TF$GE01UMR, DMR_neurospheres_TF$GE02UMR), Ratio = c(DMR_neurospheres_TF$Ratio01, DMR_neurospheres_TF$Ratio02))
DMR_neurospheres_TF_tall$Asymmetry <- NA
DMR_neurospheres_TF_tall[DMR_neurospheres_TF_tall$Ratio > 1,]$Asymmetry <- "Cortex enriched"
DMR_neurospheres_TF_tall[DMR_neurospheres_TF_tall$Ratio < 1,]$Asymmetry <- "GE enriched"
DMR_neurospheres_TF_tall[DMR_neurospheres_TF_tall$Ratio == 1,]$Asymmetry <- "Equal"
(DMR_neurospheres_TF_figure <- ggplot(DMR_neurospheres_TF_tall, aes(x = Cortex_UMR, y = GE_UMR, label = TF, color = Asymmetry)) + 
   geom_point() + 
   geom_abline(intercept=0, slope=1) + 
   geom_text(data = subset(DMR_neurospheres_TF_tall, (Ratio >= 2 | Ratio <= 0.5)), angle = 30, size = 4, hjust = 0.2, vjust = 0.2) + 
   facet_wrap(~ Donor) + 
   scale_x_log10() + 
   scale_y_log10() + 
   scale_color_hue(l = 50) + 
   theme_bw())
ggsave(DMR_neurospheres_TF_figure, file = "DMR_neurospheres_TF.pdf", height = 6)

save(Cortex01_GE01_DMR, Cortex02_GE02_DMR, 
     DMR_length_neurospheres_figure, DMR_count_neurospheres_figure, DMR_dis_neurospheres_figure, DMR_freq_neurospheres_figure, DMR_pos_neurospheres_figure, DMR_neurospheres_summary, 
     GREAT_Cortex01.UMR, GREAT_GE01.UMR, GREAT_Cortex02.UMR, GREAT_GE02.UMR, 
     genomicBreak_neurospheres, genomicBreak_neurospheres_figure, DMR_proximal_neurospheres_summary, 
     DMR_DE_Cortex01_GE01_hyper, DMR_DE_Cortex01_GE01_hypo, DMR_DE_Cortex02_GE02_hyper, DMR_DE_Cortex02_GE02_hypo, 
     DMR_proximal_neurospheres_hyper, DMR_proximal_neurospheres_hypo, venn_DMR_proximal_neurospheres_hyper, venn_DMR_proximal_neurospheres_hypo, 
     Cortex02_UMR_proximal_DAVID, GE02_UMR_proximal_DAVID, DMR_neurospheres_TF, DMR_neurospheres_TF_figure, 
     file = "/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR_neurospheres.Rdata")


