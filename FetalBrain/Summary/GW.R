# DMR analysis between 13-week and 17-week 
# DMR identification: ~/HirstLab/FetalBrain/WGBS/WGBS.DM.sh

library(ggplot2)
library(dplyr)
library(ggbio)
library(GenomicRanges)
library(labeling)
library(VennDiagram)
library(gridExtra)
library(gplots)
library(reshape)
source('~/HirstLab/Pipeline/R/enrich.R')
source("~/HirstLab/Pipeline/R/DMR.figures.R")
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
source('~/HirstLab/Pipeline/R/DMR_DE.R')
source("~/HirstLab/Pipeline/R/isoform.R")
load("~/hg19/hg19v65_genes.Rdata")

# ============ Sanity =========
setwd("/projects/epigenomics/users/lli/FetalBrain/GW/DMR/")
GW_DMR_summary <- read.delim("DMR.summary.stats", head = F, as.is = T, col.names = c("Sample", "size", "CpG.cut", "Median.length", "Median.CpG", "Total.DMR", "Hyper.DMR", "Hypo.DMR"))
GW_DMR_summary$Sample <- gsub(".m0.75.p0.005.d0.5", "", GW_DMR_summary$Sample)
GW_DMR_summary <- select(GW_DMR_summary, Sample, Total.DMR, Hyper.DMR, Hypo.DMR)
col <- c("chr", "start", "end", "ID", "DM", "CpG_count", "length") # format of DMR files
Cortex02_Cortex04_DMR <- read.delim("DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3", as.is = T, head = F, col.names = col)
GE02_GE04_DMR <- read.delim("DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3", as.is = T, head = F, col.names = col)
Cortex02_Cortex04_DMR_figures <- DMR_figures(Cortex02_Cortex04_DMR, sample1 = "Cortex-HuFNSC02", sample2 = "Cortex-HuFNSC04")
(Cortex02_Cortex04_DMR_figures$length)
(Cortex02_Cortex04_DMR_figures$count)
(Cortex02_Cortex04_DMR_figures$dis)
(Cortex02_Cortex04_DMR_figures$freq)
(Cortex02_Cortex04_DMR_figures$pos)
GE02_GE04_DMR_figures <- DMR_figures(GE02_GE04_DMR, sample1 = "GE-HuFNSC02", sample2 = "GE-HuFNSC04")
(GE02_GE04_DMR_figures$length)
(GE02_GE04_DMR_figures$count)
(GE02_GE04_DMR_figures$dis)
(GE02_GE04_DMR_figures$freq)
(GE02_GE04_DMR_figures$pos)

DMR_length_GW <- rbind(cbind(Cortex02_Cortex04_DMR_figures$length$data, cell = "Cortex"), cbind(GE02_GE04_DMR_figures$length$data, cell = "GE"))
DMR_length_GW_figure <- ggplot(DMR_length_GW, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("DMR length (bp)") + 
  facet_wrap(DM ~ cell) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_length_GW_figure, file = "DMRlength_GW.pdf", width = 10, height = 10)
(DMR_length_GW_figure + ggtitle("DMR length between 13-week and 17-week"))

DMR_count_GW <- rbind(cbind(Cortex02_Cortex04_DMR_figures$count$data, cell = "Cortex"), cbind(GE02_GE04_DMR_figures$count$data, cell = "GE"))
DMR_count_GW_figure <- ggplot(DMR_count_GW, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("No. of CpGs per DMR") + 
  facet_wrap(DM ~ cell) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_count_GW_figure, file = "DMRcount_GW.pdf", width = 10, height = 10)
(DMR_count_GW_figure + ggtitle("No. of CpGs per DMR between 13-week and 17-week"))

DMR_dis_GW <- rbind(cbind(Cortex02_Cortex04_DMR_figures$dis$data, cell = "Cortex"), cbind(GE02_GE04_DMR_figures$dis$data, cell = "GE"))
DMR_dis_GW_figure <- ggplot(DMR_dis_GW, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("Distance between adjacent DMRs (bp)") + 
  facet_wrap(DM ~ cell) + 
  guides(fill = F) + 
  theme_bw()
ggsave(DMR_dis_GW_figure, file = "DMRdis_GW.pdf", width = 10, height = 10)
(DMR_dis_GW_figure + ggtitle("Distance between adjacent DMRs between 13-week and 17-week"))

DMR_freq_GW <- rbind(cbind(Cortex02_Cortex04_DMR_figures$freq$data, cell = "Cortex"), cbind(GE02_GE04_DMR_figures$freq$data, cell = "GE"))
DMR_freq_GW_figure <- ggplot(DMR_freq_GW, aes(x = chr, y = freq, fill = DM)) + 
  geom_bar(position = "identity", stat = "identity", width = 0.8) + 
  facet_wrap(~ cell) + 
  xlab("") + 
  ylab("DMR frequency (bp/MB)") + 
  coord_flip() + 
  scale_fill_manual(values = c("red", "blue"), labels = c("hyper", "hypo"), name = "") + 
  theme_bw()
ggsave(DMR_freq_GW_figure, file = "DMRfreq_GW.pdf", width = 10, height = 10)
(DMR_freq_GW_figure + ggtitle("DMR frequency asymmetry between 13-week and 17-week"))

DMR_pos_GW <- mutate(rbind(cbind(Cortex02_Cortex04_DMR, cell = "Cortex"), cbind(GE02_GE04_DMR, cell = "GE")), pos = (start + end)/2, y = DM * 15 + 5)
DMR_pos_GW_gr <- keepSeqlevels(as(DMR_pos_GW, "GRanges"), paste0("chr", c(1:22, "X")))
data(hg19IdeogramCyto, package = "biovizBase")
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X")))
(DMR_pos_GW_figure <- ggplot(hg19) + 
   layout_karyogram(cytoband = TRUE) + 
   layout_karyogram(data = DMR_pos_GW_gr, geom = "point", aes(y = y, x = pos, color = factor(DM, levels = c("1", "-1"))), position = position_jitter(height = 2), size = 1, alpha = 0.2) + 
   ylab("") + 
   xlab("Position of DMRs on the chromosome") +
   facet_grid(seqnames ~ cell) + 
   coord_cartesian(ylim = c(-15, 25)) + 
   ggtitle(paste("DMR position - 13-week vs 17-week DMRs")) + 
   scale_color_manual(values = c("red", "blue"), labels = c("hyper", "hypo"), name = "") + 
   theme_bw() + 
   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_line(color = "transparent"), panel.grid.minor.y = element_line(color = "transparent")))
ggsave(DMR_pos_GW_figure@ggplot, file = "DMRpos_GW.pdf", width = 10, height = 10)

GW_UMR_intersect_hyper <- as.integer(system("/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed -b DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed | wc -l", intern = T))
GW_UMR_intersect_hypo <- as.integer(system("/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed -b DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed | wc -l", intern = T))
Venn_GW_UMR_hyper <- draw.pairwise.venn(area1 = GW_DMR_summary$Hyper.DMR[1], area2 = GW_DMR_summary$Hyper.DMR[2], cross.area = GW_UMR_intersect_hyper, category = c("Cortex", "GE"), fill = c("red", "blue"))
Venn_GW_UMR_hypo <- draw.pairwise.venn(area1 = GW_DMR_summary$Hypo.DMR[1], area2 = GW_DMR_summary$Hypo.DMR[2], cross.area = GW_UMR_intersect_hypo, category = c("Cortex", "GE"), fill = c("red", "blue"))
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DMR/Venn_GW_UMR.pdf", height = 5, width = 10)
grid.arrange(gTree(children = Venn_GW_UMR_hyper), gTree(children = Venn_GW_UMR_hypo), nrow = 1)
grid.text("Hyper", x = unit(0.25, "npc"), y = unit(0.9, "npc"))
grid.text("Hypo", x = unit(0.75, "npc"), y = unit(0.9, "npc"))
dev.off()

# ======= GREAT ======== 
setwd("/projects/epigenomics/users/lli/FetalBrain/GW/DMR/")
(GREAT_GW_Cortex02.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.hypo", name = "Cortex-HuFNSC02.UMRs", height = 15))
(GREAT_GW_Cortex04.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.hyper", name = "Cortex-HuFNSC04.UMRs", height = 4))
(GREAT_GW_GE02.UMR <- enrich_GREAT(file = "DMR.GE-HuFNSC02_GE-HuFNSC04.hypo", name = "GE-HuFNSC02.UMRs", height = 12))
(GREAT_GW_GE04.UMR <- enrich_GREAT(file = "DMR.GE-HuFNSC02_GE-HuFNSC04.hyper", name = "GE-HuFNSC04.UMRs", height = 15))
(GREAT_GW_17week.UMR <- enrich_GREAT(file = "DMR.HuFNSC02_HuFNSC04.hypo", name = "17-week.UMRs", height = 15))

# ======= Genomic breakdown ======
setwd("/projects/epigenomics/users/lli/FetalBrain/GW/DMR/CpG/")
genomicBreak_GW <- read.delim("genomic.breakdown.summary", head = F, as.is = T, row.names = 1, col.names = c("Name", "Total", "Intergenic", "Intron", "Exon", "Gene", "Promoter", "CGI"))
genomicBreak_GW_tall <- genomicBreak_GW[2:5, -1]
genomicBreak_GW_tall <- data.frame(Sample = rep(row.names(genomicBreak_GW), ncol(genomicBreak_GW)), Region = factor(rep(colnames(genomicBreak_GW), each = nrow(genomicBreak_GW)), levels = colnames(genomicBreak_GW)), FC = as.vector(as.matrix(genomicBreak_GW)))
genomicBreak_GW_tall$DM <- gsub(".*c3.", "", genomicBreak_GW_tall$Sample)
genomicBreak_GW_tall$Sample <- gsub(".m0.75.*", "", genomicBreak_GW_tall$Sample)
genomicBreak_GW_figure <- ggplot(genomicBreak_GW_tall, aes(x = Region, y = log2(FC), fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5) + 
  xlab("") + 
  ylab("log2 Fold enrichment") + 
  facet_wrap(~ DM) + 
  scale_fill_hue(l = 50, labels = c("Cortex", "GE")) + 
  theme_bw()
ggsave(genomicBreak_GW_figure, file = "genomicBreak_GW.pdf")
(genomicBreak_GW_figure + ggtitle("DMR breakdown between 13-week and 17-week"))

# ======= expression pattern ======
gene_FetalBrain <- read.delim("/home/lli/FetalBrain/RNAseq/rpkm_pc.txt", as.is = T)
colnames(gene_FetalBrain) <- gsub("\\.HuFNSC", "", colnames(gene_FetalBrain))
gene_summary <- data.frame(Sample = colnames(gene_FetalBrain)[-1], 
                          N_0.1 = c(filter(gene_FetalBrain, Brain01 > 0.1) %>% nrow, filter(gene_FetalBrain, Brain02 > 0.1) %>% nrow, filter(gene_FetalBrain, Cortex01 > 0.1) %>% nrow, filter(gene_FetalBrain, Cortex02 > 0.1) %>% nrow, filter(gene_FetalBrain, Cortex03 > 0.1) %>% nrow, filter(gene_FetalBrain, Cortex04 > 0.1) %>% nrow, filter(gene_FetalBrain, GE01 > 0.1) %>% nrow, filter(gene_FetalBrain, GE02 > 0.1) %>% nrow, filter(gene_FetalBrain, GE03 > 0.1) %>% nrow, filter(gene_FetalBrain, GE04 > 0.1) %>% nrow), 
                          N_1 = c(filter(gene_FetalBrain, Brain01 > 1) %>% nrow, filter(gene_FetalBrain, Brain02 > 1) %>% nrow, filter(gene_FetalBrain, Cortex01 > 1) %>% nrow, filter(gene_FetalBrain, Cortex02 > 1) %>% nrow, filter(gene_FetalBrain, Cortex03 > 1) %>% nrow, filter(gene_FetalBrain, Cortex04 > 1) %>% nrow, filter(gene_FetalBrain, GE01 > 1) %>% nrow, filter(gene_FetalBrain, GE02 > 1) %>% nrow, filter(gene_FetalBrain, GE03 > 1) %>% nrow, filter(gene_FetalBrain, GE04 > 1) %>% nrow), 
                          N_10 = c(filter(gene_FetalBrain, Brain01 > 10) %>% nrow, filter(gene_FetalBrain, Brain02 > 10) %>% nrow, filter(gene_FetalBrain, Cortex01 > 10) %>% nrow, filter(gene_FetalBrain, Cortex02 > 10) %>% nrow, filter(gene_FetalBrain, Cortex03 > 10) %>% nrow, filter(gene_FetalBrain, Cortex04 > 10) %>% nrow, filter(gene_FetalBrain, GE01 > 10) %>% nrow, filter(gene_FetalBrain, GE02 > 10) %>% nrow, filter(gene_FetalBrain, GE03 > 10) %>% nrow, filter(gene_FetalBrain, GE04 > 10) %>% nrow), 
                          N_100 = c(filter(gene_FetalBrain, Brain01 > 100) %>% nrow, filter(gene_FetalBrain, Brain02 > 100) %>% nrow, filter(gene_FetalBrain, Cortex01 > 100) %>% nrow, filter(gene_FetalBrain, Cortex02 > 100) %>% nrow, filter(gene_FetalBrain, Cortex03 > 100) %>% nrow, filter(gene_FetalBrain, Cortex04 > 100) %>% nrow, filter(gene_FetalBrain, GE01 > 100) %>% nrow, filter(gene_FetalBrain, GE02 > 100) %>% nrow, filter(gene_FetalBrain, GE03 > 100) %>% nrow, filter(gene_FetalBrain, GE04 > 100) %>% nrow))
rownames(gene_summary) <- colnames(gene_FetalBrain)[-1]
gene_summary_tall <- melt(gene_summary, id = "Sample")
(gene_summary_figure <- ggplot(gene_summary_tall, aes(x = Sample, y = value, fill = variable)) + 
   geom_bar(stat = "identity", position = "dodge") + 
   xlab("") + 
   ylab("No. of genes") + 
   scale_fill_discrete(labels=c("N > 0.1","N > 1","N > 10", "N > 100"), name = "") + 
   theme_bw())
ggsave(gene_summary_figure, file = "/projects/epigenomics/users/lli/FetalBrain/GW/gene_summary_figure.pdf")

exon_FetalBrain <- read.delim("/home/lli/FetalBrain/RNAseq/rpkm_exon.txt", as.is = T)
colnames(exon_FetalBrain) <- gsub("\\.HuFNSC", "", colnames(exon_FetalBrain))
exon_summary <- data.frame(Sample = colnames(exon_FetalBrain)[-1], 
                           N_0.1 = c(filter(exon_FetalBrain, Brain01 > 0.1) %>% nrow, filter(exon_FetalBrain, Brain02 > 0.1) %>% nrow, filter(exon_FetalBrain, Cortex01 > 0.1) %>% nrow, filter(exon_FetalBrain, Cortex02 > 0.1) %>% nrow, filter(exon_FetalBrain, Cortex03 > 0.1) %>% nrow, filter(exon_FetalBrain, Cortex04 > 0.1) %>% nrow, filter(exon_FetalBrain, GE01 > 0.1) %>% nrow, filter(exon_FetalBrain, GE02 > 0.1) %>% nrow, filter(exon_FetalBrain, GE03 > 0.1) %>% nrow, filter(exon_FetalBrain, GE04 > 0.1) %>% nrow), 
                           N_1 = c(filter(exon_FetalBrain, Brain01 > 1) %>% nrow, filter(exon_FetalBrain, Brain02 > 1) %>% nrow, filter(exon_FetalBrain, Cortex01 > 1) %>% nrow, filter(exon_FetalBrain, Cortex02 > 1) %>% nrow, filter(exon_FetalBrain, Cortex03 > 1) %>% nrow, filter(exon_FetalBrain, Cortex04 > 1) %>% nrow, filter(exon_FetalBrain, GE01 > 1) %>% nrow, filter(exon_FetalBrain, GE02 > 1) %>% nrow, filter(exon_FetalBrain, GE03 > 1) %>% nrow, filter(exon_FetalBrain, GE04 > 1) %>% nrow), 
                           N_10 = c(filter(exon_FetalBrain, Brain01 > 10) %>% nrow, filter(exon_FetalBrain, Brain02 > 10) %>% nrow, filter(exon_FetalBrain, Cortex01 > 10) %>% nrow, filter(exon_FetalBrain, Cortex02 > 10) %>% nrow, filter(exon_FetalBrain, Cortex03 > 10) %>% nrow, filter(exon_FetalBrain, Cortex04 > 10) %>% nrow, filter(exon_FetalBrain, GE01 > 10) %>% nrow, filter(exon_FetalBrain, GE02 > 10) %>% nrow, filter(exon_FetalBrain, GE03 > 10) %>% nrow, filter(exon_FetalBrain, GE04 > 10) %>% nrow), 
                           N_100 = c(filter(exon_FetalBrain, Brain01 > 100) %>% nrow, filter(exon_FetalBrain, Brain02 > 100) %>% nrow, filter(exon_FetalBrain, Cortex01 > 100) %>% nrow, filter(exon_FetalBrain, Cortex02 > 100) %>% nrow, filter(exon_FetalBrain, Cortex03 > 100) %>% nrow, filter(exon_FetalBrain, Cortex04 > 100) %>% nrow, filter(exon_FetalBrain, GE01 > 100) %>% nrow, filter(exon_FetalBrain, GE02 > 100) %>% nrow, filter(exon_FetalBrain, GE03 > 100) %>% nrow, filter(exon_FetalBrain, GE04 > 100) %>% nrow))
rownames(exon_summary) <- colnames(exon_FetalBrain)[-1]
exon_summary_tall <- melt(exon_summary, id = "Sample")
(exon_summary_figure <- ggplot(exon_summary_tall, aes(x = Sample, y = value, fill = variable)) + 
   geom_bar(stat = "identity", position = "dodge") + 
   xlab("") + 
   ylab("No. of exons") + 
   scale_fill_discrete(labels=c("N > 0.1","N > 1","N > 10", "N > 100"), name = "") + 
   theme_bw())
ggsave(exon_summary_figure, file = "/projects/epigenomics/users/lli/FetalBrain/GW/exon_summary_figure.pdf")

# ======= DE genes =======
setwd("/home/lli/FetalBrain/RNAseq/DEfine/gene/")
col <- c("ID", "rpkm1", "rpkm2", "p.value", "corrected_P.value")
# cortex
cortex01_cortex03UP <- mutate(read.delim("UP.Cortex-HuFNSC01_Cortex-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
cortex01_cortex04UP <- mutate(read.delim("UP.Cortex-HuFNSC01_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
cortex02_cortex03UP <- mutate(read.delim("UP.Cortex-HuFNSC02_Cortex-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
cortex02_cortex04UP <- mutate(read.delim("UP.Cortex-HuFNSC02_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
cortex03_cortex04UP <- mutate(read.delim("UP.Cortex-HuFNSC03_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
cortex01_cortex03DN <- mutate(read.delim("DN.Cortex-HuFNSC01_Cortex-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
cortex01_cortex04DN <- mutate(read.delim("DN.Cortex-HuFNSC01_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
cortex02_cortex03DN <- mutate(read.delim("DN.Cortex-HuFNSC02_Cortex-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
cortex02_cortex04DN <- mutate(read.delim("DN.Cortex-HuFNSC02_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
cortex03_cortex04DN <- mutate(read.delim("DN.Cortex-HuFNSC03_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
cortex01_cortex03DE <- rbind(cortex01_cortex03UP, cortex01_cortex03DN)
row.names(cortex01_cortex03DE) <- cortex01_cortex03DE$ID
cortex01_cortex04DE <- rbind(cortex01_cortex04UP, cortex01_cortex04DN)
row.names(cortex01_cortex04DE) <- cortex01_cortex04DE$ID
cortex02_cortex03DE <- rbind(cortex02_cortex03UP, cortex02_cortex03DN)
row.names(cortex02_cortex03DE) <- cortex02_cortex03DE$ID
cortex02_cortex04DE <- rbind(cortex02_cortex04UP, cortex02_cortex04DN)
row.names(cortex02_cortex04DE) <- cortex02_cortex04DE$ID
cortex03_cortex04DE <- rbind(cortex03_cortex04UP, cortex03_cortex04DN)
row.names(cortex03_cortex04DE) <- cortex03_cortex04DE$ID
GW_17_13_UP_cortex <- list(cortex01_cortex04 = cortex01_cortex04UP$ID, cortex02_cortex04 = cortex02_cortex04UP$ID)
venn_GW_17_13_UP_cortex <- venn.diagram(GW_17_13_UP_cortex, filename = NULL, fill = c("red", "blue"), main = "GW17 up-regulated genes in cortex")
GW_17_15_UP_cortex <- list(cortex01_cortex03 = cortex01_cortex03UP$ID, cortex02_cortex03 = cortex02_cortex03UP$ID)
venn_GW_17_15_UP_cortex <- venn.diagram(GW_17_15_UP_cortex, filename = NULL, fill = c("red", "blue"), main = "GW17 up-regulated genes in cortex")
GW_17_13_DN_cortex <- list(cortex01_cortex04 = cortex01_cortex04DN$ID, cortex02_cortex04 = cortex02_cortex04DN$ID)
venn_GW_17_13_DN_cortex <- venn.diagram(GW_17_13_DN_cortex, filename = NULL, fill = c("red", "blue"), main = "GW13 up-regulated genes in cortex")
GW_17_15_DN_cortex <- list(cortex01_cortex03 = cortex01_cortex03DN$ID, cortex02_cortex03 = cortex02_cortex03DN$ID)
venn_GW_17_15_DN_cortex <- venn.diagram(GW_17_15_DN_cortex, filename = NULL, fill = c("red", "blue"), main = "GW15 up-regulated genes in cortex")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_13_cortex.pdf")
grid.arrange(gTree(children = venn_GW_17_13_UP_cortex), gTree(children = venn_GW_17_13_DN_cortex), nrow = 1)
dev.off()
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_15_cortex.pdf")
grid.arrange(gTree(children = venn_GW_17_15_UP_cortex), gTree(children = venn_GW_17_15_DN_cortex), nrow = 1)
dev.off()
GW_17_13_UP_duplicated_cortex <- ensembl[intersect(cortex01_cortex04UP$ID, cortex02_cortex04UP$ID), ]
GW_17_13_DN_duplicated_cortex <- ensembl[intersect(cortex01_cortex04DN$ID, cortex02_cortex04DN$ID), ]
write.table(GW_17_13_UP_duplicated_cortex, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_UP_duplicated_cortex.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(GW_17_13_DN_duplicated_cortex, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_DN_duplicated_cortex.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_17_15_UP_duplicated_cortex <- ensembl[intersect(cortex01_cortex03UP$ID, cortex02_cortex03UP$ID), ]
GW_17_15_DN_duplicated_cortex <- ensembl[intersect(cortex01_cortex03DN$ID, cortex02_cortex03DN$ID), ]
write.table(GW_17_15_UP_duplicated_cortex, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_15_UP_duplicated_cortex.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(GW_17_15_DN_duplicated_cortex, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_15_DN_duplicated_cortex.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_cortex <- list(GW17_GW13 = c(GW_17_13_UP_duplicated_cortex$id, GW_17_13_DN_duplicated_cortex$id), GW17_GW15 = c(GW_17_15_UP_duplicated_cortex$id, GW_17_15_DN_duplicated_cortex$id), GW15_GW13 = cortex03_cortex04DE$ID)
venn_GW_cortex <- venn.diagram(GW_cortex, filename = NULL, fill = c("red", "blue", "green"), main = "GW-associated DE genes in cortex")
# GE
GE01_GE03UP <- mutate(read.delim("UP.GE-HuFNSC01_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
GE01_GE04UP <- mutate(read.delim("UP.GE-HuFNSC01_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
GE02_GE03UP <- mutate(read.delim("UP.GE-HuFNSC02_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
GE02_GE04UP <- mutate(read.delim("UP.GE-HuFNSC02_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
GE03_GE04UP <- mutate(read.delim("UP.GE-HuFNSC03_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
GE01_GE03DN <- mutate(read.delim("DN.GE-HuFNSC01_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
GE01_GE04DN <- mutate(read.delim("DN.GE-HuFNSC01_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
GE02_GE03DN <- mutate(read.delim("DN.GE-HuFNSC02_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
GE02_GE04DN <- mutate(read.delim("DN.GE-HuFNSC02_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
GE03_GE04DN <- mutate(read.delim("DN.GE-HuFNSC03_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
GE01_GE03DE <- rbind(GE01_GE03UP, GE01_GE03DN)
row.names(GE01_GE03DE) <- GE01_GE03DE$ID
GE01_GE04DE <- rbind(GE01_GE04UP, GE01_GE04DN)
row.names(GE01_GE04DE) <- GE01_GE04DE$ID
GE02_GE03DE <- rbind(GE02_GE03UP, GE02_GE03DN)
row.names(GE02_GE03DE) <- GE02_GE03DE$ID
GE02_GE04DE <- rbind(GE02_GE04UP, GE02_GE04DN)
row.names(GE02_GE04DE) <- GE02_GE04DE$ID
GE03_GE04DE <- rbind(GE03_GE04UP, GE03_GE04DN)
row.names(GE03_GE04DE) <- GE03_GE04DE$ID
GW_17_13_UP_GE <- list(GE01_GE04 = GE01_GE04UP$ID, GE02_GE04 = GE02_GE04UP$ID)
venn_GW_17_13_UP_GE <- venn.diagram(GW_17_13_UP_GE, filename = NULL, fill = c("red", "blue"), main = "GW17 up-regulated genes in GE")
GW_17_15_UP_GE <- list(GE01_GE03 = GE01_GE03UP$ID, GE02_GE03 = GE02_GE03UP$ID)
venn_GW_17_15_UP_GE <- venn.diagram(GW_17_15_UP_GE, filename = NULL, fill = c("red", "blue"), main = "GW17 up-regulated genes in GE")
GW_17_13_DN_GE <- list(GE01_GE04 = GE01_GE04DN$ID, GE02_GE04 = GE02_GE04DN$ID)
venn_GW_17_13_DN_GE <- venn.diagram(GW_17_13_DN_GE, filename = NULL, fill = c("red", "blue"), main = "GW13 up-regulated genes in GE")
GW_17_15_DN_GE <- list(GE01_GE03 = GE01_GE03DN$ID, GE02_GE03 = GE02_GE03DN$ID)
venn_GW_17_15_DN_GE <- venn.diagram(GW_17_15_DN_GE, filename = NULL, fill = c("red", "blue"), main = "GW15 up-regulated genes in GE")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_13_GE.pdf")
grid.arrange(gTree(children = venn_GW_17_13_UP_GE), gTree(children = venn_GW_17_13_DN_GE), nrow = 1)
dev.off()
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_15_GE.pdf")
grid.arrange(gTree(children = venn_GW_17_15_UP_GE), gTree(children = venn_GW_17_15_DN_GE), nrow = 1)
dev.off()
GW_17_13_UP_duplicated_GE <- ensembl[intersect(GE01_GE04UP$ID, GE02_GE04UP$ID), ]
GW_17_13_DN_duplicated_GE <- ensembl[intersect(GE01_GE04DN$ID, GE02_GE04DN$ID), ]
write.table(GW_17_13_UP_duplicated_GE, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_UP_duplicated_GE.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(GW_17_13_DN_duplicated_GE, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_DN_duplicated_GE.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_17_15_UP_duplicated_GE <- ensembl[intersect(GE01_GE03UP$ID, GE02_GE03UP$ID), ]
GW_17_15_DN_duplicated_GE <- ensembl[intersect(GE01_GE03DN$ID, GE02_GE03DN$ID), ]
write.table(GW_17_15_UP_duplicated_GE, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_15_UP_duplicated_GE.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(GW_17_15_DN_duplicated_GE, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_15_DN_duplicated_GE.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_GE <- list(GW17_GW13 = c(GW_17_13_UP_duplicated_GE$id, GW_17_13_DN_duplicated_GE$id), GW17_GW15 = c(GW_17_15_UP_duplicated_GE$id, GW_17_15_DN_duplicated_GE$id), GW15_GW13 = GE03_GE04DE$ID)
venn_GW_GE <- venn.diagram(GW_GE, filename = NULL, fill = c("red", "blue", "green"), main = "GW-associated DE genes in GE")
# Shared by two cell types
GW_17_13_UP_duplicated <- ensembl[intersect(GW_17_13_UP_duplicated_cortex$id, GW_17_13_UP_duplicated_GE$id), ] %>% filter(chr != "Y")
write.table(GW_17_13_UP_duplicated, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_UP_duplicated.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_17_13_UP <- list(Cortex = GW_17_13_UP_duplicated_cortex$id, GE = GW_17_13_UP_duplicated_GE$id)
venn_GW_17_13_UP <- venn.diagram(GW_17_13_UP, filename = NULL, fill = c("red", "blue"), main = "GW17 up-regulated genes")
GW_17_13_DN_duplicated <- ensembl[intersect(GW_17_13_DN_duplicated_cortex$id, GW_17_13_DN_duplicated_GE$id), ] %>% filter(chr != "Y")
write.table(GW_17_13_DN_duplicated, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_DN_duplicated.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_17_13_DN <- list(Cortex = GW_17_13_DN_duplicated_cortex$id, GE = GW_17_13_DN_duplicated_GE$id)
venn_GW_17_13_DN <- venn.diagram(GW_17_13_DN, filename = NULL, fill = c("red", "blue"), main = "GW13 up-regulated genes")
pcgene <- 19819 # wc -l /projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID
phyper(nrow(GW_17_13_UP_duplicated), nrow(GW_17_13_UP_duplicated_cortex), pcgene - nrow(GW_17_13_UP_duplicated_cortex), nrow(GW_17_13_UP_duplicated_GE), lower.tail = F, log = T) 
phyper(nrow(GW_17_13_DN_duplicated), nrow(GW_17_13_DN_duplicated_cortex), pcgene - nrow(GW_17_13_DN_duplicated_cortex), nrow(GW_17_13_DN_duplicated_GE), lower.tail = F, log = T) 
GW_17_15_UP_duplicated <- ensembl[intersect(GW_17_15_UP_duplicated_cortex$id, GW_17_15_UP_duplicated_GE$id), ] %>% filter(chr != "Y")
write.table(GW_17_15_UP_duplicated, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_15_UP_duplicated.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_17_15_UP <- list(Cortex = GW_17_15_UP_duplicated_cortex$id, GE = GW_17_15_UP_duplicated_GE$id)
venn_GW_17_15_UP <- venn.diagram(GW_17_15_UP, filename = NULL, fill = c("red", "blue"), main = "GW17 up-regulated genes")
GW_17_15_DN_duplicated <- ensembl[intersect(GW_17_15_DN_duplicated_cortex$id, GW_17_15_DN_duplicated_GE$id), ] %>% filter(chr != "Y")
write.table(GW_17_15_DN_duplicated, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_15_DN_duplicated.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_17_15_DN <- list(Cortex = GW_17_15_DN_duplicated_cortex$id, GE = GW_17_15_DN_duplicated_GE$id)
venn_GW_17_15_DN <- venn.diagram(GW_17_15_DN, filename = NULL, fill = c("red", "blue"), main = "GW15 up-regulated genes")
phyper(nrow(GW_17_15_UP_duplicated), nrow(GW_17_15_UP_duplicated_cortex), pcgene - nrow(GW_17_15_UP_duplicated_cortex), nrow(GW_17_15_UP_duplicated_GE), lower.tail = F, log = T) 
phyper(nrow(GW_17_15_DN_duplicated), nrow(GW_17_15_DN_duplicated_cortex), pcgene - nrow(GW_17_15_DN_duplicated_cortex), nrow(GW_17_15_DN_duplicated_GE), lower.tail = F, log = T) 
GW_15_13_UP_duplicated <- ensembl[intersect(cortex03_cortex04UP$ID, GE03_GE04UP$ID), ] %>% filter(chr != "Y")
write.table(GW_15_13_UP_duplicated, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_15_13_UP_duplicated.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_15_13_DN_duplicated <- ensembl[intersect(cortex03_cortex04DN$ID, GE03_GE04DN$ID), ] %>% filter(chr != "Y")
write.table(GW_15_13_DN_duplicated, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_15_13_DN_duplicated.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_15_13_UP <- list(Cortex = cortex03_cortex04UP$ID, GE = GE03_GE04UP$ID)
venn_GW_15_13_UP <- venn.diagram(GW_15_13_UP, filename = NULL, fill = c("red", "blue"), main = "GW15 up-regulated genes")
GW_15_13_DN <- list(Cortex = cortex03_cortex04DN$ID, GE = GE03_GE04DN$ID)
venn_GW_15_13_DN <- venn.diagram(GW_15_13_DN, filename = NULL, fill = c("red", "blue"), main = "GW13 up-regulated genes")
phyper(length(intersect(cortex03_cortex04UP$ID, GE03_GE04UP$ID)), nrow(cortex03_cortex04UP), pcgene - nrow(cortex03_cortex04UP), nrow(GE03_GE04UP), lower.tail = F, log = T) 
phyper(length(intersect(cortex03_cortex04DN$ID, GE03_GE04DN$ID)), nrow(cortex03_cortex04DN), pcgene - nrow(cortex03_cortex04DN), nrow(GE03_GE04DN), lower.tail = F, log = T) 
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_13.pdf")
grid.arrange(gTree(children = venn_GW_17_13_UP), gTree(children = venn_GW_17_13_DN), 
             gTree(children = venn_GW_17_15_UP), gTree(children = venn_GW_17_15_DN), 
             gTree(children = venn_GW_15_13_UP), gTree(children = venn_GW_15_13_DN), nrow = 3)
dev.off()
GW_DE_summary <- data.frame(UP = c(nrow(cortex01_cortex03UP), nrow(cortex01_cortex04UP), nrow(cortex02_cortex03UP), nrow(cortex02_cortex04UP), nrow(cortex03_cortex04UP), nrow(GE01_GE03UP), nrow(GE01_GE04UP), nrow(GE02_GE03UP), nrow(GE02_GE04UP), nrow(GE03_GE04UP), length(intersect(cortex01_cortex03UP$ID, GE01_GE03UP$ID)), length(intersect(cortex01_cortex04UP$ID, GE01_GE04UP$ID)), length(intersect(cortex02_cortex03UP$ID, GE02_GE03UP$ID)), length(intersect(cortex02_cortex04UP$ID, GE02_GE04UP$ID)), length(intersect(cortex03_cortex04UP$ID, GE03_GE04UP$ID))), 
                            DN = c(nrow(cortex01_cortex03DN), nrow(cortex01_cortex04DN), nrow(cortex02_cortex03DN), nrow(cortex02_cortex04DN), nrow(cortex03_cortex04DN), nrow(GE01_GE03DN), nrow(GE01_GE04DN), nrow(GE02_GE03DN), nrow(GE02_GE04DN), nrow(GE03_GE04DN), length(intersect(cortex01_cortex03DN$ID, GE01_GE03DN$ID)), length(intersect(cortex01_cortex04DN$ID, GE01_GE04DN$ID)), length(intersect(cortex02_cortex03DN$ID, GE02_GE03DN$ID)), length(intersect(cortex02_cortex04DN$ID, GE02_GE04DN$ID)), length(intersect(cortex03_cortex04DN$ID, GE03_GE04DN$ID))), 
                            DE = c(nrow(cortex01_cortex03DE), nrow(cortex01_cortex04DE), nrow(cortex02_cortex03DE), nrow(cortex02_cortex04DE), nrow(cortex03_cortex04DE), nrow(GE01_GE03DE), nrow(GE01_GE04DE), nrow(GE02_GE03DE), nrow(GE02_GE04DE), nrow(GE03_GE04DE), length(intersect(cortex01_cortex03DE$ID, GE01_GE03DE$ID)), length(intersect(cortex01_cortex04DE$ID, GE01_GE04DE$ID)), length(intersect(cortex02_cortex03DE$ID, GE02_GE03DE$ID)), length(intersect(cortex02_cortex04DE$ID, GE02_GE04DE$ID)), length(intersect(cortex03_cortex04DE$ID, GE03_GE04DE$ID))),  
                            GW = rep(c("GW17-GW15", "GW17-GW13", "GW17-GW15", "GW17-GW13", "GW15-GW13"), 3), 
                            type = rep(c("cortex", "GE", "shared"), each = 5),
                            Sample = c("GW17-GW15.1_cortex", "GW17-GW13.1_cortex", "GW17-GW15.2_cortex", "GW17-GW13.2_cortex", "GW15-GW13_cortex", "GW17-GW15.1_GE", "GW17-GW13.1_GE", "GW17-GW15.2_GE", "GW17-GW13.2_GE", "GW15-GW13_GE", "GW17-GW15.1_shared", "GW17-GW13.1_shared", "GW17-GW15.2_shared", "GW17-GW13.2_shared", "GW15-GW13_shared"))
GW_DE_summary_tall <- melt(GW_DE_summary %>% select(-DE), id = c("GW", "Sample", "type")) %>% arrange(Sample) %>% mutate(GW = factor(GW, levels = c("GW17-GW13", "GW17-GW15", "GW15-GW13")))
GW_DE_summary_tall[GW_DE_summary_tall$variable == "DN", "value"] <- - GW_DE_summary_tall[GW_DE_summary_tall$variable == "DN", "value"]
(GW_DE_summary_figure <- ggplot(GW_DE_summary_tall, aes(Sample, value, fill = type)) + 
   geom_bar(stat = "identity", position = "identity", width = 0.5) + 
   geom_hline(aes(yintercept = 0)) + 
   facet_grid(. ~ GW, scales = "free_x", space = "free_x") + 
   scale_fill_manual(values = c("red", "blue", "purple"), name = "") + 
   scale_y_continuous(breaks = c(-250, 0, 250, 500), labels = c(250, 0, 250, 500)) + 
   xlab("") + 
   ylab("No. of genes") + 
   theme_bw() + 
   theme(axis.text.x = element_text(angle = 90)))
ggsave(GW_DE_summary_figure, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_DE_summary_figure.pdf")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_DE.pdf")
grid.arrange(gTree(children = venn_GW_cortex), gTree(children = venn_GW_GE), nrow = 1)
dev.off()
# Stage specific DE genes
## DE gene trend profiles: for MZ twins, DE are taken as DE in both individual, and not DE are taken as DE in neither individual. 
## 1) GW13-15:UP, GW15-17:UP, (0,1,2); 2) GW13-15:UP, GW15-17:DN, (0,2,0); 3) GW13-15:DN, GW15-17:UP, (2,0,2); 4) GW13-15:DN, GW15-17:DN, (2,1,0);
## 5) GW13-15:UP, GW15-17:not, (0,2,2); 6) GW13-15:DN, GW15-17:not, (2,0,0); 7) GW13-15:not, GW15-17:UP, (0,0,2); 8) GW13-15:not, GW15-17:DN, (2,2,0); 
### shared by cortex and GE
GW_UP_UP <- ensembl[intersect(GW_15_13_UP_duplicated$id, GW_17_15_UP_duplicated$id), ]
GW_UP_DN <- ensembl[intersect(GW_15_13_UP_duplicated$id, GW_17_15_DN_duplicated$id), ] 
GW_DN_UP <- ensembl[intersect(GW_15_13_DN_duplicated$id, GW_17_15_UP_duplicated$id), ]
GW_DN_DN <- ensembl[intersect(GW_15_13_DN_duplicated$id, GW_17_15_DN_duplicated$id), ]
GW_UP_no <- ensembl[setdiff(GW_15_13_UP_duplicated$id, c(cortex01_cortex03DE$ID, cortex02_cortex03DE$ID, GE01_GE03DE$ID, GE02_GE03DE$ID)), ]
GW_DN_no <- ensembl[setdiff(GW_15_13_DN_duplicated$id, c(cortex01_cortex03DE$ID, cortex02_cortex03DE$ID, GE01_GE03DE$ID, GE02_GE03DE$ID)), ]
GW_no_UP <- ensembl[setdiff(GW_17_15_UP_duplicated$id, c(cortex03_cortex04DE$ID, GE03_GE04DE$ID)), ]
GW_no_DN <- ensembl[setdiff(GW_17_15_DN_duplicated$id, c(cortex03_cortex04DE$ID, GE03_GE04DE$ID)), ]
e <- 0.08
GW_DE_trend <- data.frame(x = rep(c("GW13", "GW15", "GW17"), 8), 
                          y = c(0,1,2, 0,2-e,0, 2,0+e,2, 2,1,0, 0,2,2, 2,0,0, 0,0,2, 2,2,0), 
                          trend = rep(c("UP-UP", "UP-DN", "DN-UP", "DN-DN", "UP-notDE", "DN-notDE", "notDE-UP", "notDE-DN"), each = 3), 
                          N = rep(c(nrow(GW_UP_UP), nrow(GW_UP_DN), nrow(GW_DN_UP), nrow(GW_DN_DN), nrow(GW_UP_no), nrow(GW_DN_no), nrow(GW_no_UP), nrow(GW_no_DN)), each = 3)) %>%
  mutate(type = (N == 0))
(GW_DE_trend_figure <- ggplot(GW_DE_trend, aes(x, y, group = trend, color = trend, size = N, linetype = type)) + 
   geom_line(arrow = arrow(angle = 15, type = "closed")) + 
   xlab("") + 
   ylab("Differential expression") + 
   ggtitle("Shared by Cortex and GE") + 
   scale_linetype_discrete(guide = "none") + 
   theme_bw() + 
   theme(axis.text.y = element_text(size = 0), axis.ticks.y = element_line(size = 0)))
ggsave(GW_DE_trend_figure, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_DE_trend_figure.pdf")
### cortex
GW_UP_UP_cortex <- ensembl[intersect(cortex03_cortex04UP$ID, GW_17_15_UP_duplicated_cortex$id), ]
GW_UP_DN_cortex <- ensembl[intersect(cortex03_cortex04UP$ID, GW_17_15_DN_duplicated_cortex$id), ] 
GW_DN_UP_cortex <- ensembl[intersect(cortex03_cortex04DN$ID, GW_17_15_UP_duplicated_cortex$id), ]
GW_DN_DN_cortex <- ensembl[intersect(cortex03_cortex04DN$ID, GW_17_15_DN_duplicated_cortex$id), ]
GW_UP_no_cortex <- ensembl[setdiff(cortex03_cortex04UP$ID, c(cortex01_cortex03DE$ID, cortex02_cortex03DE$ID)), ]
GW_DN_no_cortex <- ensembl[setdiff(cortex03_cortex04DN$ID, c(cortex01_cortex03DE$ID, cortex02_cortex03DE$ID)), ]
GW_no_UP_cortex <- ensembl[setdiff(GW_17_15_UP_duplicated_cortex$id, c(cortex03_cortex04DE$ID)), ]
GW_no_DN_cortex <- ensembl[setdiff(GW_17_15_DN_duplicated_cortex$id, c(cortex03_cortex04DE$ID)), ]
e <- 0.08
GW_DE_trend_cortex <- data.frame(x = rep(c("GW13", "GW15", "GW17"), 8), 
                          y = c(0,1,2, 0,2-e,0, 2,0+e,2, 2,1,0, 0,2,2, 2,0,0, 0,0,2, 2,2,0), 
                          trend = rep(c("UP-UP", "UP-DN", "DN-UP", "DN-DN", "UP-notDE", "DN-notDE", "notDE-UP", "notDE-DN"), each = 3), 
                          N = rep(c(nrow(GW_UP_UP_cortex), nrow(GW_UP_DN_cortex), nrow(GW_DN_UP_cortex), nrow(GW_DN_DN_cortex), nrow(GW_UP_no_cortex), nrow(GW_DN_no_cortex), nrow(GW_no_UP_cortex), nrow(GW_no_DN_cortex)), each = 3)) 
(GW_DE_trend_cortex_figure <- ggplot(GW_DE_trend_cortex, aes(x, y, group = trend, color = trend, size = N)) + 
   geom_line(arrow = arrow(angle = 15, type = "closed")) + 
   xlab("") + 
   ylab("Differential expression") + 
   ggtitle("Cortex") + 
   theme_bw() + 
   theme(axis.text.y = element_text(size = 0), axis.ticks.y = element_line(size = 0)))
ggsave(GW_DE_trend_cortex_figure, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_DE_trend_cortex_figure.pdf")
### GE
GW_UP_UP_GE <- ensembl[intersect(GE03_GE04UP$ID, GW_17_15_UP_duplicated_GE$id), ]
GW_UP_DN_GE <- ensembl[intersect(GE03_GE04UP$ID, GW_17_15_DN_duplicated_GE$id), ] 
GW_DN_UP_GE <- ensembl[intersect(GE03_GE04DN$ID, GW_17_15_UP_duplicated_GE$id), ]
GW_DN_DN_GE <- ensembl[intersect(GE03_GE04DN$ID, GW_17_15_DN_duplicated_GE$id), ]
GW_UP_no_GE <- ensembl[setdiff(GE03_GE04UP$ID, c(GE01_GE03DE$ID, GE02_GE03DE$ID)), ]
GW_DN_no_GE <- ensembl[setdiff(GE03_GE04DN$ID, c(GE01_GE03DE$ID, GE02_GE03DE$ID)), ]
GW_no_UP_GE <- ensembl[setdiff(GW_17_15_UP_duplicated_GE$id, c(GE03_GE04DE$ID)), ]
GW_no_DN_GE <- ensembl[setdiff(GW_17_15_DN_duplicated_GE$id, c(GE03_GE04DE$ID)), ]
e <- 0.08
GW_DE_trend_GE <- data.frame(x = rep(c("GW13", "GW15", "GW17"), 8), 
                                 y = c(0,1,2, 0,2-e,0, 2,0+e,2, 2,1,0, 0,2,2, 2,0,0, 0,0,2, 2,2,0), 
                                 trend = rep(c("UP-UP", "UP-DN", "DN-UP", "DN-DN", "UP-notDE", "DN-notDE", "notDE-UP", "notDE-DN"), each = 3), 
                                 N = rep(c(nrow(GW_UP_UP_GE), nrow(GW_UP_DN_GE), nrow(GW_DN_UP_GE), nrow(GW_DN_DN_GE), nrow(GW_UP_no_GE), nrow(GW_DN_no_GE), nrow(GW_no_UP_GE), nrow(GW_no_DN_GE)), each = 3)) 
(GW_DE_trend_GE_figure <- ggplot(GW_DE_trend_GE, aes(x, y, group = trend, color = trend, size = N)) + 
   geom_line(arrow = arrow(angle = 15, type = "closed")) + 
   xlab("") + 
   ylab("Differential expression") + 
   ggtitle("GE") + 
   theme_bw() + 
   theme(axis.text.y = element_text(size = 0), axis.ticks.y = element_line(size = 0)))
ggsave(GW_DE_trend_GE_figure, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_DE_trend_GE_figure.pdf")
## heatmap
rownames(gene_FetalBrain) <- gene_FetalBrain$Ensembl
GW_DE_rpkm <- gene_FetalBrain[c(GW_UP_UP_cortex$id, GW_UP_UP_GE$id, GW_UP_DN_cortex$id, GW_UP_DN_GE$id, GW_DN_UP_cortex$id, GW_DN_UP_GE$id, GW_DN_DN_cortex$id, GW_DN_DN_GE$id, 
                                GW_UP_no_cortex$id, GW_UP_no_GE$id, GW_no_UP_cortex$id, GW_no_UP_GE$id, GW_DN_no_cortex$id, GW_DN_no_GE$id, GW_no_DN_cortex$id, GW_no_DN_GE$id),c("Cortex04", "Cortex03", "Cortex02", "Cortex01", "GE04", "GE03", "GE02", "GE01")]
clab = rep(c(rgb(250,192,144,maxColorValue = 255), rgb(247,150,70,maxColorValue = 255), rgb(228,108,10,maxColorValue = 255), rgb(228,108,10,maxColorValue = 255)), 2)
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/heat_GW_DE.pdf")
heatmap.2(as.matrix(GW_DE_rpkm), Rowv = F, Colv = F, labRow = F, ColSideColors = clab, scale = "row", trace = "none", margins = c(11, 6), keysize = 1, density.info = "none", col = bluered(256), key.title = "", key.xlab = "")
legend(x = 0.25, y = 1, c("GW13", "GW15", "GW17"), col = c(rgb(250,192,144,maxColorValue = 255), rgb(247,150,70,maxColorValue = 255), rgb(228,108,10,maxColorValue = 255)), lwd = 8, ncol = 3)
dev.off()
## DAVID enrichment
setwd("/projects/epigenomics/users/lli/FetalBrain/GW/DE/")
(enrich_GW_17_13_UP <- enrich(name = "GW_17_13_UP", fdr = 0.05, p = "FDR", erminej = F, height = 2, width = 9))
(enrich_GW_17_13_DN <- enrich(name = "GW_17_13_DN", fdr = 0.05, p = "FDR", erminej = F, height = 3, width = 9))
(enrich_GW_17_13_UP_cortex <- enrich(name = "GW_17_13_UP_cortex", fdr = 0.05, p = "FDR", erminej = F, height = 2, width = 9))
# calcium ion binding
(enrich_GW_17_13_DN_cortex <- enrich(name = "GW_17_13_DN_cortex", fdr = 0.05, p = "FDR", erminej = F, height = 3, width = 9))
(enrich_GW_17_13_UP_GE <- enrich(name = "GW_17_13_UP_GE", fdr = 0.05, p = "FDR", erminej = F, height = 6, width = 9))
# calcium ion binding, EGF
(enrich_GW_17_13_DN_GE <- enrich(name = "GW_17_13_DN_GE", fdr = 0.05, p = "FDR", erminej = F, height = 6, width = 9))
# neurogenesis
(enrich_GW_17_15_UP <- enrich(name = "GW_17_15_UP", fdr = 0.05, p = "FDR", erminej = F, height = 2, width = 9))
(enrich_GW_17_15_DN <- enrich(name = "GW_17_15_DN", fdr = 0.05, p = "FDR", erminej = F, height = 3, width = 9))
(enrich_GW_17_15_UP_cortex <- enrich(name = "GW_17_15_UP_cortex", fdr = 0.05, p = "FDR", erminej = F, height = 2, width = 9))
(enrich_GW_17_15_DN_cortex <- enrich(name = "GW_17_15_DN_cortex", fdr = 0.05, p = "FDR", erminej = F, height = 3, width = 9))
(enrich_GW_17_15_UP_GE <- enrich(name = "GW_17_15_UP_GE", fdr = 0.05, p = "FDR", erminej = F, height = 8, width = 9))
# neurogenesis, calcium ion binding, EGF
(enrich_GW_17_15_DN_GE <- enrich(name = "GW_17_15_DN_GE", fdr = 0.05, p = "FDR", erminej = F, height = 4, width = 9))
# neurogenesis
(enrich_GW_15_13_UP <- enrich(name = "GW_15_13_UP", fdr = 0.05, p = "FDR", erminej = F, height = 2, width = 9))
(enrich_GW_15_13_DN <- enrich(name = "GW_15_13_DN", fdr = 0.05, p = "FDR", erminej = F, height = 3, width = 9))
(enrich_GW_15_13_UP_cortex <- enrich(name = "GW_15_13_UP_cortex", fdr = 0.05, p = "FDR", erminej = F, height = 2, width = 9))
# calcium ion binding
(enrich_GW_15_13_DN_cortex <- enrich(name = "GW_15_13_DN_cortex", fdr = 0.05, p = "FDR", erminej = F, height = 5, width = 9))
# neurogenesis
(enrich_GW_15_13_UP_GE <- enrich(name = "GW_15_13_UP_GE", fdr = 0.05, p = "FDR", erminej = F, height = 3, width = 9))
# amine catabolic process
(enrich_GW_15_13_DN_GE <- enrich(name = "GW_15_13_DN_GE", fdr = 0.05, p = "FDR", erminej = F, height = 2, width = 9))
# calcium ion binding

# ======= Proximal UMR and DE genes ========
setwd("/projects/epigenomics/users/lli/FetalBrain/GW/DMR/CpG/")
DMR_proximal_Cortex02_Cortex04_hyper <- read.delim("DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_Cortex02_Cortex04_hyper <- DMR_DE(DMR_gene = DMR_proximal_Cortex02_Cortex04_hyper, DE = cortex02_cortex04DE, DM = "hyper", name = "Cortex-HuFNSC02_Cortex-HuFNSC04.hyper.proximal")
DMR_proximal_Cortex02_Cortex04_hypo <- read.delim("DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_Cortex02_Cortex04_hypo <- DMR_DE(DMR_gene = DMR_proximal_Cortex02_Cortex04_hypo, DE = cortex02_cortex04DE, DM = "hypo", name = "Cortex-HuFNSC02_Cortex-HuFNSC04.hypo.proximal")
DMR_proximal_GE02_GE04_hyper <- read.delim("DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_GE02_GE04_hyper <- DMR_DE(DMR_gene = DMR_proximal_GE02_GE04_hyper, DE = GE02_GE04DE, DM = "hyper", name = "GE-HuFNSC02_GE-HuFNSC04.hyper.proximal")
DMR_proximal_GE02_GE04_hypo <- read.delim("DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.CpG_promoter_pc.bed", head = F, as.is = T)
DMR_DE_GE02_GE04_hypo <- DMR_DE(DMR_gene = DMR_proximal_GE02_GE04_hypo, DE = GE02_GE04DE, DM = "hypo", name = "GE-HuFNSC02_GE-HuFNSC04.hypo.proximal")
DMR_proximal_GW_summary <- rbind(DMR_DE_Cortex02_Cortex04_hyper$summary, DMR_DE_Cortex02_Cortex04_hypo$summary, DMR_DE_GE02_GE04_hyper$summary, DMR_DE_GE02_GE04_hypo$summary)
rownames(DMR_proximal_GW_summary) <- c("Cortex04.UMRs", "Cortex02.UMRs", "GE04.UMRs", "GE02.UMRs")

DMR_proximal_GW_hyper <- list(Cortex = DMR_DE_Cortex02_Cortex04_hyper$DMR_gene$id, GE = DMR_DE_GE02_GE04_hyper$DMR_gene$id)
venn_DMR_proximal_GW_hyper <- venn.diagram(DMR_proximal_GW_hyper, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of proximal 13-week UMRs", force.unique = T)
pdf("venn_DMR_proximal_GW_hyper.pdf")
plot.new()
grid.draw(venn_DMR_proximal_GW_hyper)
dev.off()
DMR_proximal_GW_hypo <- list(Cortex = DMR_DE_Cortex02_Cortex04_hypo$DMR_gene$id, GE = DMR_DE_GE02_GE04_hypo$DMR_gene$id)
venn_DMR_proximal_GW_hypo <- venn.diagram(DMR_proximal_GW_hypo, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of proximal 17-week UMRs", force.unique = T, na = "remove")
pdf("venn_DMR_proximal_GW_hypo.pdf")
plot.new()
grid.draw(venn_DMR_proximal_GW_hypo)
dev.off()
pcgene <- 19819 # wc -l /projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID
phyper(length(intersect(DMR_DE_Cortex02_Cortex04_hyper$DMR_gene$id, DMR_DE_GE02_GE04_hyper$DMR_gene$id)), DMR_proximal_GW_summary["Cortex04.UMRs", "unique.genes"], pcgene - DMR_proximal_GW_summary["GE04.UMRs", "unique.genes"], DMR_proximal_GW_summary["GE04.UMRs", "unique.genes"], lower.tail = F, log = T) 
phyper(length(intersect(DMR_DE_Cortex02_Cortex04_hypo$DMR_gene$id, DMR_DE_GE02_GE04_hypo$DMR_gene$id)), DMR_proximal_GW_summary["Cortex02.UMRs", "unique.genes"], pcgene - DMR_proximal_GW_summary["Cortex02.UMRs", "unique.genes"], DMR_proximal_GW_summary["GE02.UMRs", "unique.genes"], lower.tail = F, log = T) 

# =============== Isoforms ====================
setwd("/projects/epigenomics/users/lli/FetalBrain/GW/isoform/")
lib1='A03473'; cell1='Cortex'; donor1='HuFNSC01';
lib2='A04599'; cell2='Cortex'; donor2='HuFNSC03';
cortex01_cortex03_isoform <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirExon = "/home/lli/FetalBrain/RNAseq/DEfine/exon/", dirGene = "/home/lli/FetalBrain/RNAseq/DEfine/gene/")
lib1='A03473'; cell1='Cortex'; donor1='HuFNSC01';
lib2='A15298'; cell2='Cortex'; donor2='HuFNSC04';
cortex01_cortex04_isoform <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirExon = "/home/lli/FetalBrain/RNAseq/DEfine/exon/", dirGene = "/home/lli/FetalBrain/RNAseq/DEfine/gene/")
lib1='A03475'; cell1='Cortex'; donor1='HuFNSC02';
lib2='A04599'; cell2='Cortex'; donor2='HuFNSC03';
cortex02_cortex03_isoform <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirExon = "/home/lli/FetalBrain/RNAseq/DEfine/exon/", dirGene = "/home/lli/FetalBrain/RNAseq/DEfine/gene/")
lib1='A03475'; cell1='Cortex'; donor1='HuFNSC02';
lib2='A15298'; cell2='Cortex'; donor2='HuFNSC04';
cortex02_cortex04_isoform <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirExon = "/home/lli/FetalBrain/RNAseq/DEfine/exon/", dirGene = "/home/lli/FetalBrain/RNAseq/DEfine/gene/")
lib1='A04599'; cell1='Cortex'; donor1='HuFNSC03';
lib2='A15298'; cell2='Cortex'; donor2='HuFNSC04';
cortex03_cortex04_isoform <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirExon = "/home/lli/FetalBrain/RNAseq/DEfine/exon/", dirGene = "/home/lli/FetalBrain/RNAseq/DEfine/gene/")
lib1='A03474'; cell1='GE'; donor1='HuFNSC01';
lib2='A15295'; cell2='GE'; donor2='HuFNSC03';
GE01_GE03_isoform <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirExon = "/home/lli/FetalBrain/RNAseq/DEfine/exon/", dirGene = "/home/lli/FetalBrain/RNAseq/DEfine/gene/")
lib1='A03474'; cell1='GE'; donor1='HuFNSC01';
lib2='A15299'; cell2='GE'; donor2='HuFNSC04';
GE01_GE04_isoform <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirExon = "/home/lli/FetalBrain/RNAseq/DEfine/exon/", dirGene = "/home/lli/FetalBrain/RNAseq/DEfine/gene/")
lib1='A03476'; cell1='GE'; donor1='HuFNSC02';
lib2='A15295'; cell2='GE'; donor2='HuFNSC03';
GE02_GE03_isoform <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirExon = "/home/lli/FetalBrain/RNAseq/DEfine/exon/", dirGene = "/home/lli/FetalBrain/RNAseq/DEfine/gene/")
lib1='A03476'; cell1='GE'; donor1='HuFNSC02';
lib2='A15299'; cell2='GE'; donor2='HuFNSC04';
GE02_GE04_isoform <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirExon = "/home/lli/FetalBrain/RNAseq/DEfine/exon/", dirGene = "/home/lli/FetalBrain/RNAseq/DEfine/gene/")
lib1='A15295'; cell1='GE'; donor1='HuFNSC03';
lib2='A15299'; cell2='GE'; donor2='HuFNSC04';
GE03_GE04_isoform <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirExon = "/home/lli/FetalBrain/RNAseq/DEfine/exon/", dirGene = "/home/lli/FetalBrain/RNAseq/DEfine/gene/")
GW_isoform_summary <- data.frame(rbind(cortex01_cortex03_isoform$summary, cortex01_cortex04_isoform$summary, cortex02_cortex03_isoform$summary, cortex02_cortex04_isoform$summary, cortex03_cortex04_isoform$summary, 
                                       GE01_GE03_isoform$summary, GE01_GE04_isoform$summary, GE02_GE03_isoform$summary, GE02_GE04_isoform$summary, GE03_GE04_isoform$summary)) %>% 
  mutate(GW = rep(c("17 vs 15", "17 vs 13", "17 vs 15", "17 vs 13", "15 vs 13"), 2))
# Venn diagram 
GW_isoform_cortex <- list(cortex01_cortex03 = cortex01_cortex03_isoform$isoform_gene$id, cortex01_cortex04 = cortex01_cortex04_isoform$isoform_gene$id, cortex02_cortex03 = cortex02_cortex03_isoform$isoform_gene$id, cortex02_cortex04 = cortex02_cortex04_isoform$isoform_gene$id)
venn_GW_isoform_cortex <- venn.diagram(GW_isoform_cortex, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of GW isoform genes in cortex")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_isoform_cortex.pdf")
plot.new()
grid.draw(venn_GW_isoform_cortex)
dev.off()
GW_isoform_gene_cortex <- c(cortex01_cortex03_isoform$isoform_gene$id, cortex01_cortex04_isoform$isoform_gene$id, cortex02_cortex03_isoform$isoform_gene$id, cortex02_cortex04_isoform$isoform_gene$id)
GW_isoform_gene_dup_cortex <- ensembl[unique(GW_isoform_gene_cortex[duplicated(GW_isoform_gene_cortex)]), ] 
write.table(GW_isoform_gene_dup_cortex, file = "GW_isoform_gene_dup_cortex.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_isoform_GE <- list(GE01_GE03 = GE01_GE03_isoform$isoform_gene$id, GE01_GE04 = GE01_GE04_isoform$isoform_gene$id, GE02_GE03 = GE02_GE03_isoform$isoform_gene$id, GE02_GE04 = GE02_GE04_isoform$isoform_gene$id)
venn_GW_isoform_GE <- venn.diagram(GW_isoform_GE, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of GW isoform genes in GE")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_isoform_GE.pdf")
plot.new()
grid.draw(venn_GW_isoform_GE)
dev.off()
GW_isoform_gene_GE <- c(GE01_GE03_isoform$isoform_gene$id, GE01_GE04_isoform$isoform_gene$id, GE02_GE03_isoform$isoform_gene$id, GE02_GE04_isoform$isoform_gene$id)
GW_isoform_gene_dup_GE <- ensembl[unique(GW_isoform_gene_GE[duplicated(GW_isoform_gene_GE)]), ] 
write.table(GW_isoform_gene_dup_GE, file = "GW_isoform_gene_dup_GE.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_isoform <- list(cortex = GW_isoform_gene_dup_cortex$id, GE = GW_isoform_gene_dup_GE$id)
venn_GW_isoform <- venn.diagram(GW_isoform, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of GW isoform genes")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_isoform.pdf")
plot.new()
grid.draw(venn_GW_isoform)
dev.off()
GW_isoform_gene_dup <- ensembl[intersect(GW_isoform_gene_dup_cortex$id, GW_isoform_gene_dup_GE$id), ] 
write.table(GW_isoform_gene_dup, file = "GW_isoform_gene_dup.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# DAVID enrichment
enrich_GW_isoform_gene_dup <- enrich(name = "GW_isoform_gene_dup", fdr = 0.01, p = "FDR", erminej = F, height = 8, width = 9)
(enrich_GW_isoform_gene_dup)

#' ================ UMR overlap with TFBSs =====================
setwd("/projects/epigenomics/users/lli/FetalBrain/GW/DMR/TF")
DMR_GW_TF <- read.delim("DMR.HuFNSC02_HuFNSC04.m0.75.p0.005.d0.5.s300.c3.TF.summary", head = F, as.is = T, col.names = c("TF", "UMR17week_Cortex", "UMR13week_Cortex", "Ratio_Cortex", "UMR17week_GE", "UMR13week_GE", "Ratio_GE"))
DMR_GW_TF_tall <- data.frame(TF = rep(DMR_GW_TF$TF, 2), Cell = rep(c("Cortex", "GE"), each = nrow(DMR_GW_TF)), UMR17week = c(DMR_GW_TF$UMR17week_Cortex, DMR_GW_TF$UMR17week_GE), UMR13week = c(DMR_GW_TF$UMR13week_Cortex, DMR_GW_TF$UMR13week_GE), Ratio = c(DMR_GW_TF$Ratio_Cortex, DMR_GW_TF$Ratio_GE))
DMR_GW_TF_tall$Asymmetry <- NA
DMR_GW_TF_tall[DMR_GW_TF_tall$Ratio > 1,]$Asymmetry <- "17-week enriched"
DMR_GW_TF_tall[DMR_GW_TF_tall$Ratio < 1,]$Asymmetry <- "13-week enriched"
DMR_GW_TF_tall[DMR_GW_TF_tall$Ratio == 1,]$Asymmetry <- "Equal"
(DMR_GW_TF_figure <- ggplot(DMR_GW_TF_tall, aes(x = UMR17week, y = UMR13week, label = TF, color = Asymmetry)) + 
   geom_point() + 
   geom_abline(intercept=0, slope=1) + 
   geom_text(data = subset(DMR_GW_TF_tall, (Ratio > 40 | Ratio <= 1/40)), angle = 45, size = 4, hjust = 0.2, vjust = 0.2) + 
   facet_wrap(~ Cell) + 
   scale_x_log10() + 
   scale_y_log10() + 
   scale_color_hue(l = 50) + 
   theme_bw())
ggsave(DMR_GW_TF_figure, file = "DMR_GW_TF.pdf", height = 6)

# ================= chrEnd_GW enrichment ====================
setwd("/projects/epigenomics/users/lli/FetalBrain/GW/DMR/")
chrEnd_GW_Cortex17UMR <- read.delim("DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed.chrom.window", head = F, as.is = T)
chrEnd_GW_Cortex13UMR <- read.delim("DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed.chrom.window", head = F, as.is = T)
chrEnd_GW_GE17UMR <- read.delim("DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed.chrom.window", head = F, as.is = T)
chrEnd_GW_GE13UMR <- read.delim("DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed.chrom.window", head = F, as.is = T)
chrEnd_GW <- rbind(data.frame(DM = "17-week.UMRs", Cell = "Cortex", chrEnd_GW_Cortex17UMR %>% group_by(V4) %>% summarise(N = sum(V5))), 
                   data.frame(DM = "13-week.UMRs", Cell = "Cortex", chrEnd_GW_Cortex13UMR %>% group_by(V4) %>% summarise(N = sum(V5))), 
                   data.frame(DM = "17-week.UMRs", Cell = "GE", chrEnd_GW_GE17UMR %>% group_by(V4) %>% summarise(N = sum(V5))), 
                   data.frame(DM = "13-week.UMRs", Cell = "GE", chrEnd_GW_GE13UMR %>% group_by(V4) %>% summarise(N = sum(V5))))
(chrEnd_GW_figure <- ggplot(chrEnd_GW, aes(x = Cell, y = N, fill = as.factor(V4))) + 
   geom_bar(stat = "identity", width = 0.5, position = "fill") + 
   facet_grid(DM ~ .) + 
   xlab("") + 
   ylab("Fraction of total UMRs") + 
   coord_flip() + 
   scale_fill_discrete(name = "Normalized length \nalong chromosome", labels = as.character(seq(0.1, 1, by = 0.1))) + 
   theme_bw() + 
   theme(strip.text = element_text(size = 15), axis.title.y = element_text(size = 15), axis.text = element_text(size = 12), legend.text = element_text(size = 15)))
ggsave(chrEnd_GW_figure, file = "chrEnd_GW.pdf")

# =============== CGI coverage =================
setwd("/projects/epigenomics/users/lli/FetalBrain/WGBS/CGI/")
Cortex02_CGI <- read.delim("A22475.WGBS.NeurospheresCortex02.sam.bedGraph.CGI.bed", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "m", "cov")) %>% mutate(sample = "Cortex02", cell = "Cortex", GW = "17-week")
GE02_CGI <- read.delim("A17784-A13819.WGBS.NeurospheresGE02.sam.bedGraph.CGI.bed", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "m", "cov")) %>% mutate(sample = "GE02", cell = "GE", GW = "17-week")
Cortex04_CGI <- read.delim("A22477.WGBS.NeurospheresCortex04.sam.bedGraph.CGI.bed", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "m", "cov")) %>% mutate(sample = "Cortex04", cell = "Cortex", GW = "13-week")
GE04_CGI <- read.delim("A22476.WGBS.NeurospheresGE04.sam.bedGraph.CGI.bed", head = F, as.is = T, col.names = c("chr", "start", "end", "ID", "m", "cov")) %>% mutate(sample = "GE04", cell = "GE", GW = "13-week")
CGI_WGBS <- rbind(Cortex02_CGI, GE02_CGI, Cortex04_CGI, GE04_CGI)
(CGI_coverage_figure <- ggplot(CGI_WGBS, aes(x = sample, y = cov, fill = GW)) + 
   geom_boxplot(outlier.shape = NA) + 
   xlab("") + 
   ylab("Average coverage per CpG in CGIs") + 
   coord_cartesian(ylim = c(0, 100)) + 
   scale_fill_hue(l = 50) + 
   theme_bw())
ggsave(CGI_coverage_figure, file = "CGI_coverage_figure.pdf")

save(GW_DMR_summary, Cortex02_Cortex04_DMR, GE02_GE04_DMR, Venn_GW_UMR_hyper, Venn_GW_UMR_hypo, 
     DMR_length_GW_figure, DMR_count_GW_figure, DMR_dis_GW_figure, DMR_freq_GW_figure, DMR_pos_GW_figure, 
     GREAT_GW_Cortex02.UMR, GREAT_GW_Cortex04.UMR, GREAT_GW_GE02.UMR, GREAT_GW_GE04.UMR, GREAT_GW_17week.UMR,  
     genomicBreak_GW, genomicBreak_GW_figure, gene_summary_figure, exon_summary_figure, 
     cortex01_cortex03DE, cortex01_cortex04DE, cortex02_cortex03DE, cortex02_cortex04DE, cortex03_cortex04DE, 
     GW_17_13_UP_duplicated_cortex, GW_17_13_DN_duplicated_cortex, GW_17_15_UP_duplicated_cortex, GW_17_15_DN_duplicated_cortex, 
     GE01_GE03DE, GE01_GE04DE, GE02_GE03DE, GE02_GE04DE, GE03_GE04DE, GW_DE_summary, GW_DE_summary_figure, 
     GW_17_13_UP_duplicated_GE, GW_17_13_DN_duplicated_GE, GW_17_15_UP_duplicated_GE, GW_17_15_DN_duplicated_GE, 
     GW_17_13_UP_duplicated, GW_17_13_DN_duplicated, GW_17_15_UP_duplicated, GW_17_15_DN_duplicated, GW_15_13_UP_duplicated, GW_15_13_DN_duplicated, 
     venn_GW_17_13_UP, venn_GW_17_13_DN, venn_GW_17_15_UP, venn_GW_17_15_DN, venn_GW_15_13_UP, venn_GW_15_13_DN, venn_GW_cortex, venn_GW_GE, GW_DE_rpkm, 
     GW_UP_UP, GW_UP_DN, GW_DN_UP, GW_DN_DN, GW_UP_no, GW_DN_no, GW_no_UP, GW_no_DN, GW_DE_trend, GW_DE_trend_figure, 
     GW_UP_UP_cortex, GW_UP_DN_cortex, GW_DN_UP_cortex, GW_DN_DN_cortex, GW_UP_no_cortex, GW_DN_no_cortex, GW_no_UP_cortex, GW_no_DN_cortex, GW_DE_trend_cortex, GW_DE_trend_cortex_figure, 
     GW_UP_UP_GE, GW_UP_DN_GE, GW_DN_UP_GE, GW_DN_DN_GE, GW_UP_no_GE, GW_DN_no_GE, GW_no_UP_GE, GW_no_DN_GE, GW_DE_trend_GE, GW_DE_trend_GE_figure, 
     enrich_GW_17_13_UP, enrich_GW_17_13_DN, enrich_GW_17_13_UP_cortex, enrich_GW_17_13_DN_cortex, enrich_GW_17_13_UP_GE, enrich_GW_17_13_DN_GE, 
     enrich_GW_17_15_UP, enrich_GW_17_15_DN, enrich_GW_17_15_UP_cortex, enrich_GW_17_15_DN_cortex, enrich_GW_17_15_UP_GE, enrich_GW_17_15_DN_GE, 
     enrich_GW_15_13_UP, enrich_GW_15_13_DN, enrich_GW_15_13_UP_cortex, enrich_GW_15_13_DN_cortex, enrich_GW_15_13_UP_GE, enrich_GW_15_13_DN_GE, 
     DMR_DE_Cortex02_Cortex04_hyper, DMR_DE_Cortex02_Cortex04_hypo, DMR_DE_GE02_GE04_hyper, DMR_DE_GE02_GE04_hypo, 
     DMR_proximal_GW_summary, venn_DMR_proximal_GW_hyper, venn_DMR_proximal_GW_hypo, 
     cortex01_cortex03_isoform, cortex01_cortex04_isoform, cortex02_cortex03_isoform, cortex02_cortex04_isoform, cortex03_cortex04_isoform, 
     GE01_GE03_isoform, GE01_GE04_isoform, GE02_GE03_isoform, GE02_GE04_isoform, GE03_GE04_isoform, GW_isoform_summary, enrich_GW_isoform_gene_dup, 
     venn_GW_isoform_cortex, GW_isoform_gene_dup_cortex, venn_GW_isoform_GE, GW_isoform_gene_dup_GE, venn_GW_isoform, GW_isoform_gene_dup, 
     DMR_GW_TF, DMR_GW_TF_figure, chrEnd_GW_figure, CGI_coverage_figure, 
     file = "/projects/epigenomics/users/lli/FetalBrain/GW/GW.Rdata")

