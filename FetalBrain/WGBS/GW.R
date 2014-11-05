# DMR analysis between 13-week and 17-week 
# DMR identification: ~/HirstLab/FetalBrain/WGBS/WGBS.DM.sh

setwd("/projects/epigenomics/users/lli/FetalBrain/GW/DMR/")
library(ggplot2)
library(dplyr)
library(ggbio)
library(GenomicRanges)
library(labeling)
library(VennDiagram)
source('~/HirstLab/Pipeline/R/enrich.R')
source("~/HirstLab/Pipeline/R/DMR.figures.R")
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
source('~/HirstLab/Pipeline/R/DMR_DE.R')
load("~/FetalBrain/RNAseq/DEfine/gene/FetalBrain_DEgenes.Rdata")

# ============ Sanity =========
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

# ======= GREAT ======== 
(GREAT_GW_Cortex02.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.hypo", name = "Cortex-HuFNSC02.UMRs", height = 15))
(GREAT_GW_Cortex04.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.hyper", name = "Cortex-HuFNSC04.UMRs", height = 4))
(GREAT_GW_GE02.UMR <- enrich_GREAT(file = "DMR.GE-HuFNSC02_GE-HuFNSC04.hypo", name = "GE-HuFNSC02.UMRs", height = 12))
(GREAT_GW_GE04.UMR <- enrich_GREAT(file = "DMR.GE-HuFNSC02_GE-HuFNSC04.hyper", name = "GE-HuFNSC04.UMRs", height = 15))

# ======= Genomic breakdown ======
setwd("/projects/epigenomics/users/lli/FetalBrain/GW/DMR/CpG/")
genomicBreak_GW <- read.delim("genomic.breakdown.summary", head = F, as.is = T, row.names = 1, col.names = c("Name", "Total", "Intergenic", "Intron", "Exon", "Gene", "Promoter", "CGI"))
genomicBreak_GW_tall <- data.frame(Sample = rep(row.names(genomicBreak_GW), ncol(genomicBreak_GW)), Region = factor(rep(colnames(genomicBreak_GW), each = nrow(genomicBreak_GW)), levels = colnames(genomicBreak_GW)), CpG = as.vector(as.matrix(genomicBreak_GW)))
genomicBreak_GW_tall$DM <- gsub(".*c3.", "", genomicBreak_GW_tall$Sample)
genomicBreak_GW_tall$Sample <- gsub(".m0.75.*", "", genomicBreak_GW_tall$Sample)
genomicBreak_GW_tall <- genomicBreak_GW_tall[genomicBreak_GW_tall$Region != "Total", ]
genomicBreak_GW_figure <- ggplot(genomicBreak_GW_tall, aes(x = Region, y = CpG, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5) + 
  xlab("") + 
  ylab("Fraction of CpG") + 
  facet_wrap(~ DM) + 
  coord_flip() + 
  scale_fill_hue(l = 50, labels = c("Cortex", "GE")) + 
  theme_bw()
ggsave(genomicBreak_GW_figure, file = "genomicBreak_GW.pdf")
(genomicBreak_GW_figure + ggtitle("DMR breakdown between 13-week and 17-week"))

save(GW_DMR_summary, Cortex02_Cortex04_DMR, GE02_GE04_DMR, 
     DMR_length_GW_figure, DMR_count_GW_figure, DMR_dis_GW_figure, DMR_freq_GW_figure, DMR_pos_GW_figure, 
     GREAT_GW_Cortex02.UMR, GREAT_GW_Cortex04.UMR, GREAT_GW_GE02.UMR, GREAT_GW_GE04.UMR, 
     genomicBreak_GW, genomicBreak_GW_figure, 
     file = "/projects/epigenomics/users/lli/FetalBrain/GW/GW.Rdata")

