# DMR analysis between 13-week and 17-week 
# DMR identification: ~/HirstLab/FetalBrain/WGBS/WGBS.DM.sh

library(ggplot2)
library(dplyr)
library(ggbio)
library(GenomicRanges)
library(labeling)
library(VennDiagram)
library(gridExtra)
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

# ======= DE genes =======
setwd("/home/lli/FetalBrain/RNAseq/DEfine/gene/")
col <- c("ID", "rpkm1", "rpkm2", "p.value", "corrected_P.value")
# cortex
cortex01_cortex03UP <- mutate(read.delim("UP.Cortex-HuFNSC01_Cortex-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
cortex01_cortex04UP <- mutate(read.delim("UP.Cortex-HuFNSC01_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
cortex02_cortex03UP <- mutate(read.delim("UP.Cortex-HuFNSC02_Cortex-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
cortex02_cortex04UP <- mutate(read.delim("UP.Cortex-HuFNSC02_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
cortex01_cortex03DN <- mutate(read.delim("DN.Cortex-HuFNSC01_Cortex-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
cortex01_cortex04DN <- mutate(read.delim("DN.Cortex-HuFNSC01_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
cortex02_cortex03DN <- mutate(read.delim("DN.Cortex-HuFNSC02_Cortex-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
cortex02_cortex04DN <- mutate(read.delim("DN.Cortex-HuFNSC02_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
cortex01_cortex03DE <- rbind(cortex01_cortex03UP, cortex01_cortex03DN)
row.names(cortex01_cortex03DE) <- cortex01_cortex03DE$ID
cortex01_cortex04DE <- rbind(cortex01_cortex04UP, cortex01_cortex04DN)
row.names(cortex01_cortex04DE) <- cortex01_cortex04DE$ID
cortex02_cortex03DE <- rbind(cortex02_cortex03UP, cortex02_cortex03DN)
row.names(cortex02_cortex03DE) <- cortex02_cortex03DE$ID
cortex02_cortex04DE <- rbind(cortex02_cortex04UP, cortex02_cortex04DN)
row.names(cortex02_cortex04DE) <- cortex02_cortex04DE$ID
GW_17_13_UP_cortex <- list(cortex01_cortex03 = cortex01_cortex03UP$ID, cortex01_cortex04 = cortex01_cortex04UP$ID, cortex02_cortex03 = cortex02_cortex03UP$ID, cortex02_cortex04 = cortex02_cortex04UP$ID)
venn_GW_17_13_UP_cortex <- venn.diagram(GW_17_13_UP_cortex, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of 17-week up-regulated genes in cortex")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_13_UP_cortex.pdf")
plot.new()
grid.draw(venn_GW_17_13_UP_cortex)
dev.off()
GW_17_13_DN_cortex <- list(cortex01_cortex03 = cortex01_cortex03DN$ID, cortex01_cortex04 = cortex01_cortex04DN$ID, cortex02_cortex03 = cortex02_cortex03DN$ID, cortex02_cortex04 = cortex02_cortex04DN$ID)
venn_GW_17_13_DN_cortex <- venn.diagram(GW_17_13_DN_cortex, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of 13-week up-regulated genes in cortex")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_13_DN_cortex.pdf")
plot.new()
grid.draw(venn_GW_17_13_DN_cortex)
dev.off()
GW_17_13_UP_duplicated_cortex <- c(cortex01_cortex03UP$ID, cortex01_cortex04UP$ID, cortex02_cortex03UP$ID, cortex02_cortex04UP$ID)
GW_17_13_UP_duplicated_cortex <- ensembl[unique(GW_17_13_UP_duplicated_cortex[duplicated(GW_17_13_UP_duplicated_cortex)]), ]
GW_17_13_DN_duplicated_cortex <- c(cortex01_cortex03DN$ID, cortex01_cortex04DN$ID, cortex02_cortex03DN$ID, cortex02_cortex04DN$ID)
GW_17_13_DN_duplicated_cortex <- ensembl[unique(GW_17_13_DN_duplicated_cortex[duplicated(GW_17_13_DN_duplicated_cortex)]), ]
write.table(GW_17_13_UP_duplicated_cortex, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_UP_duplicated_cortex.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(GW_17_13_DN_duplicated_cortex, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_DN_duplicated_cortex.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# GE
GE01_GE03UP <- mutate(read.delim("UP.GE-HuFNSC01_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
GE01_GE04UP <- mutate(read.delim("UP.GE-HuFNSC01_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
GE02_GE03UP <- mutate(read.delim("UP.GE-HuFNSC02_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
GE02_GE04UP <- mutate(read.delim("UP.GE-HuFNSC02_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "UP")
GE01_GE03DN <- mutate(read.delim("DN.GE-HuFNSC01_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
GE01_GE04DN <- mutate(read.delim("DN.GE-HuFNSC01_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
GE02_GE03DN <- mutate(read.delim("DN.GE-HuFNSC02_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
GE02_GE04DN <- mutate(read.delim("DN.GE-HuFNSC02_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T, col.names = col), DE = "DN")
GE01_GE03DE <- rbind(GE01_GE03UP, GE01_GE03DN)
row.names(GE01_GE03DE) <- GE01_GE03DE$ID
GE01_GE04DE <- rbind(GE01_GE04UP, GE01_GE04DN)
row.names(GE01_GE04DE) <- GE01_GE04DE$ID
GE02_GE03DE <- rbind(GE02_GE03UP, GE02_GE03DN)
row.names(GE02_GE03DE) <- GE02_GE03DE$ID
GE02_GE04DE <- rbind(GE02_GE04UP, GE02_GE04DN)
row.names(GE02_GE04DE) <- GE02_GE04DE$ID
GW_17_13_UP_GE <- list(GE01_GE03 = GE01_GE03UP$ID, GE01_GE04 = GE01_GE04UP$ID, GE02_GE03 = GE02_GE03UP$ID, GE02_GE04 = GE02_GE04UP$ID)
venn_GW_17_13_UP_GE <- venn.diagram(GW_17_13_UP_GE, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of 17-week up-regulated genes in GE")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_13_UP_GE.pdf")
plot.new()
grid.draw(venn_GW_17_13_UP_GE)
dev.off()
GW_17_13_DN_GE <- list(GE01_GE03 = GE01_GE03DN$ID, GE01_GE04 = GE01_GE04DN$ID, GE02_GE03 = GE02_GE03DN$ID, GE02_GE04 = GE02_GE04DN$ID)
venn_GW_17_13_DN_GE <- venn.diagram(GW_17_13_DN_GE, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of 13-week up-regulated genes in GE")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_13_DN_GE.pdf")
plot.new()
grid.draw(venn_GW_17_13_DN_GE)
dev.off()
GW_17_13_UP_duplicated_GE <- c(GE01_GE03UP$ID, GE01_GE04UP$ID, GE02_GE03UP$ID, GE02_GE04UP$ID)
GW_17_13_UP_duplicated_GE <- ensembl[unique(GW_17_13_UP_duplicated_GE[duplicated(GW_17_13_UP_duplicated_GE)]), ]
write.table(GW_17_13_UP_duplicated_GE, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_UP_duplicated_GE.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_17_13_DN_duplicated_GE <- c(GE01_GE03DN$ID, GE01_GE04DN$ID, GE02_GE03DN$ID, GE02_GE04DN$ID)
GW_17_13_DN_duplicated_GE <- ensembl[unique(GW_17_13_DN_duplicated_GE[duplicated(GW_17_13_DN_duplicated_GE)]), ]
write.table(GW_17_13_UP_duplicated_GE, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_UP_duplicated_GE.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(GW_17_13_DN_duplicated_GE, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_DN_duplicated_GE.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# Shared by two cell types
GW_17_13_UP_duplicated <- ensembl[intersect(GW_17_13_UP_duplicated_cortex$id, GW_17_13_UP_duplicated_GE$id), ]
write.table(GW_17_13_UP_duplicated, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_UP_duplicated.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_17_13_UP <- list(Cortex = GW_17_13_UP_duplicated_cortex$id, GE = GW_17_13_UP_duplicated_GE$id)
venn_GW_17_13_UP <- venn.diagram(GW_17_13_UP, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of 17-week up-regulated genes")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_13_UP.pdf")
plot.new()
grid.draw(venn_GW_17_13_UP)
dev.off()
GW_17_13_DN_duplicated <- ensembl[intersect(GW_17_13_DN_duplicated_cortex$id, GW_17_13_DN_duplicated_GE$id), ]
write.table(GW_17_13_DN_duplicated, file = "/projects/epigenomics/users/lli/FetalBrain/GW/DE/GW_17_13_DN_duplicated.txt", sep = "\t", quote = F, row.names = F, col.names = F)
GW_17_13_DN <- list(Cortex = GW_17_13_DN_duplicated_cortex$id, GE = GW_17_13_DN_duplicated_GE$id)
venn_GW_17_13_DN <- venn.diagram(GW_17_13_DN, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of 17-week DN-regulated genes")
pdf("/projects/epigenomics/users/lli/FetalBrain/GW/DE/venn_GW_17_13_DN.pdf")
plot.new()
grid.draw(venn_GW_17_13_DN)
dev.off()
pcgene <- 19819 # wc -l /projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes.pc.EnsID
phyper(nrow(GW_17_13_UP_duplicated), nrow(GW_17_13_UP_duplicated_cortex), pcgene - nrow(GW_17_13_UP_duplicated_cortex), nrow(GW_17_13_UP_duplicated_GE), lower.tail = F, log = T) 
phyper(nrow(GW_17_13_DN_duplicated), nrow(GW_17_13_DN_duplicated_cortex), pcgene - nrow(GW_17_13_DN_duplicated_cortex), nrow(GW_17_13_DN_duplicated_GE), lower.tail = F, log = T) 
GW_DE_summary <- data.frame(UP = c(nrow(cortex01_cortex03UP), nrow(cortex01_cortex04UP), nrow(cortex02_cortex03UP), nrow(cortex02_cortex04UP), nrow(GE01_GE03UP), nrow(GE01_GE04UP), nrow(GE02_GE03UP), nrow(GE02_GE04UP)), 
                            DN = c(nrow(cortex01_cortex03DN), nrow(cortex01_cortex04DN), nrow(cortex02_cortex03DN), nrow(cortex02_cortex04DN), nrow(GE01_GE03DN), nrow(GE01_GE04DN), nrow(GE02_GE03DN), nrow(GE02_GE04DN)), 
                            DE = c(nrow(cortex01_cortex03DE), nrow(cortex01_cortex04DE), nrow(cortex02_cortex03DE), nrow(cortex02_cortex04DE), nrow(GE01_GE03DE), nrow(GE01_GE04DE), nrow(GE02_GE03DE), nrow(GE02_GE04DE)))
rownames(GW_DE_summary) <- c("cortex01_cortex03", "cortex01_cortex04", "cortex02_cortex03", "cortex02_cortex04", "GE01_GE03", "GE01_GE04", "GE02_GE03", "GE02_GE04")
# DAVID enrichment
setwd("/projects/epigenomics/users/lli/FetalBrain/GW/DE/")
enrich_GW_17_13_UP <- enrich(name = "GW_17_13_UP_duplicated", fdr = 0.05, p = "FDR", erminej = F, height = 2, width = 9)
(enrich_GW_17_13_UP)
enrich_GW_17_13_DN <- enrich(name = "GW_17_13_DN_duplicated", fdr = 0.05, p = "FDR", erminej = F, height = 3, width = 9)
(enrich_GW_17_13_DN)
enrich_GW_17_13_UP_cortex <- enrich(name = "GW_17_13_UP_duplicated_cortex", fdr = 0.05, p = "FDR", erminej = F, height = 3, width = 9)
(enrich_GW_17_13_UP_cortex)
enrich_GW_17_13_DN_cortex <- enrich(name = "GW_17_13_DN_duplicated_cortex", fdr = 0.05, p = "FDR", erminej = F, height = 3, width = 9)
(enrich_GW_17_13_DN_cortex)
enrich_GW_17_13_UP_GE <- enrich(name = "GW_17_13_UP_duplicated_GE", fdr = 0.05, p = "FDR", erminej = F, height = 8, width = 9)
(enrich_GW_17_13_UP_GE)
enrich_GW_17_13_DN_GE <- enrich(name = "GW_17_13_DN_duplicated_GE", fdr = 0.05, p = "FDR", erminej = F, height = 5, width = 9)
(enrich_GW_17_13_DN_GE)

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
GW_isoform_summary <- rbind(cortex01_cortex03_isoform$summary, cortex01_cortex04_isoform$summary, cortex02_cortex03_isoform$summary, cortex02_cortex04_isoform$summary, 
                            GE01_GE03_isoform$summary, GE01_GE04_isoform$summary, GE02_GE03_isoform$summary, GE02_GE04_isoform$summary)
rownames(GW_isoform_summary) <- c("cortex01_cortex03", "cortex01_cortex04", "cortex02_cortex03", "cortex02_cortex04", "GE01_GE03", "GE01_GE04", "GE02_GE03", "GE02_GE04")
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

save(GW_DMR_summary, Cortex02_Cortex04_DMR, GE02_GE04_DMR, Venn_GW_UMR_hyper, Venn_GW_UMR_hypo, 
     DMR_length_GW_figure, DMR_count_GW_figure, DMR_dis_GW_figure, DMR_freq_GW_figure, DMR_pos_GW_figure, 
     GREAT_GW_Cortex02.UMR, GREAT_GW_Cortex04.UMR, GREAT_GW_GE02.UMR, GREAT_GW_GE04.UMR, GREAT_GW_17week.UMR,  
     genomicBreak_GW, genomicBreak_GW_figure, 
     cortex01_cortex03DE, cortex01_cortex04DE, cortex02_cortex03DE, cortex02_cortex04DE, 
     venn_GW_17_13_UP_cortex, venn_GW_17_13_DN_cortex, GW_17_13_UP_duplicated_cortex, GW_17_13_DN_duplicated_cortex, 
     GE01_GE03DE, GE01_GE04DE, GE02_GE03DE, GE02_GE04DE, GW_DE_summary, 
     venn_GW_17_13_UP_GE, venn_GW_17_13_DN_GE, GW_17_13_UP_duplicated_GE, GW_17_13_DN_duplicated_GE, 
     GW_17_13_UP_duplicated, GW_17_13_DN_duplicated, venn_GW_17_13_UP, venn_GW_17_13_DN, 
     enrich_GW_17_13_UP, enrich_GW_17_13_DN, enrich_GW_17_13_UP_cortex, enrich_GW_17_13_DN_cortex, enrich_GW_17_13_UP_GE, enrich_GW_17_13_DN_GE, 
     DMR_DE_Cortex02_Cortex04_hyper, DMR_DE_Cortex02_Cortex04_hypo, DMR_proximal_GE02_GE04_hyper, DMR_proximal_GE02_GE04_hypo, 
     DMR_proximal_GW_summary, venn_DMR_proximal_GW_hyper, venn_DMR_proximal_GW_hypo, 
     cortex01_cortex03_isoform, cortex01_cortex04_isoform, cortex02_cortex03_isoform, cortex02_cortex04_isoform, 
     GE01_GE03_isoform, GE01_GE04_isoform, GE02_GE03_isoform, GE02_GE04_isoform, GW_isoform_summary, enrich_GW_isoform_gene_dup, 
     venn_GW_isoform_cortex, GW_isoform_gene_dup_cortex, venn_GW_isoform_GE, GW_isoform_gene_dup_GE, venn_GW_isoform, GW_isoform_gene_dup, 
     file = "/projects/epigenomics/users/lli/FetalBrain/GW/GW.Rdata")

