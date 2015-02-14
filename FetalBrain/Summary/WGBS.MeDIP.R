#' Compare WGBS and MeDIP DMRs - HuFNSC02 cortex vs GE
setwd("/projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/")
library(ggplot2)
library(labeling)
library(VennDiagram)
library(dplyr)
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
load("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR_neurospheres.Rdata")
load("/projects/epigenomics/users/lli/FetalBrain/WGBS/WGBS.DMR.Rdata")

#' genomic breakdown 
genomicBreak <- read.delim("./CpG/genomic.breakdown.summary", head = F, as.is = T, row.names = 1, col.names = c("Name", "Total", "Intergenic", "Intron", "Exon", "Gene", "Promoter", "CGI"))
genomicBreak_tall <- genomicBreak[, -1]
genomicBreak_tall <- data.frame(Intersect = rep(row.names(genomicBreak_tall), ncol(genomicBreak_tall)), Region = factor(rep(colnames(genomicBreak_tall), each = nrow(genomicBreak_tall)), levels = colnames(genomicBreak_tall)), CpG = as.vector(as.matrix(genomicBreak_tall)))
genomicBreak_tall$DM <- gsub("DM.", "", genomicBreak_tall$Intersect)
genomicBreak_tall$DM <- gsub("\\..+", "", genomicBreak_tall$DM)
genomicBreak_tall$Intersect <- gsub("DM\\.hyp[ero]+\\.", "", genomicBreak_tall$Intersect)
(genomicBreak_figure <- ggplot(genomicBreak_tall, aes(x = Region, y = CpG, fill = Intersect)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5) + 
  xlab("") + 
  ylab("Fraction of CpG") + 
  ggtitle("WGBS and MeDIP DM CpG breakdown") + 
  facet_wrap(~ DM) + 
  coord_flip() + 
  scale_fill_hue(l = 50) + 
  theme_bw())
ggsave(genomicBreak_figure, file = "WGBS_MeDIP_DM_genomicBreak.pdf")

#' significance of overlap: hypergeometric test, return log(p-value)  
total = 28217448    # Total No. of CpGs
hyper.intersect = genomicBreak["DM.hyper.intersect", "Total"]
hyper.WGBS = hyper.intersect + genomicBreak["DM.hyper.WGBS", "Total"]
hyper.MeDIP = hyper.intersect + genomicBreak["DM.hyper.MeDIP", "Total"]
hypo.intersect = genomicBreak["DM.hypo.intersect", "Total"]
hypo.WGBS = hypo.intersect + genomicBreak["DM.hypo.WGBS", "Total"]
hypo.MeDIP = hypo.intersect + genomicBreak["DM.hypo.MeDIP", "Total"]
phyper(hyper.intersect, hyper.WGBS, total - hyper.WGBS, hyper.MeDIP, lower.tail = F, log = T) 
phyper(hypo.intersect, hypo.WGBS, total - hypo.WGBS, hypo.MeDIP, lower.tail = F, log = T)
pdf("Venn_CpG_Cortex_UMR.pdf")
plot.new()
grid.text("Intersect of Cortex UMR CpGs between MeDIP and WGBS", x = unit(0.5, "npc"), y = unit(0.95, "npc"))
venn_Cortex_UMR_WGBS_MeDIP <- draw.pairwise.venn(area1 = hypo.WGBS, area2 = hypo.MeDIP, cross.area = hypo.intersect, category = c("WGBS", "MeDIP"), fill = c("red", "blue"))
dev.off()
pdf("Venn_CpG_GE_UMR.pdf")
plot.new()
grid.text("Intersect of GE UMR CpGs between MeDIP and WGBS", x = unit(0.5, "npc"), y = unit(0.95, "npc"))
venn_GE_UMR_WGBS_MeDIP <- draw.pairwise.venn(area1 = hyper.WGBS, area2 = hyper.MeDIP, cross.area = hyper.intersect, category = c("WGBS", "MeDIP"), fill = c("red", "blue"))
dev.off()

#' compare GC content 
GC.hyper.intersect <- read.delim("GC.GE_UMR.intersect.txt", head = F, as.is = T, skip = 1)
GC.hyper.WGBS <- read.delim("GC.GE_UMR.WGBS.txt", head = F, as.is = T, skip = 1)
GC.hyper.MeDIP <- read.delim("GC.GE_UMR.MeDIP.txt", head = F, as.is = T, skip = 1)
GC.hypo.intersect <- read.delim("GC.Cortex_UMR.intersect.txt", head = F, as.is = T, skip = 1)
GC.hypo.WGBS <- read.delim("GC.Cortex_UMR.WGBS.txt", head = F, as.is = T, skip = 1)
GC.hypo.MeDIP <- read.delim("GC.Cortex_UMR.MeDIP.txt", head = F, as.is = T, skip = 1)
GC <- rbind(data.frame(DM = "hyper", Intersect = "intersect", ID = GC.hyper.intersect$V4, GC = GC.hyper.intersect$V6), 
            data.frame(DM = "hyper", Intersect = "WGBS", ID = GC.hyper.WGBS$V4, GC = GC.hyper.WGBS$V6), 
            data.frame(DM = "hyper", Intersect = "MeDIP", ID = GC.hyper.MeDIP$V4, GC = GC.hyper.MeDIP$V6), 
            data.frame(DM = "hypo", Intersect = "intersect", ID = GC.hypo.intersect$V4, GC = GC.hypo.intersect$V6), 
            data.frame(DM = "hypo", Intersect = "WGBS", ID = GC.hypo.WGBS$V4, GC = GC.hypo.WGBS$V6), 
            data.frame(DM = "hypo", Intersect = "MeDIP", ID = GC.hypo.MeDIP$V4, GC = GC.hypo.MeDIP$V6))
(GC_figure <- ggplot(GC, aes(x = Intersect, y = GC, color = Intersect)) + 
  geom_boxplot() + 
  xlab("") + 
  ylab("GC content") + 
  ggtitle("GC content for WGBS and MeDIP DMRs") + 
  facet_wrap(~ DM) + 
  scale_color_hue(l = 50) + 
  theme_bw())
ggsave(GC_figure, file = "WGBS_MeDIP_DMR_GCcontent.pdf")
t.test(GC.hyper.WGBS$V6, GC.hyper.MeDIP$V6)$p.value
t.test(GC.hypo.WGBS$V6, GC.hypo.MeDIP$V6)$p.value

#' intersect closest genes
total = 52475  # wc -l /home/lli/hg19/hg19v65_genes.bed 
gene.hyper.intersect = as.integer(unlist(strsplit(system("wc -l GE_UMR_closest_gene.intersect", intern = T), " "))[1])
gene.hyper.WGBS = gene.hyper.intersect + as.integer(unlist(strsplit(system("wc -l GE_UMR_closest_gene.WGBS", intern = T), " "))[1])
gene.hyper.MeDIP = gene.hyper.intersect + as.integer(unlist(strsplit(system("wc -l GE_UMR_closest_gene.MeDIP", intern = T), " "))[1])
gene.hypo.intersect = as.integer(unlist(strsplit(system("wc -l Cortex_UMR_closest_gene.intersect", intern = T), " "))[1])
gene.hypo.WGBS = gene.hypo.intersect + as.integer(unlist(strsplit(system("wc -l Cortex_UMR_closest_gene.WGBS", intern = T), " "))[1])
gene.hypo.MeDIP = gene.hypo.intersect + as.integer(unlist(strsplit(system("wc -l Cortex_UMR_closest_gene.MeDIP", intern = T), " "))[1])
phyper(gene.hyper.intersect, gene.hyper.WGBS, total - gene.hyper.WGBS, gene.hyper.MeDIP, lower.tail = F, log = T) 
phyper(gene.hypo.intersect, gene.hypo.WGBS, total - gene.hypo.WGBS, gene.hypo.MeDIP, lower.tail = F, log = T)
pdf("Venn_genes_Cortex_UMR.pdf")
plot.new()
grid.text("Intersect of Cortex UMR closest genes between MeDIP and WGBS", x = unit(0.5, "npc"), y = unit(0.95, "npc"))
venn_Cortex_UMR_WGBS_MeDIP <- draw.pairwise.venn(area1 = gene.hypo.WGBS, area2 = gene.hypo.MeDIP, cross.area = gene.hypo.intersect, category = c("WGBS", "MeDIP"), fill = c("red", "blue"))
dev.off()
pdf("Venn_genes_GE_UMR.pdf")
plot.new()
grid.text("Intersect of GE UMR closest genes between MeDIP and WGBS", x = unit(0.5, "npc"), y = unit(0.95, "npc"))
venn_GE_UMR_WGBS_MeDIP <- draw.pairwise.venn(area1 = gene.hyper.WGBS, area2 = gene.hyper.MeDIP, cross.area = gene.hyper.intersect, category = c("WGBS", "MeDIP"), fill = c("red", "blue"))
dev.off()

#' differences in DMR length and No. of CpGs 
DMR_length_MeDIP <- select(mutate(filter(DMR_length_neurospheres_figure$data, donor == "HuFNSC02"), assay = "MeDIP"), -min, -max, -donor)
DMR_count_MeDIP <- select(mutate(filter(DMR_count_neurospheres_figure$data, donor == "HuFNSC02"), assay = "MeDIP"), -min, -max, -donor)
DMR_length_WGBS <- mutate(Cortex02_GE02_DMR_figures$length$data, assay = "WGBS")
DMR_count_WGBS <- mutate(Cortex02_GE02_DMR_figures$count$data, assay = "WGBS")
DMR_length_MeDIP_WGBS <- rbind(DMR_length_MeDIP, DMR_length_WGBS)
DMR_count_MeDIP_WGBS <- rbind(DMR_count_MeDIP, DMR_count_WGBS)
(DMR_length_figure <- ggplot(DMR_length_MeDIP_WGBS, aes(x = chr, fill = chr)) + 
  geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
  xlab("") + 
  ylab("DMR length (bp)") + 
  ggtitle(paste("DMR length - WGBS vs MeDIP")) + 
  guides(fill = F) + 
  facet_grid(DM ~ assay, scale = "free") + 
  theme_bw())
ggsave(DMR_length_figure, file = "DMRlength_WGBS_MeDIP.pdf")
(DMR_count_figure <- ggplot(DMR_count_MeDIP_WGBS, aes(x = chr, fill = chr)) + 
   geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
   xlab("") + 
   ylab("No. of CpGs per DMR") + 
   ggtitle(paste("No. of CpGs per DMR - WGBS vs MeDIP")) + 
   guides(fill = F) + 
   facet_grid(DM ~ assay, scale = "free") + 
   theme_bw())
ggsave(DMR_count_figure, file = "DMRcount_WGBS_MeDIP.pdf")

#' intersect proximal genes & DE genes 
load("~/hg19/hg19v65_genes.Rdata")
WGBS_Cortex_UMR_proximal <- read.delim("/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/CpG/DMR_pcPromoter.Cortex-HuFNSC02_GE-HuFNSC02.hypo.txt", head = F, as.is = F)
WGBS_GE_UMR_proximal <- read.delim("/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/CpG/DMR_pcPromoter.Cortex-HuFNSC02_GE-HuFNSC02.hyper.txt", head = F, as.is = F)
MeDIP_Cortex_UMR_proximal <- read.delim("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/CpG/DMR_gene.Cortex-HuFNSC02_GE-HuFNSC02.hypo.proximal.txt", head = F, as.is = F)
MeDIP_GE_UMR_proximal <- read.delim("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/CpG/DMR_gene.Cortex-HuFNSC02_GE-HuFNSC02.hyper.proximal.txt", head = F, as.is = F)
pdf("Venn_proximal_Cortex_UMR.pdf")
plot.new()
grid.text("Intersect of Cortex UMR proximal genes between MeDIP and WGBS", x = unit(0.5, "npc"), y = unit(0.95, "npc"))
venn_Cortex_UMR_WGBS_MeDIP_proximal <- draw.pairwise.venn(area1 = nrow(WGBS_Cortex_UMR_proximal), area2 = nrow(MeDIP_Cortex_UMR_proximal), cross.area = sum(as.character(WGBS_Cortex_UMR_proximal$V5) %in% as.character(MeDIP_Cortex_UMR_proximal$V2)), category = c("WGBS", "MeDIP"), fill = c("red", "blue"), ext.text = F)
dev.off()
pdf("Venn_proximal_GE_UMR.pdf")
plot.new()
grid.text("Intersect of GE UMR proximal genes between MeDIP and WGBS", x = unit(0.5, "npc"), y = unit(0.95, "npc"))
venn_GE_UMR_WGBS_MeDIP_proximal <- draw.pairwise.venn(area1 = nrow(WGBS_GE_UMR_proximal), area2 = nrow(MeDIP_GE_UMR_proximal), cross.area = sum(as.character(WGBS_GE_UMR_proximal$V5) %in% as.character(MeDIP_GE_UMR_proximal$V2)), category = c("WGBS", "MeDIP"), fill = c("red", "blue"), ext.text = F)
dev.off()

WGBS_Cortex_UMR_proximal_DE <- read.delim("/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/CpG/DMR_pcPromoter_DE.Cortex-HuFNSC02_GE-HuFNSC02.hypo.txt", head = F, as.is = F)
WGBS_GE_UMR_proximal_DE <- read.delim("/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/CpG/DMR_pcPromoter_DE.Cortex-HuFNSC02_GE-HuFNSC02.hyper.txt", head = F, as.is = F)
MeDIP_Cortex_UMR_proximal_DE <- read.delim("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/CpG/DMR_gene_DE.Cortex-HuFNSC02_GE-HuFNSC02.hypo.proximal.txt", head = F, as.is = F)
MeDIP_GE_UMR_proximal_DE <- read.delim("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/CpG/DMR_gene_DE.Cortex-HuFNSC02_GE-HuFNSC02.hyper.proximal.txt", head = F, as.is = F)
pdf("Venn_proximal_DE_Cortex_UMR.pdf")
plot.new()
grid.text("Intersect of Cortex UMR proximal DE genes between MeDIP and WGBS", x = unit(0.5, "npc"), y = unit(0.95, "npc"))
venn_Cortex_UMR_WGBS_MeDIP_proximal_DE <- draw.pairwise.venn(area1 = nrow(WGBS_Cortex_UMR_proximal_DE), area2 = nrow(MeDIP_Cortex_UMR_proximal_DE), cross.area = sum(as.character(WGBS_Cortex_UMR_proximal_DE$V6) %in% as.character(MeDIP_Cortex_UMR_proximal_DE$V2)), category = c("WGBS", "MeDIP"), fill = c("red", "blue"), ext.text = F)
dev.off()
pdf("Venn_proximal_DE_GE_UMR.pdf")
plot.new()
grid.text("Intersect of GE UMR proximal DE genes between MeDIP and WGBS", x = unit(0.5, "npc"), y = unit(0.95, "npc"))
venn_GE_UMR_WGBS_MeDIP_proximal_DE <- draw.pairwise.venn(area1 = nrow(WGBS_GE_UMR_proximal_DE), area2 = nrow(MeDIP_GE_UMR_proximal_DE), cross.area = sum(as.character(WGBS_GE_UMR_proximal_DE$V6) %in% as.character(MeDIP_GE_UMR_proximal_DE$V2)), category = c("WGBS", "MeDIP"), fill = c("red", "blue"), ext.text = F)
dev.off()
ensembl[intersect(as.character(WGBS_Cortex_UMR_proximal_DE$V6), as.character(MeDIP_Cortex_UMR_proximal_DE$V2)), ]
ensembl[intersect(as.character(WGBS_GE_UMR_proximal_DE$V6), as.character(MeDIP_GE_UMR_proximal_DE$V2)),]

#' TF 
WGBS_TF <- read.delim("/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/TF/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.TF.summary", head = F, as.is = T, col.names = c("TF", "Cortex02UMR", "GE02UMR", "Ratio02"))
MeDIP_TF <- read.delim("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/TF/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.s300.c4.TF.summary", head = F, as.is = T, col.names = c("TF", "Cortex02UMR", "GE02UMR", "Ratio02"))
cat(intersect(WGBS_TF[WGBS_TF$Ratio02 <= 1/2, "TF"], MeDIP_TF[MeDIP_TF$Ratio02 <= 1/2, "TF"]), sep = ", ")
cat(intersect(WGBS_TF[WGBS_TF$Ratio02 >= 2, "TF"], MeDIP_TF[MeDIP_TF$Ratio02 >= 2, "TF"]), sep = ", ")
combined <- merge(WGBS_TF, MeDIP_TF, by = "TF")
cor(combined$Ratio02.x, combined$Ratio02.y)

#' Distance to nearest genes
MeDIP_closestGene_hyper <- read.delim("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.s300.c4.hyper.closest.gene", head = F, as.is = T)
MeDIP_closestGene_hypo <- read.delim("/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.s300.c4.hypo.closest.gene", head = F, as.is = T)
WGBS_closestGene_hyper <- read.delim("/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.closest.gene", head = F, as.is = T)
WGBS_closestGene_hypo <- read.delim("/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.closest.gene", head = F, as.is = T)
dis_closestGene <- rbind(mutate(select(MeDIP_closestGene_hyper, V8, V9), Assay = "MeDIP", DM = "hyper"), 
                         mutate(select(MeDIP_closestGene_hypo, V8, V9), Assay = "MeDIP", DM = "hypo"), 
                         mutate(select(WGBS_closestGene_hyper, V8, V9), Assay = "WGBS", DM = "hyper"), 
                         mutate(select(WGBS_closestGene_hypo, V8, V9), Assay = "WGBS", DM = "hypo"))
(non_genebody <- nrow(filter(dis_closestGene, V9 > 0)) / nrow(dis_closestGene))
(non_genebody_MeDIP <- nrow(filter(dis_closestGene, V9 > 0, Assay == "MeDIP")) / nrow(filter(dis_closestGene, Assay == "MeDIP")))
(non_genebody_WGBS <- nrow(filter(dis_closestGene, V9 > 0, Assay == "WGBS")) / nrow(filter(dis_closestGene, Assay == "WGBS")))
dis_closestGene <- filter(dis_closestGene, V9 > 0)
(dis_closestGene_figure <- ggplot(dis_closestGene, aes(color = Assay)) + 
   geom_boxplot(aes(Assay, V9), outlier.shape = NA) + 
   # geom_density(aes(V9)) + 
   facet_wrap(~ DM) + 
   xlab("") + 
   ylab("Distance to nearest genes") + 
   coord_cartesian(ylim = c(0, 1.5e+5)) + 
   # coord_cartesian(xlim = c(0, 0.5e+5)) + 
   theme_bw())
ggsave(dis_closestGene_figure, file = "dis_closestGene_figure.pdf", height = 8, width = 10)
t.test(select(filter(dis_closestGene, DM == "hyper", Assay == "MeDIP"), V9), select(filter(dis_closestGene, DM == "hyper", Assay == "WGBS"), V9))
t.test(select(filter(dis_closestGene, DM == "hypo", Assay == "MeDIP"), V9), select(filter(dis_closestGene, DM == "hypo", Assay == "WGBS"), V9))

#' hydroxymethylation: high WGBS & low MeDIP
## WGBS-MeDIP distribution
Cortex02_diff_summary <- read.delim("Cortex02_WGBS_MeDIP.diff.summary", head = F, as.is = T)
GE02_diff_summary <- read.delim("GE02_WGBS_MeDIP.diff.summary", head = F, as.is = T)
WGBS_MeDIP_diff_summary <- rbind(data.frame(Cortex02_diff_summary, Cell = "Cortex"), data.frame(GE02_diff_summary, Cell = "GE"))
(WGBS_MeDIP_diff_summary_figure <- ggplot(WGBS_MeDIP_diff_summary, aes(x = V1/10, y = V2, fill = Cell)) + 
   geom_bar(stat = "identity", position = "dodge") + 
   xlab("WGBS-MeDIP") + 
   ylab("No. of CpGs") + 
   scale_fill_manual(values = c("blue", "red"), name = "") + 
   theme_bw())
ggsave(WGBS_MeDIP_diff_summary_figure, file = "WGBS_MeDIP_diff_summary_figure.pdf")

## intersect Cortex and GE potential hydroxy sites
Cortex_hydroxy_CpG <- as.numeric(system("less /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/Cortex02_WGBS_MeDIP.diff.bed | wc -l", intern = T))
GE_hydroxy_CpG <- as.numeric(system("less /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/GE02_WGBS_MeDIP.diff.bed | wc -l", intern = T))
Cortex_GE_hydroxy_CpG <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/Cortex02_WGBS_MeDIP.diff.bed -b /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/GE02_WGBS_MeDIP.diff.bed | wc -l", intern = T))
pdf("Venn_Cortex_GE_hydroxy_CpG.pdf")
plot.new()
grid.text("Intersect of Cortex and GE potential hydroxymethylation CpGs", x = unit(0.5, "npc"), y = unit(0.95, "npc"))
Venn_Cortex_GE_hydroxy_CpG <- draw.pairwise.venn(area1 = Cortex_hydroxy_CpG, area2 = GE_hydroxy_CpG, cross.area = Cortex_GE_hydroxy_CpG, category = c("Cortex", "GE"), fill = c("blue", "red"), ext.text = F)
dev.off()
total = 28217448    # Total No. of CpGs
phyper(Cortex_GE_hydroxy_CpG, Cortex_hydroxy_CpG, total - Cortex_hydroxy_CpG, GE_hydroxy_CpG, lower.tail = F, log = T) 

Cortex_hydroxy_region <- as.numeric(system("less /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/Cortex02_WGBS_MeDIP.region.bed | wc -l", intern = T))
GE_hydroxy_region <- as.numeric(system("less /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/GE02_WGBS_MeDIP.region.bed | wc -l", intern = T))
Cortex_GE_hydroxy_region <- as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/Cortex02_WGBS_MeDIP.region.bed -b /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/GE02_WGBS_MeDIP.region.bed | wc -l", intern = T))
pdf("Venn_Cortex_GE_hydroxy_region.pdf")
plot.new()
grid.text("Intersect of Cortex and GE potential hydroxymethylation regions", x = unit(0.5, "npc"), y = unit(0.95, "npc"))
Venn_Cortex_GE_hydroxy_region <- draw.pairwise.venn(area1 = Cortex_hydroxy_region, area2 = GE_hydroxy_region, cross.area = Cortex_GE_hydroxy_region, category = c("Cortex", "GE"), fill = c("blue", "red"), ext.text = F)
dev.off()

## genomic breakdown of potential hydroxy CpGs
genomicBreak_hydroxy <- read.delim("hydroxy.breakdown.summary", head = F, as.is = T, row.names = 1, col.names = c("Name", "Total", "Intergenic", "Intron", "Exon", "Gene", "Promoter", "CGI"))
genomicBreak_hydroxy_tall <- data.frame(Sample = rep(row.names(genomicBreak_hydroxy), ncol(genomicBreak_hydroxy[,-1])), Region = factor(rep(colnames(genomicBreak_hydroxy[,-1]), each = nrow(genomicBreak_hydroxy[,-1])), levels = colnames(genomicBreak_hydroxy[,-1])), CpG = as.vector(as.matrix(genomicBreak_hydroxy[,-1])))
genomicBreak_hydroxy_figure <- ggplot(genomicBreak_hydroxy_tall, aes(x = Region, y = CpG, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5) + 
  xlab("") + 
  ylab("Fraction of CpG") + 
  coord_flip() + 
  scale_fill_manual(values = c("blue", "red"), name = "") + 
  theme_bw()
ggsave(genomicBreak_hydroxy_figure, file = "genomicBreak_hydroxy.pdf")
(genomicBreak_hydroxy_figure + ggtitle("Genomic breakdown of potential hydroxy CpGs in Cortex02 and GE02"))

## GREAT on potential hdroxy regions
(GREAT_Cortex02.hydroxy <- enrich_GREAT(file = "hydroxy_Cortex02", name = "Cortex02.hydroxy", height = 12))
(GREAT_GE02.hydroxy <- enrich_GREAT(file = "hydroxy_GE02", name = "GE02.hydroxy", height = 12))

hydroxy.summary <- data.frame(CpG = c(Cortex_hydroxy_CpG, GE_hydroxy_CpG), Region = c(Cortex_hydroxy_region, GE_hydroxy_region), 
                      CpG.in.WGBS.UMRs = c(as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/Cortex02_WGBS_MeDIP.diff.bed -b /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed | wc -l", intern = T)), 
                                           as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/GE02_WGBS_MeDIP.diff.bed -b /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed | wc -l", intern = T))), 
                      WGBS.UMRs.with.CpG = c(as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed -b /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/Cortex02_WGBS_MeDIP.diff.bed -u | wc -l", intern = T)), 
                                             as.numeric(system("/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed -b /projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/GE02_WGBS_MeDIP.diff.bed -u | wc -l", intern = T))))
rownames(hydroxy.summary) <- c("Cortex02", "GE02")
save(hydroxy.summary, Cortex_GE_hydroxy_CpG, Cortex_GE_hydroxy_region, WGBS_MeDIP_diff_summary_figure, Venn_Cortex_GE_hydroxy_CpG, genomicBreak_hydroxy_figure, GREAT_Cortex02.hydroxy, GREAT_GE02.hydroxy, 
     file = "hydroxy.Rdata")




