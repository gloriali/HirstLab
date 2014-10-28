#' DMR analysis from WGBS fractional calls  

#' DMR identification: ~/HirstLab/FetalBrain/WGBS/WGBS.DM.sh

setwd("~/快盘/FetalBrain/WGBS/")
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(ggbio)
library(GenomicRanges)
source('~/HirstLab/Pipeline/R/DMR.figures.R')
source('~/HirstLab/Pipeline/R/enrich_GREAT.R')
source('~/HirstLab/Pipeline/R/enrich.R')
load("~/快盘/FetalBrain/RNAseq/DEfine/gene/FetalBrain_DEgenes.Rdata")

#' testing DMR script parameters
DM_test_summary <- read.delim("DM.summary.stats", head = F, as.is = T, col.names = c("Sample", "p", "delta", "Total.DM.CpGs", "Hyper.DM.CpGs", "Hypo.DM.CpGs"))
DM_test_summary$delta <- paste0("delta=", DM_test_summary$delta)
DMR_test_summary <- read.delim("DMR.summary.stats", head = F, as.is = T, col.names = c("Sample", "size", "CpG.cut", "Median.length", "Median.CpG", "Total.DMR", "Hyper.DMR", "Hypo.DMR"))
DMR_test_summary$p <- gsub(".d[0-9.]+", "", DMR_test_summary$Sample)
DMR_test_summary$p <- paste0("p=", gsub(".*.p", "", DMR_test_summary$p))
DMR_test_summary$delta <- paste0("delta=", gsub(".*.d", "", DMR_test_summary$Sample))
DMR_test_summary$Sample <- gsub(".m.*.p.*.d.*", "", DMR_test_summary$Sample)
pdf("DMR_parameters.pdf")
(ggplot(DM_test_summary, aes(x = p, y = Total.DM.CpGs)) + scale_x_log10(breaks = c(0.05, 0.005, 0.0005)) + geom_point(aes(color = Sample)) + geom_line(aes(color = Sample)) + facet_wrap(~ delta) + theme_bw() + xlab("p-value cutoff") + ylab("Total No. of DM CpGs"))
(ggplot(DM_test_summary, aes(x = p, y = log2(Hyper.DM.CpGs/Hypo.DM.CpGs), color = Sample)) + scale_x_log10(breaks = c(0.05, 0.005, 0.0005)) + geom_point() + geom_line() + geom_hline(aes(yintercept = 0)) + facet_wrap(~ delta) + theme_bw() + xlab("p-value cutoff") + ylab("log2(Hyper/Hypo)"))
(ggplot(DMR_test_summary, aes(x = size, y = Median.length)) + geom_point(aes(color = Sample)) + geom_line(aes(color = Sample)) + facet_grid(p ~ delta) + geom_hline(aes(yintercept = 250)) + theme_bw() + xlab("Max distance between adjacent DM CpGs") + ylab("Median DMR length"))
(ggplot(DMR_test_summary, aes(x = size, y = Total.DMR, color = Sample)) + geom_point() + geom_line() + facet_grid(p ~ delta) + theme_bw() + xlab("Max distance between adjacent DM CpGs") + ylab("Total No. of DMRs"))
(ggplot(DMR_test_summary, aes(x = p, y = Total.DMR, color = Sample)) + geom_point() + geom_line() + facet_grid(size ~ delta) + theme_bw() + xlab("p-value cutoff") + ylab("Total No. of DMRs"))
(ggplot(DMR_test_summary, aes(x = size, y = log2(Hyper.DMR/Hypo.DMR), color = Sample)) + geom_point() + geom_line() + geom_hline(aes(yintercept = 0)) + facet_grid(p ~ delta) + theme_bw() + xlab("Max distance between adjacent DM CpGs") + ylab("log2(Hyper/Hypo)"))
dev.off()

#' preferable parameters: p = 0.005 (two-sided), delta = 0.5, m = 0.75, size = 300
setwd("~/快盘/FetalBrain/WGBS/DMR/")
DMR_WGBS_summary <- select(filter(DMR_test_summary, size==300, p=="p=0.005", delta=="delta=0.5"), Sample, contains("DMR"))
Cortex02_GE02_DMR_WGBS <- read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3", as.is = T, head = F, col.names = c("chr", "start", "end", "ID", "DM", "count", "length"))
Cortex02_GE02_DMR_WGBS_figures <- DMR_figures(Cortex02_GE02_DMR_WGBS, sample1 = "Cortex-HuFNSC02", sample2 = "GE-HuFNSC02")
(Cortex02_GE02_DMR_WGBS_figures$length)
(Cortex02_GE02_DMR_WGBS_figures$count)
(Cortex02_GE02_DMR_WGBS_figures$dis)
(Cortex02_GE02_DMR_WGBS_figures$freq)
(Cortex02_GE02_DMR_WGBS_figures$pos)
Cortex04_GE04_DMR_WGBS <- read.delim("DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3", as.is = T, head = F, col.names = c("chr", "start", "end", "ID", "DM", "count", "length"))
Cortex04_GE04_DMR_WGBS_figures <- DMR_figures(Cortex04_GE04_DMR_WGBS, sample1 = "Cortex-HuFNSC04", sample2 = "GE-HuFNSC04")
(Cortex04_GE04_DMR_WGBS_figures$length)
(Cortex04_GE04_DMR_WGBS_figures$count)
(Cortex04_GE04_DMR_WGBS_figures$dis)
(Cortex04_GE04_DMR_WGBS_figures$freq)
(Cortex04_GE04_DMR_WGBS_figures$pos)
DMR_freq_WGBS <- rbind(cbind(Cortex02_GE02_DMR_WGBS_figures$freq$data, donor = "HuFNSC02"), cbind(Cortex04_GE04_DMR_WGBS_figures$freq$data, donor = "HuFNSC04"))
DMR_freq_WGBS_figure <- ggplot(DMR_freq_WGBS, aes(x = chr, y = freq, fill = DM)) + 
  geom_bar(position = "identity", stat = "identity", width = 0.8) + 
  facet_wrap(~ donor) + 
  xlab("") + 
  ylab("DMR frequency (bp/MB)") + 
  coord_flip() + 
  scale_fill_manual(values = c("red", "blue"), labels = c("GE UMRs","Cortex UMRs"), name = "") + 
  theme_bw()
ggsave(DMR_freq_WGBS_figure, file = "DMRfreq_neurospheres.pdf", width = 14, height = 9)
(DMR_freq_WGBS_figure + ggtitle("DMR frequency asymmetry between cortex an GE neurospheres"))

DMR_pos_WGBS <- rbind(cbind(Cortex02_GE02_DMR_WGBS, donor = "HuFNSC02"), cbind(Cortex04_GE04_DMR_WGBS, donor = "HuFNSC04")) %>% mutate(pos = (start + end)/2, y = DM * 15 + 5)
DMR_pos_WGBS_gr <- keepSeqlevels(as(DMR_pos_WGBS, "GRanges"), paste0("chr", c(1:22, "X")))
data(hg19IdeogramCyto, package = "biovizBase")
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X")))
(DMR_pos_WGBS_figure <- ggplot(hg19) + 
   layout_karyogram(cytoband = TRUE) + 
   layout_karyogram(data = DMR_pos_WGBS_gr, geom = "point", aes(y = y, x = pos, color = factor(DM, levels = c("1", "-1"))), position = position_jitter(height = 2), size = 1.5, alpha = 0.5) + 
   ylab("") + 
   xlab("Position of DMRs on the chromosome") +
   facet_grid(seqnames ~ donor) + 
   coord_cartesian(ylim = c(-15, 25)) + 
   ggtitle(paste("DMR position - Cortex vs GE DMRs")) + 
   scale_color_manual(values = c("red", "blue"), labels = c("hyper", "hypo"), name = "") + 
   theme_bw() + 
   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_line(color = "transparent"), panel.grid.minor.y = element_line(color = "transparent")))
ggsave(DMR_pos_WGBS_figure@ggplot, file = "DMRpos_neurospheres.pdf", width = 10, height = 10)

#` GREAT enrichment analysis
(GREAT_Cortex02.UMR_WGBS <- enrich_GREAT(file = "DMR.Cortex-HuFNSC02_GE-HuFNSC02.hypo", name = "HuFNSC02-Cortex.UMRs", height = 12))
(GREAT_GE02.UMR_WGBS <- enrich_GREAT(file = "DMR.Cortex-HuFNSC02_GE-HuFNSC02.hyper", name = "HuFNSC02-GE.UMRs"))
(GREAT_Cortex04.UMR_WGBS <- enrich_GREAT(file = "DMR.Cortex-HuFNSC04_GE-HuFNSC04.hypo", name = "HuFNSC04-Cortex.UMRs", height = 13))
(GREAT_GE04.UMR_WGBS <- enrich_GREAT(file = "DMR.Cortex-HuFNSC04_GE-HuFNSC04.hyper", name = "HuFNSC04-GE.UMRs", height = 3))

GOBP_Cortex02.UMR <- cbind(GREAT_Cortex02.UMR_WGBS$data, UMR = "Cortex UMRs")
GOBP_Cortex02.UMR <- GOBP_Cortex02.UMR[as.character(GOBP_Cortex02.UMR$Category) == "GOBP", ]
GOBP_Cortex02.UMR_plot <- ggplot(data = GOBP_Cortex02.UMR, aes(Term, -log10(FDR))) +
  geom_bar(fill = "blue", stat = "identity", width = .5) + 
  facet_grid(UMR ~ .) + 
  coord_flip() + 
  geom_text(aes(label = round(-log10(FDR), 2), hjust = 0)) + 
  xlab("") + 
  theme_bw() +
  scale_fill_hue(l = 40)
(GOBP_Cortex02.UMR_plot + ggtitle("GOBP GREAT enrichment for Cortex UMRs"))
ggsave(GOBP_Cortex02.UMR_plot, file = "./enrich/GOBP_HuFNSC02_Cortex.UMR_enrich.pdf", height = 3, width = 6)
GOBP_GE02.UMR <- cbind(GREAT_GE02.UMR_WGBS$data, UMR = "GE UMRs")
GOBP_GE02.UMR <- GOBP_GE02.UMR[as.character(GOBP_GE02.UMR$Category) == "GOBP", ]
GOBP_GE02.UMR_plot <- ggplot(data = GOBP_GE02.UMR, aes(Term, -log10(FDR))) +
  geom_bar(fill = "blue", stat = "identity", width = .5) + 
  facet_grid(UMR ~ .) + 
  coord_flip() + 
  geom_text(aes(label = round(-log10(FDR), 2), hjust = 0)) + 
  xlab("") + 
  theme_bw() +
  scale_fill_hue(l = 40)
(GOBP_GE02.UMR_plot + ggtitle("GOBP GREAT enrichment for GE UMRs"))
ggsave(GOBP_GE02.UMR_plot, file = "./enrich/GOBP_HuFNSC02_GE.UMR_enrich.pdf", height = 3, width = 6)

GOBP_Cortex02.UMR_WGBS <- read.delim("./enrich/ALL_GREAT_GOBP_DMR.Cortex-HuFNSC02_GE-HuFNSC02.hypo.tsv", head = F, as.is = T, skip = 4)
GOBP_GE02.UMR_WGBS <- read.delim("./enrich/ALL_GREAT_GOBP_DMR.Cortex-HuFNSC02_GE-HuFNSC02.hyper.tsv", head = F, as.is = T, skip = 4)
GOBP_Cortex02_GE02_WGBS <- rbind(select(mutate(filter(GOBP_Cortex02.UMR_WGBS, V8 >= 2 & V7 <= 0.05 & V16 <= 0.05), UMR = "Cortex.UMR", FDR = -log10(V7), Term = V3), UMR, Term, FDR),
                                 select(mutate(filter(GOBP_GE02.UMR_WGBS, V8 >= 2 & V7 <= 0.05 & V16 <= 0.05), UMR = "GE.UMR", FDR = log10(V7), Term = V3), UMR, Term, FDR))
Terms <- unique(c(arrange(filter(GOBP_Cortex02_GE02_WGBS, UMR == "Cortex.UMR"), desc(FDR))[1:20, "Term"], rev(arrange(filter(GOBP_Cortex02_GE02_WGBS, UMR == "GE.UMR"), FDR)[1:20, "Term"])))
GOBP_Cortex02_GE02_WGBS <- mutate(filter(GOBP_Cortex02_GE02_WGBS, Term %in% Terms), Term = factor(Term, levels = rev(Terms)))
GOBP_Cortex02_GE02_WGBS_figure <- ggplot(data = GOBP_Cortex02_GE02_WGBS, aes(Term, FDR, fill = UMR)) +
  geom_bar(position = "identity", stat = "identity", width = .5) + 
  coord_flip() + 
  xlab("") + 
  ylab("-log10(FDR)") + 
  scale_fill_manual(values = c("blue", "red"), name = "") + 
  theme_bw()
(GOBP_Cortex02_GE02_WGBS_figure + ggtitle("GOBP GREAT enrichment for Cortex and GE UMRs"))
ggsave(GOBP_Cortex02_GE02_WGBS_figure, file = "./enrich/GOBP_Cortex02_GE02_WGBS_figure.pdf", height = 3, width = 6)

#' genomic break down  
genomicBreak_WGBS <- read.delim("./CpG/genomic.breakdown.summary", head = F, as.is = T, row.names = 1, col.names = c("Name", "Total", "Intergenic", "Intron", "Exon", "Gene", "Promoter", "CGI"))
genomicBreak_WGBS_tall <- genomicBreak_WGBS[, -1]
genomicBreak_WGBS_tall <- data.frame(Sample = rep(gsub(".m.*", "", rownames(genomicBreak_WGBS_tall)), ncol(genomicBreak_WGBS_tall)), DM = rep(gsub(".*.c3.", "", rownames(genomicBreak_WGBS_tall)), ncol(genomicBreak_WGBS_tall)), Region = factor(rep(colnames(genomicBreak_WGBS_tall), each = nrow(genomicBreak_WGBS_tall)), levels = colnames(genomicBreak_WGBS_tall)), CpG = as.vector(as.matrix(genomicBreak_WGBS_tall)))
genomicBreak_WGBS_figure <- ggplot(genomicBreak_WGBS_tall, aes(x = Region, y = CpG, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5) + 
  xlab("") + 
  ylab("Fraction of CpG") + 
  coord_flip() + 
  scale_fill_hue(l = 50) + 
  facet_wrap(~ DM) + 
  theme_bw()
(genomicBreak_WGBS_figure + ggtitle("DM CpG breakdown between WGBS Cortex and GE"))
ggsave(genomicBreak_WGBS_figure, file = "./CpG/genomicBreak_WGBS.pdf", height = 6)

#' intersecting with genes/promoters
setwd("~/快盘/FetalBrain/WGBS/DMR/CpG")
DMR_WGBS_gene_summary <- data.frame(matrix(NA, ncol = 7, nrow = 4, dimnames = list(c("GE_UMRs-HuFNSC02", "Cortex_UMRs-HuFNSC02", "GE_UMRs-HuFNSC04", "Cortex_UMRs-HuFNSC04"), c("pc.Genes", "unique.Genes", "pc.Promoters", "unique.Promoters", "proximal.DE.Genes", "same.direction", "unique.DE.Genes"))))
#' HuFNSC02
GE02_UMR_pcGene <- read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.CpG_gene_pc.bed", head = F, as.is = T)
GE02_UMR_pcGene$V4 <- gsub("_.*", "", GE02_UMR_pcGene$V4)
GE02_UMR_pcPromoter <- read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.CpG_promoter_pc.bed", head = F, as.is = T)
GE02_UMR_pcPromoter$V4 <- gsub("_.*", "", GE02_UMR_pcPromoter$V4)
GE02_UMR_pcPromoter_DE <- GE02_UMR_pcPromoter[GE02_UMR_pcPromoter$V5 %in% cortex02_GE02DE$V1, ]
GE02_UMR_pcPromoter_DE <- GE02_UMR_pcPromoter_DE[!duplicated(GE02_UMR_pcPromoter_DE$V4), ]
GE02_UMR_pcPromoter_DE$DE <- cortex02_GE02DE[GE02_UMR_pcPromoter_DE$V5, "DE"]
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC02", "pc.Genes"] <- length(unique(GE02_UMR_pcGene$V4))
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC02", "unique.Genes"] <- length(unique(GE02_UMR_pcGene$V5))
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC02", "pc.Promoters"] <- length(unique(GE02_UMR_pcPromoter$V4))
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC02", "unique.Promoters"] <- length(unique(GE02_UMR_pcPromoter$V5))
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC02", "proximal.DE.Genes"] <- nrow(GE02_UMR_pcPromoter_DE)
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC02", "same.direction"] <- sum(GE02_UMR_pcPromoter_DE$DE == "DN")
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC02", "unique.DE.Genes"] <- length(unique(GE02_UMR_pcPromoter_DE$V5))
GE02_UMR_pcPromoter <- GE02_UMR_pcPromoter[!duplicated(GE02_UMR_pcPromoter$V5), ]
GE02_UMR_pcPromoter <- cbind(GE02_UMR_pcPromoter, ensembl[GE02_UMR_pcPromoter$V5, ])
GE02_UMR_pcPromoter$V5 <- NULL
write.table(GE02_UMR_pcPromoter, file = "DMR_pcPromoter.Cortex-HuFNSC02_GE-HuFNSC02.hyper.txt", sep = "\t", quote = F, col.names = F, row.names = F)
# (GE02_UMR_pcPromoter_DAVID <- enrich(name = "DMR_pcPromoter."GE_UMRs-HuFNSC02", erminej = F))
GE02_UMR_pcPromoter_DE <- GE02_UMR_pcPromoter_DE[!duplicated(GE02_UMR_pcPromoter_DE$V5), ]
GE02_UMR_pcPromoter_DE <- cbind(GE02_UMR_pcPromoter_DE, ensembl[GE02_UMR_pcPromoter_DE$V5, ])
GE02_UMR_pcPromoter_DE$V5 <- NULL
write.table(GE02_UMR_pcPromoter_DE, file = "DMR_pcPromoter_DE.Cortex-HuFNSC02_GE-HuFNSC02.hyper.txt", sep = "\t", quote = F, col.names = F, row.names = F)
# (GE02_UMR_pcPromoter_DE_DAVID <- enrich(name = "DMR_pcPromoter_DE.GE_UMRs-HuFNSC02", erminej = F))
Cortex02_UMR_pcGene <- read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.CpG_gene_pc.bed", head = F, as.is = T)
Cortex02_UMR_pcGene$V4 <- gsub("_.*", "", Cortex02_UMR_pcGene$V4)
Cortex02_UMR_pcPromoter <- read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.CpG_promoter_pc.bed", head = F, as.is = T)
Cortex02_UMR_pcPromoter$V4 <- gsub("_.*", "", Cortex02_UMR_pcPromoter$V4)
Cortex02_UMR_pcPromoter_DE <- Cortex02_UMR_pcPromoter[Cortex02_UMR_pcPromoter$V5 %in% cortex02_GE02DE$V1, ]
Cortex02_UMR_pcPromoter_DE <- Cortex02_UMR_pcPromoter_DE[!duplicated(Cortex02_UMR_pcPromoter_DE$V4), ]
Cortex02_UMR_pcPromoter_DE$DE <- cortex02_GE02DE[Cortex02_UMR_pcPromoter_DE$V5, "DE"]
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC02", "pc.Genes"] <- length(unique(Cortex02_UMR_pcGene$V4))
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC02", "unique.Genes"] <- length(unique(Cortex02_UMR_pcGene$V5))
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC02", "pc.Promoters"] <- length(unique(Cortex02_UMR_pcPromoter$V4))
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC02", "unique.Promoters"] <- length(unique(Cortex02_UMR_pcPromoter$V5))
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC02", "proximal.DE.Genes"] <- nrow(Cortex02_UMR_pcPromoter_DE)
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC02", "same.direction"] <- sum(Cortex02_UMR_pcPromoter_DE$DE == "UP")
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC02", "unique.DE.Genes"] <- length(unique(Cortex02_UMR_pcPromoter_DE$V5))
Cortex02_UMR_pcPromoter <- Cortex02_UMR_pcPromoter[!duplicated(Cortex02_UMR_pcPromoter$V5), ]
Cortex02_UMR_pcPromoter <- cbind(Cortex02_UMR_pcPromoter, ensembl[Cortex02_UMR_pcPromoter$V5, ])
Cortex02_UMR_pcPromoter$V5 <- NULL
write.table(Cortex02_UMR_pcPromoter, file = "DMR_pcPromoter.Cortex-HuFNSC02_GE-HuFNSC02.hypo.txt", sep = "\t", quote = F, col.names = F, row.names = F)
# (Cortex02_UMR_pcPromoter_DAVID <- enrich(name = "DMR_pcPromoter."Cortex_UMRs-HuFNSC02", erminej = F))
Cortex02_UMR_pcPromoter_DE <- Cortex02_UMR_pcPromoter_DE[!duplicated(Cortex02_UMR_pcPromoter_DE$V5), ]
Cortex02_UMR_pcPromoter_DE <- cbind(Cortex02_UMR_pcPromoter_DE, ensembl[Cortex02_UMR_pcPromoter_DE$V5, ])
Cortex02_UMR_pcPromoter_DE$V5 <- NULL
write.table(Cortex02_UMR_pcPromoter_DE, file = "DMR_pcPromoter_DE.Cortex-HuFNSC02_GE-HuFNSC02.hypo.txt", sep = "\t", quote = F, col.names = F, row.names = F)
# (Cortex02_UMR_pcPromoter_DE_DAVID <- enrich(name = "DMR_pcPromoter_DE.Cortex_UMRs-HuFNSC02", erminej = F))

#' HuFNSC04
GE04_UMR_pcGene <- read.delim("DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.CpG_gene_pc.bed", head = F, as.is = T)
GE04_UMR_pcGene$V4 <- gsub("_.*", "", GE04_UMR_pcGene$V4)
GE04_UMR_pcPromoter <- read.delim("DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.CpG_promoter_pc.bed", head = F, as.is = T)
GE04_UMR_pcPromoter$V4 <- gsub("_.*", "", GE04_UMR_pcPromoter$V4)
GE04_UMR_pcPromoter_DE <- GE04_UMR_pcPromoter[GE04_UMR_pcPromoter$V5 %in% cortex04_GE04DE$V1, ]
GE04_UMR_pcPromoter_DE <- GE04_UMR_pcPromoter_DE[!duplicated(GE04_UMR_pcPromoter_DE$V4), ]
GE04_UMR_pcPromoter_DE$DE <- cortex04_GE04DE[GE04_UMR_pcPromoter_DE$V5, "DE"]
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC04", "pc.Genes"] <- length(unique(GE04_UMR_pcGene$V4))
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC04", "unique.Genes"] <- length(unique(GE04_UMR_pcGene$V5))
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC04", "pc.Promoters"] <- length(unique(GE04_UMR_pcPromoter$V4))
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC04", "unique.Promoters"] <- length(unique(GE04_UMR_pcPromoter$V5))
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC04", "proximal.DE.Genes"] <- nrow(GE04_UMR_pcPromoter_DE)
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC04", "same.direction"] <- sum(GE04_UMR_pcPromoter_DE$DE == "DN")
DMR_WGBS_gene_summary["GE_UMRs-HuFNSC04", "unique.DE.Genes"] <- length(unique(GE04_UMR_pcPromoter_DE$V5))
GE04_UMR_pcPromoter <- GE04_UMR_pcPromoter[!duplicated(GE04_UMR_pcPromoter$V5), ]
GE04_UMR_pcPromoter <- cbind(GE04_UMR_pcPromoter, ensembl[GE04_UMR_pcPromoter$V5, ])
GE04_UMR_pcPromoter$V5 <- NULL
write.table(GE04_UMR_pcPromoter, file = "DMR_pcPromoter.Cortex-HuFNSC04_GE-HuFNSC04.hyper.txt", sep = "\t", quote = F, col.names = F, row.names = F)
# (GE04_UMR_pcPromoter_DAVID <- enrich(name = "DMR_pcPromoter."GE_UMRs-HuFNSC04", erminej = F))
GE04_UMR_pcPromoter_DE <- GE04_UMR_pcPromoter_DE[!duplicated(GE04_UMR_pcPromoter_DE$V5), ]
GE04_UMR_pcPromoter_DE <- cbind(GE04_UMR_pcPromoter_DE, ensembl[GE04_UMR_pcPromoter_DE$V5, ])
GE04_UMR_pcPromoter_DE$V5 <- NULL
write.table(GE04_UMR_pcPromoter_DE, file = "DMR_pcPromoter_DE.Cortex-HuFNSC04_GE-HuFNSC04.hyper.txt", sep = "\t", quote = F, col.names = F, row.names = F)
# (GE04_UMR_pcPromoter_DE_DAVID <- enrich(name = "DMR_pcPromoter_DE.GE_UMRs-HuFNSC04", erminej = F))
Cortex04_UMR_pcGene <- read.delim("DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.CpG_gene_pc.bed", head = F, as.is = T)
Cortex04_UMR_pcGene$V4 <- gsub("_.*", "", Cortex04_UMR_pcGene$V4)
Cortex04_UMR_pcPromoter <- read.delim("DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.CpG_promoter_pc.bed", head = F, as.is = T)
Cortex04_UMR_pcPromoter$V4 <- gsub("_.*", "", Cortex04_UMR_pcPromoter$V4)
Cortex04_UMR_pcPromoter_DE <- Cortex04_UMR_pcPromoter[Cortex04_UMR_pcPromoter$V5 %in% cortex04_GE04DE$V1, ]
Cortex04_UMR_pcPromoter_DE <- Cortex04_UMR_pcPromoter_DE[!duplicated(Cortex04_UMR_pcPromoter_DE$V4), ]
Cortex04_UMR_pcPromoter_DE$DE <- cortex04_GE04DE[Cortex04_UMR_pcPromoter_DE$V5, "DE"]
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC04", "pc.Genes"] <- length(unique(Cortex04_UMR_pcGene$V4))
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC04", "unique.Genes"] <- length(unique(Cortex04_UMR_pcGene$V5))
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC04", "pc.Promoters"] <- length(unique(Cortex04_UMR_pcPromoter$V4))
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC04", "unique.Promoters"] <- length(unique(Cortex04_UMR_pcPromoter$V5))
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC04", "proximal.DE.Genes"] <- nrow(Cortex04_UMR_pcPromoter_DE)
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC04", "same.direction"] <- sum(Cortex04_UMR_pcPromoter_DE$DE == "UP")
DMR_WGBS_gene_summary["Cortex_UMRs-HuFNSC04", "unique.DE.Genes"] <- length(unique(Cortex04_UMR_pcPromoter_DE$V5))
Cortex04_UMR_pcPromoter <- Cortex04_UMR_pcPromoter[!duplicated(Cortex04_UMR_pcPromoter$V5), ]
Cortex04_UMR_pcPromoter <- cbind(Cortex04_UMR_pcPromoter, ensembl[Cortex04_UMR_pcPromoter$V5, ])
Cortex04_UMR_pcPromoter$V5 <- NULL
write.table(Cortex04_UMR_pcPromoter, file = "DMR_pcPromoter.Cortex-HuFNSC04_GE-HuFNSC04.hypo.txt", sep = "\t", quote = F, col.names = F, row.names = F)
# (Cortex04_UMR_pcPromoter_DAVID <- enrich(name = "DMR_pcPromoter."Cortex_UMRs-HuFNSC04", erminej = F))
Cortex04_UMR_pcPromoter_DE <- Cortex04_UMR_pcPromoter_DE[!duplicated(Cortex04_UMR_pcPromoter_DE$V5), ]
Cortex04_UMR_pcPromoter_DE <- cbind(Cortex04_UMR_pcPromoter_DE, ensembl[Cortex04_UMR_pcPromoter_DE$V5, ])
Cortex04_UMR_pcPromoter_DE$V5 <- NULL
write.table(Cortex04_UMR_pcPromoter_DE, file = "DMR_pcPromoter_DE.Cortex-HuFNSC04_GE-HuFNSC04.hypo.txt", sep = "\t", quote = F, col.names = F, row.names = F)
# (Cortex04_UMR_pcPromoter_DE_DAVID <- enrich(name = "DMR_pcPromoter_DE.Cortex_UMRs-HuFNSC04", erminej = F))

#' Venn diagrams 
venn_Cortex_UMR_WGBS <- draw.pairwise.venn(area1 = nrow(filter(Cortex02_GE02_DMR_WGBS, V5 == -1)), area2 = nrow(filter(Cortex04_GE04_DMR_WGBS, V5 == -1)), cross.area = 179, category = c("HuFNSC02", "HuFNSC04"), fill = c("red", "blue"))
venn_GE_UMR_WGBS <- draw.pairwise.venn(area1 = nrow(filter(Cortex02_GE02_DMR_WGBS, V5 == 1)), area2 = nrow(filter(Cortex04_GE04_DMR_WGBS, V5 == 1)), cross.area = 10, category = c("HuFNSC02", "HuFNSC04"), fill = c("red", "blue"))
Cortex_UMR_pcGene_WGBS <- list(HuFNSC02 = Cortex02_UMR_pcGene$V5, HuFNSC04 = Cortex04_UMR_pcGene$V5)
venn_Cortex_UMR_pcGene_WGBS <- venn.diagram(Cortex_UMR_pcGene_WGBS, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of genes with Cortex UMRs")
GE_UMR_pcGene_WGBS <- list(HuFNSC02 = GE02_UMR_pcGene$V5, HuFNSC04 = GE04_UMR_pcGene$V5)
venn_GE_UMR_pcGene_WGBS <- venn.diagram(GE_UMR_pcGene_WGBS, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of genes with GE UMRs")
Cortex_UMR_pcPromoter_WGBS <- list(HuFNSC02 = Cortex02_UMR_pcPromoter$id, HuFNSC04 = Cortex04_UMR_pcPromoter$id)
venn_Cortex_UMR_pcPromoter_WGBS <- venn.diagram(Cortex_UMR_pcPromoter_WGBS, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of genes with promoter Cortex UMRs")
GE_UMR_pcPromoter_WGBS <- list(HuFNSC02 = GE02_UMR_pcPromoter$id, HuFNSC04 = GE04_UMR_pcPromoter$id)
venn_GE_UMR_pcPromoter_WGBS <- venn.diagram(GE_UMR_pcPromoter_WGBS, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of genes with promoter GE UMRs")
Cortex_UMR_pcPromoter_DE_WGBS <- list(HuFNSC02 = Cortex02_UMR_pcPromoter_DE$id, HuFNSC04 = Cortex04_UMR_pcPromoter_DE$id)
venn_Cortex_UMR_pcPromoter_DE_WGBS <- venn.diagram(Cortex_UMR_pcPromoter_DE_WGBS, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of DE genes with promoter Cortex UMRs")
GE_UMR_pcPromoter_DE_WGBS <- list(HuFNSC02 = GE02_UMR_pcPromoter_DE$id, HuFNSC04 = GE04_UMR_pcPromoter_DE$id)
venn_GE_UMR_pcPromoter_DE_WGBS <- venn.diagram(GE_UMR_pcPromoter_DE_WGBS, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of DE genes with promoter GE UMRs")
pdf("venn_Cortex_UMR_WGBS.pdf")
plot.new()
grid.text("Venn diagram of Cortex UMRs", x = unit(0.5, "npc"), y = unit(0.9, "npc"))
grid.draw(venn_Cortex_UMR_WGBS)
dev.off()
pdf("venn_GE_UMR_WGBS.pdf")
plot.new()
grid.text("Venn diagram of GE UMRs", x = unit(0.5, "npc"), y = unit(0.9, "npc"))
grid.draw(venn_GE_UMR_WGBS)
dev.off()
pdf("venn_Cortex_UMR_pcGene_WGBS.pdf")
plot.new()
grid.draw(venn_Cortex_UMR_pcGene_WGBS)
dev.off()
pdf("venn_GE_UMR_pcGene_WGBS.pdf")
plot.new()
grid.draw(venn_GE_UMR_pcGene_WGBS)
dev.off()
pdf("venn_Cortex_UMR_pcPromoter_WGBS.pdf")
plot.new()
grid.draw(venn_Cortex_UMR_pcPromoter_WGBS)
dev.off()
pdf("venn_GE_UMR_pcPromoter_WGBS.pdf")
plot.new()
grid.draw(venn_GE_UMR_pcPromoter_WGBS)
dev.off()
pdf("venn_Cortex_UMR_pcPromoter_DE_WGBS.pdf")
plot.new()
grid.draw(venn_Cortex_UMR_pcPromoter_DE_WGBS)
dev.off()
pdf("venn_GE_UMR_pcPromoter_DE_WGBS.pdf")
plot.new()
grid.draw(venn_GE_UMR_pcPromoter_DE_WGBS)
dev.off()

#' validate WGBS DMRs with MeDIP/MRE
setwd("~/快盘/FetalBrain/WGBS/DMR/valid/")
CortexUMR_MeDIP <- data.frame(UMR = "Cortex_UMRs", Assay = "MeDIP", 
                              Cortex = read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed.Cortex02_UMRs_MeDIP_cortex02.coverage", head = F, as.is = T)$V5,
                              GE = read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed.Cortex02_UMRs_MeDIP_GE02.coverage", head = F, as.is = T)$V5)
CortexUMR_MRE <- data.frame(UMR = "Cortex_UMRs", Assay = "MRE", 
                            Cortex = read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed.Cortex02_UMRs_MRE_cortex02.coverage", head = F, as.is = T)$V5,
                            GE = read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed.Cortex02_UMRs_MRE_GE02.coverage", head = F, as.is = T)$V5)
GEUMR_MeDIP <- data.frame(UMR = "GE_UMRs", Assay = "MeDIP", 
                          Cortex = read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed.GE02_UMRs_MeDIP_cortex02.coverage", head = F, as.is = T)$V5,
                          GE = read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed.GE02_UMRs_MeDIP_GE02.coverage", head = F, as.is = T)$V5)
GEUMR_MRE <- data.frame(UMR = "GE_UMRs", Assay = "MRE", 
                        Cortex = read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed.GE02_UMRs_MRE_cortex02.coverage", head = F, as.is = T)$V5,
                        GE = read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed.GE02_UMRs_MRE_GE02.coverage", head = F, as.is = T)$V5)
valid <- rbind(CortexUMR_MeDIP, CortexUMR_MRE, GEUMR_MeDIP, GEUMR_MRE)
valid <- mutate(valid, Asymmetry = (Cortex - GE)/(Cortex + GE))
(valid_boxplot <- ggplot(valid, aes(x = UMR, y = Asymmetry, color = UMR)) + 
   geom_boxplot(outlier.shape = 20, width = 0.8, size = 1.1) + 
   facet_grid(Assay ~ ., scales = "free") + 
   coord_flip() + 
   xlab("") + 
   ylab("(Cortex - GE)/(Cortex + GE)") + 
   scale_color_manual(values = c("blue", "red")) + 
   theme_bw())
ggsave(valid_boxplot, file = "valid_boxplot.pdf")

#' overlap with TFBSs
setwd("~/快盘/FetalBrain/WGBS/DMR/TF/")
DMR_TF <- read.delim("DMR.Cortex_GE.m0.75.p0.005.d0.5.s300.c3.TF.summary", head = F, as.is = T, col.names = c("TF", "Cortex02UMR", "GE02UMR", "Ratio02", "Cortex04UMR", "GE04UMR", "Ratio04"))
DMR_TF_tall <- data.frame(TF = rep(DMR_TF$TF, 2), Donor = rep(c("HuFNSC02", "HuFNSC04"), each = nrow(DMR_TF)), Cortex_UMR = c(DMR_TF$Cortex02UMR, DMR_TF$Cortex04UMR), GE_UMR = c(DMR_TF$GE02UMR, DMR_TF$GE04UMR), Ratio = c(DMR_TF$Ratio02, DMR_TF$Ratio04))
DMR_TF_tall$Asymmetry <- NA
DMR_TF_tall[DMR_TF_tall$Ratio > 1,]$Asymmetry <- "Cortex enriched"
DMR_TF_tall[DMR_TF_tall$Ratio < 1,]$Asymmetry <- "GE enriched"
DMR_TF_tall[DMR_TF_tall$Ratio == 1,]$Asymmetry <- "Equal"
(DMR_TF_figure <- ggplot(DMR_TF_tall, aes(x = Cortex_UMR, y = GE_UMR, label = TF, color = Asymmetry)) + 
   geom_point() + 
   geom_abline(intercept=0, slope=1) + 
   geom_text(data = subset(DMR_TF_tall, (Donor == "HuFNSC02" & (Ratio > 3 | Ratio <= 1/3)) | (Donor == "HuFNSC04" & (Ratio > 15 | Ratio <= 0.5))), angle = 45, size = 4, hjust = 0.2, vjust = 0.2) + 
   facet_wrap(~ Donor) + 
   scale_x_log10() + 
   scale_y_log10() + 
   scale_color_hue(l = 50) + 
   theme_bw())
ggsave(DMR_TF_figure, file = "DMR_TF_WGBS.pdf", height = 6)

save(DM_test_summary, DMR_test_summary, DMR_WGBS_summary, Cortex02_GE02_DMR_WGBS, Cortex02_GE02_DMR_WGBS_figures, Cortex04_GE04_DMR_WGBS, Cortex04_GE04_DMR_WGBS_figures, DMR_freq_WGBS_figure, DMR_pos_WGBS_figure, 
     GREAT_Cortex02.UMR_WGBS, GREAT_GE02.UMR_WGBS, GREAT_Cortex04.UMR_WGBS, GREAT_GE04.UMR_WGBS, valid_boxplot, genomicBreak_WGBS_figure, DMR_WGBS_gene_summary,
     venn_Cortex_UMR_WGBS, venn_Cortex_UMR_pcGene_WGBS, venn_Cortex_UMR_pcPromoter_WGBS, venn_Cortex_UMR_pcPromoter_DE_WGBS, 
     venn_GE_UMR_WGBS, venn_GE_UMR_pcGene_WGBS, venn_GE_UMR_pcPromoter_WGBS, venn_GE_UMR_pcPromoter_DE_WGBS, 
     Cortex02_UMR_pcPromoter, Cortex02_UMR_pcPromoter_DE, GE02_UMR_pcPromoter, GE02_UMR_pcPromoter_DE, 
     Cortex04_UMR_pcPromoter, Cortex04_UMR_pcPromoter_DE, GE04_UMR_pcPromoter, GE04_UMR_pcPromoter_DE, 
     DMR_TF, DMR_TF_figure, 
     file = "~/快盘/FetalBrain/WGBS/WGBS.DMR.Rdata")
