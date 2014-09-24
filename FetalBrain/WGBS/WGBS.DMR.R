#' DMR analysis from MeDIP fractional calls  

#' DMR identification: ~/HirstLab/FetalBrain/WGBS/WGBS.DM.sh

setwd("~/快盘/FetalBrain/WGBS/")
library(ggplot2)
source('~/HirstLab/Pipeline/DMR.figures.R')
source('~/HirstLab/Pipeline/enrich_GREAT.R')

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
DMR_WGBS_summary <- DMR_test_summary[(DMR_test_summary$size==300 & DMR_test_summary$p=="p=0.005" & DMR_test_summary$delta=="delta=0.5"), ]
Cortex02_GE02_DMR <- read.delim("DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3", as.is = T, head = F)
Cortex02_GE02_DMR_figures <- DMR_figures(Cortex02_GE02_DMR, sample1 = "Cortex-HuFNSC02", sample2 = "GE-HuFNSC02")
(Cortex02_GE02_DMR_figures$length)
(Cortex02_GE02_DMR_figures$count)
(Cortex02_GE02_DMR_figures$dis)
(Cortex02_GE02_DMR_figures$freq)
(Cortex02_GE02_DMR_figures$pos + geom_point(aes(x = (as.numeric(chr) + 0.25*DM), y = pos, color = factor(DM, levels = c("1", "-1"))), position = position_jitter(width = 0.05), size = 1.5, alpha = 0.5))
Cortex04_GE04_DMR <- read.delim("DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3", as.is = T, head = F)
Cortex04_GE04_DMR_figures <- DMR_figures(Cortex04_GE04_DMR, sample1 = "Cortex-HuFNSC04", sample2 = "GE-HuFNSC04")
(Cortex04_GE04_DMR_figures$length)
(Cortex04_GE04_DMR_figures$count)
(Cortex04_GE04_DMR_figures$dis)
(Cortex04_GE04_DMR_figures$freq)
(Cortex04_GE04_DMR_figures$pos + geom_point(aes(x = (as.numeric(chr) + 0.25*DM), y = pos, color = factor(DM, levels = c("1", "-1"))), position = position_jitter(width = 0.05), size = 1.5, alpha = 0.5))

#` GREAT enrichment analysis
(GREAT_Cortex02.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC02_GE-HuFNSC02.hypo", name = "HuFNSC02-Cortex.UMRs", height = 9))
(GREAT_GE02.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC02_GE-HuFNSC02.hyper", name = "HuFNSC02-GE.UMRs"))
(GREAT_Cortex04.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC04_GE-HuFNSC04.hypo", name = "HuFNSC04-Cortex.UMRs", height = 13))
(GREAT_GE04.UMR <- enrich_GREAT(file = "DMR.Cortex-HuFNSC04_GE-HuFNSC04.hyper", name = "HuFNSC04-GE.UMRs", height = 3))

#' genomic break down  
genomicBreak_WGBS <- read.delim("genomic.breakdown.summary", head = F, as.is = T, row.names = 1, col.names = c("Name", "Total", "Intergenic", "Intron", "Exon", "Gene", "Promoter"))
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
ggsave(genomicBreak_WGBS_figure, file = "genomicBreak_WGBS.pdf")


