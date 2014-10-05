#' Compare WGBS and MeDIP DMRs - HuFNSC02 cortex vs GE
setwd("/projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP")
require(ggplot2, lib.loc = "~/bin/R-3.0.2/")
require(labeling, lib.loc = "~/bin/R-3.0.2/")

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

