# DMR analysis from MeDIP fractional calls  

# DMR identification:  
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.signal.sh
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.fractional.m
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.DM.sh

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
# /home/lli/bin/shell/DMR.intersect.sh -d /projects/epigenomics/lli/FetalBrain/MeDIP/DMR
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

save(DMR_length_MZ_figure, DMR_count_MZ_figure, DMR_dis_MZ_figure, DMR_freq_MZ_figure, DMR_pos_MZ_figure, 
     DMR_MZ_summary, genomicBreak_MZ, genomicBreak_MZ_figure, 
     file = "DMR_MZ.Rdata")


