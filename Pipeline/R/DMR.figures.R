DMR_figures <- function(DMR, sample1, sample2, dirOut = getwd(), width = 8, height = 8, figures = c("length", "count", "adjacentDis", "frequency", "position", "circos"), chrs = c(paste0("chr", as.character(1:22)), "chrX"), colname = c("chr", "start", "end", "ID", "DM", "count", "length"), chrlength = read.csv("~/hg19/chrlen_hg19.csv", as.is = T, row.names = 1), hist_width = 100){
  library(ggplot2)
  library(dplyr)
  library(ggbio)
  library(GenomicRanges)
  library(RCircos)
  DMR_length_figure <- NULL; DMR_count_figure <- NULL; DMR_dis_figure <- NULL; DMR_freq_figure <- NULL; DMR_pos_figure <- NULL; 
  chrlength <- chrlength[chrlength$chr %in% chrs, ]
  chrlength$chr <- factor(chrlength$chr, levels = chrs[1:length(chrs)])
  colnames(DMR) <- colname
  DMR <- DMR[DMR$chr %in% chrs, ]
  DMR$chr <- factor(DMR$chr, levels = chrs)
  DMR <- DMR[order(DMR$chr, DMR$start, decreasing = F), ]
  DMR$pos <- (DMR$start + DMR$end)/2
  if("length" %in% figures){
    DMR_length_stat <- mutate(summarise(group_by(DMR, chr, DM), lower = quantile(length, 0.25), middle = as.numeric(median(length)), upper = quantile(length, 0.75), min = min(length), max = max(length)), ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)))
    DMR_length_figure <- ggplot(DMR_length_stat, aes(x = chr, fill = chr)) + 
      geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
      xlab("") + 
      ylab("DMR length (bp)") + 
      ggtitle(paste("DMR length -", sample1, "vs", sample2, "DMRs")) + 
      guides(fill = F) + 
      facet_wrap(~ DM) + 
      theme_bw() + 
    	theme(axis.text.x = element_text(angle = 90))
    ggsave(DMR_length_figure, file = paste0(dirOut, "/DMRlength_", sample1, "_", sample2, ".pdf"), width = width, height = height)
  }
  if("count" %in% figures){
    DMR_count_stat <- mutate(summarise(group_by(DMR, chr, DM), lower = quantile(count, 0.25), middle = median(count), upper = quantile(count, 0.75), min = min(count), max = max(count)), ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)))
    DMR_count_figure <- ggplot(DMR_count_stat, aes(x = chr, fill = chr)) + 
      geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
      xlab("") + 
      ylab("No. of CpGs per DMR") + 
      ggtitle(paste("No. of CpGs per DMR -", sample1, "vs", sample2, "DMRs")) + 
      guides(fill = F) + 
      facet_wrap(~ DM) + 
      theme_bw()
    ggsave(DMR_count_figure, file = paste0(dirOut, "/CpGcount_", sample1, "_", sample2, ".pdf"), width = width, height = height)
  }
  if("adjacentDis" %in% figures){
  	DMR$dis <- c(-100000, DMR[2:nrow(DMR), ]$pos - DMR[1:nrow(DMR)-1,]$pos)
  	DMR[DMR$dis < 0, ]$dis <- -100000
  	DMR_dis_stat <- mutate(summarise(group_by(DMR, chr, DM), lower = quantile(dis, 0.25), middle = median(dis), upper = quantile(dis, 0.75), min = min(dis), max = max(dis)), ymin = max(min, lower - 1.5*(upper-lower)), ymax = min(max, upper + 1.5*(upper-lower)))
    DMR_dis_figure <- ggplot(DMR_dis_stat, aes(x = chr, fill = chr)) + 
      geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
      xlab("") + 
      ylab("Distance between adjacent DMRs (bp)") + 
      ggtitle(paste("Distance between DMRs -", sample1, "vs", sample2, "DMRs")) + 
      guides(fill = F) + 
      facet_wrap(~ DM) + 
      theme_bw()
    ggsave(DMR_dis_figure, file = paste0(dirOut, "/DMRdis_", sample1, "_", sample2, ".pdf"), width = width, height = height)
  }
  if("frequency" %in% figures){
    DMR_freq <- summarise(group_by(DMR, chr, DM), DMRlen = sum(DM*length))
    DMR_freq$chrlen <- chrlength[DMR_freq$chr, ]$length
    DMR_freq$freq <- DMR_freq$DMRlen / DMR_freq$chrlen * 10^6 
    DMR_freq$DM <- factor(DMR_freq$DM, levels = c("1", "-1"))
    DMR_freq$chr <- factor(DMR_freq$chr, levels = chrs[length(chrs):1])
    DMR_freq_figure <- ggplot(DMR_freq, aes(x = chr, y = freq, fill = DM)) + 
      geom_bar(position = "identity", stat = "identity", width = 0.8) + 
      xlab("") + 
      ylab("DMR frequency (bp/MB)") + 
      coord_flip() + 
      ggtitle(paste("DMR frequency -", sample1, "vs", sample2, "DMRs")) + 
      scale_fill_manual(values = c("red", "blue"), labels = c("hyper", "hypo"), name = "") + 
      theme_bw()
    ggsave(DMR_freq_figure, file = paste0(dirOut, "/DMRfreq_", sample1, "_", sample2, ".pdf"), width = width, height = height)
  }
  if("position" %in% figures){
    DMR$y <- DMR$DM * 15 + 5
    DMR_gr <- keepSeqlevels(as(DMR, "GRanges"), paste0("chr", c(1:22, "X")))
    data(hg19IdeogramCyto, package = "biovizBase")
    hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X")))
    DMR_pos_figure <- ggplot(hg19) + 
       layout_karyogram(cytoband = TRUE) + 
       layout_karyogram(data = DMR_gr, geom = "point", aes(y = y, x = pos, color = factor(DM, levels = c("1", "-1"))), position = position_jitter(height = 2), size = 1.5, alpha = 0.5) + 
       ylab("") + 
       xlab("Position of DMRs on the chromosome") +
       coord_cartesian(ylim = c(-15, 25)) + 
       ggtitle(paste("DMR positions -", sample1, "vs", sample2, "DMRs")) + 
       scale_color_manual(values = c("red", "blue"), labels = c("hyper", "hypo"), name = "") + 
       theme_bw() + 
       theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_line(color = "transparent"), panel.grid.minor.y = element_line(color = "transparent"))
    ggsave(DMR_pos_figure@ggplot, file = paste0(dirOut, "/DMRpos_", sample1, "_", sample2, ".pdf"), width = width, height = height)    
  }
  if("circos" %in% figures){
  	DMR_hyper <- DMR %>% filter(DM == 1) %>% select(chr, start, end, ID) %>% mutate(data = 1, PlotColors = "red")
  	colnames(DMR_hyper) <- c("Chromosome", "chromStart", "chromEnd", "ID", "data", "PlotColors")
  	DMR_hypo <- DMR %>% filter(DM == -1) %>% select(chr, start, end, ID) %>% mutate(data = 1, PlotColors = "blue")
  	colnames(DMR_hypo) <- c("Chromosome", "chromStart", "chromEnd", "ID", "data", "PlotColors")
  	data(UCSC.HG19.Human.CytoBandIdeogram);
  	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
  	chr.exclude <- "chrY";
  	num.inside <- 10;
  	num.outside <- 0;
  	pdf(paste0(dirOut, "/DMRcircos_", sample1, "_", sample2, ".pdf"), width = width, height = height)
  	RCircos.Set.Core.Components(cyto.info, chr.exclude, num.inside, num.outside)
  	params <- RCircos.Get.Plot.Parameters()
  	params$hist.width <- hist_width
  	RCircos.Reset.Plot.Parameters(params)
  	RCircos.Set.Plot.Area()
  	RCircos.Chromosome.Ideogram.Plot()
  	RCircos.Histogram.Plot(DMR_hyper, data.col=5, track.num=1, side="in")
  	RCircos.Histogram.Plot(DMR_hypo, data.col=5, track.num=2, side="in")
  	legend("bottomright", c("hyper", "hypo"), col = c("red", "blue"), lwd = 8, cex = 0.8)
  	grid.text(paste0(sample1, "_", sample2), x = unit(0.5, "npc"), y = unit(0.95, "npc"), gp=gpar(fontsize=12))
  	dev.off()
  }
  return(list(length = DMR_length_figure, count = DMR_count_figure, dis = DMR_dis_figure, freq = DMR_freq_figure, pos = DMR_pos_figure))
}
# DMR analysis and visualization  
# Required input: 
#   DMR files with chr, start, end, [id], hyper(1)/hypo(-1), CpG count, length
# Parameters: 
#   DMR: data frame with DMRs   
#   sample1, sample2: sample names   
#   [dirOut]: output directory, default to current working directories       
#   [width, height]: dimensions of output figures, default to `wdith = 9, height = 7` 
#   [figures]: figures to analysis and output, default to all, `length: DMR length boxplot; count: No. of CpGs per DMR; adjacentDis: distance between adjacent DMRs; frequency: DMR frequency (bp/MB); position: position of DMRs on the chromosomes` 
#   [chrs]: chromosomes to use, default to `chr1-22, chrX`
#   [colname]: colnames of input DMR data frame, default to `c("chr", "start", "end", "ID", "DM", "count", "length")`
#   [hist_width]: histogram width for R circos plot, default to 100. 
# Output:
#   pdf figures in `dirOut/<figure>_<sample1>_<sample2>.pdf`   
#   return list of ggplot2 objects: unprocessed ones return `NULL`
#     $length: DMR length boxplot; 
#     $count: No. of CpGs per DMR; 
#     $dis: distance between adjacent DMRs; 
#     $freq: DMR frequency (bp/MB); 
#     $pos: position of DMRs on the chromosomes, DMR_pos_figure@ggplot 
#		R circos plot: pdf file only. 
