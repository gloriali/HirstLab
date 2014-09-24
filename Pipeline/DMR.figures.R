DMR_figures <- function(DMR, sample1, sample2, dirOut = getwd(), width = 9, height = 7, figures = c("length", "count", "adjacentDis", "frequency", "position"), chrs = c(paste0("chr", as.character(1:22)), "chrX"), colname = c("chr", "start", "end", "ID", "DM", "count", "length")){
  library(ggplot2)
  library(plyr)
  DMR_length_figure <- NULL; DMR_count_figure <- NULL; DMR_dis_figure <- NULL; DMR_freq_figure <- NULL; DMR_position_figure <- NULL;
  chrlength <- read.csv("~/快盘/hg19/chrlen_hg19.csv", as.is = T, row.names = 1)
  chrlength <- chrlength[chrlength$chr %in% chrs, ]
  chrlength$chr <- factor(chrlength$chr, levels = chrs[1:length(chrs)])
  colnames(DMR) <- colname
  DMR <- DMR[DMR$chr %in% chrs, ]
  DMR$chr <- factor(DMR$chr, levels = chrs)
  DMR <- DMR[order(DMR$chr, DMR$start, decreasing = F), ]
  DMR$pos <- (DMR$start + DMR$end)/2
  DMR$dis <- c(-100000, DMR[2:nrow(DMR), ]$pos - DMR[1:nrow(DMR)-1,]$pos)
  DMR[DMR$dis < 0, ]$dis <- -100000
  if("length" %in% figures){
    DMR_length_stat <- ddply(DMR, .(chr, DM), summarize, ymin = boxplot.stats(length)$stats[1], lower = boxplot.stats(length)$stats[2], middle = boxplot.stats(length)$stats[3], upper = boxplot.stats(length)$stats[4], ymax = boxplot.stats(length)$stats[5])
    DMR_length_figure <- ggplot(DMR_length_stat, aes(x = chr, fill = chr)) + 
      geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax), stat = "identity", outlier.shape = NA, width = 0.8) + 
      xlab("") + 
      ylab("DMR length (bp)") + 
      ggtitle(paste("DMR length -", sample1, "vs", sample2, "DMRs")) + 
      guides(fill = F) + 
      facet_wrap(~ DM) + 
      theme_bw()
    ggsave(DMR_length_figure, file = paste0(dirOut, "/DMRlength_", sample1, "_", sample2, ".pdf"), width = width, height = height)
  }
  if("count" %in% figures){
    DMR_count_stat <- ddply(DMR, .(chr, DM), summarize, ymin = boxplot.stats(count)$stats[1], lower = boxplot.stats(count)$stats[2], middle = boxplot.stats(count)$stats[3], upper = boxplot.stats(count)$stats[4], ymax = boxplot.stats(count)$stats[5])
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
    DMR_dis_stat <- ddply(DMR, .(chr, DM), summarize, ymin = boxplot.stats(dis)$stats[1], lower = boxplot.stats(dis)$stats[2], middle = boxplot.stats(dis)$stats[3], upper = boxplot.stats(dis)$stats[4], ymax = boxplot.stats(dis)$stats[5])
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
    DMR_freq <- ddply(DMR, .(chr, DM), summarize, DMRlen = sum(DM*length))
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
    DMR$chr <- factor(DMR$chr, levels = chrs[length(chrs):1])
    DMR_position_figure <- ggplot(DMR) + 
      geom_linerange(aes(x = factor(chr, levels = chr[length(chr):1]), ymin = 0, ymax = length), data = chrlength, alpha = 0.5) + 
      geom_point(aes(x = (as.numeric(chr) + 0.25*DM), y = pos, color = factor(DM, levels = c("1", "-1"))), position = position_jitter(width = 0.05), size = 0.5, alpha = 0.5) +  
      xlab("") + 
      ylab("Position of DMRs on the chromosome") +
      coord_flip() + 
      ggtitle(paste("DMR position -", sample1, "vs", sample2, "DMRs")) + 
      scale_color_manual(values = c("red", "blue"), labels = c("hyper", "hypo"), name = "") + 
      theme_bw()
    ggsave(DMR_position_figure, file = paste0(dirOut, "/DMRpos_", sample1, "_", sample2, ".pdf"), width = width, height = height)
  }
  return(list(length = DMR_length_figure, count = DMR_count_figure, dis = DMR_dis_figure, freq = DMR_freq_figure, pos = DMR_position_figure))
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
# Output:
#   pdf figures in `dirOut/<figure>_<sample1>_<sample2>.pdf`   
#   return list of ggplot2 objects: unprocessed ones return `NULL`
#     $length: DMR length boxplot; 
#     $count: No. of CpGs per DMR; 
#     $dis: distance between adjacent DMRs; 
#     $freq: DMR frequency (bp/MB); 
#     $pos: position of DMRs on the chromosomes.   
