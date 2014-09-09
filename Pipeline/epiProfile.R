epiProfile <- function(mark, cell1, cell2, donor1, donor2, dirIn, dirOut = "", both, neither, cell1_specific, cell2_specific, geneRPKM = F, CpG = F, 
                       color1 = rgb(200,50,0, maxColorValue = 255), color2 = rgb(50,200,50, maxColorValue = 255), color_both = "purple", color_neither = "blue"){
  library(ggplot2)
  library(grid) 
  sample1 = paste0(cell1, "-", donor1); sample2 = paste0(cell2, "-", donor2)
  if(mark == "H3K36me3"){
    cell1_mark_exons <- read.delim(paste0(dirIn, "exons/hg19v65_exons_for_genes.", cell1, donor1, "_", mark, ".coverage"), head = F, as.is = T)
    cell2_mark_exons <- read.delim(paste0(dirIn, "exons/hg19v65_exons_for_genes.", cell2, donor2, "_", mark, ".coverage"), head = F, as.is = T)
    # normalize signal
    (norm <- sum(cell2_mark_exons$V6)/sum(cell1_mark_exons$V6))
    mark_exons <- data.frame(id = c(cell1_mark_exons$V4, cell2_mark_exons$V4), Sample = c(rep(sample1, nrow(cell1_mark_exons)), rep(sample2, nrow(cell2_mark_exons))), Cell_type = c(rep(cell1, nrow(cell1_mark_exons)), rep(cell2, nrow(cell2_mark_exons))), geneRPKM = NA, Expression = NA, mark = c(cell1_mark_exons$V6 * norm, cell2_mark_exons$V6))
    mark_exons[mark_exons$id %in% both, "Expression"] <- "expressed_in_both"
    mark_exons[mark_exons$id %in% neither, "Expression"] <- "not_expressed"
    mark_exons[mark_exons$id %in% cell1_specific, "Expression"] <- paste0(sample1, "-specific")
    mark_exons[mark_exons$id %in% cell2_specific, "Expression"] <- paste0(sample2, "-specific")
    mark_exons$Expression <- factor(mark_exons$Expression, level = c("expressed_in_both", paste0(sample1, "-specific"), paste0(sample2, "-specific"), "not_expressed"))
    mark_exons$Sample <- factor(mark_exons$Sample, level = c(sample1, sample2))
    mark_exons$Cell_type <- factor(mark_exons$Cell_type, level = c(cell1, cell2))
    mark_exons$geneRPKM <- "All" 
    mark_exons$group <- interaction(mark_exons$Sample, mark_exons$Expression)
    if(geneRPKM){
      for(i in 1:length(geneRPKM)){
        mark_exons[mark_exons$id %in% unlist(geneRPKM[i]), "geneRPKM"] <- names(geneRPKM)[i]
      }
      mark_exons$geneRPKM <- factor(mark_exons$geneRPKM)
      mark_exons <- droplevels(na.omit(mark_exons))
      mark_exons$group <- interaction(interaction(mark_exons$Cell_type, mark_exons$Expression), mark_exons$geneRPKM)
    }
    library(plyr)
    mark_exons_stat <- ddply(mark_exons, ~ group, summarize, Sample = Sample[1], Cell_type = Cell_type[1], Expression = Expression[1], geneRPKM = geneRPKM[1], ymin = boxplot.stats(mark)$stats[1], lower = boxplot.stats(mark)$stats[2], middle = mean(mark), upper = boxplot.stats(mark)$stats[4], ymax = boxplot.stats(mark)$stats[5])
    mark_exons_profile <- ggplot(mark_exons_stat, aes(x = Expression, group = group)) + 
       geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax, fill = Sample), stat = "identity", position = "dodge", outlier.shape = NA, width = 0.8) + 
       facet_grid(geneRPKM ~ ., scales = "free") + 
       xlab("Exon group") + 
       ylab(paste0("Average ", mark, " signal")) + 
       scale_fill_manual(values = c(color1, color2)) + 
      theme(panel.border = element_rect(linetype = "solid", fill = "transparent"), panel.margin = unit(0.75, "lines"), axis.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 12, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black"))
    ggsave(mark_exons_profile, file = paste0(dirOut, mark, "_exons_profile_", cell1, "-", donor1, "_", cell2, "-", donor2,".pdf"), width = 9, height = 8)
    return(list(data = mark_exons, profile = mark_exons_stat, figure = mark_exons_profile))
  }
  else{
    cell1_mark_3p <- read.table(paste0(dirIn, "exons3p_200/", cell1, donor1, "_", mark, ".hg19v65_exons_for_genes.3prime_200.unique.profile"), sep = " ", head = F, as.is = T, row.names = 1)
    cell2_mark_3p <- read.table(paste0(dirIn, "exons3p_200/", cell2, donor2, "_", mark, ".hg19v65_exons_for_genes.3prime_200.unique.profile"), sep = " ", head = F, as.is = T, row.names = 1)
    cell1_mark_5p <- read.table(paste0(dirIn, "exons5p_200/", cell1, donor1, "_", mark, ".hg19v65_exons_for_genes.5prime_200.unique.profile"), sep = " ", head = F, as.is = T, row.names = 1)
    cell2_mark_5p <- read.table(paste0(dirIn, "exons5p_200/", cell2, donor2, "_", mark, ".hg19v65_exons_for_genes.5prime_200.unique.profile"), sep = " ", head = F, as.is = T, row.names = 1)
    if(mark != "WGBS"){
      cell1_mark_3p <- sum(cell2_mark_3p)/sum(cell1_mark_3p) * cell1_mark_3p
      cell1_mark_5p <- sum(cell2_mark_5p)/sum(cell1_mark_5p) * cell1_mark_5p
    }
    mark_3p <- data.frame(Sample = rep(c(sample1, sample2), each = 20*4), Cell_type = rep(c(cell1, cell2), each = 20*4), Expression = rep(c(rep(c(paste0(sample1, "-specific"), paste0(sample2, "-specific"), "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), mark = -1)
    mark_3p$group <- interaction(mark_3p$Sample, mark_3p$Expression)
    mark_3p[mark_3p$group == paste0(sample1, ".", sample1, "-specific"), "mark"] <- colMeans(cell1_mark_3p[cell1_specific,], na.rm = T)
    mark_3p[mark_3p$group == paste0(sample1, ".", sample2, "-specific"), "mark"] <- colMeans(cell1_mark_3p[cell2_specific,], na.rm = T)
    mark_3p[mark_3p$group == paste0(sample1, ".expressed_in_both"), "mark"] <- colMeans(cell1_mark_3p[both,], na.rm = T)
    mark_3p[mark_3p$group == paste0(sample1, ".not_expressed"), "mark"] <- colMeans(cell1_mark_3p[neither,], na.rm = T)
    mark_3p[mark_3p$group == paste0(sample2, ".", sample1, "-specific"), "mark"] <- colMeans(cell2_mark_3p[cell1_specific,], na.rm = T)
    mark_3p[mark_3p$group == paste0(sample2, ".", sample2, "-specific"), "mark"] <- colMeans(cell2_mark_3p[cell2_specific,], na.rm = T)
    mark_3p[mark_3p$group == paste0(sample2, ".expressed_in_both"), "mark"] <- colMeans(cell2_mark_3p[both,], na.rm = T)
    mark_3p[mark_3p$group == paste0(sample2, ".not_expressed"), "mark"] <- colMeans(cell2_mark_3p[neither,], na.rm = T)
    mark_5p <- data.frame(Sample = rep(c(sample1, sample2), each = 20*4), Cell_type = rep(c(cell1, cell2), each = 20*4), Expression = rep(c(rep(c(paste0(sample1, "-specific"), paste0(sample2, "-specific"), "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), mark = -1)
    mark_5p$group <- interaction(mark_5p$Sample, mark_5p$Expression)
    mark_5p[mark_5p$group == paste0(sample1, ".", sample1, "-specific"), "mark"] <- colMeans(cell1_mark_5p[cell1_specific,], na.rm = T)
    mark_5p[mark_5p$group == paste0(sample1, ".", sample2, "-specific"), "mark"] <- colMeans(cell1_mark_5p[cell2_specific,], na.rm = T)
    mark_5p[mark_5p$group == paste0(sample1, ".expressed_in_both"), "mark"] <- colMeans(cell1_mark_5p[both,], na.rm = T)
    mark_5p[mark_5p$group == paste0(sample1, ".not_expressed"), "mark"] <- colMeans(cell1_mark_5p[neither,], na.rm = T)
    mark_5p[mark_5p$group == paste0(sample2, ".", sample1, "-specific"), "mark"] <- colMeans(cell2_mark_5p[cell1_specific,], na.rm = T)
    mark_5p[mark_5p$group == paste0(sample2, ".", sample2, "-specific"), "mark"] <- colMeans(cell2_mark_5p[cell2_specific,], na.rm = T)
    mark_5p[mark_5p$group == paste0(sample2, ".expressed_in_both"), "mark"] <- colMeans(cell2_mark_5p[both,], na.rm = T)
    mark_5p[mark_5p$group == paste0(sample2, ".not_expressed"), "mark"] <- colMeans(cell2_mark_5p[neither,], na.rm = T)
    mark_boundaries <- data.frame(rbind(mark_3p, mark_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(mark_5p)), levels = c("5-prime", "3-prime")))
    mark_boundaries$Expression <- factor(mark_boundaries$Expression, levels = c(paste0(sample1, "-specific"), paste0(sample2, "-specific"), "expressed_in_both", "not_expressed"))
    mark_boundaries_profile <- ggplot(mark_boundaries, aes(x = Position, y = mark, group = group)) + 
       geom_line(aes(color = Expression, linetype = Sample)) + 
       geom_point(aes(color = Expression, shape = Sample)) + 
       facet_wrap(~ End) + 
       ylab(paste0("Average ", mark, " signal")) + 
       scale_color_manual(values = c(color1, color2, color_both, color_neither)) + 
       theme(panel.margin = unit(0.75, "lines"), axis.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 12, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black"))
    ggsave(mark_boundaries_profile, file = paste0(dirOut, mark, "_boundaries_profile_", cell1, "-", donor1, "_", cell2, "-", donor2,".pdf"))
    if(CpG){
      CpG_content_3p <- read.table("~/hg19/CpG.hg19v65_exons_for_genes.3prime_200.unique", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
      CpG_content_5p <- read.table("~/hg19/CpG.hg19v65_exons_for_genes.5prime_200.unique", sep = " ", head = F, as.is = T, row.names = 1, fill = T)
      mark_boundaries <- data.frame(data = c(rep(mark, nrow(mark_boundaries)), rep("CpG", ncol(CpG_content_3p)*2)), 
                                    Sample = c(as.character(mark_boundaries$Sample), rep(sample1, ncol(CpG_content_3p)*2)), 
                                    Cell_type = c(as.character(mark_boundaries$Cell_type), rep(cell1, ncol(CpG_content_3p)*2)), 
                                    Expression = c(as.character(mark_boundaries$Expression), rep("CpG density", ncol(CpG_content_3p)*2)), 
                                    Position = c(mark_boundaries$Position, rep(seq(-190, 190, by = ncol(CpG_content_3p)), times = 2)), 
                                    value = c(mark_boundaries$mark, colMeans(CpG_content_3p), colMeans(CpG_content_5p)), 
                                    End = c(as.character(mark_boundaries$End), rep(c("3-prime", "5-prime"), each = ncol(CpG_content_3p))))
      mark_boundaries$data <- factor(mark_boundaries$data, levels = c(mark, "CpG"))
      mark_boundaries$End <- factor(mark_boundaries$End, levels = c("5-prime", "3-prime"))
      mark_boundaries$Expression <- factor(mark_boundaries$Expression, levels = c(paste0(sample1, "-specific"), paste0(sample2, "-specific"), "expressed_in_both", "not_expressed", "CpG density"))
      mark_boundaries$group <- interaction(mark_boundaries$Sample, mark_boundaries$Expression)
      (mark_boundaries_profile <- ggplot(mark_boundaries, aes(x = Position, y = value, group = group)) + 
         geom_line(aes(color = Expression, linetype = Sample)) + 
         geom_point(aes(color = Expression, shape = Sample)) + 
         facet_grid(data ~ End, scales = "free_y") + 
         ylab("") + 
         scale_color_manual(values = c(color1, color2, color_both, color_neither, "black")) + 
         theme(panel.border = element_rect(linetype = "solid", fill = "transparent"), panel.margin = unit(0.75, "lines"), axis.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 12, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
      library(grid) 
      gt <- ggplot_gtable(ggplot_build(mark_boundaries_profile)) 
      # gt$layout 
      gt$heights[[4]] <- unit(3.5, "null") 
      # grid.newpage()
      pdf(paste0(dirOut, mark, "_CpG_profile_", cell1, "-", donor1, "_", cell2, "-", donor2,".pdf"), width = 9)
      grid.draw(gt) 
      grid.text("Average DNA methylation", x = unit(0.015, "npc"), y = unit(0.65, "npc"), rot = 90)
      grid.text("CpG density", x = unit(0.015, "npc"), y = unit(0.15, "npc"), rot = 90)
      dev.off()
    }
    return(list(profile = mark_boundaries, figure = mark_boundaries_profile))
  }
}
# Epigenetic profile for isoform exons
# Required input: 
#   H3K36me3: exon coverage results from RegionsCoverageFromWigCalculator.jar, name <cell><donor>_H3K36me3
#   WGBS: exon boundary profile results from RegionsProfileFromBEDCalculator.jar, name <cell><donor>_WGBS
#   other marks: exon boundary profile results from RegionsProfileFromWigCalculator.jar, name <cell><donor>_mark 
#   isoform results from isoform function between two cell types 
#   exon ID in exon profiles and isoforms should match format 
# Parameters: 
#   mark: name of the epigenetic mark
#   cell1, cell2: name of the two cell types used for epigenetic profile  
#   donor1, donor2: name of the two individuals used for epigenetic profile  
#   dirIn: input directory to exon signal coverage and exon boundary profiles 
#   dirOut: output directory, default to <current working directory> 
#   both, neither, cell1_specific, cell2_specific: exon IDs for exons expressed in both, not expressed, expressed only in cell1, and expressed only in cell2
#   geneRPKM: if to separate different gene RPKM groups, default to False, otherwise, specify list of exon IDs for different gene RPKM groups
#   CpG: logical, if True, add No. of CpG profile track, default to False
# Output: 
#   pdf figure in dirOut
#   return exon info ($data) and summary stats ($profile) and ggplot2 object ($figure) for H3K36me3
#   return exon boundaries profile summary ($profile) and ggplot2 object ($figure) for all other marks

