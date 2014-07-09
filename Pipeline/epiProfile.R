epiProfile <- function(mark, cell1, cell2, donor1, donor2, dirIn, dirOut = "", both, neither, cell1_specific, cell2_specific, geneRPKM){
  library(ggplot2)
  if(mark == "H3K36me3"){
    cell1_mark_exons <- read.delim(paste0(dirIn, "exons/hg19v65_exons_for_genes.", cell1, donor1, "_", mark, ".coverage"), head = F, as.is = T)
    cell2_mark_exons <- read.delim(paste0(dirIn, "exons/hg19v65_exons_for_genes.", cell2, donor2, "_", mark, ".coverage"), head = F, as.is = T)
    # normalize signal
    (norm <- sum(cell2_mark_exons$V6)/sum(cell1_mark_exons$V6))
    mark_exons <- data.frame(id = c(cell1_mark_exons$V4, cell2_mark_exons$V4), Cell_type = c(rep(cell1, nrow(cell1_mark_exons)), rep(cell2, nrow(cell2_mark_exons))), geneRPKM = NA, Expression = NA, mark = c(cell1_mark_exons$V6 * norm, cell2_mark_exons$V6))
    mark_exons[mark_exons$id %in% both, "Expression"] <- "expressed_in_both"
    mark_exons[mark_exons$id %in% neither, "Expression"] <- "not_expressed"
    mark_exons[mark_exons$id %in% cell1_specific, "Expression"] <- paste0(cell1, "-specific")
    mark_exons[mark_exons$id %in% cell2_specific, "Expression"] <- paste0(cell2, "-specific")
    mark_exons$Expression <- factor(mark_exons$Expression)
    for(i in 1:length(geneRPKM)){
      mark_exons[mark_exons$id %in% unlist(geneRPKM[i]), "geneRPKM"] <- names(geneRPKM)[i]
    }
    mark_exons$geneRPKM <- factor(mark_exons$geneRPKM)
    mark_exons <- droplevels(na.omit(mark_exons))
    mark_exons$group <- interaction(interaction(mark_exons$Cell_type, mark_exons$Expression), mark_exons$geneRPKM)
    library(plyr)
    mark_exons_stat <- ddply(mark_exons, ~ group, summarize, Cell_type = Cell_type[1], Expression = Expression[1], geneRPKM = geneRPKM[1], ymin = boxplot.stats(mark)$stats[1], lower = boxplot.stats(mark)$stats[2], middle = mean(mark), upper = boxplot.stats(mark)$stats[4], ymax = boxplot.stats(mark)$stats[5])
    (mark_exons_profile <- ggplot(mark_exons_stat, aes(x = Expression, group = group)) + 
       geom_boxplot(aes(lower = lower, middle = middle, upper = upper, ymin = ymin, ymax = ymax, fill = Cell_type), stat = "identity", position = "dodge", outlier.shape = NA, width = 0.8) + 
       facet_grid(geneRPKM ~ ., scales = "free") + 
       xlab("Exon group") + 
       ylab(paste0("Average ", mark, " signal")) + 
       theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 15, color = "black"), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent"), strip.text = element_text(color = "black", size = 20, hjust = 0.5, vjust = 0.5), strip.background = element_rect(color = "black")))
    ggsave(mark_exons_profile, file = paste0(dirOut, mark, "_exons_profile.pdf"), width = 9, height = 8)
    return(list(exons = mark_exons, exons_stat = mark_exons_stat))
  }
  else{
    cell1_mark_3p <- read.table(paste0(dirIn, "exons3p_200/", cell1, donor1, "_", mark, ".hg19v65_exons_for_genes.3prime_200.unique.profile"), sep = " ", head = F, as.is = T, row.names = 1)
    cell2_mark_3p <- read.table(paste0(dirIn, "exons3p_200/", cell2, donor2, "_", mark, ".hg19v65_exons_for_genes.3prime_200.unique.profile"), sep = " ", head = F, as.is = T, row.names = 1)
    cell1_mark_5p <- read.table(paste0(dirIn, "exons5p_200/", cell1, donor1, "_", mark, ".hg19v65_exons_for_genes.5prime_200.unique.profile"), sep = " ", head = F, as.is = T, row.names = 1)
    cell2_mark_5p <- read.table(paste0(dirIn, "exons5p_200/", cell2, donor2, "_", mark, ".hg19v65_exons_for_genes.5prime_200.unique.profile"), sep = " ", head = F, as.is = T, row.names = 1)
    cell1_mark_3p <- sum(cell2_mark_3p)/sum(cell1_mark_3p) * cell1_mark_3p
    cell1_mark_5p <- sum(cell2_mark_5p)/sum(cell1_mark_5p) * cell1_mark_5p
    mark_3p <- data.frame(Cell_type = rep(c(cell1, cell2), each = 20*4), Expression = rep(c(rep(c(paste0(cell1, "-specific"), paste0(cell2, "-specific"), "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), mark = -1)
    mark_3p$group <- interaction(mark_3p$Cell_type, mark_3p$Expression)
    mark_3p[mark_3p$group == paste0(cell1, ".", cell1, "-specific"), "mark"] <- colMeans(cell1_mark_3p[cell1_specific,], na.rm = T)
    mark_3p[mark_3p$group == paste0(cell1, ".", cell2, "-specific"), "mark"] <- colMeans(cell1_mark_3p[cell2_specific,], na.rm = T)
    mark_3p[mark_3p$group == paste0(cell1, ".expressed_in_both"), "mark"] <- colMeans(cell1_mark_3p[both,], na.rm = T)
    mark_3p[mark_3p$group == paste0(cell1, ".not_expressed"), "mark"] <- colMeans(cell1_mark_3p[neither,], na.rm = T)
    mark_3p[mark_3p$group == paste0(cell2, ".", cell1, "-specific"), "mark"] <- colMeans(cell2_mark_3p[cell1_specific,], na.rm = T)
    mark_3p[mark_3p$group == paste0(cell2, ".", cell2, "-specific"), "mark"] <- colMeans(cell2_mark_3p[cell2_specific,], na.rm = T)
    mark_3p[mark_3p$group == paste0(cell2, ".expressed_in_both"), "mark"] <- colMeans(cell2_mark_3p[both,], na.rm = T)
    mark_3p[mark_3p$group == paste0(cell2, ".not_expressed"), "mark"] <- colMeans(cell2_mark_3p[neither,], na.rm = T)
    mark_5p <- data.frame(Cell_type = rep(c(cell1, cell2), each = 20*4), Expression = rep(c(rep(c(paste0(cell1, "-specific"), paste0(cell2, "-specific"), "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), mark = -1)
    mark_5p$group <- interaction(mark_5p$Cell_type, mark_5p$Expression)
    mark_5p[mark_5p$group == paste0(cell1, ".", cell1, "-specific"), "mark"] <- colMeans(cell1_mark_5p[cell1_specific,], na.rm = T)
    mark_5p[mark_5p$group == paste0(cell1, ".", cell2, "-specific"), "mark"] <- colMeans(cell1_mark_5p[cell2_specific,], na.rm = T)
    mark_5p[mark_5p$group == paste0(cell1, ".expressed_in_both"), "mark"] <- colMeans(cell1_mark_5p[both,], na.rm = T)
    mark_5p[mark_5p$group == paste0(cell1, ".not_expressed"), "mark"] <- colMeans(cell1_mark_5p[neither,], na.rm = T)
    mark_5p[mark_5p$group == paste0(cell2, ".", cell1, "-specific"), "mark"] <- colMeans(cell2_mark_5p[cell1_specific,], na.rm = T)
    mark_5p[mark_5p$group == paste0(cell2, ".", cell2, "-specific"), "mark"] <- colMeans(cell2_mark_5p[cell2_specific,], na.rm = T)
    mark_5p[mark_5p$group == paste0(cell2, ".expressed_in_both"), "mark"] <- colMeans(cell2_mark_5p[both,], na.rm = T)
    mark_5p[mark_5p$group == paste0(cell2, ".not_expressed"), "mark"] <- colMeans(cell2_mark_5p[neither,], na.rm = T)
    mark_boundaries <- data.frame(rbind(mark_3p, mark_5p), End = factor(rep(c("3-prime", "5-prime"), each = nrow(mark_5p)), levels = c("5-prime", "3-prime")))
    (mark_boundaries_profile <- ggplot(mark_boundaries, aes(x = Position, y = mark, group = group)) + 
       geom_line(aes(color = Expression, linetype = Cell_type)) + 
       geom_point(aes(color = Expression, shape = Cell_type)) + 
       facet_wrap(~ End) + 
       ylab(paste0("Average ", mark, " signal")) + 
       theme_bw())
    ggsave(mark_boundaries_profile, file = paste0(dirOut, mark, "_boundaries_profile.pdf"))
    return(mark_boundaries)
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
#   dirOut: output directory, default to current working directory
#   both, neither, cell1_specific, cell2_specific: exon IDs for exons expressed in both, not expressed, expressed only in cell1, and expressed only in cell2
#   geneRPKM: required for H3K36me3, list of exon IDs for different gene RPKM groups
# Output: 
#   pdf figure in dirOut
#   return exon info and summary stats for H3K36me3
#   return exon boundaries profile summary for all other marks

