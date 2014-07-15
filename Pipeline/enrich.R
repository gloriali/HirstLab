enrich <- function(name, dirIn = "./enrich/", dirOut = "./enrich/", fdr = 0.01, p = "FDR", erminej = T, category = c("GOBP", "GOMF", "KEGG_PATHWAY", "PANTHER_PATHWAY", "REACTOME_PATHWAY", "INTERPRO"), height = 6, width = 9){
  library(ggplot2)
  DAVID <- read.delim(paste0(dirIn, name, "_DAVID.txt"), as.is=T)
  if(erminej){
    erminej_BP <- read.delim(paste0(dirIn, name, "_erminej_BP.txt"), as.is=T, skip = 34, head = F)
    if(sum(is.na(erminej_BP$V10)) == nrow(erminej_BP)){erminej_BP$V10 <- erminej_BP$V8}
    erminej_MF <- read.delim(paste0(dirIn, name, "_erminej_MF.txt"), as.is=T, skip = 34, head = F)
    if(sum(is.na(erminej_MF$V10)) == nrow(erminej_MF)){erminej_MF$V10 <- erminej_MF$V8}
    enrich <- na.omit(rbind(data.frame(Category = "GOBP", Term = erminej_BP$V2, FDR = erminej_BP$V10), data.frame(Category = "GOMF", Term = erminej_MF$V2, FDR = erminej_MF$V10), data.frame(Category = DAVID$Category, Term = DAVID$Term, FDR = DAVID[, p])))
  }
  else{
    DAVID[grep("GOTERM_MF", DAVID$Category), "Category"] <- "GOMF"
    DAVID[grep("GOTERM_BP", DAVID$Category), "Category"] <- "GOBP"
    enrich <- na.omit(data.frame(Category = DAVID$Category, Term = DAVID$Term, FDR = DAVID[, p]))
  }
  enrich <- enrich[enrich$FDR <= fdr, ]
  enrich <- enrich[as.character(enrich$Category) %in% category, ]
  enrich$Term <- gsub("IPR[0-9]+:", "", enrich$Term)
  enrich$Term <- gsub("hsa[0-9]+:", "", enrich$Term)
  enrich$Term <- gsub("P[0-9]+:", "", enrich$Term)
  enrich$Term <- gsub("REACT_[0-9]+:", "", enrich$Term)
  enrich$Term <- gsub("GO:[0-9]+~", "", enrich$Term)
  enrich$Term <- strtrim(enrich$Term, 50)
  enrich$Term <- factor(enrich$Term, levels = enrich[order(enrich$Category, enrich$FDR, decreasing = T),]$Term)
  Enrich_plot <- ggplot(data = enrich, aes(Term, -log10(FDR))) +
     geom_bar(aes(fill = Category), width = .5) + 
     coord_flip() + 
     geom_text(aes(label = round(-log10(FDR), 2), hjust = 0)) + 
     theme_bw() +
     ggtitle(paste0("Functional enrichment for ", name)) + 
     scale_fill_manual(values = c("GOBP" = "blue", "GOMF" = "purple", "KEGG_PATHWAY" = "darkgreen", "PANTHER_PATHWAY" = "darkblue", "REACTOME_PATHWAY" = "steelblue", "INTERPRO" = "lightblue"))
  ggsave(Enrich_plot, file = paste0(dirOut, name, "_enrich.pdf"), height = height, width = width)
  return(list(enrich = enrich, figure = Enrich_plot))
}
# Plot enrichment terms from ermineJ (GO) and DAVID (KEGG, PANTHER, REACTOME, INTERPRO)
# Criteria: 
#   DAVID (KEGG, PANTHER, REACTOME, INTERPRO, [GOBP], [GOMF]) result: <name>_DAVID.txt
#   [ermineJ (GO) results]: optional, <name>_erminej_BP.txt and <name>_erminej_MF.txt
# Parameters: 
#   name: file name prefix 
#   [dirIn]: directory with ermineJ and DAVID results, default to <current wd>/enrich/ 
#   [dirOut]: output directory, default to <current wd>/enrich/
#   [fdr]: FDR cutoff for enriched terms, default to 0.01
#   [p]: multi-test p-value correction method (Bonferroni, Benjamini or FDR), default to FDR
#   [erminej]: logical, whether erminej is used for GO, default to True 
#   [category]: categories of terms to include, default to `c("GOBP", "GOMF", "KEGG_PATHWAY", "PANTHER_PATHWAY", "REACTOME_PATHWAY", "INTERPRO")`
#   [height], [width]: for output pdf plot, default to 6*9
# Output: 
#   pdf figure with enrichment plot
#   return list
#     $enrich: data frame of enriched terms, category and corrected p-values 
#     $figure: ggplot2 object for enrichment plot


