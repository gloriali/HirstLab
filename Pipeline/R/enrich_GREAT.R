enrich_GREAT <- function(file, name, top = 10, dirIn = paste0(getwd(), "/enrich/"), dirOut = paste0(getwd(), "/enrich/"), categories = NULL, height = 6, width = 9){
  library(ggplot2)
  enrich <- data.frame()
  for(f in list.files(path = dirIn, pattern = paste0("GREAT.*", file))){
    category <- sub("GREAT_", "", f)
    category <- sub("_.*", "", category)
    great <- read.delim(file = paste0(dirIn, "/", f), head = F, as.is = T, skip = 2)
    enrich <- rbind(enrich, data.frame(Category = rep(category, times = min(c(top, nrow(great)))), Term = great$V1[1:min(c(top, nrow(great)))], FDR = great$V4[1:min(c(top, nrow(great)))]))
  }
  if(!is.null(categories)){
    enrich <- enrich[as.character(enrich$Category) %in% categories, ]
  }
  enrich$Term <- as.character(enrich$Term)
  for(i in 1:nrow(enrich)){
    if(nchar(enrich$Term[i]) > 120){
      enrich$Term[i] <- paste0(substr(enrich$Term[i], 1, as.integer(nchar(enrich$Term[i])/3)), "-\n", 
                              substr(enrich$Term[i], as.integer(nchar(enrich$Term[i])/3) + 1, 2*as.integer(nchar(enrich$Term[i])/3)), "-\n", 
                              substr(enrich$Term[i], 2*as.integer(nchar(enrich$Term[i])/3) + 1, nchar(enrich$Term[i])))
    }
    if(nchar(enrich$Term[i]) > 60 & nchar(enrich$Term[i]) <= 120){
      enrich$Term[i] <- paste0(substr(enrich$Term[i], 1, as.integer(nchar(enrich$Term[i])/2)), "-\n", substr(enrich$Term[i], as.integer(nchar(enrich$Term[i])/2)+1, nchar(enrich$Term[i])))
    }
  }
  enrich$Term <- factor(enrich$Term, levels = enrich[order(enrich$Category, enrich$FDR, decreasing = T),]$Term)
  Enrich_plot <- ggplot(data = enrich, aes(Term, -log10(FDR))) +
    geom_bar(aes(fill = Category), stat = "identity", width = .5) + 
    coord_flip() + 
    geom_text(aes(label = round(-log10(FDR), 2), hjust = 0)) + 
    xlab("") + 
    theme_bw() +
    ggtitle(paste0("GREAT enrichment for ", name)) + 
    scale_fill_hue(l = 40)
  ggsave(Enrich_plot, file = paste0(dirOut, name, "_enrich.pdf"), height = height, width = width)
  return(Enrich_plot)
}
# Plot enrichment terms from GREAT
# Criteria: 
#   GREAT result: GREAT_<Category>_<file>.tsv
# Parameters: 
#   file: input file name
#   name: output name prefix 
#   [top]: top terms only for plotting, default to top 10 terms   
#   [dirIn]: directory with ermineJ and DAVID results, default to `<current wd>/enrich/` 
#   [dirOut]: output directory, default to `<current wd>/enrich/`
#   [categories]: categories of terms to include, default to `NULL`
#   [height], [width]: for output pdf plot, default to 6*9
# Output: 
#   pdf figure with enrichment plot
#   return ggplot2 object for enrichment plot




