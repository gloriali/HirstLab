DE_FC <- function(data, name = paste(colnames(data), collapse = "_"), dirIn = getwd(), dirOut = paste0(dirIn, "/DE/"), FC_cut = 1, log_cut = -1){
  library(dplyr)
  library(ggplot2)
  library(grid)
  e <- 1e-6
  system(paste0("mkdir -p ", dirOut))
  print(paste("Processing", name))
  data <- data %>% mutate(geneID = rownames(data), logFC = log2((data[,1] + e) / (data[,2] + e)), Ave = (log10(data[,1] + e) + log10(data[,2] + e))/2, col = abs(logFC) >= FC_cut & Ave >= log_cut) 
  DE <- data %>% filter(col) %>% arrange(-abs(logFC)) 
  DE <- DE[, c("geneID", colnames(data)[1:2], "logFC")]
  write.table(DE, file = paste0(dirOut, "DE_", name, "_gene.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  UP <- DE %>% filter(logFC > 0)
  DN <- DE %>% filter(logFC < 0)
  DE_figure <- ggplot(data, aes(x = Ave, y = logFC, color = col)) + 
     geom_point() + 
     xlab("Average log10 expression level") + 
     ylab("log2 fold change") + 
     ggtitle(name) + 
     scale_color_manual(values = c("black", "red")) + 
     theme_bw() + 
     theme(legend.position = "none")
  gt_DE_figure <- ggplot_gtable(ggplot_build(DE_figure)) 
  grid.newpage()
  grid.draw(gt_DE_figure) 
  grid.text(paste("UP", nrow(UP)), x = unit(0.9, "npc"), y = unit(0.9, "npc"))
  grid.text(paste("DN", nrow(DN)), x = unit(0.9, "npc"), y = unit(0.87, "npc"))
  pdf(paste0(dirOut, "DE_", name, "_figure.pdf"))
  grid.newpage()
  grid.draw(gt_DE_figure) 
  grid.text(paste("UP", DE %>% filter(logFC > 0) %>% nrow), x = unit(0.9, "npc"), y = unit(0.9, "npc"))
  grid.text(paste("DN", DE %>% filter(logFC < 0) %>% nrow), x = unit(0.9, "npc"), y = unit(0.87, "npc"))
  dev.off()
  summary <- data.frame(Sample = name, UP = DE %>% filter(logFC > 0) %>% nrow, DN = DE %>% filter(logFC < 0) %>% nrow, DE = nrow(DE))
  return(list(summary = summary, UP = UP, DN = DN))
}
# Differential expression based on fold change
## Parameters: 
### data: data frame of expression values for two samples, rownames are gene IDs, colnames are sample IDs.   
### name: output name, default to `paste(colnames(data), sep = "_")`.     
### dirIn: input directory, default to `getwd()`.      
### dirOut: output directory, default to `paste0(dirIn, "/DE/")`.    
### FC_cut: cutoff for log2(fold change), default to 1 (2-fold).  
### log_cut: cutoff for average log10 expression level, default to -1 (0.1).  
## Output: 
### Figure: MA-plot pdf figure, file name `dirOut, "DE_", name, "_figure.pdf"`.   
### DE gene file: gene ID, expression level of two samples, and logFC for DE genes, file name `dirOut, "DE_", name, "_gene.txt"`.  
### return list: 
#### summary: name, UP, DN, DE
#### UP gene list: gene ID, expression level of two samples, and logFC for UP genes
#### DN gene list: gene ID, expression level of two samples, and logFC for DN genes

