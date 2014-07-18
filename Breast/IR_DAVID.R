# DAVID enrichment for IR events in lum084 myo084
library(ggplot2)
lum_IR_DAVID <- read.delim("~/快盘/REMC/epiProfile/IR/DAVID/IR_lum_spe_ENSG_DAVID_July2014_tab4Gloria.txt", as.is=T)
myo_IR_DAVID <- read.delim("~/快盘/REMC/epiProfile/IR/DAVID/IR_myo_spe_ENSG_DAVID_July2014_tab4Gloria.txt", as.is=T)
IR_DAVID <- data.frame(Category = c(lum_IR_DAVID$Category, myo_IR_DAVID$Category), Term = c(lum_IR_DAVID$Term, myo_IR_DAVID$Term), Benjamini = c(lum_IR_DAVID$Benjamini, myo_IR_DAVID$Benjamini), CellType = c(rep("lum", nrow(lum_IR_DAVID)), rep("myo", nrow(myo_IR_DAVID))))
(terms <- unique(as.character(IR_DAVID[IR_DAVID$Benjamini <= 1e-11, ]$Term))) # terms with Benjamini <= 0.01 in at least one comparisons
IR_DAVID <- IR_DAVID[IR_DAVID$Term %in% terms, ]
IR_DAVID$Term <- factor(IR_DAVID$Term, levels = unique(IR_DAVID[order(IR_DAVID$Benjamini, decreasing = T),]$Term))
IR_DAVID$CellType <- factor(IR_DAVID$CellType, levels = c("myo", "lum"))
(Enrich_IR_DAVID <- ggplot(data = IR_DAVID, aes(Term, -log10(Benjamini))) +
   geom_bar(aes(fill = CellType), width = .5, position="dodge") + 
   ylab("-log10(Benjamini P-value)") + 
   coord_flip() + 
   scale_fill_manual(values = c("lum" = rgb(200,50,0, maxColorValue = 255), "myo" = rgb(50,200,50, maxColorValue = 255))) + 
   theme(axis.title.y = element_text(size = 0), axis.title.x = element_text(size = 15), axis.text = element_text(size = 15, color = "black"), legend.text = element_text(size = 20), legend.title = element_text(size = 20), panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent")))
ggsave(Enrich_IR_DAVID, file = "~/快盘/REMC/figures/Enrich_IR_DAVID.pdf", height = 3)
