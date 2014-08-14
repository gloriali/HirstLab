# Figure S17
setwd("~/快盘/Publications/breast/revision/sup/FigureS17")
lum.GOBP <- read.delim("shown-GOBiologicalProcess.lum.UMR.tsv", head = T, as.is = T)
lum.MSig <- read.delim("shown-MSigDBGeneSetsPerturbation.lum.UMR.tsv", head = T, as.is = T)
myo.GOBP <- read.delim("shown-GOBiologicalProcess.myo.UMR.tsv", head = T, as.is = T)
myo.pathway <- read.delim("shown-pathwayCommons.myo.UMR.tsv", head = T, as.is = T)
GREAT <- data.frame(cell = rep(c("lum.UMR", "myo.UMR"), each = 20), Category = rep(c("GOBP", "MSigPerturbation", "GOBP", "PathwaysCommon"), each = 10), 
                    Term = c(lum.GOBP$Term.Name[1:10], lum.MSig$Term.Name[1:10], myo.GOBP$Term.Name[1:10], myo.pathway$Term.Name[1:10]), 
                    FDR = c(lum.GOBP$Binom.FDR.Q.Val[1:10], lum.MSig$Binom.FDR.Q.Val[1:10], myo.GOBP$Binom.FDR.Q.Val[1:10], myo.pathway$Binom.FDR.Q.Val[1:10]))
GREAT$Term <- as.character(GREAT$Term)
for(i in 1:nrow(GREAT)){
  if(nchar(GREAT$Term[i]) > 120){
    GREAT$Term[i] <- paste0(substr(GREAT$Term[i], 1, as.integer(nchar(GREAT$Term[i])/3)), "-\n", 
                            substr(GREAT$Term[i], as.integer(nchar(GREAT$Term[i])/3) + 1, 2*as.integer(nchar(GREAT$Term[i])/3)), "-\n", 
                            substr(GREAT$Term[i], 2*as.integer(nchar(GREAT$Term[i])/3) + 1, nchar(GREAT$Term[i])))
  }
  if(nchar(GREAT$Term[i]) > 60 & nchar(GREAT$Term[i]) <= 120){
    GREAT$Term[i] <- paste0(substr(GREAT$Term[i], 1, as.integer(nchar(GREAT$Term[i])/2)), "-\n", substr(GREAT$Term[i], as.integer(nchar(GREAT$Term[i])/2)+1, nchar(GREAT$Term[i])))
  }
}
GREAT_lum <- droplevels(GREAT[1:20,])
GREAT_lum$Term <- factor(GREAT_lum$Term, levels = GREAT_lum$Term[length(GREAT_lum$Term):1])
GREAT_myo <- droplevels(GREAT[21:40,])
GREAT_myo$Term <- factor(GREAT_myo$Term, levels = GREAT_myo$Term[length(GREAT_myo$Term):1])
library(ggplot2)
(GREAT_lum_plot <- ggplot(data = GREAT_lum, aes(Term, -log10(FDR))) +
  geom_bar(aes(fill = Category), width = .5) + 
  coord_flip() + 
  geom_text(aes(label = round(-log10(FDR), 2), hjust = 0)) + 
  facet_grid(cell ~ .) + 
  xlab("") + 
  ylab("") + 
  theme_bw() +
  scale_fill_manual(values = c("GOBP" = "blue", "MSigPerturbation" = "purple", "PathwaysCommon" = "darkgreen")))
ggsave(GREAT_lum_plot, file = "GREAT_lum.pdf", width = 12, height = 8)
(GREAT_myo_plot <- ggplot(data = GREAT_myo, aes(Term, -log10(FDR))) +
   geom_bar(aes(fill = Category), width = .5) + 
   coord_flip() + 
   geom_text(aes(label = round(-log10(FDR), 2), hjust = 0)) + 
   facet_grid(cell ~ .) + 
   xlab("") + 
   ylab("-log10(Binomial FDR)") + 
   theme_bw() +
   scale_fill_manual(values = c("GOBP" = "blue", "MSigPerturbation" = "purple", "PathwaysCommon" = "darkgreen")))
ggsave(GREAT_myo_plot, file = "GREAT_myo.pdf", width = 12, height = 7)


# Figure S26
setwd("~/快盘/Publications/breast/revision/sup/FigureS26")
lum.GOBP <- read.delim("shown-GOBiologicalProcess.lum.UMR.TF.tsv", head = T, as.is = T)
lum.MSig <- read.delim("shown-MSigDBGeneSetsPerturbation.lum.UMR.TF.tsv", head = T, as.is = T)
myo.GOMF <- read.delim("shown-GOMolecularFunction.myo.UMR.TF.tsv", head = T, as.is = T)
myo.GOBP <- read.delim("shown-GOBiologicalProcess.myo.UMR.TF.tsv", head = T, as.is = T)
myo.pathway <- read.delim("shown-pathway.myo.UMR.TF.tsv", head = T, as.is = T)
GREAT <- data.frame(cell = c(rep("lum.UMRs.with.TFs", 20), rep("myo.UMRs.with.TFs", 30)), Category = rep(c("GOBP", "MSigPerturbation", "GOMF", "GOBP", "PathwaysCommon"), each = 10), 
                    Term = c(lum.GOBP$Term.Name[1:10], lum.MSig$Term.Name[1:10], myo.GOMF$Term.Name[1:10], myo.GOBP$Term.Name[1:10], myo.pathway$Term.Name[1:10]), 
                    FDR = c(lum.GOBP$Binom.FDR.Q.Val[1:10], lum.MSig$Binom.FDR.Q.Val[1:10], myo.GOMF$Binom.FDR.Q.Val[1:10], myo.GOBP$Binom.FDR.Q.Val[1:10], myo.pathway$Binom.FDR.Q.Val[1:10]))
GREAT$Term <- as.character(GREAT$Term)
GREAT <- na.omit(GREAT)
for(i in 1:nrow(GREAT)){
  if(nchar(GREAT$Term[i]) > 120){
    GREAT$Term[i] <- paste0(substr(GREAT$Term[i], 1, as.integer(nchar(GREAT$Term[i])/3)), "-\n", 
                            substr(GREAT$Term[i], as.integer(nchar(GREAT$Term[i])/3) + 1, 2*as.integer(nchar(GREAT$Term[i])/3)), "-\n", 
                            substr(GREAT$Term[i], 2*as.integer(nchar(GREAT$Term[i])/3) + 1, nchar(GREAT$Term[i])))
  }
  if(nchar(GREAT$Term[i]) > 60 & nchar(GREAT$Term[i]) <= 120){
    GREAT$Term[i] <- paste0(substr(GREAT$Term[i], 1, as.integer(nchar(GREAT$Term[i])/2)), "-\n", substr(GREAT$Term[i], as.integer(nchar(GREAT$Term[i])/2)+1, nchar(GREAT$Term[i])))
  }
}
GREAT_lum <- droplevels(GREAT[1:20,])
GREAT_lum$Term <- factor(GREAT_lum$Term, levels = GREAT_lum$Term[length(GREAT_lum$Term):1])
GREAT_myo <- droplevels(GREAT[21:43,])
GREAT_myo$Term <- factor(GREAT_myo$Term, levels = GREAT_myo$Term[length(GREAT_myo$Term):1])
library(ggplot2)
(GREAT_lum_plot <- ggplot(data = GREAT_lum, aes(Term, -log10(FDR))) +
   geom_bar(aes(fill = Category), width = .5) + 
   coord_flip() + 
   geom_text(aes(label = round(-log10(FDR), 2), hjust = 0)) + 
   facet_grid(cell ~ .) + 
   xlab("") + 
   ylab("") + 
   theme_bw() +
   scale_fill_manual(values = c("GOMF" = "brown", "GOBP" = "blue", "MSigPerturbation" = "purple", "PathwaysCommon" = "darkgreen")))
ggsave(GREAT_lum_plot, file = "GREAT_lum.pdf", width = 12, height = 9)
(GREAT_myo_plot <- ggplot(data = GREAT_myo, aes(Term, -log10(FDR))) +
   geom_bar(aes(fill = Category), width = .5) + 
   coord_flip() + 
   geom_text(aes(label = round(-log10(FDR), 2), hjust = 0)) + 
   facet_grid(cell ~ .) + 
   xlab("") + 
   ylab("-log10(Binomial FDR)") + 
   theme_bw() +
   scale_fill_manual(values = c("GOMF" = "brown", "GOBP" = "blue", "MSigPerturbation" = "purple", "PathwaysCommon" = "darkgreen")))
ggsave(GREAT_myo_plot, file = "GREAT_myo.pdf", width = 12, height = 7)




