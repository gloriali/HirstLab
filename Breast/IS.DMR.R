# visualizing individual specific DMRs on chrX
source("http://bioconductor.org/biocLite.R")
biocLite("ggbio")
library(ggplot2)
library(ggbio)
(p.ideo <- Ideogram(genome = "hg19", subchr = "chrX"))
ggsave(p.ideo@ggplot, file = "ideo.X.pdf", height = 2, width = 10)
setwd("~/快盘/REMC/IS.DMR.X/lum")
# setwd("~/快盘/REMC/IS.DMR.X/myo")
chrlen <- 156040895
DMR <- read.delim("IS.DMR.bed", head = F, as.is = T, colClasses = c("character", "integer", "integer", "integer", "factor"))
DMR$V5 <- factor(DMR$V5, levels = c("1000", "0111", "0100", "1011", "0010", "1101", "0001", "1110"))
(dmr.pos <- ggplot(data = DMR, aes(x = V5, y = (V2 + V3)/2, color = V5)) + 
   geom_point(size = 1, alpha = 0.5, position = position_jitter(width = 0.1)) + 
   coord_cartesian(ylim = c(0, 156040895)) + 
   xlab("lum IS DMR Pattern") + 
#    xlab("myo IS DMR Pattern") + 
   ylab("") + 
   coord_flip() + 
   scale_color_manual(values = c("1000" = "blue4", "0111" = "blue", "0100" = "cornflowerblue", "1011" = "cyan2", "0010" = "greenyellow", "1101" = "orange", "0001" = "red", "1110" = "firebrick4")) + 
   theme(axis.title.x = element_text(size = 18), axis.text.y = element_text(size = 15, color = "black"), axis.text.x = element_text(size = 0), legend.position = "none", panel.background = element_rect(fill = "transparent", color = "black"), plot.background = element_rect(fill = "transparent")))
ggsave(dmr.pos, file = "IS_DMR_X_pos_lum.pdf", width = 10)
# ggsave(dmr.pos, file = "IS_DMR_X_pos_myo.pdf", width = 10)


# GREAT on IS DMRs on chrX
library(ggplot2)
setwd("~/快盘/REMC/IS.DMR/chrX")
lum.InterPro <- read.delim("./lum/shown-interpro.tsv", head = T, as.is = T)
lum.MSig <- read.delim("./lum/shown-MSigDBGeneSetsPerturbation.tsv", head = T, as.is = T)
myo.InterPro <- read.delim("./myo/shown-interpro.tsv", head = T, as.is = T)
myo.MSig <- read.delim("./myo/shown-MSigDBGeneSetsPerturbation.tsv", head = T, as.is = T)
GREAT <- data.frame(cell = rep(c("lum", "myo"), each = 20), Category = rep(c("InterPro", "MSigPerturbation", "InterPro", "MSigPerturbation"), each = 10), 
                    Term = c(lum.InterPro$Term.Name[1:10], lum.MSig$Term.Name[1:10], myo.InterPro$Term.Name[1:10], myo.MSig$Term.Name[1:10]), 
                    FDR = c(lum.InterPro$Binom.FDR.Q.Val[1:10], lum.MSig$Binom.FDR.Q.Val[1:10], myo.InterPro$Binom.FDR.Q.Val[1:10], myo.MSig$Binom.FDR.Q.Val[1:10]))
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
GREAT[GREAT$FDR == 0, "FDR"] <- 10^-310
GREAT_lum <- droplevels(GREAT[1:20,])
GREAT_lum$Term <- factor(GREAT_lum$Term, levels = GREAT_lum$Term[length(GREAT_lum$Term):1])
GREAT_myo <- droplevels(GREAT[21:40,])
GREAT_myo$Term <- factor(GREAT_myo$Term, levels = GREAT_myo$Term[length(GREAT_myo$Term):1])
(GREAT_lum_plot <- ggplot(data = GREAT_lum, aes(Term, -log10(FDR))) +
   geom_bar(aes(fill = Category), width = .5) + 
   coord_flip() + 
   geom_text(aes(label = round(-log10(FDR), 2), hjust = 0)) + 
   facet_grid(cell ~ .) + 
   xlab("") + 
   ylab("") + 
   theme_bw() +
   scale_fill_manual(values = c("InterPro" = "blue", "MSigPerturbation" = "purple")))
ggsave(GREAT_lum_plot, file = "GREAT_IS_DMR_chrX_lum.pdf", width = 12, height = 9)
(GREAT_myo_plot <- ggplot(data = GREAT_myo, aes(Term, -log10(FDR))) +
   geom_bar(aes(fill = Category), width = .5) + 
   coord_flip() + 
   geom_text(aes(label = round(-log10(FDR), 2), hjust = 0)) + 
   facet_grid(cell ~ .) + 
   xlab("") + 
   ylab("-log10(Binomial FDR)") + 
   theme_bw() +
   scale_fill_manual(values = c("InterPro" = "blue", "MSigPerturbation" = "purple")))
ggsave(GREAT_myo_plot, file = "GREAT_IS_DMR_chrX_myo.pdf", width = 12, height = 9)


