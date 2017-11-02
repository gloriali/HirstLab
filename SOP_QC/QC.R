# QC for optimizing MeDIP and hMeDIP SOP

## QC for MeDIP positive and negative qPCR primers
library(dplyr)
library(ggplot2)
setwd("/projects/epigenomics2/users/lli/")
meDIP_qPCR <- read.delim("meDIP_qPCR_primers_hg19.5mC.bed", as.is = T) %>%
  mutate(type = gsub("[0-9]+", "", ID), CpG_density = n/(end-start)*100)
(meDIP_qPCR_fractional_figure <- ggplot(meDIP_qPCR, aes(ID, fractional, fill = type)) + 
    geom_boxplot() + 
    scale_fill_manual(values = c("red", "blue")) + 
    xlab("") + 
    ylab("Fractional methylation") +
    theme_bw())
ggsave(meDIP_qPCR_fractional_figure, file = "meDIP_qPCR_fractional_figure.pdf", height = 4, width = 5)
meDIP_qPCR_CGdensity <- read.delim("meDIP_qPCR_primers_hg19.CpGdensity.bed", as.is = T, head = F, col.names = c("chr", "start", "end", "ID", "n", "density")) %>%
  mutate(type = gsub("[0-9]+", "", ID))
(meDIP_qPCR_CGdensity_figure <- ggplot(meDIP_qPCR_CGdensity, aes(ID, density, color = type)) + 
    geom_point(size = 3) + 
    scale_color_manual(values = c("red", "blue")) + 
    xlab("") + 
    ylab("CpG density - #CpG per 100bp") +
    theme_bw())
ggsave(meDIP_qPCR_CGdensity_figure, file = "meDIP_qPCR_CGdensity_figure.pdf", height = 4, width = 5)
