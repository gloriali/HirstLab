# DMR analysis from MeDIP fractional calls  
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.signal.sh
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.fractional.m
# ~/HirstLab/FetalBrain/MeDIP/MeDIP.DM.sh

source("~/HirstLab/Pipeline/DMR.plots.R")
col <- c("chr", "start", "end", "ID", "DM", "CpG_count", "length") # format of DMR files
# Between MZ twins 
setwd("~/快盘/FetalBrain/MeDIP/DMR/MZ/")
Brain01_Brain02_DMR <- read.delim("DMR.Brain-HuFNSC01_Brain-HuFNSC02.d0.4.s200.c3.bed", as.is = T, head = F)
Brain01_Brain02_DMR_figures <- DMR_figures(Brain01_Brain02_DMR, sample1 = "Brain-HuFNSC01", sample2 = "Brain-HuFNSC02")
(Brain01_Brain02_DMR_figures$length)
(Brain01_Brain02_DMR_figures$count)
(Brain01_Brain02_DMR_figures$dis)
(Brain01_Brain02_DMR_figures$freq)
(Brain01_Brain02_DMR_figures$pos)
Cortex01_Cortex02_DMR <- read.delim("DMR.Cortex-HuFNSC01_Cortex-HuFNSC02.d0.4.s200.c3.bed", as.is = T, head = F)
Cortex01_Cortex02_DMR_figures <- DMR_figures(Cortex01_Cortex02_DMR, sample1 = "Cortex-HuFNSC01", sample2 = "Cortex-HuFNSC02")
(Cortex01_Cortex02_DMR_figures$length)
(Cortex01_Cortex02_DMR_figures$count)
(Cortex01_Cortex02_DMR_figures$dis)
(Cortex01_Cortex02_DMR_figures$freq)
(Cortex01_Cortex02_DMR_figures$pos)
GE01_GE02_DMR <- read.delim("DMR.GE-HuFNSC01_GE-HuFNSC02.d0.4.s200.c3.bed", as.is = T, head = F)
GE01_GE02_DMR_figures <- DMR_figures(GE01_GE02_DMR, sample1 = "GE-HuFNSC01", sample2 = "GE-HuFNSC02")
(GE01_GE02_DMR_figures$length)
(GE01_GE02_DMR_figures$count)
(GE01_GE02_DMR_figures$dis)
(GE01_GE02_DMR_figures$freq)
(GE01_GE02_DMR_figures$pos)


save.image("DMR_MZ.Rdata")


