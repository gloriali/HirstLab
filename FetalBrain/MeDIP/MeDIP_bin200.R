# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R
# install.packages("plyr", lib = "/home/lli/R")
require(plyr, lib.loc = "/home/lli/R")
setwd("/home/lli/FetalBrain/MeDIP/")
load("MeDIP_FetalBrain_fractional.Rdata")
load("~/hg19/CpG_200bin.Rdata")

# MeDIP
MeDIP$bin <- CG_bin[rownames(MeDIP), ]$V2 
MeDIP_bin <- ddply(MeDIP, ~ bin, summarise, brain01 = mean(brain01), brain02 = mean(brain02), cortex01 = mean(cortex01), cortex02 = mean(cortex02), ge01 = mean(ge01), ge02 = mean(ge02))
rownames(MeDIP_bin) <- MeDIP_bin$bin
cortex02_ge02_MeDIP <- ddply(MeDIP, ~bin, summarize, p = try(t.test(cortex02, ge02, paired = T)$p.value, silent = T), cortex02_ge02 = abs(mean(cortex02)-mean(ge02)))
cortex02_ge02_MeDIP <- cortex02_ge02_MeDIP[!grepl("Error", cortex02_ge02_MeDIP$p), ]
cortex02_ge02_MeDIP$adjusted.p <- p.adjust(cortex02_ge02_MeDIP$p, method = "fdr")
rownames(cortex02_ge02_MeDIP) <- cortex02_ge02_MeDIP$bin
nrow(cortex02_ge02_MeDIP)

save(MeDIP_bin, cortex02_ge02_MeDIP, file = "MeDIP_bin200.Rdata")

