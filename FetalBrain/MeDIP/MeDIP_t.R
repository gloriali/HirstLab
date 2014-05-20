# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# install.packages("plyr", lib = "/home/lli/R/")
require(plyr, lib.loc = "/home/lli/R/")
cutoff <- 0.05
summary <- matrix(NA, nrow = 3, ncol = 3, dimnames = list(c("MeDIP", "MethylCRF", "WGBS"), c("No.bins", "valid.t.test", "No.DMRs")))

# MeDIP
setwd("~/FetalBrain/MeDIP/")
load("MeDIP_FetalBrain_fractional.Rdata")
cortex02_ge02_MeDIP <- ddply(MeDIP, ~bin, summarize, p = try(t.test(cortex02, ge02, paired = T)$p.value, silent = T), cortex02_ge02 = abs(mean(cortex02)-mean(ge02)))
cortex02_ge02_MeDIP <- na.omit(cortex02_ge02_MeDIP[!grepl("Error", cortex02_ge02_MeDIP$p), ])
cortex02_ge02_MeDIP$adjusted.p <- p.adjust(cortex02_ge02_MeDIP$p, method = "fdr")
rownames(cortex02_ge02_MeDIP) <- cortex02_ge02_MeDIP$bin
(summary["MeDIP", "valid.t.test"] <- nrow(cortex02_ge02_MeDIP))
cortex02_ge02_MeDIP_dmr <- cortex02_ge02_MeDIP[cortex02_ge02_MeDIP$adjusted.p < cutoff,]
(summary["MeDIP", "No.DMRs"] <- nrow(cortex02_ge02_MeDIP_dmr))
save(cutoff, summary, cortex02_ge02_MeDIP, cortex02_ge02_MeDIP_dmr, file = "DMR_MeDIP_bin200.Rdata")
