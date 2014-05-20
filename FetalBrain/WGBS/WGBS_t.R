# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# install.packages("plyr", lib = "/home/lli/R/")
require(plyr, lib.loc = "/home/lli/R/")
cutoff <- 0.05
summary <- matrix(NA, nrow = 3, ncol = 3, dimnames = list(c("MeDIP", "MethylCRF", "WGBS"), c("No.bins", "valid.t.test", "No.DMRs")))
# WGBS_28M
setwd("~/FetalBrain/WGBS/")
load("WGBS_28M.Rdata")
cortex02_ge02_WGBS_28M <- ddply(WGBS_28M, ~bin, summarize, p = try(t.test(cortex02, ge02, paired = T)$p.value, silent = T), cortex02_ge02 = abs(mean(cortex02, na.rm = T)-mean(ge02, na.rm = T)))
cortex02_ge02_WGBS_28M <- cortex02_ge02_WGBS_28M[!grepl("Error", cortex02_ge02_WGBS_28M$p), ]
cortex02_ge02_WGBS_28M <- cortex02_ge02_WGBS_28M[!is.na(cortex02_ge02_WGBS_28M$p), ]
cortex02_ge02_WGBS_28M$adjusted.p <- p.adjust(cortex02_ge02_WGBS_28M$p, method = "fdr")
rownames(cortex02_ge02_WGBS_28M) <- cortex02_ge02_WGBS_28M$bin
(summary["WGBS", "valid.t.test"] <- nrow(cortex02_ge02_WGBS_28M))
cortex02_ge02_WGBS_28M_dmr <- cortex02_ge02_WGBS_28M[cortex02_ge02_WGBS_28M$adjusted.p < cutoff,]
(summary["WGBS", "No.DMRs"] <- nrow(cortex02_ge02_WGBS_28M_dmr))
save(cutoff, summary, cortex02_ge02_WGBS_28M, cortex02_ge02_WGBS_28M_dmr, file = "DMR_WGBS_bin200.Rdata")
