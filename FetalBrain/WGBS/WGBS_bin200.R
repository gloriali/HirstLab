# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# install.packages("plyr", lib = "/home/lli/R/")
require(plyr, lib.loc = "/home/lli/R/")
cutoff <- 0.05
summary <- matrix(NA, nrow = 3, ncol = 3, dimnames = list(c("MeDIP", "MethylCRF", "WGBS"), c("No.bins", "valid.t.test", "No.DMRs")))

# WGBS_28M
setwd("~/FetalBrain/WGBS/")
load("WGBS_28M.Rdata")
# WGBS_28M <- WGBS_28M[WGBS_28M$chr %in% paste0("chr", c(as.character(1:22), "X")),]
# WGBS_28M$binstart <- (WGBS_28M$start %/% 200) * 200 + 1
# WGBS_28M$binend <- (WGBS_28M$start %/% 200 + 1) * 200
# WGBS_28M$bin <- factor(paste0(WGBS_28M$chr, ":", WGBS_28M$binstart, "-", WGBS_28M$binend))
# (summary["WGBS", "No.bins"] <- length(levels(WGBS_28M$bin)))
# WGBS_28M <- WGBS_28M[, c("cortex02", "ge02", "cortex04", "ge04", "bin")] 
# save(WGBS_28M, file = "WGBS_28M.Rdata")
WGBS_28M_bin <- ddply(WGBS_28M, ~ bin, summarise, cortex02 = mean(cortex02, na.rm = T), ge02 = mean(ge02, na.rm = T), cortex04 = mean(cortex04, na.rm = T), ge04 = mean(ge04, na.rm = T))
# rownames(WGBS_28M_bin) <- WGBS_28M_bin$bin
# WGBS_28M_bin$bin <- NULL
nrow(WGBS_28M_bin)
cortex02_ge02_WGBS_28M <- ddply(WGBS_28M, ~bin, summarize, p = try(t.test(cortex02, ge02, paired = T)$p.value, silent = T), cortex02_ge02 = abs(mean(cortex02, na.rm = T)-mean(ge02, na.rm = T)))
cortex02_ge02_WGBS_28M <- cortex02_ge02_WGBS_28M[!grepl("Error", cortex02_ge02_WGBS_28M$p), ]
cortex02_ge02_WGBS_28M <- cortex02_ge02_WGBS_28M[!is.na(cortex02_ge02_WGBS_28M$p), ]
cortex02_ge02_WGBS_28M$adjusted.p <- p.adjust(cortex02_ge02_WGBS_28M$p, method = "fdr")
rownames(cortex02_ge02_WGBS_28M) <- cortex02_ge02_WGBS_28M$bin
(summary["WGBS", "valid.t.test"] <- nrow(cortex02_ge02_WGBS_28M))
cortex02_ge02_WGBS_28M_dmr <- cortex02_ge02_WGBS_28M[cortex02_ge02_WGBS_28M$adjusted.p < cutoff,]
(summary["WGBS", "No.DMRs"] <- nrow(cortex02_ge02_WGBS_28M_dmr))
save(WGBS_28M_bin, cortex02_ge02_WGBS_28M, file = "WGBS_28M_bin200.Rdata")
rm(WGBS_28M, cortex02_ge02_WGBS_28M)

save(cutoff, summary, cortex02_ge02_WGBS_28M_dmr, file = "DMR_3methods_bin200.Rdata")
