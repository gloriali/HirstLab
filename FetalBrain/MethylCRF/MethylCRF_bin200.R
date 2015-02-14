# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# install.packages("plyr", lib = "/home/lli/R/")
require(plyr, lib.loc = "/home/lli/R/")
cutoff <- 0.05
summary <- matrix(NA, nrow = 3, ncol = 3, dimnames = list(c("MeDIP", "MethylCRF", "WGBS"), c("No.bins", "valid.t.test", "No.DMRs")))

# MethylCRF
setwd("~/FetalBrain/MeDIPMRE/")
load("MethylCRF.Rdata")
MethylCRF <- MethylCRF[MethylCRF$chr %in% paste0("chr", c(as.character(1:22), "X")),]
MethylCRF$binstart <- (MethylCRF$start %/% 200) * 200 + 1
MethylCRF$binend <- (MethylCRF$start %/% 200 + 1) * 200
MethylCRF$bin <- factor(paste0(MethylCRF$chr, ":", MethylCRF$binstart, "-", MethylCRF$binend))
(summary["MethylCRF", "No.bins"] <- length(levels(MethylCRF$bin)))
MethylCRF <- MethylCRF[, c("brain01", "brain02", "cortex01", "cortex02", "ge01", "ge02", "bin")] 
save(MethylCRF, file = "MethylCRF.Rdata")
MethylCRF_bin <- ddply(MethylCRF, ~ bin, summarise, brain01 = mean(brain01), brain02 = mean(brain02), cortex01 = mean(cortex01), cortex02 = mean(cortex02), ge01 = mean(ge01), ge02 = mean(ge02))
rownames(MethylCRF_bin) <- MethylCRF_bin$bin
MethylCRF_bin$bin <- NULL
nrow(MethylCRF_bin)
cortex02_ge02_MethylCRF <- ddply(MethylCRF, ~bin, summarize, p = try(t.test(cortex02, ge02, paired = T)$p.value, silent = T), cortex02_ge02 = abs(mean(cortex02)-mean(ge02)))
cortex02_ge02_MethylCRF <- na.omit(cortex02_ge02_MethylCRF[!grepl("Error", cortex02_ge02_MethylCRF$p), ])
cortex02_ge02_MethylCRF$adjusted.p <- p.adjust(cortex02_ge02_MethylCRF$p, method = "fdr")
rownames(cortex02_ge02_MethylCRF) <- cortex02_ge02_MethylCRF$bin
(summary["MethylCRF", "valid.t.test"] <- nrow(cortex02_ge02_MethylCRF))
cortex02_ge02_MethylCRF_dmr <- cortex02_ge02_MethylCRF[cortex02_ge02_MethylCRF$adjusted.p < cutoff,]
(summary["MethylCRF", "No.DMRs"] <- nrow(cortex02_ge02_MethylCRF_dmr))
save(MethylCRF_bin, cortex02_ge02_MethylCRF, file = "MethylCRF_bin200.Rdata")
save(cutoff, summary, cortex02_ge02_MethylCRF_dmr, file = "DMR_3methods_bin200.Rdata")
