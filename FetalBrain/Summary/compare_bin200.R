# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# install.packages("plyr", lib = "/home/lli/R/")
require(plyr, lib.loc = "/home/lli/R/")
cutoff <- 0.05
summary <- matrix(NA, nrow = 3, ncol = 3, dimnames = list(c("MeDIP", "MethylCRF", "WGBS"), c("No.bins", "valid.t.test", "No.DMRs")))

# MeDIP
setwd("~/FetalBrain/MeDIP/")
load("MeDIP_FetalBrain_fractional.Rdata")
MeDIP <- MeDIP[MeDIP$chr %in% paste0("chr", c(as.character(1:22), "X")),]
MeDIP$binstart <- (MeDIP$start %/% 200) * 200 + 1
MeDIP$binend <- (MeDIP$start %/% 200 + 1) * 200
MeDIP$bin <- factor(paste0(MeDIP$chr, ":", MeDIP$binstart, "-", MeDIP$binend))
(summary["MeDIP", "No.bins"] <- length(levels(MeDIP$bin)))
MeDIP <- MeDIP[, c("brain01", "brain02", "cortex01", "cortex02", "ge01", "ge02", "bin")] 
save(MeDIP, file = "MeDIP_FetalBrain_fractional.Rdata")
MeDIP_bin <- ddply(MeDIP, ~ bin, summarise, brain01 = mean(brain01), brain02 = mean(brain02), cortex01 = mean(cortex01), cortex02 = mean(cortex02), ge01 = mean(ge01), ge02 = mean(ge02))
rownames(MeDIP_bin) <- MeDIP_bin$bin
MeDIP_bin$bin <- NULL
nrow(MeDIP_bin)
cortex02_ge02_MeDIP <- ddply(MeDIP, ~bin, summarize, p = try(t.test(cortex02, ge02, paired = T)$p.value, silent = T), cortex02_ge02 = abs(mean(cortex02)-mean(ge02)))
cortex02_ge02_MeDIP <- na.omit(cortex02_ge02_MeDIP[!grepl("Error", cortex02_ge02_MeDIP$p), ])
cortex02_ge02_MeDIP$adjusted.p <- p.adjust(cortex02_ge02_MeDIP$p, method = "fdr")
rownames(cortex02_ge02_MeDIP) <- cortex02_ge02_MeDIP$bin
(summary["MeDIP", "valid.t.test"] <- nrow(cortex02_ge02_MeDIP))
cortex02_ge02_MeDIP_dmr <- cortex02_ge02_MeDIP[cortex02_ge02_MeDIP$adjusted.p < cutoff,]
(summary["MeDIP", "No.DMRs"] <- nrow(cortex02_ge02_MeDIP_dmr))
save(MeDIP_bin, cortex02_ge02_MeDIP, file = "MeDIP_bin200.Rdata")
rm(MeDIP, cortex02_ge02_MeDIP)

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
rm(MethylCRF, cortex02_ge02_MethylCRF)

# WGBS_28M
setwd("~/FetalBrain/WGBS/")
load("WGBS_28M.Rdata")
WGBS_28M <- WGBS_28M[WGBS_28M$chr %in% paste0("chr", c(as.character(1:22), "X")),]
WGBS_28M$binstart <- (WGBS_28M$start %/% 200) * 200 + 1
WGBS_28M$binend <- (WGBS_28M$start %/% 200 + 1) * 200
WGBS_28M$bin <- factor(paste0(WGBS_28M$chr, ":", WGBS_28M$binstart, "-", WGBS_28M$binend))
(summary["WGBS", "No.bins"] <- length(levels(WGBS_28M$bin)))
WGBS_28M <- WGBS_28M[, c("cortex02", "ge02", "cortex04", "ge04", "bin")] 
save(WGBS_28M, file = "WGBS_28M.Rdata")
WGBS_28M_bin <- ddply(WGBS_28M, ~ bin, summarise, cortex02 = mean(cortex02, na.rm = T), ge02 = mean(ge02, na.rm = T), cortex04 = mean(cortex04, na.rm = T), ge04 = mean(ge04, na.rm = T))
rownames(WGBS_28M_bin) <- WGBS_28M_bin$bin
WGBS_28M_bin$bin <- NULL
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

setwd("~/FetalBrain/MeDIP/")
save(cutoff, summary, cortex02_ge02_MeDIP_dmr, cortex02_ge02_MethylCRF_dmr, cortex02_ge02_WGBS_28M_dmr, file = "DMR_3methods_bin200.Rdata")
rm(summary, cortex02_ge02_MeDIP_dmr, cortex02_ge02_MethylCRF_dmr, cortex02_ge02_WGBS_28M_dmr)

# correlation between MeDIP, MethylCRF and WGBS
install.packages("hexbin", lib = "/home/lli/R/")
require(hexbin, lib.loc = "/home/lli/R")
inx <- intersect(rownames(MeDIP_bin), intersect(rownames(MethylCRF_bin), rownames(WGBS_28M_bin)))
cortex02_bin200 <- data.frame(MeDIP = MeDIP_bin[inx, "cortex02"], MethylCRF = MethylCRF_bin[inx, "cortex02"], WGBS = WGBS_28M_bin[inx, "cortex02"])
ge02_bin200 <- data.frame(MeDIP = MeDIP_bin[inx, "ge02"], MethylCRF = MethylCRF_bin[inx, "ge02"], WGBS = WGBS_28M_bin[inx, "ge02"])
rm(MeDIP_bin, MethylCRF_bin, WGBS_28M_bin)
pdf("bin200_3methods_correlation.pdf")
splom(cortex02_bin200, panel = panel.hexbinplot, main = "correlation between methods - cortex02")
splom(ge02_bin200, panel = panel.hexbinplot, main = "correlation between methods - ge02")
dev.off()
save(cortex02_bin200, ge02_bin200, file = "cortex02_ge02_3methods.Rdata")
