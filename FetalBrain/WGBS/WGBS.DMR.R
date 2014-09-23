#' DMR analysis from MeDIP fractional calls  

#' DMR identification:  
#' ~/HirstLab/FetalBrain/WGBS/WGBS.DM.sh

setwd("~/快盘/FetalBrain/WGBS/")
library(ggplot2)

#' testing DMR script parameters
DM_test_summary <- read.delim("DM.summary.stats", head = F, as.is = T, col.names = c("Sample", "p", "delta", "Total.DM.CpGs", "Hyper.DM.CpGs", "Hypo.DM.CpGs"))
DM_test_summary$delta <- paste0("delta=", DM_test_summary$delta)
DMR_test_summary <- read.delim("DMR.summary.stats", head = F, as.is = T, col.names = c("Sample", "size", "CpG.cut", "Median.length", "Median.CpG", "Total.DMR", "Hyper.DMR", "Hypo.DMR"))
DMR_test_summary$p <- gsub(".d[0-9.]+", "", DMR_test_summary$Sample)
DMR_test_summary$p <- paste0("p=", gsub(".*.p", "", DMR_test_summary$p))
DMR_test_summary$delta <- paste0("delta=", gsub(".*.d", "", DMR_test_summary$Sample))
DMR_test_summary$Sample <- gsub(".p.*.d.*", "", DMR_test_summary$Sample)
pdf("DMR_parameters.pdf")
(ggplot(DM_test_summary, aes(x = p, y = Total.DM.CpGs)) + scale_x_log10(breaks = c(0.05, 0.005, 0.0005)) + geom_point(aes(color = Sample)) + geom_line(aes(color = Sample)) + facet_wrap(~ delta) + theme_bw() + xlab("p-value cutoff") + ylab("Total No. of DM CpGs"))
(ggplot(DM_test_summary, aes(x = p, y = log2(Hyper.DM.CpGs/Hypo.DM.CpGs), color = Sample)) + scale_x_log10(breaks = c(0.05, 0.005, 0.0005)) + geom_point() + geom_line() + geom_hline(aes(yintercept = 0)) + facet_wrap(~ delta) + theme_bw() + xlab("p-value cutoff") + ylab("log2(Hyper/Hypo)"))
(ggplot(DMR_test_summary, aes(x = size, y = Median.length)) + geom_point(aes(color = Sample)) + geom_line(aes(color = Sample)) + facet_grid(p ~ delta) + geom_hline(aes(yintercept = 250)) + theme_bw() + xlab("Max distance between adjacent DM CpGs") + ylab("Median DMR length"))
(ggplot(DMR_test_summary, aes(x = size, y = Total.DMR, color = Sample)) + geom_point() + geom_line() + facet_grid(p ~ delta) + theme_bw() + xlab("Max distance between adjacent DM CpGs") + ylab("Total No. of DMRs"))
(ggplot(DMR_test_summary, aes(x = p, y = Total.DMR, color = Sample)) + geom_point() + geom_line() + facet_grid(size ~ delta) + theme_bw() + xlab("p-value cutoff") + ylab("Total No. of DMRs"))
(ggplot(DMR_test_summary, aes(x = size, y = log2(Hyper.DMR/Hypo.DMR), color = Sample)) + geom_point() + geom_line() + geom_hline(aes(yintercept = 0)) + facet_grid(p ~ delta) + theme_bw() + xlab("Max distance between adjacent DM CpGs") + ylab("log2(Hyper/Hypo)"))
dev.off()

