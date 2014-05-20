# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
setwd("~/FetalBrain/MeDIP/")
load("MeDIP_FetalBrain_fractional.Rdata")

e <- 1e-6
MeDIP$logfold_brain <- log2((MeDIP$brain01 + e) / (MeDIP$brain02 + e))
MeDIP$logfold_cortex <- log2((MeDIP$cortex01 + e) / (MeDIP$cortex02 + e))
MeDIP$logfold_ge <- log2((MeDIP$ge01 + e) / (MeDIP$ge02 + e))
# MeDIP$logfold_brain <- log2((MeDIP$brain01) / (MeDIP$brain02))
# MeDIP$logfold_cortex <- log2((MeDIP$cortex01) / (MeDIP$cortex02))
# MeDIP$logfold_ge <- log2((MeDIP$ge01) / (MeDIP$ge02))
# fold_brain <- MeDIP[is.finite(MeDIP$logfold_brain), ]
# fold_cortex <- MeDIP[is.finite(MeDIP$logfold_cortex), ]
# fold_ge <- MeDIP[is.finite(MeDIP$logfold_ge), ]
fold_brain <- MeDIP[abs(MeDIP$logfold_brain) > 1, ]
fold_cortex <- MeDIP[abs(MeDIP$logfold_cortex) > 1, ]
fold_ge <- MeDIP[abs(MeDIP$logfold_ge) > 1, ]

# xrange <- range(fold_brain$logfold_brain, fold_cortex$logfold_cortex, fold_ge$logfold_ge)
# yrange <- range(density(fold_brain$logfold_brain)$y, density(fold_cortex$logfold_cortex)$y, density(fold_ge$logfold_ge)$y)
xrange <- c(-10, 10)
yrange <- c(0, 0.05)
pdf(file="Asymmetry.pdf")
plot(xrange, yrange, type="n", main="Global DNA Methylation", xlab="log2(MeDIP01/MeDIP02)", ylab="density")
lines(density(fold_brain$logfold_brain, adjust = 0.8), col=1, lty=1, lwd=3)
lines(density(fold_cortex$logfold_cortex, adjust = 0.8), col=2, lty=1, lwd=3)
lines(density(fold_ge$logfold_ge, adjust = 0.8), col=3, lty=1, lwd=3)
legend("topright", c("brain","cortex","GE"), col=c(1:3), lty=1, lwd=5)
dev.off()
