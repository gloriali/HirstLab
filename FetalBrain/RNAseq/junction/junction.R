####################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# read junction files
setwd("~/FetalBrain/RNAseq/junction")
brain01 <- read.delim("A03484.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
brain02 <- read.delim("A07825.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
cortex01 <- read.delim("A03473.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
cortex02 <- read.delim("A03475.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
cortex03 <- read.delim("A04599.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
cortex04 <- read.delim("A15298.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
ge01 <- read.delim("A03474.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
ge02 <- read.delim("A03476.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
ge03 <- read.delim("A15295.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
ge04 <- read.delim("A15299.q5.F516.jb.bedGraph.gz", head = F, as.is = T)

pdf("Ecdf_coverage.pdf")
plot(c(0, 500), c(0, 1), main = "FetalBrain junction coverage ecdf", xlab = "coverage", ylab = "Ecdf")
lines(ecdf(brain01$V3), col = 1, lwd = 2)
lines(ecdf(brain02$V3), col = 1, lwd = 2)
lines(ecdf(cortex01$V3), col = 1, lwd = 2)
lines(ecdf(cortex02$V3), col = 1, lwd = 2)
lines(ecdf(cortex03$V3), col = 1, lwd = 2)
lines(ecdf(cortex04$V3), col = 2, lwd = 2)
lines(ecdf(ge01$V3), col = 1, lwd = 2)
lines(ecdf(ge02$V3), col = 1, lwd = 2)
lines(ecdf(ge03$V3), col = 2, lwd = 2)
lines(ecdf(ge04$V3), col = 2, lwd = 2)
legend("topleft", c("single lane", "pooled"), col = 1:2, lwd = 5)
plot(c(0, 50), c(0, 1), main = "FetalBrain junction coverage ecdf", xlab = "coverage", ylab = "Ecdf")
lines(ecdf(brain01$V3), col = 1, lwd = 2)
lines(ecdf(brain02$V3), col = 1, lwd = 2)
lines(ecdf(cortex01$V3), col = 1, lwd = 2)
lines(ecdf(cortex02$V3), col = 1, lwd = 2)
lines(ecdf(cortex03$V3), col = 1, lwd = 2)
lines(ecdf(cortex04$V3), col = 2, lwd = 2)
lines(ecdf(ge01$V3), col = 1, lwd = 2)
lines(ecdf(ge02$V3), col = 1, lwd = 2)
lines(ecdf(ge03$V3), col = 2, lwd = 2)
lines(ecdf(ge04$V3), col = 2, lwd = 2)
legend("topleft", c("single lane", "pooled"), col = 1:2, lwd = 5)
dev.off()

junction <- data.frame(id = unique(c(brain01$V1, brain02$V1, cortex01$V1, cortex02$V1, cortex03$V1, cortex04$V1, ge01$V1, ge02$V1, ge03$V1, ge04$V1)))
rownames(junction) <- junction$id
junction[, c("brain01N", "brain02N", "cortex01N", "cortex02N", "cortex03N", "cortex04N", "ge01N", "ge02N", "ge03N", "ge04N")] <- 0
junction[brain01$V1, "brain01N"] <- brain01$V3
junction[brain02$V1, "brain02N"] <- brain02$V3
junction[cortex01$V1, "cortex01N"] <- cortex01$V3
junction[cortex02$V1, "cortex02N"] <- cortex02$V3
junction[cortex03$V1, "cortex03N"] <- cortex03$V3
junction[cortex04$V1, "cortex04N"] <- cortex04$V3
junction[ge01$V1, "ge01N"] <- ge01$V3
junction[ge02$V1, "ge02N"] <- ge02$V3
junction[ge03$V1, "ge03N"] <- ge03$V3
junction[ge04$V1, "ge04N"] <- ge04$V3
nrow(junction)
summary(junction)

# calculate RPKM: 10^3*10^6 / (N*read_length)  
# check read length: /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools view /projects/epigenomics/ep50/internal/jqc.1.7.6/A03484/A03484.bam -- 75M
read_length <- 75
Nbrain01 <- 61400917   # from /projects/epigenomics/ep50/internal/jqc.1.7.6/A03484/A03484.report: Total number of exonic reads for RPKM (protein coding; no MT, no ribo proteins, top expressed 0.005 exons excluded) 
Nbrain02 <- 69798425
Ncortex01 <- 64807727
Ncortex02 <- 68340512
Ncortex03 <- 76475869
Ncortex04 <- 226802068
Nge01 <- 68292449
Nge02 <- 78885239
Nge03 <- 211087273
Nge04 <- 209694372
junction$brain01rpkm <- junction$brain01N*10^3*10^6 / (Nbrain01*read_length) 
junction$brain02rpkm <- junction$brain02N*10^3*10^6 / (Nbrain02*read_length) 
junction$cortex01rpkm <- junction$cortex01N*10^3*10^6 / (Ncortex01*read_length) 
junction$cortex02rpkm <- junction$cortex02N*10^3*10^6 / (Ncortex02*read_length) 
junction$cortex03rpkm <- junction$cortex03N*10^3*10^6 / (Ncortex03*read_length) 
junction$cortex04rpkm <- junction$cortex04N*10^3*10^6 / (Ncortex04*read_length) 
junction$ge01rpkm <- junction$ge01N*10^3*10^6 / (Nge01*read_length) 
junction$ge02rpkm <- junction$ge02N*10^3*10^6 / (Nge02*read_length) 
junction$ge03rpkm <- junction$ge03N*10^3*10^6 / (Nge03*read_length) 
junction$ge04rpkm <- junction$ge04N*10^3*10^6 / (Nge04*read_length) 
summary(junction)

pdf("Ecdf_rpkm.pdf")
plot(c(0, 50), c(0, 1), main = "FetalBrain junction RPKM ecdf", xlab = "RPKM", ylab = "Ecdf")
lines(ecdf(junction$brain01rpkm), col = 2, lwd = 2)
lines(ecdf(junction$brain02rpkm), col = 2, lwd = 2)
lines(ecdf(junction$cortex01rpkm), col = 3, lwd = 2)
lines(ecdf(junction$cortex02rpkm), col = 3, lwd = 2)
lines(ecdf(junction$cortex03rpkm), col = 3, lwd = 2)
lines(ecdf(junction$cortex04rpkm), col = 3, lwd = 2)
lines(ecdf(junction$ge01rpkm), col = 4, lwd = 2)
lines(ecdf(junction$ge02rpkm), col = 4, lwd = 2)
lines(ecdf(junction$ge03rpkm), col = 4, lwd = 2)
lines(ecdf(junction$ge04rpkm), col = 4, lwd = 2)
legend("topleft", c("brain", "cortex", "GE"), col = 2:4, lwd = 5)
plot(c(0, 5), c(0, 1), main = "FetalBrain junction RPKM ecdf", xlab = "RPKM", ylab = "Ecdf")
lines(ecdf(junction$brain01rpkm), col = 2, lwd = 2)
lines(ecdf(junction$brain02rpkm), col = 2, lwd = 2)
lines(ecdf(junction$cortex01rpkm), col = 3, lwd = 2)
lines(ecdf(junction$cortex02rpkm), col = 3, lwd = 2)
lines(ecdf(junction$cortex03rpkm), col = 3, lwd = 2)
lines(ecdf(junction$cortex04rpkm), col = 3, lwd = 2)
lines(ecdf(junction$ge01rpkm), col = 4, lwd = 2)
lines(ecdf(junction$ge02rpkm), col = 4, lwd = 2)
lines(ecdf(junction$ge03rpkm), col = 4, lwd = 2)
lines(ecdf(junction$ge04rpkm), col = 4, lwd = 2)
legend("topleft", c("brain", "cortex", "GE"), col = 2:4, lwd = 5)
dev.off()
save(junction, file = "junction.Rdata")
####################################################################################





