####################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# clustering on histone mods promoter/genebody signal
# promoter: TSS+/-2kb H3K4me3(missing cortex1, GE1), H3K4me1, H3K27me3, input
# genebody: H3K36me3
setwd("~/FetalBrain/HisMod/signal/")
samples <- c("brain1", "brain2", "cortex1", "cortex2", "GE1", "GE2")

pdf("signal_cluster.pdf")
# H3K4me1
H3K4me1_TSS2000_pc <- read.delim("~/FetalBrain/HisMod/signal/hg19v65_genes_TSS_2000.H3K4me1.pc", head = F, as.is = T, row.names = 1)
colnames(H3K4me1_TSS2000_pc) <- samples
cor_H3K4me1_TSS2000_pc <- cor(H3K4me1_TSS2000_pc, method="pearson")
dist_H3K4me1_TSS2000_pc <- as.dist(1 - cor_H3K4me1_TSS2000_pc)
hc_H3K4me1_TSS2000_pc <- as.dendrogram(hclust(dist_H3K4me1_TSS2000_pc, method="complete"))
plot(hc_H3K4me1_TSS2000_pc, main="H3K4me1 TSS_2000 clustering on pc genes", horiz=TRUE)
H3K4me1_TSS2000_nc <- read.delim("~/FetalBrain/HisMod/signal/hg19v65_genes_TSS_2000.H3K4me1.nc", head = F, as.is = T, row.names = 1)
colnames(H3K4me1_TSS2000_nc) <- samples
cor_H3K4me1_TSS2000_nc <- cor(H3K4me1_TSS2000_nc, method="pearson")
dist_H3K4me1_TSS2000_nc <- as.dist(1 - cor_H3K4me1_TSS2000_nc)
hc_H3K4me1_TSS2000_nc <- as.dendrogram(hclust(dist_H3K4me1_TSS2000_nc, method="complete"))
plot(hc_H3K4me1_TSS2000_nc, main="H3K4me1 TSS_2000 clustering on nc genes", horiz=TRUE)

# H3K4me3
H3K4me3_TSS2000_pc <- read.delim("~/FetalBrain/HisMod/signal/hg19v65_genes_TSS_2000.H3K4me3.pc", head = F, as.is = T, row.names = 1)
colnames(H3K4me3_TSS2000_pc) <- c("brain1", "brain2", "cortex2", "GE2")
cor_H3K4me3_TSS2000_pc <- cor(H3K4me3_TSS2000_pc, method="pearson")
dist_H3K4me3_TSS2000_pc <- as.dist(1 - cor_H3K4me3_TSS2000_pc)
hc_H3K4me3_TSS2000_pc <- as.dendrogram(hclust(dist_H3K4me3_TSS2000_pc, method="complete"))
plot(hc_H3K4me3_TSS2000_pc, main="H3K4me3 TSS_2000 clustering on pc genes", horiz=TRUE)
H3K4me3_TSS2000_nc <- read.delim("~/FetalBrain/HisMod/signal/hg19v65_genes_TSS_2000.H3K4me3.nc", head = F, as.is = T, row.names = 1)
colnames(H3K4me3_TSS2000_nc) <- c("brain1", "brain2", "cortex2", "GE2")
cor_H3K4me3_TSS2000_nc <- cor(H3K4me3_TSS2000_nc, method="pearson")
dist_H3K4me3_TSS2000_nc <- as.dist(1 - cor_H3K4me3_TSS2000_nc)
hc_H3K4me3_TSS2000_nc <- as.dendrogram(hclust(dist_H3K4me3_TSS2000_nc, method="complete"))
plot(hc_H3K4me3_TSS2000_nc, main="H3K4me3 TSS_2000 clustering on nc genes", horiz=TRUE)

# H3K27me3
H3K27me3_TSS2000_pc <- read.delim("~/FetalBrain/HisMod/signal/hg19v65_genes_TSS_2000.H3K27me3.pc", head = F, as.is = T, row.names = 1)
colnames(H3K27me3_TSS2000_pc) <- samples
cor_H3K27me3_TSS2000_pc <- cor(H3K27me3_TSS2000_pc, method="pearson")
dist_H3K27me3_TSS2000_pc <- as.dist(1 - cor_H3K27me3_TSS2000_pc)
hc_H3K27me3_TSS2000_pc <- as.dendrogram(hclust(dist_H3K27me3_TSS2000_pc, method="complete"))
plot(hc_H3K27me3_TSS2000_pc, main="H3K27me3 TSS_2000 clustering on pc genes", horiz=TRUE)
H3K27me3_TSS2000_nc <- read.delim("~/FetalBrain/HisMod/signal/hg19v65_genes_TSS_2000.H3K27me3.nc", head = F, as.is = T, row.names = 1)
colnames(H3K27me3_TSS2000_nc) <- samples
cor_H3K27me3_TSS2000_nc <- cor(H3K27me3_TSS2000_nc, method="pearson")
dist_H3K27me3_TSS2000_nc <- as.dist(1 - cor_H3K27me3_TSS2000_nc)
hc_H3K27me3_TSS2000_nc <- as.dendrogram(hclust(dist_H3K27me3_TSS2000_nc, method="complete"))
plot(hc_H3K27me3_TSS2000_nc, main="H3K27me3 TSS_2000 clustering on nc genes", horiz=TRUE)

# H3K36me3
H3K36me3_genebody_pc <- read.delim("~/FetalBrain/HisMod/signal/hg19v65_genes_body.H3K36me3.pc", head = F, as.is = T, row.names = 1)
colnames(H3K36me3_genebody_pc) <- samples
cor_H3K36me3_genebody_pc <- cor(H3K36me3_genebody_pc, method="pearson")
dist_H3K36me3_genebody_pc <- as.dist(1 - cor_H3K36me3_genebody_pc)
hc_H3K36me3_genebody_pc <- as.dendrogram(hclust(dist_H3K36me3_genebody_pc, method="complete"))
plot(hc_H3K36me3_genebody_pc, main="H3K36me3 genebody clustering on pc genes", horiz=TRUE)
H3K36me3_genebody_nc <- read.delim("~/FetalBrain/HisMod/signal/hg19v65_genes_body.H3K36me3.nc", head = F, as.is = T)[, 2:7] # duplicated gene IDs, cannot use as rownames
colnames(H3K36me3_genebody_nc) <- samples
cor_H3K36me3_genebody_nc <- cor(H3K36me3_genebody_nc, method="pearson")
dist_H3K36me3_genebody_nc <- as.dist(1 - cor_H3K36me3_genebody_nc)
hc_H3K36me3_genebody_nc <- as.dendrogram(hclust(dist_H3K36me3_genebody_nc, method="complete"))
plot(hc_H3K36me3_genebody_nc, main="H3K36me3 genebody clustering on nc genes", horiz=TRUE)

# input
input_TSS2000_pc <- read.delim("~/FetalBrain/HisMod/signal/hg19v65_genes_TSS_2000.input.pc", head = F, as.is = T, row.names = 1)
colnames(input_TSS2000_pc) <- samples
cor_input_TSS2000_pc <- cor(input_TSS2000_pc, method="pearson")
dist_input_TSS2000_pc <- as.dist(1 - cor_input_TSS2000_pc)
hc_input_TSS2000_pc <- as.dendrogram(hclust(dist_input_TSS2000_pc, method="complete"))
plot(hc_input_TSS2000_pc, main="input TSS_2000 clustering on pc genes", horiz=TRUE)
input_TSS2000_nc <- read.delim("~/FetalBrain/HisMod/signal/hg19v65_genes_TSS_2000.input.nc", head = F, as.is = T, row.names = 1)
colnames(input_TSS2000_nc) <- samples
cor_input_TSS2000_nc <- cor(input_TSS2000_nc, method="pearson")
dist_input_TSS2000_nc <- as.dist(1 - cor_input_TSS2000_nc)
hc_input_TSS2000_nc <- as.dendrogram(hclust(dist_input_TSS2000_nc, method="complete"))
plot(hc_input_TSS2000_nc, main="input TSS_2000 clustering on nc genes", horiz=TRUE)
dev.off()

save(cor_H3K4me1_TSS2000_pc, cor_H3K4me1_TSS2000_nc, cor_H3K4me3_TSS2000_pc, cor_H3K4me3_TSS2000_nc, cor_H3K27me3_TSS2000_pc, cor_H3K27me3_TSS2000_nc, cor_input_TSS2000_pc, cor_input_TSS2000_nc, cor_H3K36me3_genebody_pc, cor_H3K36me3_genebody_nc, file = "cor.Rdata")

