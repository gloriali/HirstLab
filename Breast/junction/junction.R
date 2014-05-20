# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R 
# read junction files
setwd("/projects/mbilenky/REMC/breast/RNA-seq/junctions/")
lum084 <- read.delim("A17918.ZJ", head = F, as.is = T)
myo084 <- read.delim("A17919.ZJ", head = F, as.is = T)
stem084 <- read.delim("A17920.ZJ", head = F, as.is = T)

setwd("~/REMC/junction")
# check strand leakage
nrow(lum084[lum084$V2*lum084$V3 != 0,])
summary(mapply(min, lum084[lum084$V2*lum084$V3 != 0,]$V2, lum084[lum084$V2*lum084$V3 != 0,]$V3))
quantile(mapply(min, lum084[lum084$V2*lum084$V3 != 0,]$V2, lum084[lum084$V2*lum084$V3 != 0,]$V3), probs = 0.9)
summary(mapply(max, lum084[lum084$V2*lum084$V3 != 0,]$V2, lum084[lum084$V2*lum084$V3 != 0,]$V3))
quantile(mapply(max, lum084[lum084$V2*lum084$V3 != 0,]$V2, lum084[lum084$V2*lum084$V3 != 0,]$V3), probs = 0.9)
nrow(myo084[myo084$V2*myo084$V3 != 0,])
summary(mapply(min, myo084[myo084$V2*myo084$V3 != 0,]$V2, myo084[myo084$V2*myo084$V3 != 0,]$V3))
quantile(mapply(min, myo084[myo084$V2*myo084$V3 != 0,]$V2, myo084[myo084$V2*myo084$V3 != 0,]$V3), probs = 0.9)
summary(mapply(max, myo084[myo084$V2*myo084$V3 != 0,]$V2, myo084[myo084$V2*myo084$V3 != 0,]$V3))
quantile(mapply(max, myo084[myo084$V2*myo084$V3 != 0,]$V2, myo084[myo084$V2*myo084$V3 != 0,]$V3), probs = 0.9)
nrow(stem084[stem084$V2*stem084$V3 != 0,])
summary(mapply(min, stem084[stem084$V2*stem084$V3 != 0,]$V2, stem084[stem084$V2*stem084$V3 != 0,]$V3))
quantile(mapply(min, stem084[stem084$V2*stem084$V3 != 0,]$V2, stem084[stem084$V2*stem084$V3 != 0,]$V3), probs = 0.9)
summary(mapply(max, stem084[stem084$V2*stem084$V3 != 0,]$V2, stem084[stem084$V2*stem084$V3 != 0,]$V3))
quantile(mapply(max, stem084[stem084$V2*stem084$V3 != 0,]$V2, stem084[stem084$V2*stem084$V3 != 0,]$V3), probs = 0.9)

# integrating pos+neg and plot distribution
junction <- data.frame(id = lum084$V1, Ensembl = lum084$V4, lum084N = lum084$V2+lum084$V3, myo084N = myo084$V2+myo084$V3, stem084N = stem084$V2+stem084$V3)
pdf("readCount_distribution.pdf")
plot(c(0, 500), c(0, 1), type = "n", main = "Ecdf of junction read count", xlab = "coverage", ylab = "Ecdf")
lines(ecdf(junction$lum084N), col = 1, lwd = 2)
lines(ecdf(junction$myo084N), col = 2, lwd = 2)
lines(ecdf(junction$stem084N), col = 3, lwd = 2)
legend("bottomright", c("RM084_lum", "RM084_myo", "RM084_stem-like"), col = 1:3, lwd = 5)
plot(c(0, 50), c(0, 1), type = "n", main = "Ecdf of junction read count", xlab = "coverage", ylab = "Ecdf")
lines(ecdf(junction$lum084N), col = 1, lwd = 2)
lines(ecdf(junction$myo084N), col = 2, lwd = 2)
lines(ecdf(junction$stem084N), col = 3, lwd = 2)
abline(v = 2) 
legend("bottomright", c("RM084_lum", "RM084_myo", "RM084_stem-like"), col = 1:3, lwd = 5)
dev.off()

# calculate RPKM: 10^3*10^6 / (N*read_length)  
Nlum084 <- 191290094   # from A17918.report: Total number of exonic reads for RPKM (protein coding; no MT, no ribo proteins, top expressed 0.005 exons excluded) 
Nmyo084 <- 215406797
Nstem084 <- 220214822
read_length <- 75
junction$lum084rpkm <- junction$lum084N*10^3*10^6 / (Nlum084*read_length) 
junction$myo084rpkm <- junction$myo084N*10^3*10^6 / (Nmyo084*read_length) 
junction$stem084rpkm <- junction$stem084N*10^3*10^6 / (Nstem084*read_length) 
junction$gene <- gsub("_ENST[0-9A-Z:-]+[,A-Z0-9_:-]*", "", junction$Ensembl)
rownames(junction) <- junction$id
save(junction, file = "junction_RM084.Rdata")

# RPKM distribution 
setwd("~/REMC/junction/")
load("junction_RM084.Rdata")
pdf("RPKM_distribution.pdf")
plot(c(0, 50), c(0, 1), type = "n", main = "Ecdf of junction RPKM", xlab = "RPKM", ylab = "Ecdf")
lines(ecdf(junction$lum084rpkm), col = 1, lwd = 2)
lines(ecdf(junction$myo084rpkm), col = 2, lwd = 2)
lines(ecdf(junction$stem084rpkm), col = 3, lwd = 2)
abline(v = 0.1) 
legend("bottomright", c("RM084_lum", "RM084_myo", "RM084_stem-like"), col = 1:3, lwd = 5)
plot(c(0, 3), c(0, 1), type = "n", main = "Ecdf of junction RPKM", xlab = "RPKM", ylab = "Ecdf")
lines(ecdf(junction$lum084rpkm), col = 1, lwd = 2)
lines(ecdf(junction$myo084rpkm), col = 2, lwd = 2)
lines(ecdf(junction$stem084rpkm), col = 3, lwd = 2)
abline(v = 0.1) 
legend("bottomright", c("RM084_lum", "RM084_myo", "RM084_stem-like"), col = 1:3, lwd = 5)
dev.off()

# bedgraph files for visualization

#############################################################################################################################################
# Isoform Approach1: fold change + RPKM cutoff
# isoform: junction RPKM < cutoff in one sample (not present) and junction RPKM >= cutoff in the other (present) and fold change > 2
setwd("~/REMC/junction")
load("junction_RM084.Rdata")
cutoff <- 0.5
junction_fold <- matrix(NA, nrow = 3, ncol = 6, dimnames = list(c("lum084_myo084", "lum084_stem084", "myo084_stem084"), c("No.junctions.2-fold", "No.junctions.RPKM.5", "exclude.DE.genes", "No.AS.genes", "previous.genes", "intersect")))

# Ensembl gene annotation 
# ensembl <- read.delim("~/hg19/hg19v65_genes", head = F, as.is = T)
# colnames(ensembl) <- c("id", "chr", "start", "end", "strand", "type", "name", "description")
# ensembl <- ensembl[!duplicated(ensembl$id),]
# ensembl <- ensembl[!is.na(ensembl$id),]
# rownames(ensembl) <- ensembl$id
# ensembl$type <- factor(ensembl$type)
# ensembl$coord <- paste0("chr", ensembl$chr, ":", ensembl$start, "-", ensembl$end)
# save(ensembl, file = "~/hg19/hg19v65_genes.Rdata")
load("~/hg19/hg19v65_genes.Rdata")
# exclude DE genes
setwd("~/REMC/breast/gene/")
lum084_myo084_DEup <- read.delim("UP.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DEdn <- read.delim("DN.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DE <- c(lum084_myo084_DEup$V1, lum084_myo084_DEdn$V1)
lum084_stem084_DEup <- read.delim("UP.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DEdn <- read.delim("DN.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DE <- c(lum084_stem084_DEup$V1, lum084_stem084_DEdn$V1)
myo084_stem084_DEup <- read.delim("UP.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DEdn <- read.delim("DN.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DE <- c(myo084_stem084_DEup$V1, myo084_stem084_DEdn$V1)

# lum RM084 vs myo RM084
lum084_myo084_fold <- junction[(junction$myo084rpkm/junction$lum084rpkm > 2)|(junction$lum084rpkm/junction$myo084rpkm > 2),]
junction_fold[1, "No.junctions.2-fold"] <- nrow(lum084_myo084_fold)
lum084_myo084_fold <- junction[(junction$lum084rpkm < cutoff & junction$myo084rpkm >= cutoff)|(junction$lum084rpkm >= cutoff & junction$myo084rpkm < cutoff),]
junction_fold[1, "No.junctions.RPKM.5"] <- nrow(lum084_myo084_fold)
lum084_myo084_fold <- lum084_myo084_fold[!(lum084_myo084_fold$gene %in% lum084_myo084_DE),]
lum084_myo084_fold <- lum084_myo084_fold[order(-(lum084_myo084_fold$lum084rpkm + lum084_myo084_fold$myo084rpkm)),]
junction_fold[1, "exclude.DE.genes"] <- nrow(lum084_myo084_fold)
lum084_myo084_fold_gene <- unique(gsub("_ENST[0-9:-]+", "", unlist(strsplit(as.character(lum084_myo084_fold$Ensembl), ","))))
junction_fold[1, "No.AS.genes"] <- length(lum084_myo084_fold_gene)
lum084_myo084_fold_transcript <- unique(gsub("ENSG[0-9]+_", "", gsub(":-*1", "", unlist(strsplit(as.character(lum084_myo084_fold$Ensembl), ",")))))
length(lum084_myo084_fold_transcript)
# compare results with previous isoform analysis with exon RPKM 
lum084_myo084_fold_gene_previous <- read.csv("~/REMC/breast/tissue/myo_lum_84_isoform_only.csv", head = T, as.is = T)
junction_fold[1, "previous.genes"] <- nrow(lum084_myo084_fold_gene_previous)
lum084_myo084_fold_gene_intersect <- intersect(lum084_myo084_fold_gene_previous$x, lum084_myo084_fold_gene)
junction_fold[1, "intersect"] <- length(lum084_myo084_fold_gene_intersect)

# lum RM084 vs stem RM084
lum084_stem084_fold <- junction[(junction$stem084rpkm/junction$lum084rpkm > 2)|(junction$lum084rpkm/junction$stem084rpkm > 2),]
junction_fold[2, "No.junctions.2-fold"] <- nrow(lum084_stem084_fold)
lum084_stem084_fold <- junction[(junction$lum084rpkm < cutoff & junction$stem084rpkm >= cutoff)|(junction$lum084rpkm >= cutoff & junction$stem084rpkm < cutoff),]
junction_fold[2, "No.junctions.RPKM.5"] <- nrow(lum084_stem084_fold)
lum084_stem084_fold <- lum084_stem084_fold[!(lum084_stem084_fold$gene %in% lum084_stem084_DE),]
lum084_stem084_fold <- lum084_stem084_fold[order(-(lum084_stem084_fold$lum084rpkm + lum084_stem084_fold$stem084rpkm)),]
junction_fold[2, "exclude.DE.genes"] <- nrow(lum084_stem084_fold)
lum084_stem084_fold_gene <- unique(gsub("_ENST[0-9:-]+", "", unlist(strsplit(as.character(lum084_stem084_fold$Ensembl), ","))))
junction_fold[2, "No.AS.genes"] <- length(lum084_stem084_fold_gene)
lum084_stem084_fold_transcript <- unique(gsub("ENSG[0-9]+_", "", gsub(":-*1", "", unlist(strsplit(as.character(lum084_stem084_fold$Ensembl), ",")))))
length(lum084_stem084_fold_transcript)
# compare results with previous isoform analysis with exon RPKM 
lum084_stem084_fold_gene_previous <- read.csv("~/REMC/breast/tissue/lum_stem_84_isoform_only.csv", head = T, as.is = T)
junction_fold[2, "previous.genes"] <- nrow(lum084_stem084_fold_gene_previous)
lum084_stem084_fold_gene_intersect <- intersect(lum084_stem084_fold_gene_previous$x, lum084_stem084_fold_gene)
junction_fold[2, "intersect"] <- length(lum084_stem084_fold_gene_intersect)

# myo RM084 vs stem RM084
myo084_stem084_fold <- junction[(junction$stem084rpkm/junction$myo084rpkm > 2)|(junction$myo084rpkm/junction$stem084rpkm > 2),]
junction_fold[3, "No.junctions.2-fold"] <- nrow(myo084_stem084_fold)
myo084_stem084_fold <- junction[(junction$myo084rpkm < cutoff & junction$stem084rpkm >= cutoff)|(junction$myo084rpkm >= cutoff & junction$stem084rpkm < cutoff),]
junction_fold[3, "No.junctions.RPKM.5"] <- nrow(myo084_stem084_fold)
myo084_stem084_fold <- myo084_stem084_fold[!(myo084_stem084_fold$gene %in% myo084_stem084_DE),]
myo084_stem084_fold <- myo084_stem084_fold[order(-(myo084_stem084_fold$myo084rpkm + myo084_stem084_fold$stem084rpkm)),]
junction_fold[3, "exclude.DE.genes"] <- nrow(myo084_stem084_fold)
myo084_stem084_fold_gene <- unique(gsub("_ENST[0-9:-]+", "", unlist(strsplit(as.character(myo084_stem084_fold$Ensembl), ","))))
junction_fold[3, "No.AS.genes"] <- length(myo084_stem084_fold_gene)
myo084_stem084_fold_transcript <- unique(gsub("ENSG[0-9]+_", "", gsub(":-*1", "", unlist(strsplit(as.character(myo084_stem084_fold$Ensembl), ",")))))
length(myo084_stem084_fold_transcript)
# compare results with previous isoform analysis with exon RPKM 
myo084_stem084_fold_gene_previous <- read.csv("~/REMC/breast/tissue/myo_stem_84_isoform_only.csv", head = T, as.is = T)
junction_fold[3, "previous.genes"] <- nrow(myo084_stem084_fold_gene_previous)
myo084_stem084_fold_gene_intersect <- intersect(myo084_stem084_fold_gene_previous$x, myo084_stem084_fold_gene)
junction_fold[3, "intersect"] <- length(myo084_stem084_fold_gene_intersect)

#############################################################################################################################################
# Isoform Approach2: DEfine + RPKM = 0.5 cutoff
# apply RPKM = 0.5 cutoff to DEfine output
cutoff <- 0.5
junction_DEfine <- matrix(NA, nrow = 3, ncol = 6, dimnames = list(c("lum084_myo084", "lum084_stem084", "myo084_stem084"), c("No.junctions.DEfine", "No.junctions.RPKM.5", "exclude.DE.genes", "No.AS.genes", "previous.genes", "intersect")))

# write DEfine input files
# setwd("~/REMC/junction/")
# load("junction_RM084.Rdata")
# write.table(junction[,c("id", "lum084N", "lum084rpkm")], file = "lum_RM084_DEfine.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# write.table(junction[,c("id", "myo084N", "myo084rpkm")], file = "myo_RM084_DEfine.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# write.table(junction[,c("id", "stem084N", "stem084rpkm")], file = "stem_RM084_DEfine.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# exclude DE genes
load("~/hg19/hg19v65_genes.Rdata")
setwd("~/REMC/breast/gene/")
lum084_myo084_DEup <- read.delim("UP.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DEdn <- read.delim("DN.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DE <- c(lum084_myo084_DEup$V1, lum084_myo084_DEdn$V1)
lum084_stem084_DEup <- read.delim("UP.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DEdn <- read.delim("DN.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DE <- c(lum084_stem084_DEup$V1, lum084_stem084_DEdn$V1)
myo084_stem084_DEup <- read.delim("UP.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DEdn <- read.delim("DN.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DE <- c(myo084_stem084_DEup$V1, myo084_stem084_DEdn$V1)

setwd("~/REMC/junction/")
# lum084 vs myo084
lum084_myo084_DEfineup <- read.delim("./DEfine/UP.lum-RM084_myo-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DEfinedn <- read.delim("./DEfine/DN.lum-RM084_myo-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DEfine <- rbind(lum084_myo084_DEfineup, lum084_myo084_DEfinedn)
colnames(lum084_myo084_DEfine) <- c("id", "lum084rpkm", "myo084rpkm", "p_value", "corrected_p")
rownames(lum084_myo084_DEfine) <- lum084_myo084_DEfine$id
junction_DEfine[1, "No.junctions.DEfine"] <- nrow(lum084_myo084_DEfine)
lum084_myo084_DEfine <- lum084_myo084_DEfine[(lum084_myo084_DEfine$lum084rpkm < cutoff & lum084_myo084_DEfine$myo084rpkm >= cutoff)|(lum084_myo084_DEfine$lum084rpkm >= cutoff & lum084_myo084_DEfine$myo084rpkm < cutoff),]
junction_DEfine[1, "No.junctions.RPKM.5"] <- nrow(lum084_myo084_DEfine)
lum084_myo084_DEfine$Ensembl <- junction[lum084_myo084_DEfine$id,]$Ensembl
lum084_myo084_DEfine$gene <- junction[lum084_myo084_DEfine$id,]$gene
lum084_myo084_DEfine <- lum084_myo084_DEfine[!(lum084_myo084_DEfine$gene %in% lum084_myo084_DE),]
junction_DEfine[1, "exclude.DE.genes"] <- nrow(lum084_myo084_DEfine)
lum084_myo084_DEfine <- lum084_myo084_DEfine[order(lum084_myo084_DEfine$corrected_p, -(lum084_myo084_DEfine$lum084rpkm + lum084_myo084_DEfine$myo084rpkm)),]
lum084_myo084_DEfine$rank <- 1:nrow(lum084_myo084_DEfine)
lum084_myo084_DEfine_gene <- unique(gsub("_ENST[0-9:-]+", "", unlist(strsplit(as.character(lum084_myo084_DEfine$Ensembl), ","))))
length(lum084_myo084_DEfine_gene)
lum084_myo084_DEfine_gene <- data.frame(Ensemble = lum084_myo084_DEfine_gene, minRank = sapply(lum084_myo084_DEfine_gene, function(x) min(c(lum084_myo084_DEfine[grepl(x, lum084_myo084_DEfine$Ensembl),]$rank))))
lum084_myo084_DEfine_gene <- lum084_myo084_DEfine_gene[order(lum084_myo084_DEfine_gene$minRank),]
lum084_myo084_DEfine_gene[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[lum084_myo084_DEfine_gene$Ensembl,c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
junction_DEfine[1, "No.AS.genes"] <- nrow(lum084_myo084_DEfine_gene)
lum084_myo084_DEfine_transcript <- unique(gsub("ENSG[0-9]+_", "", gsub(":-*1", "", unlist(strsplit(as.character(lum084_myo084_DEfine$Ensembl), ",")))))
length(lum084_myo084_DEfine_transcript)
write.table(lum084_myo084_DEfine, file = "lum084_myo084_DEfine_isoforms.txt", sep = "\t", quote = F, row.names = F)
write.table(lum084_myo084_DEfine_gene, file = "lum084_myo084_DEfine_isoforms_gene.txt", sep = "\t", quote = F, row.names = F)
write.table(lum084_myo084_DEfine_transcript, file = "lum084_myo084_DEfine_isoforms_transcript.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# compare results with previous isoform analysis with exon RPKM 
lum084_myo084_DEfine_gene_previous <- read.csv("~/REMC/breast/tissue/myo_lum_84_isoform_only.csv", head = T, as.is = T)
junction_DEfine[1, "previous.genes"] <- nrow(lum084_myo084_DEfine_gene_previous)
lum084_myo084_DEfine_gene_intersect <- intersect(lum084_myo084_DEfine_gene_previous$x, lum084_myo084_DEfine_gene$Ensemble)
junction_DEfine[1, "intersect"] <- length(lum084_myo084_DEfine_gene_intersect)

# lum084 vs stem084
lum084_stem084_DEfineup <- read.delim("./DEfine/UP.lum-RM084_stem-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DEfinedn <- read.delim("./DEfine/DN.lum-RM084_stem-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DEfine <- rbind(lum084_stem084_DEfineup, lum084_stem084_DEfinedn)
colnames(lum084_stem084_DEfine) <- c("id", "lum084rpkm", "stem084rpkm", "p_value", "corrected_p")
rownames(lum084_stem084_DEfine) <- lum084_stem084_DEfine$id
junction_DEfine[2, "No.junctions.DEfine"] <- nrow(lum084_stem084_DEfine)
lum084_stem084_DEfine <- lum084_stem084_DEfine[(lum084_stem084_DEfine$lum084rpkm < cutoff & lum084_stem084_DEfine$stem084rpkm >= cutoff)|(lum084_stem084_DEfine$lum084rpkm >= cutoff & lum084_stem084_DEfine$stem084rpkm < cutoff),]
junction_DEfine[2, "No.junctions.RPKM.5"] <- nrow(lum084_stem084_DEfine)
lum084_stem084_DEfine$Ensembl <- junction[lum084_stem084_DEfine$id,]$Ensembl
lum084_stem084_DEfine$gene <- junction[lum084_stem084_DEfine$id,]$gene
lum084_stem084_DEfine <- lum084_stem084_DEfine[!(lum084_stem084_DEfine$gene %in% lum084_stem084_DE),]
junction_DEfine[2, "exclude.DE.genes"] <- nrow(lum084_stem084_DEfine)
lum084_stem084_DEfine <- lum084_stem084_DEfine[order(lum084_stem084_DEfine$corrected_p, -(lum084_stem084_DEfine$lum084rpkm + lum084_stem084_DEfine$stem084rpkm)),]
lum084_stem084_DEfine$rank <- 1:nrow(lum084_stem084_DEfine)
lum084_stem084_DEfine_gene <- unique(gsub("_ENST[0-9:-]+", "", unlist(strsplit(as.character(lum084_stem084_DEfine$Ensembl), ","))))
length(lum084_stem084_DEfine_gene)
lum084_stem084_DEfine_gene <- data.frame(Ensemble = lum084_stem084_DEfine_gene, minRank = sapply(lum084_stem084_DEfine_gene, function(x) min(c(lum084_stem084_DEfine[grepl(x, lum084_stem084_DEfine$Ensembl),]$rank))))
lum084_stem084_DEfine_gene <- lum084_stem084_DEfine_gene[order(lum084_stem084_DEfine_gene$minRank),]
lum084_stem084_DEfine_gene[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[lum084_stem084_DEfine_gene$Ensembl,c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
junction_DEfine[2, "No.AS.genes"] <- nrow(lum084_stem084_DEfine_gene)
lum084_stem084_DEfine_transcript <- unique(gsub("ENSG[0-9]+_", "", gsub(":-*1", "", unlist(strsplit(as.character(lum084_stem084_DEfine$Ensembl), ",")))))
length(lum084_stem084_DEfine_transcript)
write.table(lum084_stem084_DEfine, file = "lum084_stem084_DEfine_isoforms.txt", sep = "\t", quote = F, row.names = F)
write.table(lum084_stem084_DEfine_gene, file = "lum084_stem084_DEfine_isoforms_gene.txt", sep = "\t", quote = F, row.names = F)
write.table(lum084_stem084_DEfine_transcript, file = "lum084_stem084_DEfine_isoforms_transcript.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# compare results with previous isoform analysis with exon RPKM 
lum084_stem084_DEfine_gene_previous <- read.csv("~/REMC/breast/tissue/lum_stem_84_isoform_only.csv", head = T, as.is = T)
junction_DEfine[2, "previous.genes"] <- nrow(lum084_stem084_DEfine_gene_previous)
lum084_stem084_DEfine_gene_intersect <- intersect(lum084_stem084_DEfine_gene_previous$x, lum084_stem084_DEfine_gene$Ensemble)
junction_DEfine[2, "intersect"] <- length(lum084_stem084_DEfine_gene_intersect)

# myo084 vs stem084
myo084_stem084_DEfineup <- read.delim("./DEfine/UP.myo-RM084_stem-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DEfinedn <- read.delim("./DEfine/DN.myo-RM084_stem-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DEfine <- rbind(myo084_stem084_DEfineup, myo084_stem084_DEfinedn)
colnames(myo084_stem084_DEfine) <- c("id", "myo084rpkm", "stem084rpkm", "p_value", "corrected_p")
rownames(myo084_stem084_DEfine) <- myo084_stem084_DEfine$id
junction_DEfine[3, "No.junctions.DEfine"] <- nrow(myo084_stem084_DEfine)
myo084_stem084_DEfine <- myo084_stem084_DEfine[(myo084_stem084_DEfine$myo084rpkm < cutoff & myo084_stem084_DEfine$stem084rpkm >= cutoff)|(myo084_stem084_DEfine$myo084rpkm >= cutoff & myo084_stem084_DEfine$stem084rpkm < cutoff),]
junction_DEfine[3, "No.junctions.RPKM.5"] <- nrow(myo084_stem084_DEfine)
myo084_stem084_DEfine$Ensembl <- junction[myo084_stem084_DEfine$id,]$Ensembl
myo084_stem084_DEfine$gene <- junction[myo084_stem084_DEfine$id,]$gene
myo084_stem084_DEfine <- myo084_stem084_DEfine[!(myo084_stem084_DEfine$gene %in% myo084_stem084_DE),]
junction_DEfine[3, "exclude.DE.genes"] <- nrow(myo084_stem084_DEfine)
myo084_stem084_DEfine <- myo084_stem084_DEfine[order(myo084_stem084_DEfine$corrected_p, -(myo084_stem084_DEfine$myo084rpkm + myo084_stem084_DEfine$stem084rpkm)),]
myo084_stem084_DEfine$rank <- 1:nrow(myo084_stem084_DEfine)
myo084_stem084_DEfine_gene <- unique(gsub("_ENST[0-9:-]+", "", unlist(strsplit(as.character(myo084_stem084_DEfine$Ensembl), ","))))
length(myo084_stem084_DEfine_gene)
myo084_stem084_DEfine_gene <- data.frame(Ensemble = myo084_stem084_DEfine_gene, minRank = sapply(myo084_stem084_DEfine_gene, function(x) min(c(myo084_stem084_DEfine[grepl(x, myo084_stem084_DEfine$Ensembl),]$rank))))
myo084_stem084_DEfine_gene <- myo084_stem084_DEfine_gene[order(myo084_stem084_DEfine_gene$minRank),]
myo084_stem084_DEfine_gene[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[myo084_stem084_DEfine_gene$Ensembl,c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
junction_DEfine[3, "No.AS.genes"] <- nrow(myo084_stem084_DEfine_gene)
myo084_stem084_DEfine_transcript <- unique(gsub("ENSG[0-9]+_", "", gsub(":-*1", "", unlist(strsplit(as.character(myo084_stem084_DEfine$Ensembl), ",")))))
length(myo084_stem084_DEfine_transcript)
write.table(myo084_stem084_DEfine, file = "myo084_stem084_DEfine_isoforms.txt", sep = "\t", quote = F, row.names = F)
write.table(myo084_stem084_DEfine_gene, file = "myo084_stem084_DEfine_isoforms_gene.txt", sep = "\t", quote = F, row.names = F)
write.table(myo084_stem084_DEfine_transcript, file = "myo084_stem084_DEfine_isoforms_transcript.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# compare results with previous isoform analysis with exon RPKM 
myo084_stem084_DEfine_gene_previous <- read.csv("~/REMC/breast/tissue/myo_stem_84_isoform_only.csv", head = T, as.is = T)
junction_DEfine[3, "previous.genes"] <- nrow(myo084_stem084_DEfine_gene_previous)
myo084_stem084_DEfine_gene_intersect <- intersect(myo084_stem084_DEfine_gene_previous$x, myo084_stem084_DEfine_gene$Ensemble)
junction_DEfine[3, "intersect"] <- length(myo084_stem084_DEfine_gene_intersect)

# intersecting across cell types
length(intersect(lum084_myo084_gene$Ensembl, lum084_stem084_gene$Ensembl))
length(intersect(lum084_myo084_gene$Ensembl, myo084_stem084_gene$Ensembl))
length(intersect(myo084_stem084_gene$Ensembl, lum084_stem084_gene$Ensembl))
lum_myo_stem_084 <- intersect(lum084_myo084_gene$Ensembl, intersect(myo084_stem084_gene$Ensembl, lum084_stem084_gene$Ensembl))
length(lum_myo_stem_084)

# compare DEfine+cutoff with foldchange+cutoff
length(intersect(lum084_myo084_fold$id, lum084_myo084$id))
length(intersect(lum084_myo084_fold_gene, lum084_myo084_gene$Ensemble))
lum084_myo084_foldnotDE <- lum084_myo084_fold[!(lum084_myo084_fold$id %in% lum084_myo084$id),]
lum084_myo084_foldnotDE <- lum084_myo084_foldnotDE[order(-(lum084_myo084_foldnotDE$lum084rpkm+lum084_myo084_foldnotDE$myo084rpkm)),]
summary(lum084_myo084_foldnotDE$lum084rpkm+lum084_myo084_foldnotDE$myo084rpkm)
length(intersect(lum084_stem084_fold$id, lum084_stem084$id))
length(intersect(lum084_stem084_fold_gene, lum084_stem084_gene$Ensemble))
lum084_stem084_foldnotDE <- lum084_stem084_fold[!(lum084_stem084_fold$id %in% lum084_stem084$id),]
lum084_stem084_foldnotDE <- lum084_stem084_foldnotDE[order(-(lum084_stem084_foldnotDE$lum084rpkm+lum084_stem084_foldnotDE$stem084rpkm)),]
summary(lum084_stem084_foldnotDE$lum084rpkm+lum084_stem084_foldnotDE$stem084rpkm)
length(intersect(myo084_stem084_fold$id, myo084_stem084$id))
length(intersect(myo084_stem084_fold_gene, myo084_stem084_gene$Ensemble))
myo084_stem084_foldnotDE <- myo084_stem084_fold[!(myo084_stem084_fold$id %in% myo084_stem084$id),]
myo084_stem084_foldnotDE <- myo084_stem084_foldnotDE[order(-(myo084_stem084_foldnotDE$myo084rpkm+myo084_stem084_foldnotDE$stem084rpkm)),]
summary(myo084_stem084_foldnotDE$myo084rpkm+myo084_stem084_foldnotDE$stem084rpkm)

###############################################################################################################################################
# estimate gene RPKM with average junction RPKM of the gene
# read gene RPKM
# /gsc/software/linux-x86_64-centos5/R-2.15.2/bin/R
setwd("~/REMC/junction/")
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'
libs=c("A01029",
       "A01030",
       "A01031",
       "A17918",
       "A17919",
       "A17920",
       "A18760",
       "HS1187",
       "HS1188",
       "HS2263",
       "A18472",
       "A18761")

pc_lib <- read.table(paste0(dirIn, "A01029", "/coverage/", "A01029", ".G.A.rpkm.pc"))
pc <- matrix(data = NA, ncol = length(libs), nrow = nrow(pc_lib))
colnames(pc) <- c("lum080", "myo080", "stem080", "lum084", "myo084", "stem084", "fibr070", "lum035", "myo035", "vHMEC035", "fibr071", "vHMEC071")
rownames(pc) <- pc_lib$V1
i <- 1
for(lib in libs){
  pc_lib <- read.table(paste0(dirIn, lib, "/coverage/", lib, ".G.A.rpkm.pc"))
  pc[,i] <- pc_lib$V3
  i=i+1
}
nc_lib <- read.table(paste0(dirIn, "A01029", "/coverage/", "A01029", ".G.A.rpkm.nc"))
nc <- matrix(data = NA, ncol = length(libs), nrow = nrow(nc_lib))
colnames(nc) <- c("lum080", "myo080", "stem080", "lum084", "myo084", "stem084", "fibr070", "lum035", "myo035", "vHMEC035", "fibr071", "vHMEC071")
rownames(nc) <- nc_lib$V1
i <- 1
for(lib in libs){
  nc_lib <- read.table(paste0(dirIn, lib, "/coverage/", lib, ".G.A.rpkm.nc"))
  nc[,i] <- nc_lib$V3
  i=i+1
}
geneRPKM <- as.data.frame(rbind(pc, nc))
save(geneRPKM, file = "geneRPKM.Rdata")

# load data
load("geneRPKM.Rdata")
junction$gene <- factor(junction$gene)
junction[,c("lum084gene", "myo084gene", "stem084gene")] <- geneRPKM[as.character(junction$gene),c("lum084", "myo084", "stem084")]
save(junction, file = "junction_RM084.Rdata")

# gene RPKM vs junction estimated RPKM
lum084gene_junction_rpkm <- aggregate(lum084rpkm ~ gene, data = junction, mean, na.action = na.omit)
lum084gene_junction_rpkm$geneRPKM <- geneRPKM[as.character(lum084gene_junction_rpkm$gene),]$lum084
summary(lum084gene_junction_rpkm)
summary(junction$lum084rpkm)
myo084gene_junction_rpkm <- aggregate(myo084rpkm ~ gene, data = junction, mean, na.action = na.omit)
myo084gene_junction_rpkm$geneRPKM <- geneRPKM[as.character(myo084gene_junction_rpkm$gene),]$myo084
summary(myo084gene_junction_rpkm)
summary(junction$myo084rpkm)
stem084gene_junction_rpkm <- aggregate(stem084rpkm ~ gene, data = junction, mean, na.action = na.omit)
stem084gene_junction_rpkm$geneRPKM <- geneRPKM[as.character(stem084gene_junction_rpkm$gene),]$stem084
summary(stem084gene_junction_rpkm)
summary(junction$stem084rpkm)

pdf("geneRPKM_junctionEstimated.pdf")
smoothScatter(log10(lum084gene_junction_rpkm$geneRPKM) ~ log10(lum084gene_junction_rpkm$lum084rpkm), main = "gene RPKM vs junction estimated RPKM - lum084", xlab = "log10(junction estimated RPKM)", ylab = "log10(gene RPKM)")
smoothScatter(log10(myo084gene_junction_rpkm$geneRPKM) ~ log10(myo084gene_junction_rpkm$myo084rpkm), main = "gene RPKM vs junction estimated RPKM - myo084", xlab = "log10(junction estimated RPKM)", ylab = "log10(gene RPKM)")
smoothScatter(log10(stem084gene_junction_rpkm$geneRPKM) ~ log10(stem084gene_junction_rpkm$stem084rpkm), main = "gene RPKM vs junction estimated RPKM - stem084", xlab = "log10(junction estimated RPKM)", ylab = "log10(gene RPKM)")
plot(c(0, 5), c(0, 1), type = "n", main = "Ecdf of junction RPKM / gene RPKM (all junctions)", xlab = "junctionRPKM/geneRPKM", ylab = "Ecdf")
lines(ecdf(na.omit((junction$lum084rpkm)/(junction$lum084gene))), col = 1, lwd = 2)
lines(ecdf(na.omit((junction$myo084rpkm)/(junction$myo084gene))), col = 2, lwd = 2)
lines(ecdf(na.omit((junction$stem084rpkm)/(junction$stem084gene))), col = 3, lwd = 2)
legend("bottomright", c("lum084", "myo084", "stem084"), col = 1:3, lwd = 5)
dev.off()

################################################################################################################################################
# Isoform Approach3: use relative junction expression: = junctionRPKM/geneRPKM
# fold change relative junction expression > 2 & at least one sample relative junction expression > 0.5 & both gene RPKM > 0.1
setwd("~/REMC/junction/")
load("junction_RM084.Rdata")
junction_relative <- matrix(NA, nrow = 3, ncol = 7, dimnames = list(c("lum084_myo084", "lum084_stem084", "myo084_stem084"), c("relative.4-fold", "junctions.relative.5", "geneRPKM.1", "exclude.DE.genes", "No.AS.genes", "previous.genes", "intersect")))
rcutoff <- 0.5
genecutoff <- 0.1

e <- 10^-6
junction$lum084relative <- (junction$lum084rpkm+e)/(junction$lum084gene+e)
junction$myo084relative <- (junction$myo084rpkm+e)/(junction$myo084gene+e)
junction$stem084relative <- (junction$stem084rpkm+e)/(junction$stem084gene+e)

pdf("relative_junctionExpression.pdf")
plot(c(0, 0.002), c(0, 0.5), type = "n", main = "Ecdf of junction expression relative to the gene", xlab = "junctionRPKM/geneRPKM", ylab = "Ecdf")
lines(ecdf(na.omit(junction$lum084relative)), col = 1, lwd = 2)
lines(ecdf(na.omit(junction$myo084relative)), col = 2, lwd = 2)
lines(ecdf(na.omit(junction$stem084relative)), col = 3, lwd = 2)
legend("bottomright", c("lum084", "myo084", "stem084"), col = 1:3, lwd = 5)
plot(c(-3, 3), c(0, 1), type = "n", main = "Ecdf of fold change in relative junction expression", xlab = "log2(fold_change(junctionRPKM/geneRPKM))", ylab = "Ecdf")
lines(ecdf(na.omit(log2(junction$lum084relative/junction$myo084relative))), col = 1, lwd = 2)
lines(ecdf(na.omit(log2(junction$lum084relative/junction$stem084relative))), col = 2, lwd = 2)
lines(ecdf(na.omit(log2(junction$myo084relative/junction$stem084relative))), col = 3, lwd = 2)
legend("bottomright", c("lum084/myo084", "lum084/stem084", "myo084/stem084"), col = 1:3, lwd = 5)
abline(v = 0)
dev.off()

# exclude DE genes
load("~/hg19/hg19v65_genes.Rdata")
setwd("~/REMC/breast/gene/")
lum084_myo084_DEup <- read.delim("UP.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DEdn <- read.delim("DN.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DE <- c(lum084_myo084_DEup$V1, lum084_myo084_DEdn$V1)
lum084_stem084_DEup <- read.delim("UP.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DEdn <- read.delim("DN.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DE <- c(lum084_stem084_DEup$V1, lum084_stem084_DEdn$V1)
myo084_stem084_DEup <- read.delim("UP.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DEdn <- read.delim("DN.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DE <- c(myo084_stem084_DEup$V1, myo084_stem084_DEdn$V1)

setwd("~/REMC/junction/")
# lum084 vs myo084
lum084_myo084_relative <- junction[abs(log2(junction$lum084relative/junction$myo084relative)) > 2,]
junction_relative[1, "relative.4-fold"] <- nrow(lum084_myo084_relative)
lum084_myo084_relative <- lum084_myo084_relative[lum084_myo084_relative$lum084relative > rcutoff | lum084_myo084_relative$myo084relative > rcutoff, ]
junction_relative[1, "junctions.relative.5"] <- nrow(lum084_myo084_relative)
lum084_myo084_relative <- lum084_myo084_relative[lum084_myo084_relative$lum084gene > genecutoff & lum084_myo084_relative$myo084gene > genecutoff, ]
junction_relative[1, "geneRPKM.1"] <- nrow(lum084_myo084_relative)
lum084_myo084_relative <- lum084_myo084_relative[!(lum084_myo084_relative$gene %in% lum084_myo084_DE),]
lum084_myo084_relative <- lum084_myo084_relative[order(abs(log2(lum084_myo084_relative$lum084relative/lum084_myo084_relative$myo084relative)), decreasing = T),]
junction_relative[1, "exclude.DE.genes"] <- nrow(lum084_myo084_relative)
lum084_myo084_relative_gene <- lum084_myo084_relative[!duplicated(lum084_myo084_relative$gene), c("gene", "lum084gene", "myo084gene", "lum084relative", "myo084relative", "id")]
lum084_myo084_relative_gene[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum084_myo084_relative_gene$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
junction_relative[1, "No.AS.genes"] <- nrow(lum084_myo084_relative_gene)
write.table(lum084_myo084_relative_gene, file = "lum084_myo084_relative_gene_relative.txt", sep = "\t", quote = F, row.names = F)
lum084_myo084_relative_gene_previous <- read.csv("~/REMC/breast/tissue/myo_lum_84_isoform_only.csv", head = T, as.is = T)
junction_relative[1, "previous.genes"] <- nrow(lum084_myo084_relative_gene_previous)
lum084_myo084_relative_gene_intersect <- intersect(lum084_myo084_relative_gene_previous$x, lum084_myo084_relative_gene$gene)
junction_relative[1, "intersect"] <- length(lum084_myo084_relative_gene_intersect)

# lum084 vs stem084
lum084_stem084_relative <- junction[abs(log2(junction$lum084relative/junction$stem084relative)) > 2,]
junction_relative[2, "relative.4-fold"] <- nrow(lum084_stem084_relative)
lum084_stem084_relative <- lum084_stem084_relative[lum084_stem084_relative$lum084relative > rcutoff | lum084_stem084_relative$stem084relative > rcutoff, ]
junction_relative[2, "junctions.relative.5"] <- nrow(lum084_stem084_relative)
lum084_stem084_relative <- lum084_stem084_relative[lum084_stem084_relative$lum084gene > genecutoff & lum084_stem084_relative$stem084gene > genecutoff, ]
junction_relative[2, "geneRPKM.1"] <- nrow(lum084_stem084_relative)
lum084_stem084_relative <- lum084_stem084_relative[!(lum084_stem084_relative$gene %in% lum084_stem084_DE),]
lum084_stem084_relative <- lum084_stem084_relative[order(abs(log2(lum084_stem084_relative$lum084relative/lum084_stem084_relative$stem084relative)), decreasing = T),]
junction_relative[2, "exclude.DE.genes"] <- nrow(lum084_stem084_relative)
lum084_stem084_relative_gene <- lum084_stem084_relative[!duplicated(lum084_stem084_relative$gene), c("gene", "lum084gene", "stem084gene", "lum084relative", "stem084relative", "id")]
lum084_stem084_relative_gene[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum084_stem084_relative_gene$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
junction_relative[2, "No.AS.genes"] <- nrow(lum084_stem084_relative_gene)
write.table(lum084_stem084_relative_gene, file = "lum084_stem084_relative_gene_relative.txt", sep = "\t", quote = F, row.names = F)
lum084_stem084_relative_gene_previous <- read.csv("~/REMC/breast/tissue/lum_stem_84_isoform_only.csv", head = T, as.is = T)
junction_relative[2, "previous.genes"] <- nrow(lum084_stem084_relative_gene_previous)
lum084_stem084_relative_gene_intersect <- intersect(lum084_stem084_relative_gene_previous$x, lum084_stem084_relative_gene$gene)
junction_relative[2, "intersect"] <- length(lum084_stem084_relative_gene_intersect)

# myo084 vs stem084
myo084_stem084_relative <- junction[abs(log2(junction$myo084relative/junction$stem084relative)) > 2,]
junction_relative[3, "relative.4-fold"] <- nrow(myo084_stem084_relative)
myo084_stem084_relative <- myo084_stem084_relative[myo084_stem084_relative$myo084relative > rcutoff | myo084_stem084_relative$stem084relative > rcutoff, ]
junction_relative[3, "junctions.relative.5"] <- nrow(myo084_stem084_relative)
myo084_stem084_relative <- myo084_stem084_relative[myo084_stem084_relative$myo084gene > genecutoff & myo084_stem084_relative$stem084gene > genecutoff, ]
junction_relative[3, "geneRPKM.1"] <- nrow(myo084_stem084_relative)
myo084_stem084_relative <- myo084_stem084_relative[!(myo084_stem084_relative$gene %in% myo084_stem084_DE),]
myo084_stem084_relative <- myo084_stem084_relative[order(abs(log2(myo084_stem084_relative$myo084relative/myo084_stem084_relative$stem084relative)), decreasing = T),]
junction_relative[3, "exclude.DE.genes"] <- nrow(myo084_stem084_relative)
myo084_stem084_relative_gene <- myo084_stem084_relative[!duplicated(myo084_stem084_relative$gene), c("gene", "myo084gene", "stem084gene", "myo084relative", "stem084relative", "id")]
myo084_stem084_relative_gene[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(myo084_stem084_relative_gene$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
junction_relative[3, "No.AS.genes"] <- nrow(myo084_stem084_relative_gene)
write.table(myo084_stem084_relative_gene, file = "myo084_stem084_relative_gene_relative.txt", sep = "\t", quote = F, row.names = F)
myo084_stem084_relative_gene_previous <- read.csv("~/REMC/breast/tissue/myo_stem_84_isoform_only.csv", head = T, as.is = T)
junction_relative[3, "previous.genes"] <- nrow(myo084_stem084_relative_gene_previous)
myo084_stem084_relative_gene_intersect <- intersect(myo084_stem084_relative_gene_previous$x, myo084_stem084_relative_gene$gene)
junction_relative[3, "intersect"] <- length(myo084_stem084_relative_gene_intersect)

################################################################################################################################################
# Isoform Approach4: combined approach: exon DE & at least one junction in the exon DE & gene NOT DE
junction_combined <- matrix(NA, nrow = 3, ncol = 5, dimnames = list(c("lum084_myo084", "lum084_stem084", "myo084_stem084"), c("DEfine.exons", "exclude.DE.genes", "with.annotated.junction", "with.DEfine.junctions", "exon.junction.correlate")))
# sed 's/\|/   /g' hg19v65_junctions_vs_exons_for_genes_pos.relations > hg19v65_junctions_vs_exons_for_genes_pos.relations1
# sed 's/\|/   /g' hg19v65_junctions_vs_exons_for_genes_neg.relations > hg19v65_junctions_vs_exons_for_genes_neg.relations1
junction_exon_pos <- read.table("~/hg19/hg19v65_junctions_vs_exons_for_genes_pos.relations1", head = F, as.is = T, sep = " ")
junction_exon_neg <- read.table("~/hg19/hg19v65_junctions_vs_exons_for_genes_neg.relations1", head = F, as.is = T, sep = " ")
junction_exon <- rbind(junction_exon_pos, junction_exon_neg)
junction_exon$chr <- paste0("chr", gsub(":[0-9-]+", "", junction_exon$V1))
junction_exon$junctionID <- paste0(junction_exon$V1, "-", junction_exon$V4)
junction_exon$exonID <- paste0(junction_exon$chr, ":", junction_exon$V13, "<", junction_exon$V16, "_", junction_exon$V19)
junction_exon$ENSG <- junction_exon$V19
junction_exon$junction <- gsub("\t[0-9A-Z.]+", "", junction_exon$V10)
junction_exon <- junction_exon[, c("chr", "junctionID", "exonID", "ENSG", "junction")]
junction_exon <- junction_exon[order(junction_exon$junctionID),]
write.table(junction_exon, file = "~/hg19/hg19v65_junctions_vs_exons_for_genes.relations", sep = "\t", row.names = F, quote = F)
junction_exon$junction <- NULL
junction_exon <- unique(junction_exon)
# junction_exon$count <- 1
# count <- aggregate(count ~ as.factor(junctionID), junction_exon, sum)
# count <- na.omit(count)
# summary(as.factor(count$count))
write.table(junction_exon, file = "~/hg19/hg19v65_junctions_vs_exons_for_genes.relations.unique", sep = "\t", row.names = F, quote = F)

junction_exon <- read.delim("~/hg19/hg19v65_junctions_vs_exons_for_genes.relations", head = T, as.is = T)
load("~/hg19/hg19v65_genes.Rdata")

# DE genes FDR = 0.015
setwd("~/REMC/breast/gene/")
lum084_myo084_DEup <- read.delim("UP.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DEdn <- read.delim("DN.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DE <- c(lum084_myo084_DEup$V1, lum084_myo084_DEdn$V1)
lum084_stem084_DEup <- read.delim("UP.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DEdn <- read.delim("DN.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DE <- c(lum084_stem084_DEup$V1, lum084_stem084_DEdn$V1)
myo084_stem084_DEup <- read.delim("UP.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DEdn <- read.delim("DN.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DE <- c(myo084_stem084_DEup$V1, myo084_stem084_DEdn$V1)

# DE exons FDR = 0.015
setwd("~/REMC/breast/exon/")
lum084_myo084_exonup <- read.delim("UP.myo-RM084_lum-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_exondn <- read.delim("DN.myo-RM084_lum-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_exon <- rbind(cbind(lum084_myo084_exonup, DE = "up"), cbind(lum084_myo084_exondn, DE = "dn"))
lum084_stem084_exonup <- read.delim("UP.lum-RM084_stem-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_exondn <- read.delim("DN.lum-RM084_stem-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_exon <- rbind(cbind(lum084_stem084_exonup, DE = "up"), cbind(lum084_stem084_exondn, DE = "dn"))
myo084_stem084_exonup <- read.delim("UP.myo-RM084_stem-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_exondn <- read.delim("DN.myo-RM084_stem-RM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_exon <- rbind(cbind(myo084_stem084_exonup, DE = "up"), cbind(myo084_stem084_exondn, DE = "dn"))

# DE junctions FDR = 0.03
setwd("~/REMC/junction/DEfine/")
lum084_myo084_junctionup <- read.delim("UP.lum-RM084_myo-RM084.FDR_0.03.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_junctiondn <- read.delim("DN.lum-RM084_myo-RM084.FDR_0.03.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_junction <- rbind(cbind(lum084_myo084_junctionup, DE = "up"), cbind(lum084_myo084_junctiondn, DE = "dn"))
lum084_stem084_junctionup <- read.delim("UP.lum-RM084_stem-RM084.FDR_0.03.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_junctiondn <- read.delim("DN.lum-RM084_stem-RM084.FDR_0.03.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_junction <- rbind(cbind(lum084_stem084_junctionup, DE = "up"), cbind(lum084_stem084_junctiondn, DE = "dn"))
myo084_stem084_junctionup <- read.delim("UP.myo-RM084_stem-RM084.FDR_0.03.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_junctiondn <- read.delim("DN.myo-RM084_stem-RM084.FDR_0.03.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_junction <- rbind(cbind(myo084_stem084_junctionup, DE = "up"), cbind(myo084_stem084_junctiondn, DE = "dn"))

# exon DE & at least one junction in the exon DE & gene NOT DE
setwd("~/REMC/junction/")
support <- function(exonID, DE, junctionDE){
  junctions <- junction_exon[junction_exon$exonID == exonID,]$junctionID
  support <- as.character(junctionDE[junctions, ]$DE)
  return(sum(support == DE)) # No. of junction supports exon DE
}
Vsupport <- Vectorize(support, vectorize.args = c("exonID", "DE"))

# lum084 vs myo084
lum084_myo084_exon$exonID <- gsub("<[-1_A-Z0-9]+", "", lum084_myo084_exon$V1)
lum084_myo084_exon$ENSG <- unlist(strsplit(lum084_myo084_exon$V1, "_"))[seq(2, 2*nrow(lum084_myo084_exon), by = 2)]
junction_combined[1, "DEfine.exons"] <- nrow(lum084_myo084_exon)
head(lum084_myo084_exon)
# exclude DE genes
lum084_myo084_exon <- lum084_myo084_exon[!(lum084_myo084_exon$ENSG %in% lum084_myo084_DE),]
junction_combined[1, "exclude.DE.genes"] <- nrow(lum084_myo084_exon)
# with annotated junctions 
lum084_myo084_exon <- lum084_myo084_exon[lum084_myo084_exon$exonID %in% junction_exon$exonID, ]
junction_combined[1, "with.annotated.junction"] <- nrow(lum084_myo084_exon)
rownames(lum084_myo084_junction) <- lum084_myo084_junction$V1
lum084_myo084_exon <- lum084_myo084_exon[order(lum084_myo084_exon$V5), ]
# with junction DE
lum084_myo084_exon$support <- Vsupport(lum084_myo084_exon$exonID, lum084_myo084_exon$DE, lum084_myo084_junction)
lum084_myo084_exon <- na.omit(lum084_myo084_exon)
junction_combined[1, "with.DEfine.junctions"] <- nrow(lum084_myo084_exon)
# with junction DE in the same direction
lum084_myo084_exon <- lum084_myo084_exon[lum084_myo084_exon$support > 0, ]
junction_combined[1, "exon.junction.correlate"] <- nrow(lum084_myo084_exon)

# lum084 vs stem084
lum084_stem084_exon$exonID <- gsub("<[-1_A-Z0-9]+", "", lum084_stem084_exon$V1)
lum084_stem084_exon$ENSG <- unlist(strsplit(lum084_stem084_exon$V1, "_"))[seq(2, 2*nrow(lum084_stem084_exon), by = 2)]
junction_combined[2, "DEfine.exons"] <- nrow(lum084_stem084_exon)
head(lum084_stem084_exon)
# exclude DE genes
lum084_stem084_exon <- lum084_stem084_exon[!(lum084_stem084_exon$ENSG %in% lum084_stem084_DE),]
junction_combined[2, "exclude.DE.genes"] <- nrow(lum084_stem084_exon)
# with annotated junctions 
lum084_stem084_exon <- lum084_stem084_exon[lum084_stem084_exon$exonID %in% junction_exon$exonID, ]
junction_combined[2, "with.annotated.junction"] <- nrow(lum084_stem084_exon)
rownames(lum084_stem084_junction) <- lum084_stem084_junction$V1
lum084_stem084_exon <- lum084_stem084_exon[order(lum084_stem084_exon$V5), ]
# with junction DE
lum084_stem084_exon$support <- Vsupport(lum084_stem084_exon$exonID, lum084_stem084_exon$DE, lum084_stem084_junction)
lum084_stem084_exon <- na.omit(lum084_stem084_exon)
junction_combined[2, "with.DEfine.junctions"] <- nrow(lum084_stem084_exon)
# with junction DE in the same direction
lum084_stem084_exon <- lum084_stem084_exon[lum084_stem084_exon$support > 0, ]
junction_combined[2, "exon.junction.correlate"] <- nrow(lum084_stem084_exon)

# myo084 vs stem084
myo084_stem084_exon$exonID <- gsub("<[-1_A-Z0-9]+", "", myo084_stem084_exon$V1)
myo084_stem084_exon$ENSG <- unlist(strsplit(myo084_stem084_exon$V1, "_"))[seq(2, 2*nrow(myo084_stem084_exon), by = 2)]
junction_combined[3, "DEfine.exons"] <- nrow(myo084_stem084_exon)
head(myo084_stem084_exon)
# exclude DE genes
myo084_stem084_exon <- myo084_stem084_exon[!(myo084_stem084_exon$ENSG %in% myo084_stem084_DE),]
junction_combined[3, "exclude.DE.genes"] <- nrow(myo084_stem084_exon)
# with annotated junctions 
myo084_stem084_exon <- myo084_stem084_exon[myo084_stem084_exon$exonID %in% junction_exon$exonID, ]
junction_combined[3, "with.annotated.junction"] <- nrow(myo084_stem084_exon)
rownames(myo084_stem084_junction) <- myo084_stem084_junction$V1
myo084_stem084_exon <- myo084_stem084_exon[order(myo084_stem084_exon$V5), ]
# with junction DE
myo084_stem084_exon$support <- Vsupport(myo084_stem084_exon$exonID, myo084_stem084_exon$DE, myo084_stem084_junction)
myo084_stem084_exon <- na.omit(myo084_stem084_exon)
junction_combined[3, "with.DEfine.junctions"] <- nrow(myo084_stem084_exon)
# with junction DE in the same direction
myo084_stem084_exon <- myo084_stem084_exon[myo084_stem084_exon$support > 0, ]
junction_combined[3, "exon.junction.correlate"] <- nrow(myo084_stem084_exon)

save(junction_fold, junction_DEfine, junction_relative, junction_combined, junction_valid, file = "junction_summary_tables.Rdata")


# junction DE & at least one exon with the junction DE & gene NOT DE
setwd("~/REMC/junction/")
load("junction_RM084.Rdata")
support2 <- function(junctionID, DE, exonDE){
  exons <- junction_exon[junction_exon$junctionID == junctionID,]$exonID
  support <- as.character(exonDE[exons, ]$DE)
  return(sum(support == DE)) # No. of junction supports exon DE
}
Vsupport2 <- Vectorize(support2, vectorize.args = c("junctionID", "DE"))

# lum084 vs myo084
lum084_myo084_junction$ENSG <- junction[lum084_myo084_junction$V1,]$gene
nrow(lum084_myo084_junction)
head(lum084_myo084_junction)
# exclude DE genes
lum084_myo084_junction <- lum084_myo084_junction[!(lum084_myo084_junction$ENSG %in% lum084_myo084_DE),]
nrow(lum084_myo084_junction)
# with annotated exon 
lum084_myo084_junction <- lum084_myo084_junction[lum084_myo084_junction$V1 %in% junction_exon$junctionID, ]
nrow(lum084_myo084_junction)
rownames(lum084_myo084_junction) <- lum084_myo084_junction$V1
lum084_myo084_junction <- lum084_myo084_junction[order(lum084_myo084_junction$V5), ]
# with exon DE
lum084_myo084_junction$support <- Vsupport2(lum084_myo084_junction$V1, lum084_myo084_junction$DE, lum084_myo084_exon)
lum084_myo084_junction <- na.omit(lum084_myo084_junction)
nrow(lum084_myo084_junction)
# with exon DE in the same direction
lum084_myo084_junction <- lum084_myo084_junction[lum084_myo084_junction$support > 0, ]
nrow(lum084_myo084_junction)

# lum084 vs stem084
lum084_stem084_junction$ENSG <- junction[lum084_stem084_junction$V1,]$gene
nrow(lum084_stem084_junction)
head(lum084_stem084_junction)
# exclude DE genes
lum084_stem084_junction <- lum084_stem084_junction[!(lum084_stem084_junction$ENSG %in% lum084_stem084_DE),]
nrow(lum084_stem084_junction)
# with annotated exon 
lum084_stem084_junction <- lum084_stem084_junction[lum084_stem084_junction$V1 %in% junction_exon$junctionID, ]
nrow(lum084_stem084_junction)
rownames(lum084_stem084_junction) <- lum084_stem084_junction$V1
lum084_stem084_junction <- lum084_stem084_junction[order(lum084_stem084_junction$V5), ]
# with exon DE
lum084_stem084_junction$support <- Vsupport2(lum084_stem084_junction$V1, lum084_stem084_junction$DE, lum084_stem084_exon)
lum084_stem084_junction <- na.omit(lum084_stem084_junction)
nrow(lum084_stem084_junction)
# with exon DE in the same direction
lum084_stem084_junction <- lum084_stem084_junction[lum084_stem084_junction$support > 0, ]
nrow(lum084_stem084_junction)

# myo084 vs stem084
myo084_stem084_junction$ENSG <- junction[myo084_stem084_junction$V1,]$gene
nrow(myo084_stem084_junction)
head(myo084_stem084_junction)
# exclude DE genes
myo084_stem084_junction <- myo084_stem084_junction[!(myo084_stem084_junction$ENSG %in% myo084_stem084_DE),]
nrow(myo084_stem084_junction)
# with annotated exon 
myo084_stem084_junction <- myo084_stem084_junction[myo084_stem084_junction$V1 %in% junction_exon$junctionID, ]
nrow(myo084_stem084_junction)
rownames(myo084_stem084_junction) <- myo084_stem084_junction$V1
myo084_stem084_junction <- myo084_stem084_junction[order(myo084_stem084_junction$V5), ]
# with exon DE
myo084_stem084_junction$support <- Vsupport2(myo084_stem084_junction$V1, myo084_stem084_junction$DE, myo084_stem084_exon)
myo084_stem084_junction <- na.omit(myo084_stem084_junction)
nrow(myo084_stem084_junction)
# with exon DE in the same direction
myo084_stem084_junction <- myo084_stem084_junction[myo084_stem084_junction$support > 0, ]
nrow(myo084_stem084_junction)

################################################################################################################################################
# validate previously identified AS genes with junction RPKM 
junction_valid <- matrix(NA, nrow = 3, ncol = 3, dimnames = list(c("lum084_myo084", "lum084_stem084", "myo084_stem084"), c("No.isoforms", "with.junction.cov", "with.junction.support")))

# previously identified AS genes & junction RPKM 
setwd("~/REMC/junction/")
load("~/hg19/hg19v65_genes.Rdata")
load("junction_RM084.Rdata")
lum084_myo084_gene_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84_isoform_only.csv", head = T, as.is = T, col.names = c("gene"))
lum084_stem084_gene_valid <- read.csv("~/REMC/breast/tissue/lum_stem_84_isoform_only.csv", head = T, as.is = T, col.names = c("gene"))
myo084_stem084_gene_valid <- read.csv("~/REMC/breast/tissue/myo_stem_84_isoform_only.csv", head = T, as.is = T, col.names = c("gene"))

# AS gene with junction support: highest fold change between junction RPKM > 2 
high_fold <- function(gene, sample1, sample2, junctions = junction){
  e <- 1e-6
  covcut <- 3 # cutoff for sum junction coverage of two samples
  junctions <- junctions[as.character(junctions$gene) == gene, ]
  junctions <- junctions[(junctions[, paste0(sample1, "N")] + junctions[, paste0(sample2, "N")]) > covcut,]
  if(nrow(junctions) == 0){
    return(c(NA, NA, NA, NA, NA)) # not enough coverage for junctionss of this gene
  }else{
    junctions$fold <- abs(log2((junctions[, paste0(sample1, "rpkm")] + e) / (junctions[, paste0(sample2, "rpkm")] + e)))
    highest <- junctions[which.max(junctions$fold),]
    return(c(highest[, paste0(sample1, "rpkm")], highest[, paste0(sample2, "rpkm")], highest[, paste0(sample1, "gene")], highest[, paste0(sample2, "gene")], highest[, "fold"]))
  }
}
Vhigh_fold <- Vectorize(high_fold, vectorize.args = "gene", SIMPLIFY = T)

logfoldcut <- 1 # cutoff for abs(log2(fold change junction rpkm))

# lum084 vs myo084
nrow(lum084_myo084_gene_valid)
lum084_myo084_gene_valid <- data.frame(gene = lum084_myo084_gene_valid$gene, t(Vhigh_fold(lum084_myo084_gene_valid$gene, "lum084", "myo084")))
colnames(lum084_myo084_gene_valid) <- c("gene", "lum084rpkm", "myo084rpkm", "lum084gene", "myo084gene", "logfold_abs")
junction_valid[1, 1] <- nrow(lum084_myo084_gene_valid)
# No. of genes with enough junction coverage 
lum084_myo084_gene_valid <- na.omit(lum084_myo084_gene_valid) 
junction_valid[1, 2] <- nrow(lum084_myo084_gene_valid)
# No. of genes with highest fold change > 2 i.e. abs(log2(fold change)) > 1 
lum084_myo084_gene_valid <- lum084_myo084_gene_valid[lum084_myo084_gene_valid$logfold_abs > logfoldcut,]
lum084_myo084_gene_valid <- lum084_myo084_gene_valid[order(lum084_myo084_gene_valid$logfold_abs, decreasing = T), ]
lum084_myo084_gene_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum084_myo084_gene_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
junction_valid[1, 3] <- nrow(lum084_myo084_gene_valid)
write.table(lum084_myo084_gene_valid, file = "lum084_myo084_gene_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# lum084 vs stem084
nrow(lum084_stem084_gene_valid)
lum084_stem084_gene_valid <- data.frame(gene = lum084_stem084_gene_valid$gene, t(Vhigh_fold(lum084_stem084_gene_valid$gene, "lum084", "stem084")))
colnames(lum084_stem084_gene_valid) <- c("gene", "lum084rpkm", "stem084rpkm", "lum084gene", "stem084gene", "logfold_abs")
junction_valid[2, 1] <- nrow(lum084_stem084_gene_valid)
# No. of genes with enough junction coverage 
lum084_stem084_gene_valid <- na.omit(lum084_stem084_gene_valid) 
junction_valid[2, 2] <- nrow(lum084_stem084_gene_valid)
# No. of genes with highest fold change > 2 i.e. abs(log2(fold change)) > 1 
lum084_stem084_gene_valid <- lum084_stem084_gene_valid[lum084_stem084_gene_valid$logfold_abs > logfoldcut,]
lum084_stem084_gene_valid <- lum084_stem084_gene_valid[order(lum084_stem084_gene_valid$logfold_abs, decreasing = T), ]
lum084_stem084_gene_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum084_stem084_gene_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
junction_valid[2, 3] <- nrow(lum084_stem084_gene_valid)
write.table(lum084_stem084_gene_valid, file = "lum084_stem084_gene_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# myo084 vs stem084
nrow(myo084_stem084_gene_valid)
myo084_stem084_gene_valid <- data.frame(gene = myo084_stem084_gene_valid$gene, t(Vhigh_fold(myo084_stem084_gene_valid$gene, "myo084", "stem084")))
colnames(myo084_stem084_gene_valid) <- c("gene", "myo084rpkm", "stem084rpkm", "myo084gene", "stem084gene", "logfold_abs")
junction_valid[3, 1] <- nrow(myo084_stem084_gene_valid)
# No. of genes with enough junction coverage 
myo084_stem084_gene_valid <- na.omit(myo084_stem084_gene_valid) 
junction_valid[3, 2] <- nrow(myo084_stem084_gene_valid)
# No. of genes with highest fold change > 2 i.e. abs(log2(fold change)) > 1 
myo084_stem084_gene_valid <- myo084_stem084_gene_valid[myo084_stem084_gene_valid$logfold_abs > logfoldcut,]
myo084_stem084_gene_valid <- myo084_stem084_gene_valid[order(myo084_stem084_gene_valid$logfold_abs, decreasing = T), ]
myo084_stem084_gene_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(myo084_stem084_gene_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
junction_valid[3, 3] <- nrow(myo084_stem084_gene_valid)
write.table(myo084_stem084_gene_valid, file = "myo084_stem084_gene_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

junction_valid_common <- merge(lum084_myo084_gene_valid, merge(lum084_stem084_gene_valid, myo084_stem084_gene_valid, by = c("gene", "coord", "chr", "start", "end", "strand", "type", "name", "description")), by = c("gene", "coord", "chr", "start", "end", "strand", "type", "name", "description"))
junction_valid_common$lum084gene <- junction_valid_common$lum084gene.x
junction_valid_common$myo084gene <- junction_valid_common$myo084gene.x
junction_valid_common$stem084gene <- junction_valid_common$stem084gene.x
junction_valid_common$lum_myo <- junction_valid_common$logfold_abs
junction_valid_common$lum_stem <- junction_valid_common$logfold_abs.x
junction_valid_common$myo_stem <- junction_valid_common$logfold_abs.y
junction_valid_common <- junction_valid_common[names(junction_valid_common) %in% c("gene", "coord", "chr", "start", "end", "strand", "type", "name", "description", "lum084gene", "myo084gene", "stem084gene", "lum_myo", "lum_stem", "myo_stem")]
junction_valid_common <- junction_valid_common[order(junction_valid_common$lum_myo + junction_valid_common$lum_stem + junction_valid_common$myo_stem, decreasing = T),]
rownames(junction_valid_common) <- junction_valid_common$gene
write.table(junction_valid_common, file = "junction_valid_common.txt", sep = "\t", quote = F, row.names = F, col.names = T)
save(junction_valid, junction_valid_common, lum084_myo084_gene_valid, lum084_stem084_gene_valid, myo084_stem084_gene_valid, file = "lum_myo_stem_084_validate.Rdata")

save(junction_fold, junction_DEfine, junction_relative, junction_combined, junction_valid, file = "junction_summary_tables.Rdata")


################################################################################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# validate previously identified isoform exons with junction RPKM 
junction_valid <- matrix(NA, nrow = 3, ncol = 6, dimnames = list(c("lum084_myo084", "lum084_stem084", "myo084_stem084"), c("No.isoform.exons", "No.isoform.genes", "No.exons.with.junction.cov", "No.genes.with.junction.cov", "No.exons.with.junction.support", "No.genes.with.junction.support")))

# DE genes FDR = 0.015
setwd("~/REMC/breast/gene/")
lum084_myo084_DEup <- read.delim("UP.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DEdn <- read.delim("DN.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_myo084_DE <- c(lum084_myo084_DEup$V1, lum084_myo084_DEdn$V1)
lum084_stem084_DEup <- read.delim("UP.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DEdn <- read.delim("DN.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum084_stem084_DE <- c(lum084_stem084_DEup$V1, lum084_stem084_DEdn$V1)
myo084_stem084_DEup <- read.delim("UP.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DEdn <- read.delim("DN.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo084_stem084_DE <- c(myo084_stem084_DEup$V1, myo084_stem084_DEdn$V1)
rm(lum084_myo084_DEup, lum084_myo084_DEdn, lum084_stem084_DEup, lum084_stem084_DEdn, myo084_stem084_DEup, myo084_stem084_DEdn)

# previously identified isoform exons & junction RPKM 
setwd("~/REMC/junction/")
load("~/hg19/hg19v65_genes.Rdata")
load("junction_RM084.Rdata")
junction_exon <- read.delim("~/hg19/hg19v65_junctions_vs_exons_for_genes.relations.unique", head = T, as.is = T)

# isoform exons with junction support: junction RPKM in one sample < 0.1(cov <= 1), > 0.1 in the other(cov >= 2) 
high_fold <- function(exon, up, sample1, sample2, junctions = junction, junction_exons = junction_exon){
  e <- 1e-6
  covcut <- 1 # cutoff for sum junction coverage of two samples
  junctions <- junctions[junction_exons[junction_exons$exonID == exon, "junctionID"], ]
  junctions <- junctions[(junctions[, paste0(sample1, "N")] + junctions[, paste0(sample2, "N")]) > covcut,]
  if(nrow(junctions) == 0){
    return(c(NA, NA, NA)) # not enough coverage for junctionss of this gene
  }else{
    junctions$fold <- log2((junctions[, paste0(sample1, "rpkm")] + e) / (junctions[, paste0(sample2, "rpkm")] + e))
    if(up > 0){
      highest <- junctions[which.max(junctions$fold),]
    }else{
      highest <- junctions[which.min(junctions$fold),]
    }
    return(c(highest[, paste0(sample1, "rpkm")], highest[, paste0(sample2, "rpkm")], highest[, "fold"]))
  }
}
Vhigh_fold <- Vectorize(high_fold, vectorize.args = c("exon", "up"), SIMPLIFY = T)

logfoldcut <- 1 # cutoff for abs(log2(fold change junction rpkm))
jcut <- 0.1 # cutoff for junction RPKM, i.e.one sample > 0.1, the other < 0.1

# lum084 vs myo084
lum084_myo084up_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84up_isoform.csv", head = T, as.is = T)
lum084_myo084dn_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84dn_isoform.csv", head = T, as.is = T)
lum084_myo084_isoform_all <- rbind(lum084_myo084up_isoform_valid, lum084_myo084dn_isoform_valid)
lum084_myo084_isoform_all <- lum084_myo084_isoform_all[!is.element(lum084_myo084_isoform_all$id, lum084_myo084_DE),]
write.csv(lum084_myo084_isoform_all, file = "~/REMC/breast/tissue/myo_lum_84up_isoform.csv", row.names = F, quote = F)
rm(lum084_myo084up_isoform_valid, lum084_myo084dn_isoform_valid)
lum084_myo084_isoform_valid <- lum084_myo084_isoform_all
(junction_valid["lum084_myo084", "No.isoform.exons"] <- nrow(lum084_myo084_isoform_valid))
(junction_valid["lum084_myo084", "No.isoform.genes"] <- length(unique(lum084_myo084_isoform_valid$id)))
sum(is.element(lum084_myo084_isoform_valid$V1, junction_exon$exonID))
lum084_myo084_isoform_valid <- data.frame(exon = lum084_myo084_isoform_valid$V1, gene = lum084_myo084_isoform_valid$id, DEfine.p = lum084_myo084_isoform_valid$V5, 
                                          lum084gene = lum084_myo084_isoform_valid$ave3, myo084gene = lum084_myo084_isoform_valid$ave2, 
                                          lum084exon = lum084_myo084_isoform_valid$V3, myo084exon = lum084_myo084_isoform_valid$V2, 
                                          t(Vhigh_fold(exon = lum084_myo084_isoform_valid$V1, up = lum084_myo084_isoform_valid$V3 - lum084_myo084_isoform_valid$V2, sample1 = "lum084", sample2 = "myo084")))
colnames(lum084_myo084_isoform_valid)[(ncol(lum084_myo084_isoform_valid)-2):ncol(lum084_myo084_isoform_valid)] <- c("lum084junction", "myo084junction", "junction_logfold")
(junction_valid["lum084_myo084", "No.isoform.exons"] <- nrow(lum084_myo084_isoform_valid))
(junction_valid["lum084_myo084", "No.isoform.genes"] <- length(unique(lum084_myo084_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
lum084_myo084_isoform_valid <- na.omit(lum084_myo084_isoform_valid)
(junction_valid["lum084_myo084", "No.exons.with.junction.cov"] <- nrow(lum084_myo084_isoform_valid))
(junction_valid["lum084_myo084", "No.genes.with.junction.cov"] <- length(unique(lum084_myo084_isoform_valid$gene)))
# No. of exons/genes with junction support
#lum084_myo084_isoform_valid <- lum084_myo084_isoform_valid[abs(lum084_myo084_isoform_valid$junction_logfold) > logfoldcut,]
lum084_myo084_isoform_valid <- lum084_myo084_isoform_valid[(lum084_myo084_isoform_valid$lum084junction >= jcut & lum084_myo084_isoform_valid$myo084junction < jcut) | (lum084_myo084_isoform_valid$lum084junction < jcut & lum084_myo084_isoform_valid$myo084junction >= jcut),]
lum084_myo084_isoform_valid$junction_logfold <- NULL
lum084_myo084_isoform_valid <- lum084_myo084_isoform_valid[order((lum084_myo084_isoform_valid$lum084gene + lum084_myo084_isoform_valid$myo084gene), decreasing = T), ]
lum084_myo084_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum084_myo084_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["lum084_myo084", "No.exons.with.junction.support"] <- nrow(lum084_myo084_isoform_valid))
(junction_valid["lum084_myo084", "No.genes.with.junction.support"] <- length(unique(lum084_myo084_isoform_valid$gene)))
write.table(lum084_myo084_isoform_valid, file = "lum084_myo084_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# lum084 vs stem084
lum084_stem084up_isoform_valid <- read.csv("~/REMC/breast/tissue/lum_stem_84up_isoform.csv", head = T, as.is = T)
lum084_stem084dn_isoform_valid <- read.csv("~/REMC/breast/tissue/lum_stem_84dn_isoform.csv", head = T, as.is = T)
lum084_stem084_isoform_all <- rbind(lum084_stem084up_isoform_valid, lum084_stem084dn_isoform_valid)
lum084_stem084_isoform_all <- lum084_stem084_isoform_all[!is.element(lum084_stem084_isoform_all$id, lum084_stem084_DE),]
write.csv(lum084_stem084_isoform_all, file = "~/REMC/breast/tissue/stem_lum_84up_isoform.csv", row.names = F, quote = F)
rm(lum084_stem084up_isoform_valid, lum084_stem084dn_isoform_valid)
lum084_stem084_isoform_valid <- lum084_stem084_isoform_all
(junction_valid["lum084_stem084", "No.isoform.exons"] <- nrow(lum084_stem084_isoform_valid))
(junction_valid["lum084_stem084", "No.isoform.genes"] <- length(unique(lum084_stem084_isoform_valid$id)))
sum(is.element(lum084_stem084_isoform_valid$V1, junction_exon$exonID))
lum084_stem084_isoform_valid <- data.frame(exon = lum084_stem084_isoform_valid$V1, gene = lum084_stem084_isoform_valid$id, DEfine.p = lum084_stem084_isoform_valid$V5, 
                                           lum084gene = lum084_stem084_isoform_valid$ave2, stem084gene = lum084_stem084_isoform_valid$ave3, 
                                           lum084exon = lum084_stem084_isoform_valid$V2, stem084exon = lum084_stem084_isoform_valid$V3, 
                                           t(Vhigh_fold(exon = lum084_stem084_isoform_valid$V1, up = lum084_stem084_isoform_valid$V2 - lum084_stem084_isoform_valid$V3, sample1 = "lum084", sample2 = "stem084")))
colnames(lum084_stem084_isoform_valid)[(ncol(lum084_stem084_isoform_valid)-2):ncol(lum084_stem084_isoform_valid)] <- c("lum084junction", "stem084junction", "junction_logfold")
(junction_valid["lum084_stem084", "No.isoform.exons"] <- nrow(lum084_stem084_isoform_valid))
(junction_valid["lum084_stem084", "No.isoform.genes"] <- length(unique(lum084_stem084_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
lum084_stem084_isoform_valid <- na.omit(lum084_stem084_isoform_valid)
(junction_valid["lum084_stem084", "No.exons.with.junction.cov"] <- nrow(lum084_stem084_isoform_valid))
(junction_valid["lum084_stem084", "No.genes.with.junction.cov"] <- length(unique(lum084_stem084_isoform_valid$gene)))
# No. of exons/genes with junction support
#lum084_stem084_isoform_valid <- lum084_stem084_isoform_valid[abs(lum084_stem084_isoform_valid$junction_logfold) > logfoldcut,]
lum084_stem084_isoform_valid <- lum084_stem084_isoform_valid[(lum084_stem084_isoform_valid$lum084junction >= jcut & lum084_stem084_isoform_valid$stem084junction < jcut) | (lum084_stem084_isoform_valid$lum084junction < jcut & lum084_stem084_isoform_valid$stem084junction >= jcut),]
lum084_stem084_isoform_valid$junction_logfold <- NULL
lum084_stem084_isoform_valid <- lum084_stem084_isoform_valid[order((lum084_stem084_isoform_valid$lum084gene + lum084_stem084_isoform_valid$stem084gene), decreasing = T), ]
lum084_stem084_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum084_stem084_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["lum084_stem084", "No.exons.with.junction.support"] <- nrow(lum084_stem084_isoform_valid))
(junction_valid["lum084_stem084", "No.genes.with.junction.support"] <- length(unique(lum084_stem084_isoform_valid$gene)))
write.table(lum084_stem084_isoform_valid, file = "lum084_stem084_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# myo084 vs stem084
myo084_stem084up_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_stem_84up_isoform.csv", head = T, as.is = T)
myo084_stem084dn_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_stem_84dn_isoform.csv", head = T, as.is = T)
myo084_stem084_isoform_all <- rbind(myo084_stem084up_isoform_valid, myo084_stem084dn_isoform_valid)
myo084_stem084_isoform_all <- myo084_stem084_isoform_all[!is.element(myo084_stem084_isoform_all$id, myo084_stem084_DE),]
write.csv(myo084_stem084_isoform_all, file = "~/REMC/breast/tissue/stem_myo_84up_isoform.csv", row.names = F, quote = F)
rm(myo084_stem084up_isoform_valid, myo084_stem084dn_isoform_valid)
myo084_stem084_isoform_valid <- myo084_stem084_isoform_all
(junction_valid["myo084_stem084", "No.isoform.exons"] <- nrow(myo084_stem084_isoform_valid))
(junction_valid["myo084_stem084", "No.isoform.genes"] <- length(unique(myo084_stem084_isoform_valid$id)))
sum(is.element(myo084_stem084_isoform_valid$V1, junction_exon$exonID))
myo084_stem084_isoform_valid <- data.frame(exon = myo084_stem084_isoform_valid$V1, gene = myo084_stem084_isoform_valid$id, DEfine.p = myo084_stem084_isoform_valid$V5, 
                                           myo084gene = myo084_stem084_isoform_valid$ave2, stem084gene = myo084_stem084_isoform_valid$ave3, 
                                           myo084exon = myo084_stem084_isoform_valid$V2, stem084exon = myo084_stem084_isoform_valid$V3, 
                                           t(Vhigh_fold(exon = myo084_stem084_isoform_valid$V1, up = myo084_stem084_isoform_valid$V2 - myo084_stem084_isoform_valid$V3, sample1 = "myo084", sample2 = "stem084")))
colnames(myo084_stem084_isoform_valid)[(ncol(myo084_stem084_isoform_valid)-2):ncol(myo084_stem084_isoform_valid)] <- c("myo084junction", "stem084junction", "junction_logfold")
(junction_valid["myo084_stem084", "No.isoform.exons"] <- nrow(myo084_stem084_isoform_valid))
(junction_valid["myo084_stem084", "No.isoform.genes"] <- length(unique(myo084_stem084_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
myo084_stem084_isoform_valid <- na.omit(myo084_stem084_isoform_valid)
(junction_valid["myo084_stem084", "No.exons.with.junction.cov"] <- nrow(myo084_stem084_isoform_valid))
(junction_valid["myo084_stem084", "No.genes.with.junction.cov"] <- length(unique(myo084_stem084_isoform_valid$gene)))
# No. of exons/genes with junction support
#myo084_stem084_isoform_valid <- myo084_stem084_isoform_valid[abs(myo084_stem084_isoform_valid$junction_logfold) > logfoldcut,]
myo084_stem084_isoform_valid <- myo084_stem084_isoform_valid[(myo084_stem084_isoform_valid$myo084junction >= jcut & myo084_stem084_isoform_valid$stem084junction < jcut) | (myo084_stem084_isoform_valid$myo084junction < jcut & myo084_stem084_isoform_valid$stem084junction >= jcut),]
myo084_stem084_isoform_valid$junction_logfold <- NULL
myo084_stem084_isoform_valid <- myo084_stem084_isoform_valid[order((myo084_stem084_isoform_valid$myo084gene + myo084_stem084_isoform_valid$stem084gene), decreasing = T), ]
myo084_stem084_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(myo084_stem084_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["myo084_stem084", "No.exons.with.junction.support"] <- nrow(myo084_stem084_isoform_valid))
(junction_valid["myo084_stem084", "No.genes.with.junction.support"] <- length(unique(myo084_stem084_isoform_valid$gene)))
write.table(myo084_stem084_isoform_valid, file = "myo084_stem084_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rm(junction_exon)

# common isoforms across cell types
junction_valid_common <- merge(lum084_myo084_isoform_valid_gene, merge(lum084_stem084_isoform_valid_gene, myo084_stem084_isoform_valid_gene, by = c("gene", "coord", "chr", "start", "end", "strand", "type", "name", "description")), by = c("gene", "coord", "chr", "start", "end", "strand", "type", "name", "description"))
junction_valid_common$lum084gene <- junction_valid_common$lum084gene.x
junction_valid_common$myo084gene <- junction_valid_common$myo084gene.x
junction_valid_common$stem084gene <- junction_valid_common$stem084gene.x
junction_valid_common <- junction_valid_common[names(junction_valid_common) %in% c("gene", "coord", "chr", "start", "end", "strand", "type", "name", "description", "lum084gene", "myo084gene", "stem084gene")]
junction_valid_common <- junction_valid_common[order((junction_valid_common$lum084gene + junction_valid_common$myo084gene + junction_valid_common$stem084gene), decreasing = T),]
rownames(junction_valid_common) <- junction_valid_common$gene
write.table(junction_valid_common, file = "junction_valid_common.txt", sep = "\t", quote = F, row.names = F, col.names = T)
save.image(file = "junction_valid.Rdata")

setwd("~/快盘/REMC/junction/")
load("junction_valid.Rdata")
# No. of exons for DE genes / isoform genes: use collapsed exons instead of Ensembl exons
exon <- read.delim("~/hg19/hg19v65_exons_collapsed.txt", head = F, as.is = T)
Nexon <- data.frame(id = levels(as.factor(exon$V2)), Nexon = sapply(levels(as.factor(exon$V2)), function(x) sum(exon$V2 == x)))
lum084_myo084_DE <- data.frame(id = lum084_myo084_DE, Nexon = Nexon[lum084_myo084_DE, ]$Nexon)
lum084_myo084_isoform_all_gene <- lum084_myo084_isoform_all[!duplicated(lum084_myo084_isoform_all$id), ]
lum084_myo084_isoform_all_gene$Nexon <- Nexon[lum084_myo084_isoform_all_gene$id, "Nexon"]
lum084_myo084_isoform_valid_gene$Nexon <- Nexon[as.character(lum084_myo084_isoform_valid_gene$gene), "Nexon"] 
lum084_stem084_DE <- data.frame(id = lum084_stem084_DE, Nexon = Nexon[lum084_stem084_DE, ]$Nexon)
lum084_stem084_isoform_all_gene <- lum084_stem084_isoform_all[!duplicated(lum084_stem084_isoform_all$id), ]
lum084_stem084_isoform_all_gene$Nexon <- Nexon[lum084_stem084_isoform_all_gene$id, "Nexon"]
lum084_stem084_isoform_valid_gene$Nexon <- Nexon[as.character(lum084_stem084_isoform_valid_gene$gene), "Nexon"] 
myo084_stem084_DE <- data.frame(id = myo084_stem084_DE, Nexon = Nexon[myo084_stem084_DE, ]$Nexon)
myo084_stem084_isoform_all_gene <- myo084_stem084_isoform_all[!duplicated(myo084_stem084_isoform_all$id), ]
myo084_stem084_isoform_all_gene$Nexon <- Nexon[myo084_stem084_isoform_all_gene$id, "Nexon"]
myo084_stem084_isoform_valid_gene$Nexon <- Nexon[as.character(myo084_stem084_isoform_valid_gene$gene), "Nexon"] 
rm(exon)
pdf("Nexon_DEvsIsoforms.pdf")
plot(c(0, 50), c(0, 0.1), type = "n", main = "No. of exons in DE genes and isoforms", xlab = "No. of exons", ylab = "density")
lines(density(na.omit(Nexon$Nexon)), col = 1, lty = 1, lwd = 5)
lines(density(na.omit(lum084_myo084_DE$Nexon)), col = 2, lty = 1, lwd = 3)
lines(density(na.omit(lum084_stem084_DE$Nexon)), col = 3, lty = 1, lwd = 3)
lines(density(na.omit(myo084_stem084_DE$Nexon)), col = 4, lty = 1, lwd = 3)
lines(density(na.omit(lum084_myo084_isoform_all_gene$Nexon)), col = 2, lty = 3, lwd = 3)
lines(density(na.omit(lum084_stem084_isoform_all_gene$Nexon)), col = 3, lty = 3, lwd = 3)
lines(density(na.omit(myo084_stem084_isoform_all_gene$Nexon)), col = 4, lty = 3, lwd = 3)
# lines(density(na.omit(lum084_myo084_isoform_valid_gene$Nexon)), col = 2, lty = 4, lwd = 3)
# lines(density(na.omit(lum084_stem084_isoform_valid_gene$Nexon)), col = 3, lty = 4, lwd = 3)
# lines(density(na.omit(myo084_stem084_isoform_valid_gene$Nexon)), col = 4, lty = 4, lwd = 3)
# legend("bottomright", c("DE genes", "all isoforms", "validated isoforms", "all genes"), col = 1, lty = c(1, 3, 4, 1), lwd = 5, cex = 0.8)
legend("bottomright", c("DE genes", "isoforms", "all genes"), col = 1, lty = c(1, 3, 1), lwd = 5, cex = 0.8)
legend("topright", c("all genes", "lum084 vs myo084", "lum084 vs stem084", "myo084 vs stem084"), col = 1:4, lty = 1, lwd = 5, cex = 0.8)
dev.off()

# Isoform genes in DE genes with small No. of exons
# all DE exons from DEfine
lum084_myo084_exon <- read.delim("~/REMC/breast/exon/myo-RM084_lum-RM084.RPKM.corrected", head = F, as.is = T)
lum084_stem084_exon <- read.delim("~/REMC/breast/exon/lum-RM084_stem-RM084.RPKM.corrected", head = F, as.is = T)
myo084_stem084_exon <- read.delim("~/REMC/breast/exon/myo-RM084_stem-RM084.RPKM.corrected", head = F, as.is = T)
# DE genes with <= 5 exons
lum084_myo084_DE5 <- lum084_myo084_DE[lum084_myo084_DE$Nexon <= 5, ]
lum084_stem084_DE5 <- lum084_stem084_DE[lum084_stem084_DE$Nexon <= 5, ]
myo084_stem084_DE5 <- myo084_stem084_DE[myo084_stem084_DE$Nexon <= 5, ]
# No. of DE exons for each gene 
DEexons <- function(id, exons){
  cutoff <- 1
  e <- 1e-6
  id_exon <- exons[grep(id, exons$V1), ]
  id_exon$logfold <- log2((id_exon$V3 + e) / (id_exon$V5 + e))
  return(sum(abs(id_exon$logfold) >= cutoff))
}
vDEexons <- Vectorize(DEexons, vectorize.args = "id")
lum084_myo084_DE5$DEexons <- vDEexons(as.character(lum084_myo084_DE5$id), lum084_myo084_exon)
lum084_stem084_DE5$DEexons <- vDEexons(as.character(lum084_stem084_DE5$id), lum084_stem084_exon)
myo084_stem084_DE5$DEexons <- vDEexons(as.character(myo084_stem084_DE5$id), myo084_stem084_exon)
pdf("PercentDEexons.pdf")
plot(c(0, 1), c(0, 30), type = "n", main = "DE exons in DE genes", xlab = "Proportion of DE exons in DE genes", ylab = "density")
lines(density(lum084_myo084_DE5$DEexons / lum084_myo084_DE5$Nexon), col = 2, lty = 1, lwd = 3)
lines(density(lum084_stem084_DE5$DEexons / lum084_stem084_DE5$Nexon), col = 3, lty = 1, lwd = 3)
lines(density(myo084_stem084_DE5$DEexons / myo084_stem084_DE5$Nexon), col = 4, lty = 1, lwd = 3)
legend("topleft", c("lum084 vs myo084", "lum084 vs stem084", "myo084 vs stem084"), col = 2:4, lty = 1, lwd = 5, cex = 0.8)
dev.off()
rm(lum084_myo084_exon, lum084_stem084_exon, myo084_stem084_exon)

# Position of isoform exons on the gene
# lum 084 vs myo 084
# start / end point of the gene 
lum084_myo084_isoform_all$start <- ensembl[lum084_myo084_isoform_all$id, "start"]
lum084_myo084_isoform_all$end <- ensembl[lum084_myo084_isoform_all$id, "end"]
lum084_myo084_isoform_all$strand <- ensembl[lum084_myo084_isoform_all$id, "strand"]
# exon position: mid point of the exon
lum084_myo084_isoform_all$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", lum084_myo084_isoform_all$V1))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", lum084_myo084_isoform_all$V1))))/2
# exon postition on the gene: calculate + strand & - strand separately
lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == 1, ]$exon_pos <- (lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == 1, ]$exon_pos - lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == 1, ]$start) / (lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == 1, ]$end - lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == 1, ]$start) * 100
lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == -1, ]$exon_pos <- 100 - (lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == -1, ]$exon_pos - lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == -1, ]$start) / (lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == -1, ]$end - lum084_myo084_isoform_all[lum084_myo084_isoform_all$strand == -1, ]$start) * 100
lum084_myo084_isoform_valid$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", lum084_myo084_isoform_valid$exon))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", lum084_myo084_isoform_valid$exon))))/2
lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == 1, ]$exon_pos <- (lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == 1, ]$exon_pos - lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == 1, ]$start) / (lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == 1, ]$end - lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == 1, ]$start) * 100
lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == -1, ]$exon_pos <- 100 - (lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == -1, ]$exon_pos - lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == -1, ]$start) / (lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == -1, ]$end - lum084_myo084_isoform_valid[lum084_myo084_isoform_valid$strand == -1, ]$start) * 100
# lum 084 vs stem 084
lum084_stem084_isoform_all$start <- ensembl[lum084_stem084_isoform_all$id, "start"]
lum084_stem084_isoform_all$end <- ensembl[lum084_stem084_isoform_all$id, "end"]
lum084_stem084_isoform_all$strand <- ensembl[lum084_stem084_isoform_all$id, "strand"]
lum084_stem084_isoform_all$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", lum084_stem084_isoform_all$V1))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", lum084_stem084_isoform_all$V1))))/2
lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == 1, ]$exon_pos <- (lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == 1, ]$exon_pos - lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == 1, ]$start) / (lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == 1, ]$end - lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == 1, ]$start) * 100
lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == -1, ]$exon_pos <- 100 - (lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == -1, ]$exon_pos - lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == -1, ]$start) / (lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == -1, ]$end - lum084_stem084_isoform_all[lum084_stem084_isoform_all$strand == -1, ]$start) * 100
lum084_stem084_isoform_valid$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", lum084_stem084_isoform_valid$exon))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", lum084_stem084_isoform_valid$exon))))/2
lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == 1, ]$exon_pos <- (lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == 1, ]$exon_pos - lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == 1, ]$start) / (lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == 1, ]$end - lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == 1, ]$start) * 100
lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == -1, ]$exon_pos <- 100 - (lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == -1, ]$exon_pos - lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == -1, ]$start) / (lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == -1, ]$end - lum084_stem084_isoform_valid[lum084_stem084_isoform_valid$strand == -1, ]$start) * 100
# myo 084 vs stem 084
myo084_stem084_isoform_all$start <- ensembl[myo084_stem084_isoform_all$id, "start"]
myo084_stem084_isoform_all$end <- ensembl[myo084_stem084_isoform_all$id, "end"]
myo084_stem084_isoform_all$strand <- ensembl[myo084_stem084_isoform_all$id, "strand"]
myo084_stem084_isoform_all$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", myo084_stem084_isoform_all$V1))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", myo084_stem084_isoform_all$V1))))/2
myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == 1, ]$exon_pos <- (myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == 1, ]$exon_pos - myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == 1, ]$start) / (myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == 1, ]$end - myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == 1, ]$start) * 100
myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == -1, ]$exon_pos <- 100 - (myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == -1, ]$exon_pos - myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == -1, ]$start) / (myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == -1, ]$end - myo084_stem084_isoform_all[myo084_stem084_isoform_all$strand == -1, ]$start) * 100
myo084_stem084_isoform_valid$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", myo084_stem084_isoform_valid$exon))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", myo084_stem084_isoform_valid$exon))))/2
myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == 1, ]$exon_pos <- (myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == 1, ]$exon_pos - myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == 1, ]$start) / (myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == 1, ]$end - myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == 1, ]$start) * 100
myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == -1, ]$exon_pos <- 100 - (myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == -1, ]$exon_pos - myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == -1, ]$start) / (myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == -1, ]$end - myo084_stem084_isoform_valid[myo084_stem084_isoform_valid$strand == -1, ]$start) * 100
pdf("exon_position.pdf")
plot(c(0, 100), c(0, 0.015), type = "n", main = "Distribution of isoform exons along genes", xlab = "gene position", ylab = "density")
lines(density(na.omit(lum084_myo084_isoform_all$exon_pos)), col = 1, lty = 1, lwd = 3)
lines(density(na.omit(lum084_stem084_isoform_all$exon_pos)), col = 2, lty = 1, lwd = 3)
lines(density(na.omit(myo084_stem084_isoform_all$exon_pos)), col = 3, lty = 1, lwd = 3)
# lines(density(na.omit(lum084_myo084_isoform_valid$exon_pos)), col = 1, lty = 3, lwd = 3)
# lines(density(na.omit(lum084_stem084_isoform_valid$exon_pos)), col = 2, lty = 3, lwd = 3)
# lines(density(na.omit(myo084_stem084_isoform_valid$exon_pos)), col = 3, lty = 3, lwd = 3)
# legend("topleft", c("all", "validated"), col = 1, lty = c(1, 3), lwd = 5, cex = 0.8)
legend("bottomright", c("lum084 vs myo084", "lum084 vs stem084", "myo084 vs stem084"), col = 1:3, lty = 1, lwd = 5)
dev.off()

# Venn Diagram and average expression level   
lum084 <- read.delim("~/REMC/gene/A17918.G.A.rpkm.pc", head = F, as.is =T)
myo084 <- read.delim("~/REMC/gene/A17919.G.A.rpkm.pc", head = F, as.is =T)
stem084 <- read.delim("~/REMC/gene/A17920.G.A.rpkm.pc", head = F, as.is =T)
RM084gene <- data.frame(id = lum084$V1, lum084 = lum084$V3, myo084 = myo084$V3, stem084 = stem084$V3)
rm(lum084, myo084, stem084)
# all isoforms
isoform_all <- list(lum084_myo084 = lum084_myo084_isoform_all_gene$id, lum084_stem084 = lum084_stem084_isoform_all_gene$id, myo084_stem084 = myo084_stem084_isoform_all_gene$id)
# gene list of different sections of Venn diagram (lm: lum084_myo084, ls: lum084_stem084, ms: myo084_stem084)
all_lm_only <- RM084gene[RM084gene$id %in% setdiff(lum084_myo084_isoform_all_gene$id, union(lum084_stem084_isoform_all_gene$id, myo084_stem084_isoform_all_gene$id)), ]
all_ls_only <- RM084gene[RM084gene$id %in% setdiff(lum084_stem084_isoform_all_gene$id, union(lum084_myo084_isoform_all_gene$id, myo084_stem084_isoform_all_gene$id)), ]
all_ms_only <- RM084gene[RM084gene$id %in% setdiff(myo084_stem084_isoform_all_gene$id, union(lum084_myo084_isoform_all_gene$id, lum084_stem084_isoform_all_gene$id)), ]
all_lm_ls_not_ms <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_myo084_isoform_all_gene$id, lum084_stem084_isoform_all_gene$id), myo084_stem084_isoform_all_gene$id), ]
all_lm_ms_not_ls <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_myo084_isoform_all_gene$id, myo084_stem084_isoform_all_gene$id), lum084_stem084_isoform_all_gene$id), ]
all_ls_ms_not_lm <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_stem084_isoform_all_gene$id, myo084_stem084_isoform_all_gene$id), lum084_myo084_isoform_all_gene$id), ]
all_lm_ls_ms <- RM084gene[RM084gene$id %in% intersect(intersect(lum084_myo084_isoform_all_gene$id, lum084_stem084_isoform_all_gene$id), myo084_stem084_isoform_all_gene$id), ]
write.table(all_lm_only, file = "./Venn/all_lm_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_ls_only, file = "./Venn/all_ls_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_ms_only, file = "./Venn/all_ms_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_lm_ls_not_ms, file = "./Venn/all_lm_ls_not_ms.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_lm_ms_not_ls, file = "./Venn/all_lm_ms_not_ls.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_ls_ms_not_lm, file = "./Venn/all_ls_ms_not_lm.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(all_lm_ls_ms, file = "./Venn/all_lm_ls_ms.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# validated isoforms
isoform_valid <- list(lum084_myo084 = as.character(lum084_myo084_isoform_valid_gene$gene), lum084_stem084 = as.character(lum084_stem084_isoform_valid_gene$gene), myo084_stem084 = as.character(myo084_stem084_isoform_valid_gene$gene))
# gene list of different sections of Venn diagram (lm: lum084_myo084, ls: lum084_stem084, ms: myo084_stem084)
valid_lm_only <- RM084gene[RM084gene$id %in% setdiff(lum084_myo084_isoform_valid_gene$gene, union(lum084_stem084_isoform_valid_gene$gene, myo084_stem084_isoform_valid_gene$gene)), ]
valid_ls_only <- RM084gene[RM084gene$id %in% setdiff(lum084_stem084_isoform_valid_gene$gene, union(lum084_myo084_isoform_valid_gene$gene, myo084_stem084_isoform_valid_gene$gene)), ]
valid_ms_only <- RM084gene[RM084gene$id %in% setdiff(myo084_stem084_isoform_valid_gene$gene, union(lum084_myo084_isoform_valid_gene$gene, lum084_stem084_isoform_valid_gene$gene)), ]
valid_lm_ls_not_ms <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_myo084_isoform_valid_gene$gene, lum084_stem084_isoform_valid_gene$gene), myo084_stem084_isoform_valid_gene$gene), ]
valid_lm_ms_not_ls <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_myo084_isoform_valid_gene$gene, myo084_stem084_isoform_valid_gene$gene), lum084_stem084_isoform_valid_gene$gene), ]
valid_ls_ms_not_lm <- RM084gene[RM084gene$id %in% setdiff(intersect(lum084_stem084_isoform_valid_gene$gene, myo084_stem084_isoform_valid_gene$gene), lum084_myo084_isoform_valid_gene$gene), ]
valid_lm_ls_ms <- RM084gene[RM084gene$id %in% intersect(intersect(lum084_myo084_isoform_valid_gene$gene, lum084_stem084_isoform_valid_gene$gene), myo084_stem084_isoform_valid_gene$gene), ]
write.table(valid_lm_only, file = "./Venn/valid_lm_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_ls_only, file = "./Venn/valid_ls_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_ms_only, file = "./Venn/valid_ms_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_lm_ls_not_ms, file = "./Venn/valid_lm_ls_not_ms.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_lm_ms_not_ls, file = "./Venn/valid_lm_ms_not_ls.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_ls_ms_not_lm, file = "./Venn/valid_ls_ms_not_lm.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_lm_ls_ms, file = "./Venn/valid_lm_ls_ms.txt", sep = "\t", quote = F, row.names = F, col.names = F)
library(VennDiagram)
pdf("venn_isoform_all.pdf")
plot.new()
venn_all <- venn.diagram(isoform_all, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of previously identified isoforms")
grid.draw(venn_all)
grid.text(as.character(round(mean(c(all_lm_only$lum084, all_lm_only$myo084)), 2)), 0.2, 0.65, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(all_ls_only$lum084, all_ls_only$stem084)), 2)), 0.8, 0.65, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(all_ms_only$myo084, all_ms_only$stem084)), 2)), 0.5, 0.25, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(all_lm_ls_not_ms$lum084, all_lm_ls_not_ms$myo084, all_lm_ls_not_ms$stem084)), 2)), 0.5, 0.7, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(all_lm_ms_not_ls$lum084, all_lm_ms_not_ls$myo084, all_lm_ms_not_ls$stem084)), 2)), 0.3, 0.45, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(all_ls_ms_not_lm$lum084, all_ls_ms_not_lm$myo084, all_ls_ms_not_lm$stem084)), 2)), 0.7, 0.45, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(all_lm_ls_ms$lum084, all_lm_ls_ms$myo084, all_lm_ls_ms$stem084)), 2)), 0.5, 0.45, gp=gpar(fontsize=20, col="white"))
dev.off()
pdf("venn_isoform_valid.pdf")
plot.new()
venn_valid <- venn.diagram(isoform_valid, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of validated isoforms")
grid.draw(venn_valid)
grid.text(as.character(round(mean(c(valid_lm_only$lum084, valid_lm_only$myo084)), 2)), 0.2, 0.65, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(valid_ls_only$lum084, valid_ls_only$stem084)), 2)), 0.8, 0.65, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(valid_ms_only$myo084, valid_ms_only$stem084)), 2)), 0.5, 0.25, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(valid_lm_ls_not_ms$lum084, valid_lm_ls_not_ms$myo084, valid_lm_ls_not_ms$stem084)), 2)), 0.5, 0.7, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(valid_lm_ms_not_ls$lum084, valid_lm_ms_not_ls$myo084, valid_lm_ms_not_ls$stem084)), 2)), 0.3, 0.45, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(valid_ls_ms_not_lm$lum084, valid_ls_ms_not_lm$myo084, valid_ls_ms_not_lm$stem084)), 2)), 0.7, 0.45, gp=gpar(fontsize=20, col="white"))
grid.text(as.character(round(mean(c(valid_lm_ls_ms$lum084, valid_lm_ls_ms$myo084, valid_lm_ls_ms$stem084)), 2)), 0.5, 0.45, gp=gpar(fontsize=20, col="white"))
dev.off()

save.image(file = "junction_valid.Rdata")


