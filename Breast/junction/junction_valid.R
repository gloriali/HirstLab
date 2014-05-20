# ###########################################################################################################
# # /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# # get junction RPKM 
# setwd("/projects/mbilenky/REMC/breast/RNA-seq/junctions_new/")
# lum084 <- read.delim("A17918.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
# myo084 <- read.delim("A17919.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
# stem084 <- read.delim("A17920.q5.F516.jb.bedGraph.gz", head = F, as.is = T)
# junction <- data.frame(id = unique(c(lum084$V1, myo084$V1, stem084$V1)))
# rownames(junction) <- junction$id
# junction[, c("lum084N", "myo084N", "stem084N")] <- 0
# junction[lum084$V1, "lum084N"] <- lum084$V3
# junction[myo084$V1, "myo084N"] <- myo084$V3
# junction[stem084$V1, "stem084N"] <- stem084$V3
# Nlum084 <- 191290094   # from /projects/epigenomics/ep50/internal/jqc.1.7.6/A17918/A17918.report: Total number of exonic reads for RPKM (protein coding; no MT, no ribo proteins, top expressed 0.005 exons excluded) 
# Nmyo084 <- 215406797
# Nstem084 <- 220214822
# read_length <- 75
# junction$lum084rpkm <- junction$lum084N*10^3*10^6 / (Nlum084*read_length) 
# junction$myo084rpkm <- junction$myo084N*10^3*10^6 / (Nmyo084*read_length) 
# junction$stem084rpkm <- junction$stem084N*10^3*10^6 / (Nstem084*read_length) 
# save(junction, file = "~/REMC/junction/all/junction_RM084_new.Rdata")
# 
###########################################################################################################
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
# junction_exon$junctionID <- paste0("chr", junction_exon$junctionID)

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
write.table(junction_valid_common, file = "junction_valid_common_RM084.txt", sep = "\t", quote = F, row.names = F, col.names = T)
save.image(file = "junction_valid_RM084.Rdata")

###########################################################################################################
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
plot(c(0, 50), c(0, 0.08), type = "n", main = "No. of exons in DE genes and isoforms", xlab = "No. of exons", ylab = "density")
lines(density(na.omit(Nexon[Nexon$id %in% RM084gene[(RM084gene$lum084 + RM084gene$myo084 + RM084gene$stem084) > 0.005,]$id,]$Nexon)), col = 1, lty = 1, lwd = 5)
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
legend("topleft", c("all expressed genes", "DE genes", "isoforms"), col = 1, lty = c(1, 1, 3), lwd = 5, cex = 0.8)
legend("topright", c("all expressed genes", "lum084 vs myo084", "lum084 vs stem084", "myo084 vs stem084"), col = 1:4, lty = 1, lwd = 5, cex = 0.8)
dev.off()

# gene RPKM
lum084 <- read.delim("~/REMC/gene/A17918.G.A.rpkm.pc", head = F, as.is =T)
myo084 <- read.delim("~/REMC/gene/A17919.G.A.rpkm.pc", head = F, as.is =T)
stem084 <- read.delim("~/REMC/gene/A17920.G.A.rpkm.pc", head = F, as.is =T)
RM084gene <- data.frame(id = lum084$V1, lum084 = lum084$V3, myo084 = myo084$V3, stem084 = stem084$V3)
rm(lum084, myo084, stem084)

# Venn Diagram with average expression level and No. of exons              
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
all_lm_only$Nexon <- Nexon[all_lm_only$id, ]$Nexon
all_ls_only$Nexon <- Nexon[all_ls_only$id, ]$Nexon
all_ms_only$Nexon <- Nexon[all_ms_only$id, ]$Nexon
all_lm_ls_not_ms$Nexon <- Nexon[all_lm_ls_not_ms$id, ]$Nexon
all_lm_ms_not_ls$Nexon <- Nexon[all_lm_ms_not_ls$id, ]$Nexon
all_ls_ms_not_lm$Nexon <- Nexon[all_ls_ms_not_lm$id, ]$Nexon
all_lm_ls_ms$Nexon <- Nexon[all_lm_ls_ms$id, ]$Nexon
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
valid_lm_only$Nexon <- Nexon[valid_lm_only$id, ]$Nexon
valid_ls_only$Nexon <- Nexon[valid_ls_only$id, ]$Nexon
valid_ms_only$Nexon <- Nexon[valid_ms_only$id, ]$Nexon
valid_lm_ls_not_ms$Nexon <- Nexon[valid_lm_ls_not_ms$id, ]$Nexon
valid_lm_ms_not_ls$Nexon <- Nexon[valid_lm_ms_not_ls$id, ]$Nexon
valid_ls_ms_not_lm$Nexon <- Nexon[valid_ls_ms_not_lm$id, ]$Nexon
valid_lm_ls_ms$Nexon <- Nexon[valid_lm_ls_ms$id, ]$Nexon
write.table(valid_lm_only, file = "./Venn/valid_lm_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_ls_only, file = "./Venn/valid_ls_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_ms_only, file = "./Venn/valid_ms_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_lm_ls_not_ms, file = "./Venn/valid_lm_ls_not_ms.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_lm_ms_not_ls, file = "./Venn/valid_lm_ms_not_ls.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_ls_ms_not_lm, file = "./Venn/valid_ls_ms_not_lm.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(valid_lm_ls_ms, file = "./Venn/valid_lm_ls_ms.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# library(VennDiagram)
# pdf("venn_isoform_all.pdf")
# plot.new()
# venn_all <- venn.diagram(isoform_all, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of previously identified isoforms")
# grid.draw(venn_all)
# grid.text(as.character(round(mean(c(all_lm_only$lum084, all_lm_only$myo084)), 2)), 0.2, 0.65, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(all_ls_only$lum084, all_ls_only$stem084)), 2)), 0.8, 0.65, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(all_ms_only$myo084, all_ms_only$stem084)), 2)), 0.5, 0.25, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(all_lm_ls_not_ms$lum084, all_lm_ls_not_ms$myo084, all_lm_ls_not_ms$stem084)), 2)), 0.5, 0.7, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(all_lm_ms_not_ls$lum084, all_lm_ms_not_ls$myo084, all_lm_ms_not_ls$stem084)), 2)), 0.3, 0.45, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(all_ls_ms_not_lm$lum084, all_ls_ms_not_lm$myo084, all_ls_ms_not_lm$stem084)), 2)), 0.7, 0.45, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(all_lm_ls_ms$lum084, all_lm_ls_ms$myo084, all_lm_ls_ms$stem084)), 2)), 0.5, 0.45, gp=gpar(fontsize=20, col="white"))
# dev.off()
# pdf("venn_isoform_valid.pdf")
# plot.new()
# venn_valid <- venn.diagram(isoform_valid, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of validated isoforms")
# grid.draw(venn_valid)
# grid.text(as.character(round(mean(c(valid_lm_only$lum084, valid_lm_only$myo084)), 2)), 0.2, 0.65, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(valid_ls_only$lum084, valid_ls_only$stem084)), 2)), 0.8, 0.65, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(valid_ms_only$myo084, valid_ms_only$stem084)), 2)), 0.5, 0.25, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(valid_lm_ls_not_ms$lum084, valid_lm_ls_not_ms$myo084, valid_lm_ls_not_ms$stem084)), 2)), 0.5, 0.7, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(valid_lm_ms_not_ls$lum084, valid_lm_ms_not_ls$myo084, valid_lm_ms_not_ls$stem084)), 2)), 0.3, 0.45, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(valid_ls_ms_not_lm$lum084, valid_ls_ms_not_lm$myo084, valid_ls_ms_not_lm$stem084)), 2)), 0.7, 0.45, gp=gpar(fontsize=20, col="white"))
# grid.text(as.character(round(mean(c(valid_lm_ls_ms$lum084, valid_lm_ls_ms$myo084, valid_lm_ls_ms$stem084)), 2)), 0.5, 0.45, gp=gpar(fontsize=20, col="white"))
# dev.off()

library(venneuler)
pdf("venn_Nexons.pdf")
venn_all_Nexon <- venneuler(c(lum084_myo084 = 45.35, lum084_stem084 = 44.36, myo084_stem084 = 43.62, 
                              "lum084_myo084&lum084_stem084" = 22.56, "lum084_stem084&myo084_stem084" = 21.46, 
                              "lum084_myo084&myo084_stem084" = 22.66, 
                              "lum084_myo084&lum084_stem084&myo084_stem084" = 11.72))
plot(venn_all_Nexon)
mtext("Average No. of exons for all isoforms", cex = 2, side = 3)
text(0.55, 0.3, labels = as.character(round(mean(all_lm_only$Nexon), 2)), col="black", cex = 2)
text(0.3, 0.65, labels = as.character(round(mean(all_ls_only$Nexon), 2)), col="black", cex = 2)
text(0.7, 0.65, labels = as.character(round(mean(all_ms_only$Nexon), 2)), col="black", cex = 2)
text(0.38, 0.42, labels = as.character(round(mean(all_lm_ls_not_ms$Nexon), 2)), col="black", cex = 2)
text(0.65, 0.45, labels = as.character(round(mean(all_lm_ms_not_ls$Nexon), 2)), col="black", cex = 2)
text(0.5, 0.65, labels = as.character(round(mean(all_ls_ms_not_lm$Nexon), 2)), col="black", cex = 2)
text(0.5, 0.52, labels = as.character(round(mean(all_lm_ls_ms$Nexon), 2)), col="black", cex = 2)
venn_valid_Nexon <- venneuler(c(lum084_myo084 = 45.35, lum084_stem084 = 44.36, myo084_stem084 = 43.62, 
                              "lum084_myo084&lum084_stem084" = 22.56, "lum084_stem084&myo084_stem084" = 21.46, 
                              "lum084_myo084&myo084_stem084" = 22.66, 
                              "lum084_myo084&lum084_stem084&myo084_stem084" = 11.72))
plot(venn_valid_Nexon)
mtext("Average No. of exons for validated isoforms", cex = 2, side = 3)
text(0.55, 0.3, labels = as.character(round(mean(valid_lm_only$Nexon), 2)), col="black", cex = 2)
text(0.3, 0.65, labels = as.character(round(mean(valid_ls_only$Nexon), 2)), col="black", cex = 2)
text(0.7, 0.65, labels = as.character(round(mean(valid_ms_only$Nexon), 2)), col="black", cex = 2)
text(0.38, 0.42, labels = as.character(round(mean(valid_lm_ls_not_ms$Nexon), 2)), col="black", cex = 2)
text(0.65, 0.45, labels = as.character(round(mean(valid_lm_ms_not_ls$Nexon), 2)), col="black", cex = 2)
text(0.5, 0.65, labels = as.character(round(mean(valid_ls_ms_not_lm$Nexon), 2)), col="black", cex = 2)
text(0.5, 0.52, labels = as.character(round(mean(valid_lm_ls_ms$Nexon), 2)), col="black", cex = 2)
dev.off()

pdf("venn_rpkm.pdf")
venn_all_rpkm <- venneuler(c(lum084_myo084 = 2.27, lum084_stem084 = 2.86, myo084_stem084 = 2.24, 
                              "lum084_myo084&lum084_stem084" = 0.81, "lum084_stem084&myo084_stem084" = 0.84, 
                              "lum084_myo084&myo084_stem084" = 0.65, 
                              "lum084_myo084&lum084_stem084&myo084_stem084" = 0.18))
plot(venn_all_rpkm)
mtext("Average gene RPKM for all isoforms", cex = 2, side = 3)
text(0.7, 0.65, labels = as.character(round(mean(c(all_lm_only$lum084, all_lm_only$myo084)), 2)), col="black", cex = 2)
text(0.5, 0.25, labels = as.character(round(mean(c(all_ls_only$lum084, all_ls_only$stem084)), 2)), col="black", cex = 2)
text(0.35, 0.7, labels = as.character(round(mean(c(all_ms_only$myo084, all_ms_only$stem084)), 2)), col="black", cex = 2)
text(0.6, 0.45, labels = as.character(round(mean(c(all_lm_ls_not_ms$lum084, all_lm_ls_not_ms$myo084, all_lm_ls_not_ms$stem084)), 2)), col="black", cex = 2)
text(0.52, 0.65, labels = as.character(round(mean(c(all_lm_ms_not_ls$lum084, all_lm_ms_not_ls$myo084, all_lm_ms_not_ls$stem084)), 2)), col="black", cex = 2)
text(0.38, 0.48, labels = as.character(round(mean(c(all_ls_ms_not_lm$lum084, all_ls_ms_not_lm$myo084, all_ls_ms_not_lm$stem084)), 2)), col="black", cex = 2)
text(0.5, 0.52, labels = as.character(round(mean(c(all_lm_ls_ms$lum084, all_lm_ls_ms$myo084, all_lm_ls_ms$stem084)), 2)), col="black", cex = 2)
venn_valid_rpkm <- venneuler(c(lum084_myo084 = 1.76, lum084_stem084 = 1.93, myo084_stem084 = 2.58, 
                             "lum084_myo084&lum084_stem084" = 0.72, "lum084_stem084&myo084_stem084" = 0.71, 
                             "lum084_myo084&myo084_stem084" = 0.51, 
                             "lum084_myo084&lum084_stem084&myo084_stem084" = 0.19))
plot(venn_valid_rpkm)
mtext("Average gene RPKM for validated isoforms", cex = 2, side = 3)
text(0.7, 0.65, labels = as.character(round(mean(c(valid_lm_only$lum084, valid_lm_only$myo084)), 2)), col="black", cex = 2)
text(0.55, 0.3, labels = as.character(round(mean(c(valid_ls_only$lum084, valid_ls_only$stem084)), 2)), col="black", cex = 2)
text(0.3, 0.65, labels = as.character(round(mean(c(valid_ms_only$myo084, valid_ms_only$stem084)), 2)), col="black", cex = 2)
text(0.62, 0.48, labels = as.character(round(mean(c(valid_lm_ls_not_ms$lum084, valid_lm_ls_not_ms$myo084, valid_lm_ls_not_ms$stem084)), 2)), col="black", cex = 2)
text(0.48, 0.6, labels = as.character(round(mean(c(valid_lm_ms_not_ls$lum084, valid_lm_ms_not_ls$myo084, valid_lm_ms_not_ls$stem084)), 2)), col="black", cex = 2)
text(0.42, 0.42, labels = as.character(round(mean(c(valid_ls_ms_not_lm$lum084, valid_ls_ms_not_lm$myo084, valid_ls_ms_not_lm$stem084)), 2)), col="black", cex = 2)
text(0.5, 0.52, labels = as.character(round(mean(c(valid_lm_ls_ms$lum084, valid_lm_ls_ms$myo084, valid_lm_ls_ms$stem084)), 2)), col="black", cex = 2)
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

# average exon length for isoforms on Venn diagram
lum084_myo084_isoform_all$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", lum084_myo084_isoform_all$V1)))
lum084_myo084_isoform_all$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", lum084_myo084_isoform_all$V1)))
lum084_myo084_isoform_all$exon_length <- abs(lum084_myo084_isoform_all$exon_end - lum084_myo084_isoform_all$exon_start)
lum084_stem084_isoform_all$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", lum084_stem084_isoform_all$V1)))
lum084_stem084_isoform_all$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", lum084_stem084_isoform_all$V1)))
lum084_stem084_isoform_all$exon_length <- abs(lum084_stem084_isoform_all$exon_end - lum084_stem084_isoform_all$exon_start)
myo084_stem084_isoform_all$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", myo084_stem084_isoform_all$V1)))
myo084_stem084_isoform_all$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", myo084_stem084_isoform_all$V1)))
myo084_stem084_isoform_all$exon_length <- abs(myo084_stem084_isoform_all$exon_end - myo084_stem084_isoform_all$exon_start)
all_lm_only$exon_length <- sapply(all_lm_only$id, function(x) mean(lum084_myo084_isoform_all[lum084_myo084_isoform_all$id == x, "exon_length"]))
all_ls_only$exon_length <- sapply(all_ls_only$id, function(x) mean(lum084_stem084_isoform_all[lum084_stem084_isoform_all$id == x, "exon_length"]))
all_ms_only$exon_length <- sapply(all_ms_only$id, function(x) mean(myo084_stem084_isoform_all[myo084_stem084_isoform_all$id == x, "exon_length"]))
all_lm_ls_not_ms$exon_length <- sapply(all_lm_ls_not_ms$id, function(x) mean(c(lum084_myo084_isoform_all[lum084_myo084_isoform_all$id == x, "exon_length"], lum084_stem084_isoform_all[lum084_stem084_isoform_all$id == x, "exon_length"])))
all_lm_ms_not_ls$exon_length <- sapply(all_lm_ms_not_ls$id, function(x) mean(c(lum084_myo084_isoform_all[lum084_myo084_isoform_all$id == x, "exon_length"], myo084_stem084_isoform_all[myo084_stem084_isoform_all$id == x, "exon_length"])))
all_ls_ms_not_lm$exon_length <- sapply(all_ls_ms_not_lm$id, function(x) mean(c(lum084_stem084_isoform_all[lum084_stem084_isoform_all$id == x, "exon_length"], myo084_stem084_isoform_all[myo084_stem084_isoform_all$id == x, "exon_length"])))
all_lm_ls_ms$exon_length <- sapply(all_lm_ls_ms$id, function(x) mean(c(lum084_myo084_isoform_all[lum084_myo084_isoform_all$id == x, "exon_length"], lum084_stem084_isoform_all[lum084_stem084_isoform_all$id == x, "exon_length"], myo084_stem084_isoform_all[myo084_stem084_isoform_all$id == x, "exon_length"])))

lum084_myo084_isoform_valid$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", lum084_myo084_isoform_valid$exon)))
lum084_myo084_isoform_valid$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", lum084_myo084_isoform_valid$exon)))
lum084_myo084_isoform_valid$exon_length <- abs(lum084_myo084_isoform_valid$exon_end - lum084_myo084_isoform_valid$exon_start)
lum084_stem084_isoform_valid$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", lum084_stem084_isoform_valid$exon)))
lum084_stem084_isoform_valid$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", lum084_stem084_isoform_valid$exon)))
lum084_stem084_isoform_valid$exon_length <- abs(lum084_stem084_isoform_valid$exon_end - lum084_stem084_isoform_valid$exon_start)
myo084_stem084_isoform_valid$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", myo084_stem084_isoform_valid$exon)))
myo084_stem084_isoform_valid$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", myo084_stem084_isoform_valid$exon)))
myo084_stem084_isoform_valid$exon_length <- abs(myo084_stem084_isoform_valid$exon_end - myo084_stem084_isoform_valid$exon_start)
valid_lm_only$exon_length <- sapply(valid_lm_only$id, function(x) mean(lum084_myo084_isoform_valid[as.character(lum084_myo084_isoform_valid$gene) == x, "exon_length"]))
valid_ls_only$exon_length <- sapply(valid_ls_only$id, function(x) mean(lum084_stem084_isoform_valid[as.character(lum084_stem084_isoform_valid$gene) == x, "exon_length"]))
valid_ms_only$exon_length <- sapply(valid_ms_only$id, function(x) mean(myo084_stem084_isoform_valid[as.character(myo084_stem084_isoform_valid$gene) == x, "exon_length"]))
valid_lm_ls_not_ms$exon_length <- sapply(valid_lm_ls_not_ms$id, function(x) mean(c(lum084_myo084_isoform_valid[as.character(lum084_myo084_isoform_valid$gene) == x, "exon_length"], lum084_stem084_isoform_valid[as.character(lum084_stem084_isoform_valid$gene) == x, "exon_length"])))
valid_lm_ms_not_ls$exon_length <- sapply(valid_lm_ms_not_ls$id, function(x) mean(c(lum084_myo084_isoform_valid[as.character(lum084_myo084_isoform_valid$gene) == x, "exon_length"], myo084_stem084_isoform_valid[as.character(myo084_stem084_isoform_valid$gene) == x, "exon_length"])))
valid_ls_ms_not_lm$exon_length <- sapply(valid_ls_ms_not_lm$id, function(x) mean(c(lum084_stem084_isoform_valid[as.character(lum084_stem084_isoform_valid$gene) == x, "exon_length"], myo084_stem084_isoform_valid[as.character(myo084_stem084_isoform_valid$gene) == x, "exon_length"])))
valid_lm_ls_ms$exon_length <- sapply(valid_lm_ls_ms$id, function(x) mean(c(lum084_myo084_isoform_valid[as.character(lum084_myo084_isoform_valid$gene) == x, "exon_length"], lum084_stem084_isoform_valid[as.character(lum084_stem084_isoform_valid$gene) == x, "exon_length"], myo084_stem084_isoform_valid[as.character(myo084_stem084_isoform_valid$gene) == x, "exon_length"])))

write.table(lum084_myo084_isoform_all_gene, file = "~/快盘/REMC/junction/Venn/lum084_myo084_isoform_all_gene.txt", sep = "\t", col.names = F, row.names = F)
write.table(lum084_stem084_isoform_all_gene, file = "~/快盘/REMC/junction/Venn/lum084_stem084_isoform_all_gene.txt", sep = "\t", col.names = F, row.names = F)
write.table(myo084_stem084_isoform_all_gene, file = "~/快盘/REMC/junction/Venn/myo084_stem084_isoform_all_gene.txt", sep = "\t", col.names = F, row.names = F)
write.table(lum084_myo084_isoform_valid_gene, file = "~/快盘/REMC/junction/Venn/lum084_myo084_isoform_valid_gene.txt", sep = "\t", col.names = F, row.names = F)
write.table(lum084_stem084_isoform_valid_gene, file = "~/快盘/REMC/junction/Venn/lum084_stem084_isoform_valid_gene.txt", sep = "\t", col.names = F, row.names = F)
write.table(myo084_stem084_isoform_valid_gene, file = "~/快盘/REMC/junction/Venn/myo084_stem084_isoform_valid_gene.txt", sep = "\t", col.names = F, row.names = F)
save.image(file = "~/快盘/REMC/junction/junction_valid.Rdata")

###########################################################################################################








