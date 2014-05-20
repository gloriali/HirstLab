################################################################################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# validate previously identified isoform exons with junction RPKM 
junction_valid <- matrix(NA, nrow = 7, ncol = 6, dimnames = list(c("lum084_myo084", "lum084_stem084", "myo084_stem084", "lum080_myo080", "lum080_stem080", "myo080_stem080", "lum035_myo035"), c("No.isoform.exons", "No.isoform.genes", "No.exons.with.junction.cov", "No.genes.with.junction.cov", "No.exons.with.junction.support", "No.genes.with.junction.support")))

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
lum080_myo080_DEup <- read.delim("UP.lumRM080_myoRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum080_myo080_DEdn <- read.delim("DN.lumRM080_myoRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum080_myo080_DE <- c(lum080_myo080_DEup$V1, lum080_myo080_DEdn$V1)
lum080_stem080_DEup <- read.delim("UP.lumRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum080_stem080_DEdn <- read.delim("DN.lumRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum080_stem080_DE <- c(lum080_stem080_DEup$V1, lum080_stem080_DEdn$V1)
myo080_stem080_DEup <- read.delim("UP.myoRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo080_stem080_DEdn <- read.delim("DN.myoRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
myo080_stem080_DE <- c(myo080_stem080_DEup$V1, myo080_stem080_DEdn$V1)
rm(lum080_myo080_DEup, lum080_myo080_DEdn, lum080_stem080_DEup, lum080_stem080_DEdn, myo080_stem080_DEup, myo080_stem080_DEdn)
lum035_myo035_DEup <- read.delim("UP.lumRM035_myoRM035.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum035_myo035_DEdn <- read.delim("DN.lumRM035_myoRM035.FDR_0.015.rmin_0.005.Nmin_25", head = F, as.is = T)
lum035_myo035_DE <- c(lum035_myo035_DEup$V1, lum035_myo035_DEdn$V1)
rm(lum035_myo035_DEup, lum035_myo035_DEdn)

# previously identified isoform exons & junction RPKM 
setwd("~/REMC/junction/all/")
# load("junction_RM084_new.Rdata")
load("junction_all.Rdata")
load("~/hg19/hg19v65_genes.Rdata")
junction_exon <- read.delim("~/hg19/hg19v65_junctions_vs_exons_for_genes.relations.unique", head = T, as.is = T)
junction_exon$junctionID <- paste0("chr", junction_exon$junctionID)

# isoform exons with junction support: junction RPKM in one sample < 0.1(cov <= 1), > 0.1 in the other(cov >= 2) 
high_fold <- function(exon, up, sample1, sample2, junctions = junction, junction_exons = junction_exon){
  e <- 1e-6
  covcut <- 1 # cutoff for sum junction coverage of two samples
  junctions <- na.omit(junctions[junction_exons[junction_exons$exonID == exon, "junctionID"], ])
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

load("RM084_isoform_all.Rdata")
# lum084 vs myo084
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
# save.image(file = "junction_valid_RM084_new.Rdata")

# lum080 vs myo080
lum080_myo080up_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84up_isoform.csv", head = T, as.is = T)
lum080_myo080dn_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84dn_isoform.csv", head = T, as.is = T)
lum080_myo080_isoform_all <- rbind(lum080_myo080up_isoform_valid, lum080_myo080dn_isoform_valid)
lum080_myo080_isoform_all <- lum080_myo080_isoform_all[!is.element(lum080_myo080_isoform_all$id, lum080_myo080_DE),]
rm(lum080_myo080up_isoform_valid, lum080_myo080dn_isoform_valid)
lum080_myo080_isoform_valid <- lum080_myo080_isoform_all
(junction_valid["lum080_myo080", "No.isoform.exons"] <- nrow(lum080_myo080_isoform_valid))
(junction_valid["lum080_myo080", "No.isoform.genes"] <- length(unique(lum080_myo080_isoform_valid$id)))
sum(is.element(lum080_myo080_isoform_valid$V1, junction_exon$exonID))
lum080_myo080_isoform_valid <- data.frame(exon = lum080_myo080_isoform_valid$V1, gene = lum080_myo080_isoform_valid$id, DEfine.p = lum080_myo080_isoform_valid$V5, 
                                          lum080gene = lum080_myo080_isoform_valid$ave3, myo080gene = lum080_myo080_isoform_valid$ave2, 
                                          lum080exon = lum080_myo080_isoform_valid$V3, myo080exon = lum080_myo080_isoform_valid$V2, 
                                          t(Vhigh_fold(exon = lum080_myo080_isoform_valid$V1, up = lum080_myo080_isoform_valid$V3 - lum080_myo080_isoform_valid$V2, sample1 = "lum080", sample2 = "myo080")))
colnames(lum080_myo080_isoform_valid)[(ncol(lum080_myo080_isoform_valid)-2):ncol(lum080_myo080_isoform_valid)] <- c("lum080junction", "myo080junction", "junction_logfold")
(junction_valid["lum080_myo080", "No.isoform.exons"] <- nrow(lum080_myo080_isoform_valid))
(junction_valid["lum080_myo080", "No.isoform.genes"] <- length(unique(lum080_myo080_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
lum080_myo080_isoform_valid <- na.omit(lum080_myo080_isoform_valid)
(junction_valid["lum080_myo080", "No.exons.with.junction.cov"] <- nrow(lum080_myo080_isoform_valid))
(junction_valid["lum080_myo080", "No.genes.with.junction.cov"] <- length(unique(lum080_myo080_isoform_valid$gene)))
# No. of exons/genes with junction support
#lum080_myo080_isoform_valid <- lum080_myo080_isoform_valid[abs(lum080_myo080_isoform_valid$junction_logfold) > logfoldcut,]
lum080_myo080_isoform_valid <- lum080_myo080_isoform_valid[(lum080_myo080_isoform_valid$lum080junction >= jcut & lum080_myo080_isoform_valid$myo080junction < jcut) | (lum080_myo080_isoform_valid$lum080junction < jcut & lum080_myo080_isoform_valid$myo080junction >= jcut),]
lum080_myo080_isoform_valid$junction_logfold <- NULL
lum080_myo080_isoform_valid <- lum080_myo080_isoform_valid[order((lum080_myo080_isoform_valid$lum080gene + lum080_myo080_isoform_valid$myo080gene), decreasing = T), ]
lum080_myo080_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum080_myo080_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["lum080_myo080", "No.exons.with.junction.support"] <- nrow(lum080_myo080_isoform_valid))
(junction_valid["lum080_myo080", "No.genes.with.junction.support"] <- length(unique(lum080_myo080_isoform_valid$gene)))
write.table(lum080_myo080_isoform_valid, file = "lum080_myo080_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# lum080 vs stem080
lum080_stem080up_isoform_valid <- read.csv("~/REMC/breast/tissue/lum_stem_84up_isoform.csv", head = T, as.is = T)
lum080_stem080dn_isoform_valid <- read.csv("~/REMC/breast/tissue/lum_stem_84dn_isoform.csv", head = T, as.is = T)
lum080_stem080_isoform_all <- rbind(lum080_stem080up_isoform_valid, lum080_stem080dn_isoform_valid)
lum080_stem080_isoform_all <- lum080_stem080_isoform_all[!is.element(lum080_stem080_isoform_all$id, lum080_stem080_DE),]
rm(lum080_stem080up_isoform_valid, lum080_stem080dn_isoform_valid)
lum080_stem080_isoform_valid <- lum080_stem080_isoform_all
(junction_valid["lum080_stem080", "No.isoform.exons"] <- nrow(lum080_stem080_isoform_valid))
(junction_valid["lum080_stem080", "No.isoform.genes"] <- length(unique(lum080_stem080_isoform_valid$id)))
sum(is.element(lum080_stem080_isoform_valid$V1, junction_exon$exonID))
lum080_stem080_isoform_valid <- data.frame(exon = lum080_stem080_isoform_valid$V1, gene = lum080_stem080_isoform_valid$id, DEfine.p = lum080_stem080_isoform_valid$V5, 
                                           lum080gene = lum080_stem080_isoform_valid$ave2, stem080gene = lum080_stem080_isoform_valid$ave3, 
                                           lum080exon = lum080_stem080_isoform_valid$V2, stem080exon = lum080_stem080_isoform_valid$V3, 
                                           t(Vhigh_fold(exon = lum080_stem080_isoform_valid$V1, up = lum080_stem080_isoform_valid$V2 - lum080_stem080_isoform_valid$V3, sample1 = "lum080", sample2 = "stem080")))
colnames(lum080_stem080_isoform_valid)[(ncol(lum080_stem080_isoform_valid)-2):ncol(lum080_stem080_isoform_valid)] <- c("lum080junction", "stem080junction", "junction_logfold")
(junction_valid["lum080_stem080", "No.isoform.exons"] <- nrow(lum080_stem080_isoform_valid))
(junction_valid["lum080_stem080", "No.isoform.genes"] <- length(unique(lum080_stem080_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
lum080_stem080_isoform_valid <- na.omit(lum080_stem080_isoform_valid)
(junction_valid["lum080_stem080", "No.exons.with.junction.cov"] <- nrow(lum080_stem080_isoform_valid))
(junction_valid["lum080_stem080", "No.genes.with.junction.cov"] <- length(unique(lum080_stem080_isoform_valid$gene)))
# No. of exons/genes with junction support
#lum080_stem080_isoform_valid <- lum080_stem080_isoform_valid[abs(lum080_stem080_isoform_valid$junction_logfold) > logfoldcut,]
lum080_stem080_isoform_valid <- lum080_stem080_isoform_valid[(lum080_stem080_isoform_valid$lum080junction >= jcut & lum080_stem080_isoform_valid$stem080junction < jcut) | (lum080_stem080_isoform_valid$lum080junction < jcut & lum080_stem080_isoform_valid$stem080junction >= jcut),]
lum080_stem080_isoform_valid$junction_logfold <- NULL
lum080_stem080_isoform_valid <- lum080_stem080_isoform_valid[order((lum080_stem080_isoform_valid$lum080gene + lum080_stem080_isoform_valid$stem080gene), decreasing = T), ]
lum080_stem080_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum080_stem080_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["lum080_stem080", "No.exons.with.junction.support"] <- nrow(lum080_stem080_isoform_valid))
(junction_valid["lum080_stem080", "No.genes.with.junction.support"] <- length(unique(lum080_stem080_isoform_valid$gene)))
write.table(lum080_stem080_isoform_valid, file = "lum080_stem080_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# myo080 vs stem080
myo080_stem080up_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_stem_84up_isoform.csv", head = T, as.is = T)
myo080_stem080dn_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_stem_84dn_isoform.csv", head = T, as.is = T)
myo080_stem080_isoform_all <- rbind(myo080_stem080up_isoform_valid, myo080_stem080dn_isoform_valid)
myo080_stem080_isoform_all <- myo080_stem080_isoform_all[!is.element(myo080_stem080_isoform_all$id, myo080_stem080_DE),]
rm(myo080_stem080up_isoform_valid, myo080_stem080dn_isoform_valid)
myo080_stem080_isoform_valid <- myo080_stem080_isoform_all
(junction_valid["myo080_stem080", "No.isoform.exons"] <- nrow(myo080_stem080_isoform_valid))
(junction_valid["myo080_stem080", "No.isoform.genes"] <- length(unique(myo080_stem080_isoform_valid$id)))
sum(is.element(myo080_stem080_isoform_valid$V1, junction_exon$exonID))
myo080_stem080_isoform_valid <- data.frame(exon = myo080_stem080_isoform_valid$V1, gene = myo080_stem080_isoform_valid$id, DEfine.p = myo080_stem080_isoform_valid$V5, 
                                           myo080gene = myo080_stem080_isoform_valid$ave2, stem080gene = myo080_stem080_isoform_valid$ave3, 
                                           myo080exon = myo080_stem080_isoform_valid$V2, stem080exon = myo080_stem080_isoform_valid$V3, 
                                           t(Vhigh_fold(exon = myo080_stem080_isoform_valid$V1, up = myo080_stem080_isoform_valid$V2 - myo080_stem080_isoform_valid$V3, sample1 = "myo080", sample2 = "stem080")))
colnames(myo080_stem080_isoform_valid)[(ncol(myo080_stem080_isoform_valid)-2):ncol(myo080_stem080_isoform_valid)] <- c("myo080junction", "stem080junction", "junction_logfold")
(junction_valid["myo080_stem080", "No.isoform.exons"] <- nrow(myo080_stem080_isoform_valid))
(junction_valid["myo080_stem080", "No.isoform.genes"] <- length(unique(myo080_stem080_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
myo080_stem080_isoform_valid <- na.omit(myo080_stem080_isoform_valid)
(junction_valid["myo080_stem080", "No.exons.with.junction.cov"] <- nrow(myo080_stem080_isoform_valid))
(junction_valid["myo080_stem080", "No.genes.with.junction.cov"] <- length(unique(myo080_stem080_isoform_valid$gene)))
# No. of exons/genes with junction support
#myo080_stem080_isoform_valid <- myo080_stem080_isoform_valid[abs(myo080_stem080_isoform_valid$junction_logfold) > logfoldcut,]
myo080_stem080_isoform_valid <- myo080_stem080_isoform_valid[(myo080_stem080_isoform_valid$myo080junction >= jcut & myo080_stem080_isoform_valid$stem080junction < jcut) | (myo080_stem080_isoform_valid$myo080junction < jcut & myo080_stem080_isoform_valid$stem080junction >= jcut),]
myo080_stem080_isoform_valid$junction_logfold <- NULL
myo080_stem080_isoform_valid <- myo080_stem080_isoform_valid[order((myo080_stem080_isoform_valid$myo080gene + myo080_stem080_isoform_valid$stem080gene), decreasing = T), ]
myo080_stem080_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(myo080_stem080_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["myo080_stem080", "No.exons.with.junction.support"] <- nrow(myo080_stem080_isoform_valid))
(junction_valid["myo080_stem080", "No.genes.with.junction.support"] <- length(unique(myo080_stem080_isoform_valid$gene)))
write.table(myo080_stem080_isoform_valid, file = "myo080_stem080_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# lum035 vs myo035
lum035_myo035up_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84up_isoform.csv", head = T, as.is = T)
lum035_myo035dn_isoform_valid <- read.csv("~/REMC/breast/tissue/myo_lum_84dn_isoform.csv", head = T, as.is = T)
lum035_myo035_isoform_all <- rbind(lum035_myo035up_isoform_valid, lum035_myo035dn_isoform_valid)
lum035_myo035_isoform_all <- lum035_myo035_isoform_all[!is.element(lum035_myo035_isoform_all$id, lum035_myo035_DE),]
rm(lum035_myo035up_isoform_valid, lum035_myo035dn_isoform_valid)
lum035_myo035_isoform_valid <- lum035_myo035_isoform_all
(junction_valid["lum035_myo035", "No.isoform.exons"] <- nrow(lum035_myo035_isoform_valid))
(junction_valid["lum035_myo035", "No.isoform.genes"] <- length(unique(lum035_myo035_isoform_valid$id)))
sum(is.element(lum035_myo035_isoform_valid$V1, junction_exon$exonID))
lum035_myo035_isoform_valid <- data.frame(exon = lum035_myo035_isoform_valid$V1, gene = lum035_myo035_isoform_valid$id, DEfine.p = lum035_myo035_isoform_valid$V5, 
                                          lum035gene = lum035_myo035_isoform_valid$ave3, myo035gene = lum035_myo035_isoform_valid$ave2, 
                                          lum035exon = lum035_myo035_isoform_valid$V3, myo035exon = lum035_myo035_isoform_valid$V2, 
                                          t(Vhigh_fold(exon = lum035_myo035_isoform_valid$V1, up = lum035_myo035_isoform_valid$V3 - lum035_myo035_isoform_valid$V2, sample1 = "lum035", sample2 = "myo035")))
colnames(lum035_myo035_isoform_valid)[(ncol(lum035_myo035_isoform_valid)-2):ncol(lum035_myo035_isoform_valid)] <- c("lum035junction", "myo035junction", "junction_logfold")
(junction_valid["lum035_myo035", "No.isoform.exons"] <- nrow(lum035_myo035_isoform_valid))
(junction_valid["lum035_myo035", "No.isoform.genes"] <- length(unique(lum035_myo035_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
lum035_myo035_isoform_valid <- na.omit(lum035_myo035_isoform_valid)
(junction_valid["lum035_myo035", "No.exons.with.junction.cov"] <- nrow(lum035_myo035_isoform_valid))
(junction_valid["lum035_myo035", "No.genes.with.junction.cov"] <- length(unique(lum035_myo035_isoform_valid$gene)))
# No. of exons/genes with junction support
#lum035_myo035_isoform_valid <- lum035_myo035_isoform_valid[abs(lum035_myo035_isoform_valid$junction_logfold) > logfoldcut,]
lum035_myo035_isoform_valid <- lum035_myo035_isoform_valid[(lum035_myo035_isoform_valid$lum035junction >= jcut & lum035_myo035_isoform_valid$myo035junction < jcut) | (lum035_myo035_isoform_valid$lum035junction < jcut & lum035_myo035_isoform_valid$myo035junction >= jcut),]
lum035_myo035_isoform_valid$junction_logfold <- NULL
lum035_myo035_isoform_valid <- lum035_myo035_isoform_valid[order((lum035_myo035_isoform_valid$lum035gene + lum035_myo035_isoform_valid$myo035gene), decreasing = T), ]
lum035_myo035_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(lum035_myo035_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid["lum035_myo035", "No.exons.with.junction.support"] <- nrow(lum035_myo035_isoform_valid))
(junction_valid["lum035_myo035", "No.genes.with.junction.support"] <- length(unique(lum035_myo035_isoform_valid$gene)))
write.table(lum035_myo035_isoform_valid, file = "lum035_myo035_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rm(junction_exon)

save.image(file = "junction_valid.Rdata")

(junction_valid)

###########################################################################################################
