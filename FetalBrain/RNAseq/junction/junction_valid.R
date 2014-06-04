################################################################################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# validate previously identified isoform exons with junction RPKM 

# previously identified isoform exons & junction RPKM 
setwd("~/FetalBrain/RNAseq/junction")
load("junction.Rdata")
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

## cortex vs GE
junction_valid_cortexge <- matrix(NA, nrow = 4, ncol = 6, dimnames = list(c("HuFNSC01", "HuFNSC02", "HuFNSC03", "HuFNSC04"), c("No.isoform.exons", "No.isoform.genes", "No.exons.with.junction.cov", "No.genes.with.junction.cov", "No.exons.with.junction.support", "No.genes.with.junction.support")))

# HuFNSC01
HuFNSC01_isoform_all <- read.delim("~/FetalBrain/RNAseq/isoform/cortex_ge_01_isoform.txt", as.is =T)
HuFNSC01_isoform_valid <- HuFNSC01_isoform_all
(junction_valid_cortexge["HuFNSC01", "No.isoform.exons"] <- nrow(HuFNSC01_isoform_valid))
(junction_valid_cortexge["HuFNSC01", "No.isoform.genes"] <- length(unique(HuFNSC01_isoform_valid$id)))
sum(is.element(HuFNSC01_isoform_valid$V1, junction_exon$exonID))
HuFNSC01_isoform_valid <- data.frame(exon = HuFNSC01_isoform_valid$V1, gene = HuFNSC01_isoform_valid$id, DEfine.p = HuFNSC01_isoform_valid$V5, 
                                          cortex01gene = HuFNSC01_isoform_valid$ave2, ge01gene = HuFNSC01_isoform_valid$ave3, 
                                          cortex01exon = HuFNSC01_isoform_valid$V2, ge01exon = HuFNSC01_isoform_valid$V3, 
                                          t(Vhigh_fold(exon = HuFNSC01_isoform_valid$V1, up = HuFNSC01_isoform_valid$V2 - HuFNSC01_isoform_valid$V3, sample1 = "cortex01", sample2 = "ge01")))
colnames(HuFNSC01_isoform_valid)[(ncol(HuFNSC01_isoform_valid) - 2):ncol(HuFNSC01_isoform_valid)] <- c("cortex01junction", "ge01junction", "junction_logfold")
(junction_valid_cortexge["HuFNSC01", "No.isoform.exons"] <- nrow(HuFNSC01_isoform_valid))
(junction_valid_cortexge["HuFNSC01", "No.isoform.genes"] <- length(unique(HuFNSC01_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
HuFNSC01_isoform_valid <- na.omit(HuFNSC01_isoform_valid)
(junction_valid_cortexge["HuFNSC01", "No.exons.with.junction.cov"] <- nrow(HuFNSC01_isoform_valid))
(junction_valid_cortexge["HuFNSC01", "No.genes.with.junction.cov"] <- length(unique(HuFNSC01_isoform_valid$gene)))
# No. of exons/genes with junction support
HuFNSC01_isoform_valid <- HuFNSC01_isoform_valid[(HuFNSC01_isoform_valid$cortex01junction >= jcut & HuFNSC01_isoform_valid$ge01junction < jcut) | (HuFNSC01_isoform_valid$cortex01junction < jcut & HuFNSC01_isoform_valid$ge01junction >= jcut),]
HuFNSC01_isoform_valid$junction_logfold <- NULL
HuFNSC01_isoform_valid <- HuFNSC01_isoform_valid[order((HuFNSC01_isoform_valid$cortex01gene + HuFNSC01_isoform_valid$ge01gene), decreasing = T), ]
HuFNSC01_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(HuFNSC01_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid_cortexge["HuFNSC01", "No.exons.with.junction.support"] <- nrow(HuFNSC01_isoform_valid))
(junction_valid_cortexge["HuFNSC01", "No.genes.with.junction.support"] <- length(unique(HuFNSC01_isoform_valid$gene)))
write.table(HuFNSC01_isoform_valid, file = "cortex_ge_01_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
HuFNSC01_isoform_valid_gene <- HuFNSC01_isoform_valid[!duplicated(HuFNSC01_isoform_valid$gene), ]
write.table(HuFNSC01_isoform_valid_gene, file = "cortex_ge_01_isoform_valid_gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# HuFNSC02
HuFNSC02_isoform_all <- read.delim("~/FetalBrain/RNAseq/isoform/cortex_ge_02_isoform.txt", as.is =T)
HuFNSC02_isoform_valid <- HuFNSC02_isoform_all
(junction_valid_cortexge["HuFNSC02", "No.isoform.exons"] <- nrow(HuFNSC02_isoform_valid))
(junction_valid_cortexge["HuFNSC02", "No.isoform.genes"] <- length(unique(HuFNSC02_isoform_valid$id)))
sum(is.element(HuFNSC02_isoform_valid$V1, junction_exon$exonID))
HuFNSC02_isoform_valid <- data.frame(exon = HuFNSC02_isoform_valid$V1, gene = HuFNSC02_isoform_valid$id, DEfine.p = HuFNSC02_isoform_valid$V5, 
                                     cortex02gene = HuFNSC02_isoform_valid$ave2, ge02gene = HuFNSC02_isoform_valid$ave3, 
                                     cortex02exon = HuFNSC02_isoform_valid$V2, ge02exon = HuFNSC02_isoform_valid$V3, 
                                     t(Vhigh_fold(exon = HuFNSC02_isoform_valid$V1, up = HuFNSC02_isoform_valid$V2 - HuFNSC02_isoform_valid$V3, sample1 = "cortex02", sample2 = "ge02")))
colnames(HuFNSC02_isoform_valid)[(ncol(HuFNSC02_isoform_valid) - 2):ncol(HuFNSC02_isoform_valid)] <- c("cortex02junction", "ge02junction", "junction_logfold")
(junction_valid_cortexge["HuFNSC02", "No.isoform.exons"] <- nrow(HuFNSC02_isoform_valid))
(junction_valid_cortexge["HuFNSC02", "No.isoform.genes"] <- length(unique(HuFNSC02_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
HuFNSC02_isoform_valid <- na.omit(HuFNSC02_isoform_valid)
(junction_valid_cortexge["HuFNSC02", "No.exons.with.junction.cov"] <- nrow(HuFNSC02_isoform_valid))
(junction_valid_cortexge["HuFNSC02", "No.genes.with.junction.cov"] <- length(unique(HuFNSC02_isoform_valid$gene)))
# No. of exons/genes with junction support
HuFNSC02_isoform_valid <- HuFNSC02_isoform_valid[(HuFNSC02_isoform_valid$cortex02junction >= jcut & HuFNSC02_isoform_valid$ge02junction < jcut) | (HuFNSC02_isoform_valid$cortex02junction < jcut & HuFNSC02_isoform_valid$ge02junction >= jcut),]
HuFNSC02_isoform_valid$junction_logfold <- NULL
HuFNSC02_isoform_valid <- HuFNSC02_isoform_valid[order((HuFNSC02_isoform_valid$cortex02gene + HuFNSC02_isoform_valid$ge02gene), decreasing = T), ]
HuFNSC02_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(HuFNSC02_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid_cortexge["HuFNSC02", "No.exons.with.junction.support"] <- nrow(HuFNSC02_isoform_valid))
(junction_valid_cortexge["HuFNSC02", "No.genes.with.junction.support"] <- length(unique(HuFNSC02_isoform_valid$gene)))
write.table(HuFNSC02_isoform_valid, file = "cortex_ge_02_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
HuFNSC02_isoform_valid_gene <- HuFNSC02_isoform_valid[!duplicated(HuFNSC02_isoform_valid$gene), ]
write.table(HuFNSC02_isoform_valid_gene, file = "cortex_ge_02_isoform_valid_gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# HuFNSC03
HuFNSC03_isoform_all <- read.delim("~/FetalBrain/RNAseq/isoform/cortex_ge_03_isoform.txt", as.is =T)
HuFNSC03_isoform_valid <- HuFNSC03_isoform_all
(junction_valid_cortexge["HuFNSC03", "No.isoform.exons"] <- nrow(HuFNSC03_isoform_valid))
(junction_valid_cortexge["HuFNSC03", "No.isoform.genes"] <- length(unique(HuFNSC03_isoform_valid$id)))
sum(is.element(HuFNSC03_isoform_valid$V1, junction_exon$exonID))
HuFNSC03_isoform_valid <- data.frame(exon = HuFNSC03_isoform_valid$V1, gene = HuFNSC03_isoform_valid$id, DEfine.p = HuFNSC03_isoform_valid$V5, 
                                     cortex03gene = HuFNSC03_isoform_valid$ave2, ge03gene = HuFNSC03_isoform_valid$ave3, 
                                     cortex03exon = HuFNSC03_isoform_valid$V2, ge03exon = HuFNSC03_isoform_valid$V3, 
                                     t(Vhigh_fold(exon = HuFNSC03_isoform_valid$V1, up = HuFNSC03_isoform_valid$V2 - HuFNSC03_isoform_valid$V3, sample1 = "cortex03", sample2 = "ge03")))
colnames(HuFNSC03_isoform_valid)[(ncol(HuFNSC03_isoform_valid) - 2):ncol(HuFNSC03_isoform_valid)] <- c("cortex03junction", "ge03junction", "junction_logfold")
(junction_valid_cortexge["HuFNSC03", "No.isoform.exons"] <- nrow(HuFNSC03_isoform_valid))
(junction_valid_cortexge["HuFNSC03", "No.isoform.genes"] <- length(unique(HuFNSC03_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
HuFNSC03_isoform_valid <- na.omit(HuFNSC03_isoform_valid)
(junction_valid_cortexge["HuFNSC03", "No.exons.with.junction.cov"] <- nrow(HuFNSC03_isoform_valid))
(junction_valid_cortexge["HuFNSC03", "No.genes.with.junction.cov"] <- length(unique(HuFNSC03_isoform_valid$gene)))
# No. of exons/genes with junction support
HuFNSC03_isoform_valid <- HuFNSC03_isoform_valid[(HuFNSC03_isoform_valid$cortex03junction >= jcut & HuFNSC03_isoform_valid$ge03junction < jcut) | (HuFNSC03_isoform_valid$cortex03junction < jcut & HuFNSC03_isoform_valid$ge03junction >= jcut),]
HuFNSC03_isoform_valid$junction_logfold <- NULL
HuFNSC03_isoform_valid <- HuFNSC03_isoform_valid[order((HuFNSC03_isoform_valid$cortex03gene + HuFNSC03_isoform_valid$ge03gene), decreasing = T), ]
HuFNSC03_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(HuFNSC03_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid_cortexge["HuFNSC03", "No.exons.with.junction.support"] <- nrow(HuFNSC03_isoform_valid))
(junction_valid_cortexge["HuFNSC03", "No.genes.with.junction.support"] <- length(unique(HuFNSC03_isoform_valid$gene)))
write.table(HuFNSC03_isoform_valid, file = "cortex_ge_03_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
HuFNSC03_isoform_valid_gene <- HuFNSC03_isoform_valid[!duplicated(HuFNSC03_isoform_valid$gene), ]
write.table(HuFNSC03_isoform_valid_gene, file = "cortex_ge_03_isoform_valid_gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# HuFNSC04
HuFNSC04_isoform_all <- read.delim("~/FetalBrain/RNAseq/isoform/cortex_ge_04_isoform.txt", as.is =T)
HuFNSC04_isoform_valid <- HuFNSC04_isoform_all
(junction_valid_cortexge["HuFNSC04", "No.isoform.exons"] <- nrow(HuFNSC04_isoform_valid))
(junction_valid_cortexge["HuFNSC04", "No.isoform.genes"] <- length(unique(HuFNSC04_isoform_valid$id)))
sum(is.element(HuFNSC04_isoform_valid$V1, junction_exon$exonID))
HuFNSC04_isoform_valid <- data.frame(exon = HuFNSC04_isoform_valid$V1, gene = HuFNSC04_isoform_valid$id, DEfine.p = HuFNSC04_isoform_valid$V5, 
                                     cortex04gene = HuFNSC04_isoform_valid$ave2, ge04gene = HuFNSC04_isoform_valid$ave3, 
                                     cortex04exon = HuFNSC04_isoform_valid$V2, ge04exon = HuFNSC04_isoform_valid$V3, 
                                     t(Vhigh_fold(exon = HuFNSC04_isoform_valid$V1, up = HuFNSC04_isoform_valid$V2 - HuFNSC04_isoform_valid$V3, sample1 = "cortex04", sample2 = "ge04")))
colnames(HuFNSC04_isoform_valid)[(ncol(HuFNSC04_isoform_valid) - 2):ncol(HuFNSC04_isoform_valid)] <- c("cortex04junction", "ge04junction", "junction_logfold")
(junction_valid_cortexge["HuFNSC04", "No.isoform.exons"] <- nrow(HuFNSC04_isoform_valid))
(junction_valid_cortexge["HuFNSC04", "No.isoform.genes"] <- length(unique(HuFNSC04_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
HuFNSC04_isoform_valid <- na.omit(HuFNSC04_isoform_valid)
(junction_valid_cortexge["HuFNSC04", "No.exons.with.junction.cov"] <- nrow(HuFNSC04_isoform_valid))
(junction_valid_cortexge["HuFNSC04", "No.genes.with.junction.cov"] <- length(unique(HuFNSC04_isoform_valid$gene)))
# No. of exons/genes with junction support
HuFNSC04_isoform_valid <- HuFNSC04_isoform_valid[(HuFNSC04_isoform_valid$cortex04junction >= jcut & HuFNSC04_isoform_valid$ge04junction < jcut) | (HuFNSC04_isoform_valid$cortex04junction < jcut & HuFNSC04_isoform_valid$ge04junction >= jcut),]
HuFNSC04_isoform_valid$junction_logfold <- NULL
HuFNSC04_isoform_valid <- HuFNSC04_isoform_valid[order((HuFNSC04_isoform_valid$cortex04gene + HuFNSC04_isoform_valid$ge04gene), decreasing = T), ]
HuFNSC04_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(HuFNSC04_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid_cortexge["HuFNSC04", "No.exons.with.junction.support"] <- nrow(HuFNSC04_isoform_valid))
(junction_valid_cortexge["HuFNSC04", "No.genes.with.junction.support"] <- length(unique(HuFNSC04_isoform_valid$gene)))
write.table(HuFNSC04_isoform_valid, file = "cortex_ge_04_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
HuFNSC04_isoform_valid_gene <- HuFNSC04_isoform_valid[!duplicated(HuFNSC04_isoform_valid$gene), ]
write.table(HuFNSC04_isoform_valid_gene, file = "cortex_ge_04_isoform_valid_gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)
(junction_valid_cortexge)

## individual: HuFNSC01 vs HuFNSC02
junction_valid_individual <- matrix(NA, nrow = 5, ncol = 6, dimnames = list(c("brain01_02", "cortex01_02", "ge01_02", "cortex03_04", "ge03_04"), c("No.isoform.exons", "No.isoform.genes", "No.exons.with.junction.cov", "No.genes.with.junction.cov", "No.exons.with.junction.support", "No.genes.with.junction.support")))

# brain01_02
brain01_02_isoform_all <- read.delim("~/FetalBrain/RNAseq/isoform/brain01_02_isoform.txt", as.is =T)
brain01_02_isoform_valid <- brain01_02_isoform_all
(junction_valid_individual["brain01_02", "No.isoform.exons"] <- nrow(brain01_02_isoform_valid))
(junction_valid_individual["brain01_02", "No.isoform.genes"] <- length(unique(brain01_02_isoform_valid$id)))
sum(is.element(brain01_02_isoform_valid$V1, junction_exon$exonID))
brain01_02_isoform_valid <- data.frame(exon = brain01_02_isoform_valid$V1, gene = brain01_02_isoform_valid$id, DEfine.p = brain01_02_isoform_valid$V5, 
                                     brain01gene = brain01_02_isoform_valid$ave2, brain02gene = brain01_02_isoform_valid$ave3, 
                                     brain01exon = brain01_02_isoform_valid$V2, brain02exon = brain01_02_isoform_valid$V3, 
                                     t(Vhigh_fold(exon = brain01_02_isoform_valid$V1, up = brain01_02_isoform_valid$V2 - brain01_02_isoform_valid$V3, sample1 = "brain01", sample2 = "brain02")))
colnames(brain01_02_isoform_valid)[(ncol(brain01_02_isoform_valid) - 2):ncol(brain01_02_isoform_valid)] <- c("brain01junction", "brain02junction", "junction_logfold")
(junction_valid_individual["brain01_02", "No.isoform.exons"] <- nrow(brain01_02_isoform_valid))
(junction_valid_individual["brain01_02", "No.isoform.genes"] <- length(unique(brain01_02_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
brain01_02_isoform_valid <- na.omit(brain01_02_isoform_valid)
(junction_valid_individual["brain01_02", "No.exons.with.junction.cov"] <- nrow(brain01_02_isoform_valid))
(junction_valid_individual["brain01_02", "No.genes.with.junction.cov"] <- length(unique(brain01_02_isoform_valid$gene)))
# No. of exons/genes with junction support
brain01_02_isoform_valid <- brain01_02_isoform_valid[(brain01_02_isoform_valid$brain01junction >= jcut & brain01_02_isoform_valid$brain02junction < jcut) | (brain01_02_isoform_valid$brain01junction < jcut & brain01_02_isoform_valid$brain02junction >= jcut),]
brain01_02_isoform_valid$junction_logfold <- NULL
brain01_02_isoform_valid <- brain01_02_isoform_valid[order((brain01_02_isoform_valid$brain01gene + brain01_02_isoform_valid$brain02gene), decreasing = T), ]
brain01_02_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(brain01_02_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid_individual["brain01_02", "No.exons.with.junction.support"] <- nrow(brain01_02_isoform_valid))
(junction_valid_individual["brain01_02", "No.genes.with.junction.support"] <- length(unique(brain01_02_isoform_valid$gene)))
write.table(brain01_02_isoform_valid, file = "brain_01_02_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
brain01_02_isoform_valid_gene <- brain01_02_isoform_valid[!duplicated(brain01_02_isoform_valid$gene), ]
write.table(brain01_02_isoform_valid_gene, file = "brain_01_02_isoform_valid_gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# cortex01_02
cortex01_02_isoform_all <- read.delim("~/FetalBrain/RNAseq/isoform/cortex01_02_isoform.txt", as.is =T)
cortex01_02_isoform_valid <- cortex01_02_isoform_all
(junction_valid_individual["cortex01_02", "No.isoform.exons"] <- nrow(cortex01_02_isoform_valid))
(junction_valid_individual["cortex01_02", "No.isoform.genes"] <- length(unique(cortex01_02_isoform_valid$id)))
sum(is.element(cortex01_02_isoform_valid$V1, junction_exon$exonID))
cortex01_02_isoform_valid <- data.frame(exon = cortex01_02_isoform_valid$V1, gene = cortex01_02_isoform_valid$id, DEfine.p = cortex01_02_isoform_valid$V5, 
                                       cortex01gene = cortex01_02_isoform_valid$ave2, cortex02gene = cortex01_02_isoform_valid$ave3, 
                                       cortex01exon = cortex01_02_isoform_valid$V2, cortex02exon = cortex01_02_isoform_valid$V3, 
                                       t(Vhigh_fold(exon = cortex01_02_isoform_valid$V1, up = cortex01_02_isoform_valid$V2 - cortex01_02_isoform_valid$V3, sample1 = "cortex01", sample2 = "cortex02")))
colnames(cortex01_02_isoform_valid)[(ncol(cortex01_02_isoform_valid) - 2):ncol(cortex01_02_isoform_valid)] <- c("cortex01junction", "cortex02junction", "junction_logfold")
(junction_valid_individual["cortex01_02", "No.isoform.exons"] <- nrow(cortex01_02_isoform_valid))
(junction_valid_individual["cortex01_02", "No.isoform.genes"] <- length(unique(cortex01_02_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
cortex01_02_isoform_valid <- na.omit(cortex01_02_isoform_valid)
(junction_valid_individual["cortex01_02", "No.exons.with.junction.cov"] <- nrow(cortex01_02_isoform_valid))
(junction_valid_individual["cortex01_02", "No.genes.with.junction.cov"] <- length(unique(cortex01_02_isoform_valid$gene)))
# No. of exons/genes with junction support
cortex01_02_isoform_valid <- cortex01_02_isoform_valid[(cortex01_02_isoform_valid$cortex01junction >= jcut & cortex01_02_isoform_valid$cortex02junction < jcut) | (cortex01_02_isoform_valid$cortex01junction < jcut & cortex01_02_isoform_valid$cortex02junction >= jcut),]
cortex01_02_isoform_valid$junction_logfold <- NULL
cortex01_02_isoform_valid <- cortex01_02_isoform_valid[order((cortex01_02_isoform_valid$cortex01gene + cortex01_02_isoform_valid$cortex02gene), decreasing = T), ]
cortex01_02_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(cortex01_02_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid_individual["cortex01_02", "No.exons.with.junction.support"] <- nrow(cortex01_02_isoform_valid))
(junction_valid_individual["cortex01_02", "No.genes.with.junction.support"] <- length(unique(cortex01_02_isoform_valid$gene)))
write.table(cortex01_02_isoform_valid, file = "cortex01_02_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
cortex01_02_isoform_valid_gene <- cortex01_02_isoform_valid[!duplicated(cortex01_02_isoform_valid$gene), ]
write.table(cortex01_02_isoform_valid_gene, file = "cortex_01_02_isoform_valid_gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# ge01_02
ge01_02_isoform_all <- read.delim("~/FetalBrain/RNAseq/isoform/ge01_02_isoform.txt", as.is =T)
ge01_02_isoform_valid <- ge01_02_isoform_all
(junction_valid_individual["ge01_02", "No.isoform.exons"] <- nrow(ge01_02_isoform_valid))
(junction_valid_individual["ge01_02", "No.isoform.genes"] <- length(unique(ge01_02_isoform_valid$id)))
sum(is.element(ge01_02_isoform_valid$V1, junction_exon$exonID))
ge01_02_isoform_valid <- data.frame(exon = ge01_02_isoform_valid$V1, gene = ge01_02_isoform_valid$id, DEfine.p = ge01_02_isoform_valid$V5, 
                                       ge01gene = ge01_02_isoform_valid$ave2, ge02gene = ge01_02_isoform_valid$ave3, 
                                       ge01exon = ge01_02_isoform_valid$V2, ge02exon = ge01_02_isoform_valid$V3, 
                                       t(Vhigh_fold(exon = ge01_02_isoform_valid$V1, up = ge01_02_isoform_valid$V2 - ge01_02_isoform_valid$V3, sample1 = "ge01", sample2 = "ge02")))
colnames(ge01_02_isoform_valid)[(ncol(ge01_02_isoform_valid) - 2):ncol(ge01_02_isoform_valid)] <- c("ge01junction", "ge02junction", "junction_logfold")
(junction_valid_individual["ge01_02", "No.isoform.exons"] <- nrow(ge01_02_isoform_valid))
(junction_valid_individual["ge01_02", "No.isoform.genes"] <- length(unique(ge01_02_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
ge01_02_isoform_valid <- na.omit(ge01_02_isoform_valid)
(junction_valid_individual["ge01_02", "No.exons.with.junction.cov"] <- nrow(ge01_02_isoform_valid))
(junction_valid_individual["ge01_02", "No.genes.with.junction.cov"] <- length(unique(ge01_02_isoform_valid$gene)))
# No. of exons/genes with junction support
ge01_02_isoform_valid <- ge01_02_isoform_valid[(ge01_02_isoform_valid$ge01junction >= jcut & ge01_02_isoform_valid$ge02junction < jcut) | (ge01_02_isoform_valid$ge01junction < jcut & ge01_02_isoform_valid$ge02junction >= jcut),]
ge01_02_isoform_valid$junction_logfold <- NULL
ge01_02_isoform_valid <- ge01_02_isoform_valid[order((ge01_02_isoform_valid$ge01gene + ge01_02_isoform_valid$ge02gene), decreasing = T), ]
ge01_02_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(ge01_02_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid_individual["ge01_02", "No.exons.with.junction.support"] <- nrow(ge01_02_isoform_valid))
(junction_valid_individual["ge01_02", "No.genes.with.junction.support"] <- length(unique(ge01_02_isoform_valid$gene)))
write.table(ge01_02_isoform_valid, file = "ge_01_02_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
ge01_02_isoform_valid_gene <- ge01_02_isoform_valid[!duplicated(ge01_02_isoform_valid$gene), ]
write.table(ge01_02_isoform_valid_gene, file = "ge_01_02_isoform_valid_gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# cortex03_04
cortex03_04_isoform_all <- read.delim("~/FetalBrain/RNAseq/isoform/cortex03_04_isoform.txt", as.is =T)
cortex03_04_isoform_valid <- cortex03_04_isoform_all
(junction_valid_individual["cortex03_04", "No.isoform.exons"] <- nrow(cortex03_04_isoform_valid))
(junction_valid_individual["cortex03_04", "No.isoform.genes"] <- length(unique(cortex03_04_isoform_valid$id)))
sum(is.element(cortex03_04_isoform_valid$V1, junction_exon$exonID))
cortex03_04_isoform_valid <- data.frame(exon = cortex03_04_isoform_valid$V1, gene = cortex03_04_isoform_valid$id, DEfine.p = cortex03_04_isoform_valid$V5, 
                                        cortex03gene = cortex03_04_isoform_valid$ave2, cortex04gene = cortex03_04_isoform_valid$ave3, 
                                        cortex03exon = cortex03_04_isoform_valid$V2, cortex04exon = cortex03_04_isoform_valid$V3, 
                                        t(Vhigh_fold(exon = cortex03_04_isoform_valid$V1, up = cortex03_04_isoform_valid$V2 - cortex03_04_isoform_valid$V3, sample1 = "cortex03", sample2 = "cortex04")))
colnames(cortex03_04_isoform_valid)[(ncol(cortex03_04_isoform_valid) - 2):ncol(cortex03_04_isoform_valid)] <- c("cortex03junction", "cortex04junction", "junction_logfold")
(junction_valid_individual["cortex03_04", "No.isoform.exons"] <- nrow(cortex03_04_isoform_valid))
(junction_valid_individual["cortex03_04", "No.isoform.genes"] <- length(unique(cortex03_04_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
cortex03_04_isoform_valid <- na.omit(cortex03_04_isoform_valid)
(junction_valid_individual["cortex03_04", "No.exons.with.junction.cov"] <- nrow(cortex03_04_isoform_valid))
(junction_valid_individual["cortex03_04", "No.genes.with.junction.cov"] <- length(unique(cortex03_04_isoform_valid$gene)))
# No. of exons/genes with junction support
cortex03_04_isoform_valid <- cortex03_04_isoform_valid[(cortex03_04_isoform_valid$cortex03junction >= jcut & cortex03_04_isoform_valid$cortex04junction < jcut) | (cortex03_04_isoform_valid$cortex03junction < jcut & cortex03_04_isoform_valid$cortex04junction >= jcut),]
cortex03_04_isoform_valid$junction_logfold <- NULL
cortex03_04_isoform_valid <- cortex03_04_isoform_valid[order((cortex03_04_isoform_valid$cortex03gene + cortex03_04_isoform_valid$cortex04gene), decreasing = T), ]
cortex03_04_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(cortex03_04_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid_individual["cortex03_04", "No.exons.with.junction.support"] <- nrow(cortex03_04_isoform_valid))
(junction_valid_individual["cortex03_04", "No.genes.with.junction.support"] <- length(unique(cortex03_04_isoform_valid$gene)))
write.table(cortex03_04_isoform_valid, file = "cortex_03_04_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
cortex03_04_isoform_valid_gene <- cortex03_04_isoform_valid[!duplicated(cortex03_04_isoform_valid$gene), ]
write.table(cortex03_04_isoform_valid_gene, file = "cortex_03_04_isoform_valid_gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# ge03_04
ge03_04_isoform_all <- read.delim("~/FetalBrain/RNAseq/isoform/ge03_04_isoform.txt", as.is =T)
ge03_04_isoform_valid <- ge03_04_isoform_all
(junction_valid_individual["ge03_04", "No.isoform.exons"] <- nrow(ge03_04_isoform_valid))
(junction_valid_individual["ge03_04", "No.isoform.genes"] <- length(unique(ge03_04_isoform_valid$id)))
sum(is.element(ge03_04_isoform_valid$V1, junction_exon$exonID))
ge03_04_isoform_valid <- data.frame(exon = ge03_04_isoform_valid$V1, gene = ge03_04_isoform_valid$id, DEfine.p = ge03_04_isoform_valid$V5, 
                                    ge03gene = ge03_04_isoform_valid$ave2, ge04gene = ge03_04_isoform_valid$ave3, 
                                    ge03exon = ge03_04_isoform_valid$V2, ge04exon = ge03_04_isoform_valid$V3, 
                                    t(Vhigh_fold(exon = ge03_04_isoform_valid$V1, up = ge03_04_isoform_valid$V2 - ge03_04_isoform_valid$V3, sample1 = "ge03", sample2 = "ge04")))
colnames(ge03_04_isoform_valid)[(ncol(ge03_04_isoform_valid) - 2):ncol(ge03_04_isoform_valid)] <- c("ge03junction", "ge04junction", "junction_logfold")
(junction_valid_individual["ge03_04", "No.isoform.exons"] <- nrow(ge03_04_isoform_valid))
(junction_valid_individual["ge03_04", "No.isoform.genes"] <- length(unique(ge03_04_isoform_valid$gene)))
# No. of exons/genes with enough junction coverage 
ge03_04_isoform_valid <- na.omit(ge03_04_isoform_valid)
(junction_valid_individual["ge03_04", "No.exons.with.junction.cov"] <- nrow(ge03_04_isoform_valid))
(junction_valid_individual["ge03_04", "No.genes.with.junction.cov"] <- length(unique(ge03_04_isoform_valid$gene)))
# No. of exons/genes with junction support
ge03_04_isoform_valid <- ge03_04_isoform_valid[(ge03_04_isoform_valid$ge03junction >= jcut & ge03_04_isoform_valid$ge04junction < jcut) | (ge03_04_isoform_valid$ge03junction < jcut & ge03_04_isoform_valid$ge04junction >= jcut),]
ge03_04_isoform_valid$junction_logfold <- NULL
ge03_04_isoform_valid <- ge03_04_isoform_valid[order((ge03_04_isoform_valid$ge03gene + ge03_04_isoform_valid$ge04gene), decreasing = T), ]
ge03_04_isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(ge03_04_isoform_valid$gene), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
(junction_valid_individual["ge03_04", "No.exons.with.junction.support"] <- nrow(ge03_04_isoform_valid))
(junction_valid_individual["ge03_04", "No.genes.with.junction.support"] <- length(unique(ge03_04_isoform_valid$gene)))
write.table(ge03_04_isoform_valid, file = "ge_03_04_isoform_valid.txt", sep = "\t", quote = F, row.names = F, col.names = T)
ge03_04_isoform_valid_gene <- ge03_04_isoform_valid[!duplicated(ge03_04_isoform_valid$gene), ]
write.table(ge03_04_isoform_valid_gene, file = "ge_03_04_isoform_valid_gene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

(junction_valid_individual)
rm(junction_exon)
save.image(file = "junction_valid.Rdata")
################################################################################################################################################
setwd("~/快盘/FetalBrain/RNAseq/junction/")
load("junction_valid.Rdata")

# DE gene list  
DE_HuFNSC01up <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/UP.Cortex-HuFNSC01_GE-HuFNSC01.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
DE_HuFNSC01dn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/DN.Cortex-HuFNSC01_GE-HuFNSC01.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
HuFNSC01_DE <- unique(c(DE_HuFNSC01up$V1, DE_HuFNSC01dn$V1))
rm(DE_HuFNSC01up, DE_HuFNSC01dn)
DE_HuFNSC02up <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/UP.Cortex-HuFNSC02_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
DE_HuFNSC02dn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/DN.Cortex-HuFNSC02_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
HuFNSC02_DE <- unique(c(DE_HuFNSC02up$V1, DE_HuFNSC02dn$V1))
rm(DE_HuFNSC02up, DE_HuFNSC02dn)
DE_HuFNSC03up <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/UP.Cortex-HuFNSC03_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
DE_HuFNSC03dn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/DN.Cortex-HuFNSC03_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
HuFNSC03_DE <- unique(c(DE_HuFNSC03up$V1, DE_HuFNSC03dn$V1))
rm(DE_HuFNSC03up, DE_HuFNSC03dn)
DE_HuFNSC04up <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/UP.Cortex-HuFNSC04_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
DE_HuFNSC04dn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/DN.Cortex-HuFNSC04_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
HuFNSC04_DE <- unique(c(DE_HuFNSC04up$V1, DE_HuFNSC04dn$V1))
rm(DE_HuFNSC04up, DE_HuFNSC04dn)

# No. of exons for DE genes / isoform genes: use collapsed exons instead of Ensembl exons
exon <- read.delim("~/hg19/hg19v65_exons_collapsed.txt", head = F, as.is = T)
Nexon <- data.frame(id = levels(as.factor(exon$V2)), Nexon = sapply(levels(as.factor(exon$V2)), function(x) sum(exon$V2 == x)))
HuFNSC01_DE <- data.frame(id = HuFNSC01_DE, Nexon = Nexon[HuFNSC01_DE, ]$Nexon)
HuFNSC01_isoform_all_gene <- HuFNSC01_isoform_all[!duplicated(HuFNSC01_isoform_all$id), ]
HuFNSC01_isoform_all_gene$Nexon <- Nexon[HuFNSC01_isoform_all_gene$id, "Nexon"]
HuFNSC01_isoform_valid_gene$Nexon <- Nexon[as.character(HuFNSC01_isoform_valid_gene$gene), "Nexon"] 

rm(exon)

