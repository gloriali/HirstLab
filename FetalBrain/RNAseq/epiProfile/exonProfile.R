setwd("~/快盘/FetalBrain/RNAseq/epiProfile/")
source("~/HirstLab/Pipeline/epiProfile.R")
library(ggplot2)
load("exonProfile_cell.Rdata")
load("exonProfile_twins.Rdata")

#################################################################################################
# cortex vs GE: HuFNSC02  
# grouping exons  
dirRPKM <- "~/FetalBrain/RNAseq/rpkm/"
dirIsoform <- "../isoform/"
dirIn <- "~/FetalBrain/RNAseq/epiProfile/"
lib1='A03475'; cell1='Cortex'; donor1='HuFNSC02';
lib2='A03476'; cell2='GE'; donor2='HuFNSC02';
cortex02 <- read.delim(paste0(dirRPKM, lib1, ".G.exn.A.rpkm"), head = F, as.is = T)
cortex02$coord <- gsub("<[-]*1", "", cortex02$V1)
cortex02$coord  <- gsub("chr", "", cortex02$coord)
cortex02$coord  <- gsub(":", "_", cortex02$coord)
cortex02$coord  <- gsub("-", "_", cortex02$coord)
cortex02$coord  <- paste0(cortex02$V2, "_", cortex02$coord)
rownames(cortex02) <- cortex02$coord
ge02 <- read.delim(paste0(dirRPKM, lib2, ".G.exn.A.rpkm"), head = F, as.is = T)
ge02$coord <- gsub("<[-]*1", "", ge02$V1)
ge02$coord  <- gsub("chr", "", ge02$coord)
ge02$coord  <- gsub(":", "_", ge02$coord)
ge02$coord  <- gsub("-", "_", ge02$coord)
ge02$coord  <- paste0(ge02$V2, "_", ge02$coord)
rownames(ge02) <- ge02$coord
cortex02gene <- read.delim(paste0(dirRPKM, lib1, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(cortex02gene) <- cortex02gene$V1
ge02gene <- read.delim(paste0(dirRPKM, lib2, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(ge02gene) <- ge02gene$V1
exons_geneRPKM <- data.frame(exon = cortex02$coord)
exons_geneRPKM$gene <- gsub("_[0-9XY]+_[0-9_]+", "", exons_geneRPKM$exon)
exons_geneRPKM$cortex02gene <- cortex02gene[exons_geneRPKM$gene, "V3"]
exons_geneRPKM$ge02gene <- ge02gene[exons_geneRPKM$gene, "V3"]
exons_geneRPKM$geneRPKM <- (exons_geneRPKM$cortex02gene + exons_geneRPKM$ge02gene)/2
exons_geneRPKM$cortex02exon <- cortex02[as.character(exons_geneRPKM$exon), "V4"]
exons_geneRPKM$ge02exon <- ge02[as.character(exons_geneRPKM$exon), "V4"]
exons_geneRPKM <- na.omit(exons_geneRPKM)
exons_geneRPKM$exon <- as.character(exons_geneRPKM$exon)
isoform <- read.delim(paste0(dirIsoform, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform.txt"), head = T, as.is = T)
isoform$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform$V1)
isoform$coord  <- gsub("chr", "", isoform$coord)
isoform$coord  <- gsub(":", "_", isoform$coord)
isoform$coord  <- gsub("-", "_", isoform$coord)
isoform$coord  <- paste0(isoform$id, "_", isoform$coord)
cortex_specific <- isoform[isoform$V2 < isoform$V3, "coord"]  # isoform exons
ge_specific <- isoform[isoform$V2 > isoform$V3, "coord"]  # isoform exons
both <- setdiff(exons_geneRPKM[exons_geneRPKM$cortex02exon > 0.1 & exons_geneRPKM$ge02exon > 0.1, "exon"], c(cortex_specific, ge_specific))
neither <- setdiff(cortex02$coord, c(cortex_specific, ge_specific, both))
cortex_exon_0 <- exons_geneRPKM[exons_geneRPKM$cortex02exon == 0 & exons_geneRPKM$cortex02gene > 0.1, "exon"]
ge_exon_0 <- exons_geneRPKM[exons_geneRPKM$ge02exon == 0 & exons_geneRPKM$ge02gene > 0.1, "exon"]
cortex_gene_0 <- exons_geneRPKM[exons_geneRPKM$cortex02gene == 0, "exon"]
ge_gene_0 <- exons_geneRPKM[exons_geneRPKM$ge02gene == 0, "exon"]
exons_1 <- as.character(exons_geneRPKM[exons_geneRPKM$geneRPKM <= 1, "exon"])
exons_1_10 <- as.character(exons_geneRPKM[exons_geneRPKM$geneRPKM > 1 & exons_geneRPKM$geneRPKM <= 10, "exon"])
exons_10_100 <- as.character(exons_geneRPKM[exons_geneRPKM$geneRPKM > 10 & exons_geneRPKM$geneRPKM <= 100, "exon"])
exons_100 <- as.character(exons_geneRPKM[exons_geneRPKM$geneRPKM > 100, "exon"])

WGBS <- epiProfile(mark = "WGBS", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both, neither = neither, cell1_specific = cortex_specific, cell2_specific = ge_specific, CpG = T)
WGBS_profile <- WGBS$profile
WGBS_figure <- WGBS$figure
library(grid) 
gt <- ggplot_gtable(ggplot_build(WGBS_figure)) 
gt$heights[[4]] <- unit(3, "null") 
grid.newpage()
grid.draw(gt) 
grid.text("Average DNA methylation", x = unit(0.015, "npc"), y = unit(0.65, "npc"), rot = 90)
grid.text("No. of CpGs", x = unit(0.015, "npc"), y = unit(0.18, "npc"), rot = 90)

MeDIP <- epiProfile(mark = "MeDIP", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both, neither = neither, cell1_specific = cortex_specific, cell2_specific = ge_specific)
MeDIP_profile <- MeDIP$profile
(MeDIP_figure <- MeDIP$figure)

H3K4me1 <- epiProfile(mark = "H3K4me1", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both, neither = neither, cell1_specific = cortex_specific, cell2_specific = ge_specific)
H3K4me1_profile <- H3K4me1$profile
(H3K4me1_figure <- H3K4me1$figure)

H3K4me3 <- epiProfile(mark = "H3K4me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both, neither = neither, cell1_specific = cortex_specific, cell2_specific = ge_specific)
H3K4me3_profile <- H3K4me3$profile
(H3K4me3_figure <- H3K4me3$figure)

H3K9me3 <- epiProfile(mark = "H3K9me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both, neither = neither, cell1_specific = cortex_specific, cell2_specific = ge_specific)
H3K9me3_profile <- H3K9me3$profile
(H3K9me3_figure <- H3K9me3$figure)

H3K27me3 <- epiProfile(mark = "H3K27me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both, neither = neither, cell1_specific = cortex_specific, cell2_specific = ge_specific)
H3K27me3_profile <- H3K27me3$profile
(H3K27me3_figure <- H3K27me3$figure)

H3K36me3 <- epiProfile(mark = "H3K36me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both, neither = neither, cell1_specific = cortex_specific, cell2_specific = ge_specific, 
                       geneRPKM = list("RPKM < 1" = exons_1, "RPKM 1-10" = exons_1_10))
H3K36me3_data <- H3K36me3$data
H3K36me3_profile <- H3K36me3$profile
(H3K36me3_figure <- H3K36me3$figure)

save(exons_geneRPKM, isoform, cortex_specific, ge_specific, both, neither, 
     WGBS_profile, WGBS_figure, MeDIP_profile, MeDIP_figure, H3K4me1_profile, H3K4me1_figure, H3K4me3_profile, H3K4me3_figure, H3K9me3_profile, H3K9me3_figure, H3K27me3_profile, H3K27me3_figure, H3K36me3_data, H3K36me3_profile, H3K36me3_figure, 
     file = "exonProfile_cell.Rdata")
#################################################################################################

