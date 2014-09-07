setwd("~/快盘/FetalBrain/RNAseq/epiProfile/")
source("~/HirstLab/Pipeline/epiProfile.R")
library(ggplot2)
library(grid) 
# load("exonProfile_cell.Rdata")
# load("exonProfile_twins.Rdata")

#################################################################################################
# cortex vs GE
dirRPKM <- "~/FetalBrain/RNAseq/rpkm/"
dirisoform <- "../isoform/"
dirIn <- "~/FetalBrain/RNAseq/epiProfile/"

# HuFNSC02  
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
exons_geneRPKM_HuFNSC02 <- data.frame(exon = cortex02$coord)
exons_geneRPKM_HuFNSC02$exon <- as.character(exons_geneRPKM_HuFNSC02$exon)
exons_geneRPKM_HuFNSC02$gene <- gsub("_[0-9XY]+_[0-9_]+", "", exons_geneRPKM_HuFNSC02$exon)
exons_geneRPKM_HuFNSC02$cortex02gene <- cortex02gene[exons_geneRPKM_HuFNSC02$gene, "V3"]
exons_geneRPKM_HuFNSC02$ge02gene <- ge02gene[exons_geneRPKM_HuFNSC02$gene, "V3"]
exons_geneRPKM_HuFNSC02$geneRPKM <- (exons_geneRPKM_HuFNSC02$cortex02gene + exons_geneRPKM_HuFNSC02$ge02gene)/2
exons_geneRPKM_HuFNSC02$cortex02exon <- cortex02[as.character(exons_geneRPKM_HuFNSC02$exon), "V4"]
exons_geneRPKM_HuFNSC02$ge02exon <- ge02[as.character(exons_geneRPKM_HuFNSC02$exon), "V4"]
exons_geneRPKM_HuFNSC02 <- na.omit(exons_geneRPKM_HuFNSC02)
isoform_HuFNSC02 <- read.delim(paste0(dirisoform, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform.txt"), head = T, as.is = T)
isoform_HuFNSC02$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform_HuFNSC02$V1)
isoform_HuFNSC02$coord  <- gsub("chr", "", isoform_HuFNSC02$coord)
isoform_HuFNSC02$coord  <- gsub(":", "_", isoform_HuFNSC02$coord)
isoform_HuFNSC02$coord  <- gsub("-", "_", isoform_HuFNSC02$coord)
isoform_HuFNSC02$coord  <- paste0(isoform_HuFNSC02$id, "_", isoform_HuFNSC02$coord)
cortex_specific_HuFNSC02 <- isoform_HuFNSC02[isoform_HuFNSC02$V2 < isoform_HuFNSC02$V3, "coord"]  # isoform_HuFNSC02 exons
ge_specific_HuFNSC02 <- isoform_HuFNSC02[isoform_HuFNSC02$V2 > isoform_HuFNSC02$V3, "coord"]  # isoform_HuFNSC02 exons
both_HuFNSC02 <- setdiff(exons_geneRPKM_HuFNSC02[exons_geneRPKM_HuFNSC02$cortex02exon > 0.1 & exons_geneRPKM_HuFNSC02$ge02exon > 0.1, "exon"], c(cortex_specific_HuFNSC02, ge_specific_HuFNSC02))
neither_HuFNSC02 <- setdiff(cortex02$coord, c(cortex_specific_HuFNSC02, ge_specific_HuFNSC02, both_HuFNSC02))
cortex_exon_0 <- exons_geneRPKM_HuFNSC02[exons_geneRPKM_HuFNSC02$cortex02exon == 0 & exons_geneRPKM_HuFNSC02$cortex02gene > 0.1, "exon"]
ge_exon_0 <- exons_geneRPKM_HuFNSC02[exons_geneRPKM_HuFNSC02$ge02exon == 0 & exons_geneRPKM_HuFNSC02$ge02gene > 0.1, "exon"]
cortex_gene_0 <- exons_geneRPKM_HuFNSC02[exons_geneRPKM_HuFNSC02$cortex02gene == 0, "exon"]
ge_gene_0 <- exons_geneRPKM_HuFNSC02[exons_geneRPKM_HuFNSC02$ge02gene == 0, "exon"]
# exons_1 <- as.character(exons_geneRPKM_HuFNSC02[exons_geneRPKM_HuFNSC02$geneRPKM <= 1, "exon"])
# exons_1_10 <- as.character(exons_geneRPKM_HuFNSC02[exons_geneRPKM_HuFNSC02$geneRPKM > 1 & exons_geneRPKM_HuFNSC02$geneRPKM <= 10, "exon"])
# exons_10_100 <- as.character(exons_geneRPKM_HuFNSC02[exons_geneRPKM_HuFNSC02$geneRPKM > 10 & exons_geneRPKM_HuFNSC02$geneRPKM <= 100, "exon"])
# exons_100 <- as.character(exons_geneRPKM_HuFNSC02[exons_geneRPKM_HuFNSC02$geneRPKM > 100, "exon"])

WGBS_HuFNSC02 <- epiProfile(mark = "WGBS", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_HuFNSC02, neither = neither_HuFNSC02, cell1_specific = cortex_specific_HuFNSC02, cell2_specific = ge_specific_HuFNSC02, CpG = T)
WGBS_profile_HuFNSC02 <- WGBS_HuFNSC02$profile
WGBS_figure_HuFNSC02 <- WGBS_HuFNSC02$figure
gt <- ggplot_gtable(ggplot_build(WGBS_figure_HuFNSC02)) 
gt$heights[[4]] <- unit(3.5, "null") 
grid.newpage()
grid.draw(gt) 
grid.text("Average DNA methylation", x = unit(0.015, "npc"), y = unit(0.65, "npc"), rot = 90)
grid.text("CpG density", x = unit(0.015, "npc"), y = unit(0.15, "npc"), rot = 90)

MeDIP_HuFNSC02 <- epiProfile(mark = "MeDIP", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_HuFNSC02, neither = neither_HuFNSC02, cell1_specific = cortex_specific_HuFNSC02, cell2_specific = ge_specific_HuFNSC02)
MeDIP_profile_HuFNSC02 <- MeDIP_HuFNSC02$profile
(MeDIP_figure_HuFNSC02 <- MeDIP_HuFNSC02$figure)

H3K36me3_HuFNSC02 <- epiProfile(mark = "H3K36me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_HuFNSC02, neither = neither_HuFNSC02, cell1_specific = cortex_specific_HuFNSC02, cell2_specific = ge_specific_HuFNSC02)
H3K36me3_data_HuFNSC02 <- H3K36me3_HuFNSC02$data
H3K36me3_profile_HuFNSC02 <- H3K36me3_HuFNSC02$profile
(H3K36me3_figure_HuFNSC02 <- H3K36me3_HuFNSC02$figure)

H3K4me1_HuFNSC02 <- epiProfile(mark = "H3K4me1", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_HuFNSC02, neither = neither_HuFNSC02, cell1_specific = cortex_specific_HuFNSC02, cell2_specific = ge_specific_HuFNSC02)
H3K4me1_profile_HuFNSC02 <- H3K4me1_HuFNSC02$profile
(H3K4me1_figure_HuFNSC02 <- H3K4me1_HuFNSC02$figure)

H3K4me3_HuFNSC02 <- epiProfile(mark = "H3K4me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_HuFNSC02, neither = neither_HuFNSC02, cell1_specific = cortex_specific_HuFNSC02, cell2_specific = ge_specific_HuFNSC02)
H3K4me3_profile_HuFNSC02 <- H3K4me3_HuFNSC02$profile
(H3K4me3_figure_HuFNSC02 <- H3K4me3_HuFNSC02$figure)

H3K9me3_HuFNSC02 <- epiProfile(mark = "H3K9me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_HuFNSC02, neither = neither_HuFNSC02, cell1_specific = cortex_specific_HuFNSC02, cell2_specific = ge_specific_HuFNSC02)
H3K9me3_profile_HuFNSC02 <- H3K9me3_HuFNSC02$profile
(H3K9me3_figure_HuFNSC02 <- H3K9me3_HuFNSC02$figure)

H3K27me3_HuFNSC02 <- epiProfile(mark = "H3K27me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_HuFNSC02, neither = neither_HuFNSC02, cell1_specific = cortex_specific_HuFNSC02, cell2_specific = ge_specific_HuFNSC02)
H3K27me3_profile_HuFNSC02 <- H3K27me3_HuFNSC02$profile
(H3K27me3_figure_HuFNSC02 <- H3K27me3_HuFNSC02$figure)

# HuFNSC04  
lib1='A15298'; cell1='Cortex'; donor1='HuFNSC04';
lib2='A15299'; cell2='GE'; donor2='HuFNSC04';
cortex04 <- read.delim(paste0(dirRPKM, lib1, ".G.exn.A.rpkm"), head = F, as.is = T)
cortex04$coord <- gsub("<[-]*1", "", cortex04$V1)
cortex04$coord  <- gsub("chr", "", cortex04$coord)
cortex04$coord  <- gsub(":", "_", cortex04$coord)
cortex04$coord  <- gsub("-", "_", cortex04$coord)
cortex04$coord  <- paste0(cortex04$V2, "_", cortex04$coord)
rownames(cortex04) <- cortex04$coord
ge04 <- read.delim(paste0(dirRPKM, lib2, ".G.exn.A.rpkm"), head = F, as.is = T)
ge04$coord <- gsub("<[-]*1", "", ge04$V1)
ge04$coord  <- gsub("chr", "", ge04$coord)
ge04$coord  <- gsub(":", "_", ge04$coord)
ge04$coord  <- gsub("-", "_", ge04$coord)
ge04$coord  <- paste0(ge04$V2, "_", ge04$coord)
rownames(ge04) <- ge04$coord
cortex04gene <- read.delim(paste0(dirRPKM, lib1, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(cortex04gene) <- cortex04gene$V1
ge04gene <- read.delim(paste0(dirRPKM, lib2, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(ge04gene) <- ge04gene$V1
exons_geneRPKM_HuFNSC04 <- data.frame(exon = cortex04$coord)
exons_geneRPKM_HuFNSC04$exon <- as.character(exons_geneRPKM_HuFNSC04$exon)
exons_geneRPKM_HuFNSC04$gene <- gsub("_[0-9XY]+_[0-9_]+", "", exons_geneRPKM_HuFNSC04$exon)
exons_geneRPKM_HuFNSC04$cortex04gene <- cortex04gene[exons_geneRPKM_HuFNSC04$gene, "V3"]
exons_geneRPKM_HuFNSC04$ge04gene <- ge04gene[exons_geneRPKM_HuFNSC04$gene, "V3"]
exons_geneRPKM_HuFNSC04$geneRPKM <- (exons_geneRPKM_HuFNSC04$cortex04gene + exons_geneRPKM_HuFNSC04$ge04gene)/2
exons_geneRPKM_HuFNSC04$cortex04exon <- cortex04[as.character(exons_geneRPKM_HuFNSC04$exon), "V4"]
exons_geneRPKM_HuFNSC04$ge04exon <- ge04[as.character(exons_geneRPKM_HuFNSC04$exon), "V4"]
exons_geneRPKM_HuFNSC04 <- na.omit(exons_geneRPKM_HuFNSC04)
isoform_HuFNSC04 <- read.delim(paste0(dirisoform, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform.txt"), head = T, as.is = T)
isoform_HuFNSC04$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform_HuFNSC04$V1)
isoform_HuFNSC04$coord  <- gsub("chr", "", isoform_HuFNSC04$coord)
isoform_HuFNSC04$coord  <- gsub(":", "_", isoform_HuFNSC04$coord)
isoform_HuFNSC04$coord  <- gsub("-", "_", isoform_HuFNSC04$coord)
isoform_HuFNSC04$coord  <- paste0(isoform_HuFNSC04$id, "_", isoform_HuFNSC04$coord)
cortex_specific_HuFNSC04 <- isoform_HuFNSC04[isoform_HuFNSC04$V2 < isoform_HuFNSC04$V3, "coord"]  # isoform_HuFNSC04 exons
ge_specific_HuFNSC04 <- isoform_HuFNSC04[isoform_HuFNSC04$V2 > isoform_HuFNSC04$V3, "coord"]  # isoform_HuFNSC04 exons
both_HuFNSC04 <- setdiff(exons_geneRPKM_HuFNSC04[exons_geneRPKM_HuFNSC04$cortex04exon > 0.1 & exons_geneRPKM_HuFNSC04$ge04exon > 0.1, "exon"], c(cortex_specific_HuFNSC04, ge_specific_HuFNSC04))
neither_HuFNSC04 <- setdiff(cortex04$coord, c(cortex_specific_HuFNSC04, ge_specific_HuFNSC04, both_HuFNSC04))
cortex_exon_0 <- exons_geneRPKM_HuFNSC04[exons_geneRPKM_HuFNSC04$cortex04exon == 0 & exons_geneRPKM_HuFNSC04$cortex04gene > 0.1, "exon"]
ge_exon_0 <- exons_geneRPKM_HuFNSC04[exons_geneRPKM_HuFNSC04$ge04exon == 0 & exons_geneRPKM_HuFNSC04$ge04gene > 0.1, "exon"]
cortex_gene_0 <- exons_geneRPKM_HuFNSC04[exons_geneRPKM_HuFNSC04$cortex04gene == 0, "exon"]
ge_gene_0 <- exons_geneRPKM_HuFNSC04[exons_geneRPKM_HuFNSC04$ge04gene == 0, "exon"]
# exons_1 <- as.character(exons_geneRPKM_HuFNSC04[exons_geneRPKM_HuFNSC04$geneRPKM <= 1, "exon"])
# exons_1_10 <- as.character(exons_geneRPKM_HuFNSC04[exons_geneRPKM_HuFNSC04$geneRPKM > 1 & exons_geneRPKM_HuFNSC04$geneRPKM <= 10, "exon"])
# exons_10_100 <- as.character(exons_geneRPKM_HuFNSC04[exons_geneRPKM_HuFNSC04$geneRPKM > 10 & exons_geneRPKM_HuFNSC04$geneRPKM <= 100, "exon"])
# exons_100 <- as.character(exons_geneRPKM_HuFNSC04[exons_geneRPKM_HuFNSC04$geneRPKM > 100, "exon"])

WGBS_HuFNSC04 <- epiProfile(mark = "WGBS", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_HuFNSC04, neither = neither_HuFNSC04, cell1_specific = cortex_specific_HuFNSC04, cell2_specific = ge_specific_HuFNSC04, CpG = T)
WGBS_profile_HuFNSC04 <- WGBS_HuFNSC04$profile
WGBS_figure_HuFNSC04 <- WGBS_HuFNSC04$figure
gt <- ggplot_gtable(ggplot_build(WGBS_figure_HuFNSC04)) 
gt$heights[[4]] <- unit(3.5, "null") 
grid.newpage()
grid.draw(gt) 
grid.text("Average DNA methylation", x = unit(0.015, "npc"), y = unit(0.65, "npc"), rot = 90)
grid.text("CpG density", x = unit(0.015, "npc"), y = unit(0.15, "npc"), rot = 90)

# HuFNSC01  
lib1='A03473'; cell1='Cortex'; donor1='HuFNSC01';
lib2='A03474'; cell2='GE'; donor2='HuFNSC01';
cortex01 <- read.delim(paste0(dirRPKM, lib1, ".G.exn.A.rpkm"), head = F, as.is = T)
cortex01$coord <- gsub("<[-]*1", "", cortex01$V1)
cortex01$coord  <- gsub("chr", "", cortex01$coord)
cortex01$coord  <- gsub(":", "_", cortex01$coord)
cortex01$coord  <- gsub("-", "_", cortex01$coord)
cortex01$coord  <- paste0(cortex01$V2, "_", cortex01$coord)
rownames(cortex01) <- cortex01$coord
ge01 <- read.delim(paste0(dirRPKM, lib2, ".G.exn.A.rpkm"), head = F, as.is = T)
ge01$coord <- gsub("<[-]*1", "", ge01$V1)
ge01$coord  <- gsub("chr", "", ge01$coord)
ge01$coord  <- gsub(":", "_", ge01$coord)
ge01$coord  <- gsub("-", "_", ge01$coord)
ge01$coord  <- paste0(ge01$V2, "_", ge01$coord)
rownames(ge01) <- ge01$coord
cortex01gene <- read.delim(paste0(dirRPKM, lib1, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(cortex01gene) <- cortex01gene$V1
ge01gene <- read.delim(paste0(dirRPKM, lib2, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(ge01gene) <- ge01gene$V1
exons_geneRPKM_HuFNSC01 <- data.frame(exon = cortex01$coord)
exons_geneRPKM_HuFNSC01$exon <- as.character(exons_geneRPKM_HuFNSC01$exon)
exons_geneRPKM_HuFNSC01$gene <- gsub("_[0-9XY]+_[0-9_]+", "", exons_geneRPKM_HuFNSC01$exon)
exons_geneRPKM_HuFNSC01$cortex01gene <- cortex01gene[exons_geneRPKM_HuFNSC01$gene, "V3"]
exons_geneRPKM_HuFNSC01$ge01gene <- ge01gene[exons_geneRPKM_HuFNSC01$gene, "V3"]
exons_geneRPKM_HuFNSC01$geneRPKM <- (exons_geneRPKM_HuFNSC01$cortex01gene + exons_geneRPKM_HuFNSC01$ge01gene)/2
exons_geneRPKM_HuFNSC01$cortex01exon <- cortex01[as.character(exons_geneRPKM_HuFNSC01$exon), "V4"]
exons_geneRPKM_HuFNSC01$ge01exon <- ge01[as.character(exons_geneRPKM_HuFNSC01$exon), "V4"]
exons_geneRPKM_HuFNSC01 <- na.omit(exons_geneRPKM_HuFNSC01)
isoform_HuFNSC01 <- read.delim(paste0(dirisoform, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform.txt"), head = T, as.is = T)
isoform_HuFNSC01$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform_HuFNSC01$V1)
isoform_HuFNSC01$coord  <- gsub("chr", "", isoform_HuFNSC01$coord)
isoform_HuFNSC01$coord  <- gsub(":", "_", isoform_HuFNSC01$coord)
isoform_HuFNSC01$coord  <- gsub("-", "_", isoform_HuFNSC01$coord)
isoform_HuFNSC01$coord  <- paste0(isoform_HuFNSC01$id, "_", isoform_HuFNSC01$coord)
cortex_specific_HuFNSC01 <- isoform_HuFNSC01[isoform_HuFNSC01$V2 < isoform_HuFNSC01$V3, "coord"]  # isoform_HuFNSC01 exons
ge_specific_HuFNSC01 <- isoform_HuFNSC01[isoform_HuFNSC01$V2 > isoform_HuFNSC01$V3, "coord"]  # isoform_HuFNSC01 exons
both_HuFNSC01 <- setdiff(exons_geneRPKM_HuFNSC01[exons_geneRPKM_HuFNSC01$cortex01exon > 0.1 & exons_geneRPKM_HuFNSC01$ge01exon > 0.1, "exon"], c(cortex_specific_HuFNSC01, ge_specific_HuFNSC01))
neither_HuFNSC01 <- setdiff(cortex01$coord, c(cortex_specific_HuFNSC01, ge_specific_HuFNSC01, both_HuFNSC01))
cortex_exon_0 <- exons_geneRPKM_HuFNSC01[exons_geneRPKM_HuFNSC01$cortex01exon == 0 & exons_geneRPKM_HuFNSC01$cortex01gene > 0.1, "exon"]
ge_exon_0 <- exons_geneRPKM_HuFNSC01[exons_geneRPKM_HuFNSC01$ge01exon == 0 & exons_geneRPKM_HuFNSC01$ge01gene > 0.1, "exon"]
cortex_gene_0 <- exons_geneRPKM_HuFNSC01[exons_geneRPKM_HuFNSC01$cortex01gene == 0, "exon"]
ge_gene_0 <- exons_geneRPKM_HuFNSC01[exons_geneRPKM_HuFNSC01$ge01gene == 0, "exon"]
# exons_1 <- as.character(exons_geneRPKM_HuFNSC01[exons_geneRPKM_HuFNSC01$geneRPKM <= 1, "exon"])
# exons_1_10 <- as.character(exons_geneRPKM_HuFNSC01[exons_geneRPKM_HuFNSC01$geneRPKM > 1 & exons_geneRPKM_HuFNSC01$geneRPKM <= 10, "exon"])
# exons_10_100 <- as.character(exons_geneRPKM_HuFNSC01[exons_geneRPKM_HuFNSC01$geneRPKM > 10 & exons_geneRPKM_HuFNSC01$geneRPKM <= 100, "exon"])
# exons_100 <- as.character(exons_geneRPKM_HuFNSC01[exons_geneRPKM_HuFNSC01$geneRPKM > 100, "exon"])

MeDIP_HuFNSC01 <- epiProfile(mark = "MeDIP", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_HuFNSC01, neither = neither_HuFNSC01, cell1_specific = cortex_specific_HuFNSC01, cell2_specific = ge_specific_HuFNSC01)
MeDIP_profile_HuFNSC01 <- MeDIP_HuFNSC01$profile
(MeDIP_figure_HuFNSC01 <- MeDIP_HuFNSC01$figure)

H3K36me3_HuFNSC01 <- epiProfile(mark = "H3K36me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_HuFNSC01, neither = neither_HuFNSC01, cell1_specific = cortex_specific_HuFNSC01, cell2_specific = ge_specific_HuFNSC01)
H3K36me3_data_HuFNSC01 <- H3K36me3_HuFNSC01$data
H3K36me3_profile_HuFNSC01 <- H3K36me3_HuFNSC01$profile
(H3K36me3_figure_HuFNSC01 <- H3K36me3_HuFNSC01$figure)

save(exons_geneRPKM_HuFNSC01, isoform_HuFNSC01, cortex_specific_HuFNSC01, ge_specific_HuFNSC01, both_HuFNSC01, neither_HuFNSC01, MeDIP_profile_HuFNSC01, MeDIP_figure_HuFNSC01, H3K36me3_data_HuFNSC01, H3K36me3_profile_HuFNSC01, H3K36me3_figure_HuFNSC01, 
     exons_geneRPKM_HuFNSC02, isoform_HuFNSC02, cortex_specific_HuFNSC02, ge_specific_HuFNSC02, both_HuFNSC02, neither_HuFNSC02, WGBS_profile_HuFNSC02, WGBS_figure_HuFNSC02, MeDIP_profile_HuFNSC02, MeDIP_figure_HuFNSC02, H3K4me1_profile_HuFNSC02, H3K4me1_figure_HuFNSC02, H3K4me3_profile_HuFNSC02, H3K4me3_figure_HuFNSC02, H3K9me3_profile_HuFNSC02, H3K9me3_figure_HuFNSC02, H3K27me3_profile_HuFNSC02, H3K27me3_figure_HuFNSC02, H3K36me3_data_HuFNSC02, H3K36me3_profile_HuFNSC02, H3K36me3_figure_HuFNSC02, 
     exons_geneRPKM_HuFNSC04, isoform_HuFNSC04, cortex_specific_HuFNSC04, ge_specific_HuFNSC04, both_HuFNSC04, neither_HuFNSC04, WGBS_profile_HuFNSC04, WGBS_figure_HuFNSC04, 
     file = "exonProfile_cell.Rdata")

#################################################################################################
#################################################################################################
# MZ twins 
dirRPKM <- "~/FetalBrain/RNAseq/rpkm/"
dirisoform <- "../isoform/"
dirIn <- "~/FetalBrain/RNAseq/epiProfile/"

# brain 
lib1='A03484'; cell1='brain'; donor1='HuFNSC01';
lib2='A07825'; cell2='brain'; donor2='HuFNSC02';
brain01 <- read.delim(paste0(dirRPKM, lib1, ".G.exn.A.rpkm"), head = F, as.is = T)
brain01$coord <- gsub("<[-]*1", "", brain01$V1)
brain01$coord  <- gsub("chr", "", brain01$coord)
brain01$coord  <- gsub(":", "_", brain01$coord)
brain01$coord  <- gsub("-", "_", brain01$coord)
brain01$coord  <- paste0(brain01$V2, "_", brain01$coord)
rownames(brain01) <- brain01$coord
brain02 <- read.delim(paste0(dirRPKM, lib2, ".G.exn.A.rpkm"), head = F, as.is = T)
brain02$coord <- gsub("<[-]*1", "", brain02$V1)
brain02$coord  <- gsub("chr", "", brain02$coord)
brain02$coord  <- gsub(":", "_", brain02$coord)
brain02$coord  <- gsub("-", "_", brain02$coord)
brain02$coord  <- paste0(brain02$V2, "_", brain02$coord)
rownames(brain02) <- brain02$coord
brain01gene <- read.delim(paste0(dirRPKM, lib1, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(brain01gene) <- brain01gene$V1
brain02gene <- read.delim(paste0(dirRPKM, lib2, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(brain02gene) <- brain02gene$V1
exons_geneRPKM_brain <- data.frame(exon = brain01$coord)
exons_geneRPKM_brain$exon <- as.character(exons_geneRPKM_brain$exon)
exons_geneRPKM_brain$gene <- gsub("_[0-9XY]+_[0-9_]+", "", exons_geneRPKM_brain$exon)
exons_geneRPKM_brain$brain01gene <- brain01gene[exons_geneRPKM_brain$gene, "V3"]
exons_geneRPKM_brain$brain02gene <- brain02gene[exons_geneRPKM_brain$gene, "V3"]
exons_geneRPKM_brain$geneRPKM <- (exons_geneRPKM_brain$brain01gene + exons_geneRPKM_brain$brain02gene)/2
exons_geneRPKM_brain$brain01exon <- brain01[as.character(exons_geneRPKM_brain$exon), "V4"]
exons_geneRPKM_brain$brain02exon <- brain02[as.character(exons_geneRPKM_brain$exon), "V4"]
exons_geneRPKM_brain <- na.omit(exons_geneRPKM_brain)
isoform_brain <- read.delim(paste0(dirisoform, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform.txt"), head = T, as.is = T)
isoform_brain$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform_brain$V1)
isoform_brain$coord  <- gsub("chr", "", isoform_brain$coord)
isoform_brain$coord  <- gsub(":", "_", isoform_brain$coord)
isoform_brain$coord  <- gsub("-", "_", isoform_brain$coord)
isoform_brain$coord  <- paste0(isoform_brain$id, "_", isoform_brain$coord)
brain01_specific_brain <- isoform_brain[isoform_brain$V2 < isoform_brain$V3, "coord"]  # isoform_brain exons
brain02_specific_brain <- isoform_brain[isoform_brain$V2 > isoform_brain$V3, "coord"]  # isoform_brain exons
both_brain <- setdiff(exons_geneRPKM_brain[exons_geneRPKM_brain$brain01exon > 0.1 & exons_geneRPKM_brain$brain02exon > 0.1, "exon"], c(brain01_specific_brain, brain02_specific_brain))
neither_brain <- setdiff(brain01$coord, c(brain01_specific_brain, brain02_specific_brain, both_brain))
brain01_exon_0 <- exons_geneRPKM_brain[exons_geneRPKM_brain$brain01exon == 0 & exons_geneRPKM_brain$brain01gene > 0.1, "exon"]
brain02_exon_0 <- exons_geneRPKM_brain[exons_geneRPKM_brain$brain02exon == 0 & exons_geneRPKM_brain$brain02gene > 0.1, "exon"]
brain01_gene_0 <- exons_geneRPKM_brain[exons_geneRPKM_brain$brain01gene == 0, "exon"]
brain02_gene_0 <- exons_geneRPKM_brain[exons_geneRPKM_brain$brain02gene == 0, "exon"]
# exons_1 <- as.character(exons_geneRPKM_brain[exons_geneRPKM_brain$geneRPKM <= 1, "exon"])
# exons_1_10 <- as.character(exons_geneRPKM_brain[exons_geneRPKM_brain$geneRPKM > 1 & exons_geneRPKM_brain$geneRPKM <= 10, "exon"])
# exons_10_100 <- as.character(exons_geneRPKM_brain[exons_geneRPKM_brain$geneRPKM > 10 & exons_geneRPKM_brain$geneRPKM <= 100, "exon"])
# exons_100 <- as.character(exons_geneRPKM_brain[exons_geneRPKM_brain$geneRPKM > 100, "exon"])

MeDIP_brain <- epiProfile(mark = "MeDIP", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_brain, neither = neither_brain, cell1_specific = brain01_specific_brain, cell2_specific = brain02_specific_brain)
MeDIP_profile_brain <- MeDIP_brain$profile
(MeDIP_figure_brain <- MeDIP_brain$figure)

H3K36me3_brain <- epiProfile(mark = "H3K36me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_brain, neither = neither_brain, cell1_specific = brain01_specific_brain, cell2_specific = brain02_specific_brain)
H3K36me3_data_brain <- H3K36me3_brain$data
H3K36me3_profile_brain <- H3K36me3_brain$profile
(H3K36me3_figure_brain <- H3K36me3_brain$figure)

# cortex 
lib1='A03473'; cell1='cortex'; donor1='HuFNSC01';
lib2='A03475'; cell2='cortex'; donor2='HuFNSC02';
cortex01 <- read.delim(paste0(dirRPKM, lib1, ".G.exn.A.rpkm"), head = F, as.is = T)
cortex01$coord <- gsub("<[-]*1", "", cortex01$V1)
cortex01$coord  <- gsub("chr", "", cortex01$coord)
cortex01$coord  <- gsub(":", "_", cortex01$coord)
cortex01$coord  <- gsub("-", "_", cortex01$coord)
cortex01$coord  <- paste0(cortex01$V2, "_", cortex01$coord)
rownames(cortex01) <- cortex01$coord
cortex02 <- read.delim(paste0(dirRPKM, lib2, ".G.exn.A.rpkm"), head = F, as.is = T)
cortex02$coord <- gsub("<[-]*1", "", cortex02$V1)
cortex02$coord  <- gsub("chr", "", cortex02$coord)
cortex02$coord  <- gsub(":", "_", cortex02$coord)
cortex02$coord  <- gsub("-", "_", cortex02$coord)
cortex02$coord  <- paste0(cortex02$V2, "_", cortex02$coord)
rownames(cortex02) <- cortex02$coord
cortex01gene <- read.delim(paste0(dirRPKM, lib1, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(cortex01gene) <- cortex01gene$V1
cortex02gene <- read.delim(paste0(dirRPKM, lib2, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(cortex02gene) <- cortex02gene$V1
exons_geneRPKM_cortex <- data.frame(exon = cortex01$coord)
exons_geneRPKM_cortex$exon <- as.character(exons_geneRPKM_cortex$exon)
exons_geneRPKM_cortex$gene <- gsub("_[0-9XY]+_[0-9_]+", "", exons_geneRPKM_cortex$exon)
exons_geneRPKM_cortex$cortex01gene <- cortex01gene[exons_geneRPKM_cortex$gene, "V3"]
exons_geneRPKM_cortex$cortex02gene <- cortex02gene[exons_geneRPKM_cortex$gene, "V3"]
exons_geneRPKM_cortex$geneRPKM <- (exons_geneRPKM_cortex$cortex01gene + exons_geneRPKM_cortex$cortex02gene)/2
exons_geneRPKM_cortex$cortex01exon <- cortex01[as.character(exons_geneRPKM_cortex$exon), "V4"]
exons_geneRPKM_cortex$cortex02exon <- cortex02[as.character(exons_geneRPKM_cortex$exon), "V4"]
exons_geneRPKM_cortex <- na.omit(exons_geneRPKM_cortex)
isoform_cortex <- read.delim(paste0(dirisoform, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform.txt"), head = T, as.is = T)
isoform_cortex$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform_cortex$V1)
isoform_cortex$coord  <- gsub("chr", "", isoform_cortex$coord)
isoform_cortex$coord  <- gsub(":", "_", isoform_cortex$coord)
isoform_cortex$coord  <- gsub("-", "_", isoform_cortex$coord)
isoform_cortex$coord  <- paste0(isoform_cortex$id, "_", isoform_cortex$coord)
cortex01_specific_cortex <- isoform_cortex[isoform_cortex$V2 < isoform_cortex$V3, "coord"]  # isoform_cortex exons
cortex02_specific_cortex <- isoform_cortex[isoform_cortex$V2 > isoform_cortex$V3, "coord"]  # isoform_cortex exons
both_cortex <- setdiff(exons_geneRPKM_cortex[exons_geneRPKM_cortex$cortex01exon > 0.1 & exons_geneRPKM_cortex$cortex02exon > 0.1, "exon"], c(cortex01_specific_cortex, cortex02_specific_cortex))
neither_cortex <- setdiff(cortex01$coord, c(cortex01_specific_cortex, cortex02_specific_cortex, both_cortex))
cortex01_exon_0 <- exons_geneRPKM_cortex[exons_geneRPKM_cortex$cortex01exon == 0 & exons_geneRPKM_cortex$cortex01gene > 0.1, "exon"]
cortex02_exon_0 <- exons_geneRPKM_cortex[exons_geneRPKM_cortex$cortex02exon == 0 & exons_geneRPKM_cortex$cortex02gene > 0.1, "exon"]
cortex01_gene_0 <- exons_geneRPKM_cortex[exons_geneRPKM_cortex$cortex01gene == 0, "exon"]
cortex02_gene_0 <- exons_geneRPKM_cortex[exons_geneRPKM_cortex$cortex02gene == 0, "exon"]
# exons_1 <- as.character(exons_geneRPKM_cortex[exons_geneRPKM_cortex$geneRPKM <= 1, "exon"])
# exons_1_10 <- as.character(exons_geneRPKM_cortex[exons_geneRPKM_cortex$geneRPKM > 1 & exons_geneRPKM_cortex$geneRPKM <= 10, "exon"])
# exons_10_100 <- as.character(exons_geneRPKM_cortex[exons_geneRPKM_cortex$geneRPKM > 10 & exons_geneRPKM_cortex$geneRPKM <= 100, "exon"])
# exons_100 <- as.character(exons_geneRPKM_cortex[exons_geneRPKM_cortex$geneRPKM > 100, "exon"])

MeDIP_cortex <- epiProfile(mark = "MeDIP", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_cortex, neither = neither_cortex, cell1_specific = cortex01_specific_cortex, cell2_specific = cortex02_specific_cortex)
MeDIP_profile_cortex <- MeDIP_cortex$profile
(MeDIP_figure_cortex <- MeDIP_cortex$figure)

H3K36me3_cortex <- epiProfile(mark = "H3K36me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_cortex, neither = neither_cortex, cell1_specific = cortex01_specific_cortex, cell2_specific = cortex02_specific_cortex)
H3K36me3_data_cortex <- H3K36me3_cortex$data
H3K36me3_profile_cortex <- H3K36me3_cortex$profile
(H3K36me3_figure_cortex <- H3K36me3_cortex$figure)

# ge 
lib1='A03474'; cell1='ge'; donor1='HuFNSC01';
lib2='A03476'; cell2='ge'; donor2='HuFNSC02';
ge01 <- read.delim(paste0(dirRPKM, lib1, ".G.exn.A.rpkm"), head = F, as.is = T)
ge01$coord <- gsub("<[-]*1", "", ge01$V1)
ge01$coord  <- gsub("chr", "", ge01$coord)
ge01$coord  <- gsub(":", "_", ge01$coord)
ge01$coord  <- gsub("-", "_", ge01$coord)
ge01$coord  <- paste0(ge01$V2, "_", ge01$coord)
rownames(ge01) <- ge01$coord
ge02 <- read.delim(paste0(dirRPKM, lib2, ".G.exn.A.rpkm"), head = F, as.is = T)
ge02$coord <- gsub("<[-]*1", "", ge02$V1)
ge02$coord  <- gsub("chr", "", ge02$coord)
ge02$coord  <- gsub(":", "_", ge02$coord)
ge02$coord  <- gsub("-", "_", ge02$coord)
ge02$coord  <- paste0(ge02$V2, "_", ge02$coord)
rownames(ge02) <- ge02$coord
ge01gene <- read.delim(paste0(dirRPKM, lib1, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(ge01gene) <- ge01gene$V1
ge02gene <- read.delim(paste0(dirRPKM, lib2, ".G.A.rpkm.pc"), head = F, as.is = T)
rownames(ge02gene) <- ge02gene$V1
exons_geneRPKM_ge <- data.frame(exon = ge01$coord)
exons_geneRPKM_ge$exon <- as.character(exons_geneRPKM_ge$exon)
exons_geneRPKM_ge$gene <- gsub("_[0-9XY]+_[0-9_]+", "", exons_geneRPKM_ge$exon)
exons_geneRPKM_ge$ge01gene <- ge01gene[exons_geneRPKM_ge$gene, "V3"]
exons_geneRPKM_ge$ge02gene <- ge02gene[exons_geneRPKM_ge$gene, "V3"]
exons_geneRPKM_ge$geneRPKM <- (exons_geneRPKM_ge$ge01gene + exons_geneRPKM_ge$ge02gene)/2
exons_geneRPKM_ge$ge01exon <- ge01[as.character(exons_geneRPKM_ge$exon), "V4"]
exons_geneRPKM_ge$ge02exon <- ge02[as.character(exons_geneRPKM_ge$exon), "V4"]
exons_geneRPKM_ge <- na.omit(exons_geneRPKM_ge)
isoform_ge <- read.delim(paste0(dirisoform, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform.txt"), head = T, as.is = T)
isoform_ge$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform_ge$V1)
isoform_ge$coord  <- gsub("chr", "", isoform_ge$coord)
isoform_ge$coord  <- gsub(":", "_", isoform_ge$coord)
isoform_ge$coord  <- gsub("-", "_", isoform_ge$coord)
isoform_ge$coord  <- paste0(isoform_ge$id, "_", isoform_ge$coord)
ge01_specific_ge <- isoform_ge[isoform_ge$V2 < isoform_ge$V3, "coord"]  # isoform_ge exons
ge02_specific_ge <- isoform_ge[isoform_ge$V2 > isoform_ge$V3, "coord"]  # isoform_ge exons
both_ge <- setdiff(exons_geneRPKM_ge[exons_geneRPKM_ge$ge01exon > 0.1 & exons_geneRPKM_ge$ge02exon > 0.1, "exon"], c(ge01_specific_ge, ge02_specific_ge))
neither_ge <- setdiff(ge01$coord, c(ge01_specific_ge, ge02_specific_ge, both_ge))
ge01_exon_0 <- exons_geneRPKM_ge[exons_geneRPKM_ge$ge01exon == 0 & exons_geneRPKM_ge$ge01gene > 0.1, "exon"]
ge02_exon_0 <- exons_geneRPKM_ge[exons_geneRPKM_ge$ge02exon == 0 & exons_geneRPKM_ge$ge02gene > 0.1, "exon"]
ge01_gene_0 <- exons_geneRPKM_ge[exons_geneRPKM_ge$ge01gene == 0, "exon"]
ge02_gene_0 <- exons_geneRPKM_ge[exons_geneRPKM_ge$ge02gene == 0, "exon"]
# exons_1 <- as.character(exons_geneRPKM_ge[exons_geneRPKM_ge$geneRPKM <= 1, "exon"])
# exons_1_10 <- as.character(exons_geneRPKM_ge[exons_geneRPKM_ge$geneRPKM > 1 & exons_geneRPKM_ge$geneRPKM <= 10, "exon"])
# exons_10_100 <- as.character(exons_geneRPKM_ge[exons_geneRPKM_ge$geneRPKM > 10 & exons_geneRPKM_ge$geneRPKM <= 100, "exon"])
# exons_100 <- as.character(exons_geneRPKM_ge[exons_geneRPKM_ge$geneRPKM > 100, "exon"])

MeDIP_ge <- epiProfile(mark = "MeDIP", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_ge, neither = neither_ge, cell1_specific = ge01_specific_ge, cell2_specific = ge02_specific_ge)
MeDIP_profile_ge <- MeDIP_ge$profile
(MeDIP_figure_ge <- MeDIP_ge$figure)

H3K36me3_ge <- epiProfile(mark = "H3K36me3", cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, dirIn = dirIn, both = both_ge, neither = neither_ge, cell1_specific = ge01_specific_ge, cell2_specific = ge02_specific_ge)
H3K36me3_data_ge <- H3K36me3_ge$data
H3K36me3_profile_ge <- H3K36me3_ge$profile
(H3K36me3_figure_ge <- H3K36me3_ge$figure)

save(exons_geneRPKM_brain, isoform_brain, brain01_specific_brain, brain02_specific_brain, both_brain, neither_brain, MeDIP_profile_brain, MeDIP_figure_brain, H3K36me3_data_brain, H3K36me3_profile_brain, H3K36me3_figure_brain, 
     exons_geneRPKM_cortex, isoform_cortex, cortex01_specific_cortex, cortex02_specific_cortex, both_cortex, neither_cortex, MeDIP_profile_cortex, MeDIP_figure_cortex, H3K36me3_data_cortex, H3K36me3_profile_cortex, H3K36me3_figure_cortex, 
     exons_geneRPKM_ge, isoform_ge, ge01_specific_ge, ge02_specific_ge, both_ge, neither_ge, MeDIP_profile_ge, MeDIP_figure_ge, H3K36me3_data_ge, H3K36me3_profile_ge, H3K36me3_figure_ge, 
     file = "exonProfile_twins.Rdata")

save(WGBS_figure_HuFNSC02, MeDIP_figure_HuFNSC02, WGBS_figure_HuFNSC04, MeDIP_figure_HuFNSC01, H3K36me3_figure_HuFNSC01, H3K36me3_figure_HuFNSC02, 
     MeDIP_figure_brain, MeDIP_figure_cortex, MeDIP_figure_ge, H3K36me3_figure_brain, H3K36me3_figure_cortex, H3K36me3_figure_ge, 
     file = "exonProfile_Rmd.Rdata")

