# intersect DMRs with promoter regions using bedtools 
# /home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_brain.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRbrain_promoter.txt
# /home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_cortex.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRcortex_promoter.txt
# /home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_ge.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRge_promoter.txt

setwd("~/快盘/FetalBrain/MeDIPMRE/Individual/")
load("../crfaverage.Rdata")

# promoter DMRs and DE
brain_promoter <- read.delim("./DMRs/intersect/DMRbrain_promoter.txt", head = F, as.is = T)
nrow(brain_promoter)  # No. of promoter DMRs
brain_promoter$Ensembl <- gsub("_[a-zA-Z0-9_]+", "", brain_promoter$V8)
brainDEup <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/individual/UP.Brain-HuFNSC01_Brain-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head=F, as.is=T)
brainDEdn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/individual/DN.Brain-HuFNSC01_Brain-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head=F, as.is=T)
brain_promoter$DE <- NA
brain_promoter[is.element(brain_promoter$Ensembl, brainDEup$V1),]$DE <- "up"
brain_promoter[is.element(brain_promoter$Ensembl, brainDEdn$V1),]$DE <- "dn"
brain_promoterDE <- na.omit(brain_promoter)
brain_promoterDE$logfold <- crfaverage[brain_promoterDE$V4,]$logfold_brain
brain_promoterDE <- brain_promoterDE[(brain_promoterDE$DE == "dn" & brain_promoterDE$logfold > 0) | (brain_promoterDE$DE == "up" & brain_promoterDE$logfold < 0),]
nrow(brain_promoterDE)  # No. of promoter DMRs correlate to DE from DEfine
write.table(brain_promoterDE, file = "brain_promoterDMR_DE.txt", quote=F, row.names=F, sep="\t")

cortex_promoter <- read.delim("./DMRs/intersect/DMRcortex_promoter.txt", head = F, as.is = T)
nrow(cortex_promoter)  # No. of promoter DMRs
cortex_promoter$Ensembl <- gsub("_[a-zA-Z0-9_]+", "", cortex_promoter$V8)
cortexDEup <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/individual/UP.Cortex-HuFNSC01_Cortex-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head=F, as.is=T)
cortexDEdn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/individual/DN.Cortex-HuFNSC01_Cortex-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head=F, as.is=T)
cortex_promoter$DE <- NA
cortex_promoter[is.element(cortex_promoter$Ensembl, cortexDEup$V1),]$DE <- "up"
cortex_promoter[is.element(cortex_promoter$Ensembl, cortexDEdn$V1),]$DE <- "dn"
cortex_promoterDE <- na.omit(cortex_promoter)
cortex_promoterDE$logfold <- crfaverage[cortex_promoterDE$V4,]$logfold_cortex
cortex_promoterDE <- cortex_promoterDE[(cortex_promoterDE$DE == "dn" & cortex_promoterDE$logfold > 0) | (cortex_promoterDE$DE == "up" & cortex_promoterDE$logfold < 0),]
nrow(cortex_promoterDE)  # No. of promoter DMRs correlate to DE from DEfine
write.table(cortex_promoterDE, file = "cortex_promoterDMR_DE.txt", quote=F, row.names=F, sep="\t")

ge_promoter <- read.delim("./DMRs/intersect/DMRge_promoter.txt", head = F, as.is = T)
nrow(ge_promoter)  # No. of promoter DMRs
ge_promoter$Ensembl <- gsub("_[a-zA-Z0-9_]+", "", ge_promoter$V8)
geDEup <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/individual/UP.GE-HuFNSC01_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head=F, as.is=T)
geDEdn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/individual/DN.GE-HuFNSC01_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head=F, as.is=T)
ge_promoter$DE <- NA
ge_promoter[is.element(ge_promoter$Ensembl, geDEup$V1),]$DE <- "up"
ge_promoter[is.element(ge_promoter$Ensembl, geDEdn$V1),]$DE <- "dn"
ge_promoterDE <- na.omit(ge_promoter)
ge_promoterDE$logfold <- crfaverage[ge_promoterDE$V4,]$logfold_ge
ge_promoterDE <- ge_promoterDE[(ge_promoterDE$DE == "dn" & ge_promoterDE$logfold > 0) | (ge_promoterDE$DE == "up" & ge_promoterDE$logfold < 0),]
nrow(ge_promoterDE)  # No. of promoter DMRs correlate to DE from DEfine
write.table(ge_promoterDE, file = "ge_promoterDMR_DE.txt", quote=F, row.names=F, sep="\t")

