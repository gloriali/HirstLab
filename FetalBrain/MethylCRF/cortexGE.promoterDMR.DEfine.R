# intersect DMRs with promoter regions 
# /Applications/bedtools-2.17.0/bin/intersectBed -a ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE.bed -b ~/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_promoter_cortexGE.txt
# /Applications/bedtools-2.17.0/bin/intersectBed -a ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE01.bed -b ~/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_promoter_cortexGE01.txt
# /Applications/bedtools-2.17.0/bin/intersectBed -a ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE02.bed -b ~/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_promoter_cortexGE02.txt

setwd("~/快盘/FetalBrain/MeDIPMRE/cortexVsGE/")
load("../crfaverage.Rdata")

# promoter DMRs and DE
cortexge01_promoter <- read.delim("./DMR_promoter_cortexGE01.txt", head = F, as.is = T)
nrow(cortexge01_promoter)  # No. of promoter DMRs
cortexge01_promoter$Ensembl <- gsub("_[a-zA-Z0-9_]+", "", cortexge01_promoter$V8)
cortexge01DEup <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/UP.Cortex-HuFNSC01_GE-HuFNSC01.FDR_0.01.rmin_0.005.Nmin_25", head=F, as.is=T)
cortexge01DEdn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/DN.Cortex-HuFNSC01_GE-HuFNSC01.FDR_0.01.rmin_0.005.Nmin_25", head=F, as.is=T)
cortexge01_promoter$DE <- NA
cortexge01_promoter[is.element(cortexge01_promoter$Ensembl, cortexge01DEup$V1),]$DE <- "up"
cortexge01_promoter[is.element(cortexge01_promoter$Ensembl, cortexge01DEdn$V1),]$DE <- "dn"
cortexge01_promoterDE <- na.omit(cortexge01_promoter)
cortexge01_promoterDE$logfold <- crfaverage[cortexge01_promoterDE$V4,]$logfold_cortexge01
cortexge01_promoterDE <- cortexge01_promoterDE[(cortexge01_promoterDE$DE == "dn" & cortexge01_promoterDE$logfold > 0) | (cortexge01_promoterDE$DE == "up" & cortexge01_promoterDE$logfold < 0),]
nrow(cortexge01_promoterDE)  # No. of promoter DMRs correlate to DE from DEfine
write.table(cortexge01_promoterDE, file = "cortexge01_promoterDMR_DE.txt", quote=F, row.names=F, sep="\t")

cortexge02_promoter <- read.delim("./DMR_promoter_cortexge02.txt", head = F, as.is = T)
nrow(cortexge02_promoter)  # No. of promoter DMRs
cortexge02_promoter$Ensembl <- gsub("_[a-zA-Z0-9_]+", "", cortexge02_promoter$V8)
cortexge02DEup <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/UP.Cortex-HuFNSC02_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head=F, as.is=T)
cortexge02DEdn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/DN.Cortex-HuFNSC02_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head=F, as.is=T)
cortexge02_promoter$DE <- NA
cortexge02_promoter[is.element(cortexge02_promoter$Ensembl, cortexge02DEup$V1),]$DE <- "up"
cortexge02_promoter[is.element(cortexge02_promoter$Ensembl, cortexge02DEdn$V1),]$DE <- "dn"
cortexge02_promoterDE <- na.omit(cortexge02_promoter)
cortexge02_promoterDE$logfold <- crfaverage[cortexge02_promoterDE$V4,]$logfold_cortexge02
cortexge02_promoterDE <- cortexge02_promoterDE[(cortexge02_promoterDE$DE == "dn" & cortexge02_promoterDE$logfold > 0) | (cortexge02_promoterDE$DE == "up" & cortexge02_promoterDE$logfold < 0),]
nrow(cortexge02_promoterDE)  # No. of promoter DMRs correlate to DE from DEfine
write.table(cortexge02_promoterDE, file = "cortexge02_promoterDMR_DE.txt", quote=F, row.names=F, sep="\t")

# differential methylated at promoters and expressed genes in both individuals
(cortexge_promoterDE <- intersect(cortexge01_promoterDE$Ensembl, cortexge02_promoterDE$Ensembl))
write.table(cortexge_promoterDE, file = "cortexge_promoterDMR_DE.txt", quote=F, row.names=F, sep="\t")

