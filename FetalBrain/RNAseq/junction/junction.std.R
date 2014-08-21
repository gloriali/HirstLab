##################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# junction validation: enough junction coverage: sum of junction coverage in two libraries > 1
# junction RPKM changes in the same direction as exon RPKM  
# junction RPKM cutoff: junction RPKM > 0.1 in one sample and < 0.1 in the other   

##################################################################################
setwd("~/FetalBrain/RNAseq/junction/")
source("/home/lli/bin/R-3.0.2/junction.R")

Nbrain01 <- 61400917   
Nbrain02 <- 69798425
Ncortex01 <- 64807727
Ncortex02 <- 68340512
Ncortex03 <- 76475869
Ncortex04 <- 226802068
Nge01 <- 68292449
Nge02 <- 78885239
Nge03 <- 211087273
Nge04 <- 209694372

# isoform_valid between cortex and GE
lib1='A03473'; cell1='Cortex'; donor1='HuFNSC01'; Nread1=Ncortex01;
lib2='A03474'; cell2='GE'; donor2='HuFNSC01'; Nread2=Nge01;
cortex01_GE01 <- junction(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, Nread1 = Nread1, Nread2 = Nread2)
cortex01_GE01_summary <- cortex01_GE01$summary
cortex01_GE01_isoform_valid_exon <- cortex01_GE01$isoform_valid_exon
cortex01_GE01_isoform_valid_gene <- cortex01_GE01$isoform_valid_gene
rm(cortex01_GE01, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A03475'; cell1='Cortex'; donor1='HuFNSC02'; Nread1=Ncortex02;
lib2='A03476'; cell2='GE'; donor2='HuFNSC02'; Nread2=Nge02;
cortex02_GE02 <- junction(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, Nread1 = Nread1, Nread2 = Nread2)
cortex02_GE02_summary <- cortex02_GE02$summary
cortex02_GE02_isoform_valid_exon <- cortex02_GE02$isoform_valid_exon
cortex02_GE02_isoform_valid_gene <- cortex02_GE02$isoform_valid_gene
rm(cortex02_GE02, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A04599'; cell1='Cortex'; donor1='HuFNSC03'; Nread1=Ncortex03;
lib2='A15295'; cell2='GE'; donor2='HuFNSC03'; Nread2=Nge03;
cortex03_GE03 <- junction(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, Nread1 = Nread1, Nread2 = Nread2)
cortex03_GE03_summary <- cortex03_GE03$summary
cortex03_GE03_isoform_valid_exon <- cortex03_GE03$isoform_valid_exon
cortex03_GE03_isoform_valid_gene <- cortex03_GE03$isoform_valid_gene
rm(cortex03_GE03, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A15298'; cell1='Cortex'; donor1='HuFNSC04'; Nread1=Ncortex04; 
lib2='A15299'; cell2='GE'; donor2='HuFNSC04'; Nread2=Nge04;
cortex04_GE04 <- junction(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, Nread1 = Nread1, Nread2 = Nread2)
cortex04_GE04_summary <- cortex04_GE04$summary
cortex04_GE04_isoform_valid_exon <- cortex04_GE04$isoform_valid_exon
cortex04_GE04_isoform_valid_gene <- cortex04_GE04$isoform_valid_gene
rm(cortex04_GE04, lib1, lib2, cell1, cell2, donor1, donor2)

cortex_GE_valid_summary <- rbind(cortex01_GE01_summary, cortex02_GE02_summary, cortex03_GE03_summary, cortex04_GE04_summary)
rm(cortex01_GE01_summary, cortex02_GE02_summary, cortex03_GE03_summary, cortex04_GE04_summary)
write.table(cortex_GE_valid_summary, file = "cortex_GE_isoform_valid_summary.txt", sep = "\t", quote = F)

# isoform_valid between individuals
lib1='A03484'; cell1='Brain'; donor1='HuFNSC01'; Nread1=Nbrain01;
lib2='A07825'; cell2='Brain'; donor2='HuFNSC02'; Nread2=Nbrain02; 
brain01_brain02 <- junction(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, Nread1 = Nread1, Nread2 = Nread2)
brain01_brain02_summary <- brain01_brain02$summary
brain01_brain02_isoform_valid_exon <- brain01_brain02$isoform_valid_exon
brain01_brain02_isoform_valid_gene <- brain01_brain02$isoform_valid_gene
rm(brain01_brain02, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A03473'; cell1='Cortex'; donor1='HuFNSC01'; Nread1=Ncortex01; 
lib2='A03475'; cell2='Cortex'; donor2='HuFNSC02'; Nread2=Ncortex02; 
cortex01_cortex02 <- junction(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, Nread1 = Nread1, Nread2 = Nread2)
cortex01_cortex02_summary <- cortex01_cortex02$summary
cortex01_cortex02_isoform_valid_exon <- cortex01_cortex02$isoform_valid_exon
cortex01_cortex02_isoform_valid_gene <- cortex01_cortex02$isoform_valid_gene
rm(cortex01_cortex02, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A03474'; cell1='GE'; donor1='HuFNSC01'; Nread1=Nge01; 
lib2='A03476'; cell2='GE'; donor2='HuFNSC02'; Nread2=Nge02; 
GE01_GE02 <- junction(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, Nread1 = Nread1, Nread2 = Nread2)
GE01_GE02_summary <- GE01_GE02$summary
GE01_GE02_isoform_valid_exon <- GE01_GE02$isoform_valid_exon
GE01_GE02_isoform_valid_gene <- GE01_GE02$isoform_valid_gene
rm(GE01_GE02, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A04599'; cell1='Cortex'; donor1='HuFNSC03'; Nread1=Ncortex03; 
lib2='A15298'; cell2='Cortex'; donor2='HuFNSC04'; Nread2=Ncortex04; 
cortex03_cortex04 <- junction(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, Nread1 = Nread1, Nread2 = Nread2)
cortex03_cortex04_summary <- cortex03_cortex04$summary
cortex03_cortex04_isoform_valid_exon <- cortex03_cortex04$isoform_valid_exon
cortex03_cortex04_isoform_valid_gene <- cortex03_cortex04$isoform_valid_gene
rm(cortex03_cortex04, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A15295'; cell1='GE'; donor1='HuFNSC03'; Nread1=Nge03; 
lib2='A15299'; cell2='GE'; donor2='HuFNSC04'; Nread2=Nge04; 
GE03_GE04 <- junction(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, Nread1 = Nread1, Nread2 = Nread2)
GE03_GE04_summary <- GE03_GE04$summary
GE03_GE04_isoform_valid_exon <- GE03_GE04$isoform_valid_exon
GE03_GE04_isoform_valid_gene <- GE03_GE04$isoform_valid_gene
rm(GE03_GE04, lib1, lib2, cell1, cell2, donor1, donor2)

individual_valid_summary <- rbind(brain01_brain02_summary, cortex01_cortex02_summary, GE01_GE02_summary, cortex03_cortex04_summary, GE03_GE04_summary)
rm(brain01_brain02_summary, cortex01_cortex02_summary, GE01_GE02_summary, cortex03_cortex04_summary, GE03_GE04_summary)
write.table(individual_valid_summary, file = "individual_isoform_valid_summary.txt", sep = "\t", quote = F)
rm(list = ls(pattern = "^N"))
save.image("FetalBrain_junction.Rdata")
rm(list = ls())

##################################################################################
setwd("~/快盘/FetalBrain/RNAseq/junction/")
load("FetalBrain_junction.Rdata")
source("~/HirstLab/Pipeline/enrich.R")
library(VennDiagram)
library(ggplot2)

# cortex vs GE
isoform_valid_cortex_GE <- list(HuFNSC01 = cortex01_GE01_isoform_valid_gene$geneID, HuFNSC02 = cortex02_GE02_isoform_valid_gene$geneID, HuFNSC03 = cortex03_GE03_isoform_valid_gene$geneID, HuFNSC04 = cortex04_GE04_isoform_valid_gene$geneID)
venn_cortex_GE_valid <- venn.diagram(isoform_valid_cortex_GE, filename = NULL, fill = c("red", "blue", "green", "purple"), main = "Venn diagram of cortex vs GE isoform_valids", main.cex = 2)
pdf("venn_cortex_GE_isoform_valid.pdf")
plot.new()
grid.draw(venn_cortex_GE_valid)
dev.off()

# HuFNSC01 vs HuFNSC02
isoform_valid_HuFNSC01_HuFNSC02 <- list(brain = brain01_brain02_isoform_valid_gene$geneID, cortex = cortex01_cortex02_isoform_valid_gene$geneID, GE = GE01_GE02_isoform_valid_gene$geneID)
venn_HuFNSC01_HuFNSC02_valid <- venn.diagram(isoform_valid_HuFNSC01_HuFNSC02, filename = NULL, fill = c("green", "red", "blue"), main = "Venn diagram of HuFNSC01 vs HuFNSC02 isoform_valids", main.cex = 2)
pdf("venn_HuFNSC01_HuFNSC02_isoform_valid.pdf")
plot.new()
grid.draw(venn_HuFNSC01_HuFNSC02_valid)
dev.off()

# HuFNSC03 vs HuFNSC04
isoform_valid_HuFNSC03_HuFNSC04 <- list(cortex = cortex03_cortex04_isoform_valid_gene$geneID, GE = GE03_GE04_isoform_valid_gene$geneID)
venn_HuFNSC03_HuFNSC04_valid <- venn.diagram(isoform_valid_HuFNSC03_HuFNSC04, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of HuFNSC03 vs HuFNSC04 isoform_valids", main.cex = 2)
pdf("venn_HuFNSC03_HuFNSC04_isoform_valid.pdf")
plot.new()
grid.draw(venn_HuFNSC03_HuFNSC04_valid)
dev.off()

save.image("FetalBrain_isoform_valid.Rdata")

##################################################################################
setwd("~/快盘/FetalBrain/RNAseq/junction/")
load("FetalBrain_isoform_valid.Rdata")
load("~/快盘/FetalBrain/RNAseq/isoform/FetalBrain_isoform.Rdata")
source("~/HirstLab/Pipeline/enrich.R")
library(VennDiagram)
library(ggplot2)

# gene RPKM
brain01 <- read.delim("~/FetalBrain/RNAseq/rpkm/A03484.G.A.rpkm.pc", head = F, as.is = T)
brain02 <- read.delim("~/FetalBrain/RNAseq/rpkm/A07825.G.A.rpkm.pc", head = F, as.is = T)
cortex01 <- read.delim("~/FetalBrain/RNAseq/rpkm/A03473.G.A.rpkm.pc", head = F, as.is = T)
cortex02 <- read.delim("~/FetalBrain/RNAseq/rpkm/A03475.G.A.rpkm.pc", head = F, as.is = T)
cortex03 <- read.delim("~/FetalBrain/RNAseq/rpkm/A04599.G.A.rpkm.pc", head = F, as.is = T)
cortex04 <- read.delim("~/FetalBrain/RNAseq/rpkm/A15298.G.A.rpkm.pc", head = F, as.is = T)
ge01 <- read.delim("~/FetalBrain/RNAseq/rpkm/A03474.G.A.rpkm.pc", head = F, as.is = T)
ge02 <- read.delim("~/FetalBrain/RNAseq/rpkm/A03476.G.A.rpkm.pc", head = F, as.is = T)
ge03 <- read.delim("~/FetalBrain/RNAseq/rpkm/A15295.G.A.rpkm.pc", head = F, as.is = T)
ge04 <- read.delim("~/FetalBrain/RNAseq/rpkm/A15299.G.A.rpkm.pc", head = F, as.is = T)
cortex01 <- cortex01[cortex01$V1 %in% brain01$V1,]
cortex04 <- cortex04[cortex04$V1 %in% brain01$V1,]
ge01 <- ge01[ge01$V1 %in% brain01$V1,]
ge03 <- ge03[ge03$V1 %in% brain01$V1,]
ge04 <- ge04[ge04$V1 %in% brain01$V1,]
geneRPKM <- data.frame(id = brain01$V1, brain01 = brain01$V3, brain02 = brain02$V3, cortex01 = cortex01$V3, cortex02 = cortex02$V3, cortex03 = cortex03$V3, cortex04 = cortex04$V3, ge01 = ge01$V3, ge02 = ge02$V3, ge03 = ge03$V3, ge04 = ge04$V3)
rm(brain01, brain02, cortex01, cortex02, cortex03, cortex04, ge01, ge02, ge03, ge04)

# DE gene list  
DE_cortex01_GE01up <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/UP.Cortex-HuFNSC01_GE-HuFNSC01.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
DE_cortex01_GE01dn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/DN.Cortex-HuFNSC01_GE-HuFNSC01.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
cortex01_GE01_DE <- unique(c(DE_cortex01_GE01up$V1, DE_cortex01_GE01dn$V1))
rm(DE_cortex01_GE01up, DE_cortex01_GE01dn)
DE_cortex02_GE02up <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/UP.Cortex-HuFNSC02_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
DE_cortex02_GE02dn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/DN.Cortex-HuFNSC02_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
cortex02_GE02_DE <- unique(c(DE_cortex02_GE02up$V1, DE_cortex02_GE02dn$V1))
rm(DE_cortex02_GE02up, DE_cortex02_GE02dn)
DE_cortex03_GE03up <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/UP.Cortex-HuFNSC03_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
DE_cortex03_GE03dn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/DN.Cortex-HuFNSC03_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
cortex03_GE03_DE <- unique(c(DE_cortex03_GE03up$V1, DE_cortex03_GE03dn$V1))
rm(DE_cortex03_GE03up, DE_cortex03_GE03dn)
DE_cortex04_GE04up <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/UP.Cortex-HuFNSC04_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
DE_cortex04_GE04dn <- read.delim("~/快盘/FetalBrain/RNAseq/DEfine/gene/cortexge/DN.Cortex-HuFNSC04_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25", head = F, as.is = T)
cortex04_GE04_DE <- unique(c(DE_cortex04_GE04up$V1, DE_cortex04_GE04dn$V1))
rm(DE_cortex04_GE04up, DE_cortex04_GE04dn)

# No. of exons for DE genes / isoform genes: use collapsed exons instead of Ensembl exons
# exon <- read.delim("~/hg19/hg19v65_exons_collapsed.txt", head = F, as.is = T)
# Nexon <- data.frame(id = levels(as.factor(exon$V2)), Nexon = sapply(levels(as.factor(exon$V2)), function(x) sum(exon$V2 == x)))
# write.table(Nexon, file = "~/hg19/hg19v65_Nexons_per_gene.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# rm(exon)
Nexon <- read.delim("~/hg19/hg19v65_Nexons_per_gene.txt", head = F, as.is = T, row.names = 1)
colnames(Nexon) <- "Nexon"
cortex01_GE01_DE <- data.frame(id = cortex01_GE01_DE, Nexon = Nexon[cortex01_GE01_DE, "Nexon"])
cortex01_GE01_isoform_gene$Nexon <- Nexon[cortex01_GE01_isoform_gene$id, "Nexon"]
cortex01_GE01_isoform_valid_gene$Nexon <- Nexon[as.character(cortex01_GE01_isoform_valid_gene$geneID), "Nexon"] 
cortex02_GE02_DE <- data.frame(id = cortex02_GE02_DE, Nexon = Nexon[cortex02_GE02_DE, "Nexon"])
cortex02_GE02_isoform_gene$Nexon <- Nexon[cortex02_GE02_isoform_gene$id, "Nexon"]
cortex02_GE02_isoform_valid_gene$Nexon <- Nexon[as.character(cortex02_GE02_isoform_valid_gene$geneID), "Nexon"] 
cortex03_GE03_DE <- data.frame(id = cortex03_GE03_DE, Nexon = Nexon[cortex03_GE03_DE, "Nexon"])
cortex03_GE03_isoform_gene$Nexon <- Nexon[cortex03_GE03_isoform_gene$id, "Nexon"]
cortex03_GE03_isoform_valid_gene$Nexon <- Nexon[as.character(cortex03_GE03_isoform_valid_gene$geneID), "Nexon"] 
cortex04_GE04_DE <- data.frame(id = cortex04_GE04_DE, Nexon = Nexon[cortex04_GE04_DE, "Nexon"])
cortex04_GE04_isoform_gene$Nexon <- Nexon[cortex04_GE04_isoform_gene$id, "Nexon"]
cortex04_GE04_isoform_valid_gene$Nexon <- Nexon[as.character(cortex04_GE04_isoform_valid_gene$geneID), "Nexon"] 
pdf("Nexon_DEvsIsoforms.pdf")
plot(c(0, 50), c(0, 0.08), type = "n", main = "No. of exons in DE genes and isoforms", xlab = "No. of exons", ylab = "density")
lines(density(na.omit(Nexon[rownames(Nexon) %in% as.character(geneRPKM[(geneRPKM$cortex01 + geneRPKM$cortex02 + geneRPKM$cortex03 + geneRPKM$cortex04 + geneRPKM$ge01 + geneRPKM$ge02 + geneRPKM$ge03 + geneRPKM$ge04) > 0.08, "id"]), "Nexon"])), col = 1, lty = 1, lwd = 5)
lines(density(na.omit(cortex01_GE01_DE$Nexon)), col = "red", lty = 1, lwd = 3)
lines(density(na.omit(cortex02_GE02_DE$Nexon)), col = "blue", lty = 1, lwd = 3)
lines(density(na.omit(cortex03_GE03_DE$Nexon)), col = "green", lty = 1, lwd = 3)
lines(density(na.omit(cortex04_GE04_DE$Nexon)), col = "purple", lty = 1, lwd = 3)
lines(density(na.omit(cortex01_GE01_isoform_gene$Nexon)), col = "red", lty = 3, lwd = 3)
lines(density(na.omit(cortex02_GE02_isoform_gene$Nexon)), col = "blue", lty = 3, lwd = 3)
lines(density(na.omit(cortex03_GE03_isoform_gene$Nexon)), col = "green", lty = 3, lwd = 3)
lines(density(na.omit(cortex04_GE04_isoform_gene$Nexon)), col = "purple", lty = 3, lwd = 3)
legend("topleft", c("all expressed genes", "DE genes", "isoforms"), col = 1, lty = c(1, 1, 3), lwd = 5, cex = 0.8)
legend("topright", c("all expressed genes", "cortex01 vs GE01", "cortex02 vs GE02", "cortex03 vs GE03", "cortex04 vs GE04"), col = c("black", "red", "blue", "green", "purple"), lty = 1, lwd = 5, cex = 0.8)
dev.off()

# Position of isoform exons (mid point of the exon) on the gene
# cortex01 vs GE01
cortex01_GE01_isoform_exon$start <- ensembl[cortex01_GE01_isoform_exon$id, "start"]
cortex01_GE01_isoform_exon$end <- ensembl[cortex01_GE01_isoform_exon$id, "end"]
cortex01_GE01_isoform_exon$strand <- ensembl[cortex01_GE01_isoform_exon$id, "strand"]
cortex01_GE01_isoform_exon$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", cortex01_GE01_isoform_exon$V1))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", cortex01_GE01_isoform_exon$V1))))/2
cortex01_GE01_isoform_exon[cortex01_GE01_isoform_exon$strand == 1, ]$exon_pos <- (cortex01_GE01_isoform_exon[cortex01_GE01_isoform_exon$strand == 1, ]$exon_pos - cortex01_GE01_isoform_exon[cortex01_GE01_isoform_exon$strand == 1, ]$start) / (cortex01_GE01_isoform_exon[cortex01_GE01_isoform_exon$strand == 1, ]$end - cortex01_GE01_isoform_exon[cortex01_GE01_isoform_exon$strand == 1, ]$start) * 100
cortex01_GE01_isoform_exon[cortex01_GE01_isoform_exon$strand == -1, ]$exon_pos <- 100 - (cortex01_GE01_isoform_exon[cortex01_GE01_isoform_exon$strand == -1, ]$exon_pos - cortex01_GE01_isoform_exon[cortex01_GE01_isoform_exon$strand == -1, ]$start) / (cortex01_GE01_isoform_exon[cortex01_GE01_isoform_exon$strand == -1, ]$end - cortex01_GE01_isoform_exon[cortex01_GE01_isoform_exon$strand == -1, ]$start) * 100
# cortex02 vs GE02
cortex02_GE02_isoform_exon$start <- ensembl[cortex02_GE02_isoform_exon$id, "start"]
cortex02_GE02_isoform_exon$end <- ensembl[cortex02_GE02_isoform_exon$id, "end"]
cortex02_GE02_isoform_exon$strand <- ensembl[cortex02_GE02_isoform_exon$id, "strand"]
cortex02_GE02_isoform_exon$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", cortex02_GE02_isoform_exon$V1))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", cortex02_GE02_isoform_exon$V1))))/2
cortex02_GE02_isoform_exon[cortex02_GE02_isoform_exon$strand == 1, ]$exon_pos <- (cortex02_GE02_isoform_exon[cortex02_GE02_isoform_exon$strand == 1, ]$exon_pos - cortex02_GE02_isoform_exon[cortex02_GE02_isoform_exon$strand == 1, ]$start) / (cortex02_GE02_isoform_exon[cortex02_GE02_isoform_exon$strand == 1, ]$end - cortex02_GE02_isoform_exon[cortex02_GE02_isoform_exon$strand == 1, ]$start) * 100
cortex02_GE02_isoform_exon[cortex02_GE02_isoform_exon$strand == -1, ]$exon_pos <- 100 - (cortex02_GE02_isoform_exon[cortex02_GE02_isoform_exon$strand == -1, ]$exon_pos - cortex02_GE02_isoform_exon[cortex02_GE02_isoform_exon$strand == -1, ]$start) / (cortex02_GE02_isoform_exon[cortex02_GE02_isoform_exon$strand == -1, ]$end - cortex02_GE02_isoform_exon[cortex02_GE02_isoform_exon$strand == -1, ]$start) * 100
# cortex03 vs GE03
cortex03_GE03_isoform_exon$start <- ensembl[cortex03_GE03_isoform_exon$id, "start"]
cortex03_GE03_isoform_exon$end <- ensembl[cortex03_GE03_isoform_exon$id, "end"]
cortex03_GE03_isoform_exon$strand <- ensembl[cortex03_GE03_isoform_exon$id, "strand"]
cortex03_GE03_isoform_exon$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", cortex03_GE03_isoform_exon$V1))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", cortex03_GE03_isoform_exon$V1))))/2
cortex03_GE03_isoform_exon[cortex03_GE03_isoform_exon$strand == 1, ]$exon_pos <- (cortex03_GE03_isoform_exon[cortex03_GE03_isoform_exon$strand == 1, ]$exon_pos - cortex03_GE03_isoform_exon[cortex03_GE03_isoform_exon$strand == 1, ]$start) / (cortex03_GE03_isoform_exon[cortex03_GE03_isoform_exon$strand == 1, ]$end - cortex03_GE03_isoform_exon[cortex03_GE03_isoform_exon$strand == 1, ]$start) * 100
cortex03_GE03_isoform_exon[cortex03_GE03_isoform_exon$strand == -1, ]$exon_pos <- 100 - (cortex03_GE03_isoform_exon[cortex03_GE03_isoform_exon$strand == -1, ]$exon_pos - cortex03_GE03_isoform_exon[cortex03_GE03_isoform_exon$strand == -1, ]$start) / (cortex03_GE03_isoform_exon[cortex03_GE03_isoform_exon$strand == -1, ]$end - cortex03_GE03_isoform_exon[cortex03_GE03_isoform_exon$strand == -1, ]$start) * 100
# cortex04 vs GE04
cortex04_GE04_isoform_exon$start <- ensembl[cortex04_GE04_isoform_exon$id, "start"]
cortex04_GE04_isoform_exon$end <- ensembl[cortex04_GE04_isoform_exon$id, "end"]
cortex04_GE04_isoform_exon$strand <- ensembl[cortex04_GE04_isoform_exon$id, "strand"]
cortex04_GE04_isoform_exon$exon_pos <- (as.numeric(gsub("chr[0-9X]+:", "", gsub("-[0-9]+<[-]*1_ENSG[0-9]+", "", cortex04_GE04_isoform_exon$V1))) + as.numeric(gsub("chr[0-9X]+:[0-9]+-", "", gsub("<[-]*1_ENSG[0-9]+", "", cortex04_GE04_isoform_exon$V1))))/2
cortex04_GE04_isoform_exon[cortex04_GE04_isoform_exon$strand == 1, ]$exon_pos <- (cortex04_GE04_isoform_exon[cortex04_GE04_isoform_exon$strand == 1, ]$exon_pos - cortex04_GE04_isoform_exon[cortex04_GE04_isoform_exon$strand == 1, ]$start) / (cortex04_GE04_isoform_exon[cortex04_GE04_isoform_exon$strand == 1, ]$end - cortex04_GE04_isoform_exon[cortex04_GE04_isoform_exon$strand == 1, ]$start) * 100
cortex04_GE04_isoform_exon[cortex04_GE04_isoform_exon$strand == -1, ]$exon_pos <- 100 - (cortex04_GE04_isoform_exon[cortex04_GE04_isoform_exon$strand == -1, ]$exon_pos - cortex04_GE04_isoform_exon[cortex04_GE04_isoform_exon$strand == -1, ]$start) / (cortex04_GE04_isoform_exon[cortex04_GE04_isoform_exon$strand == -1, ]$end - cortex04_GE04_isoform_exon[cortex04_GE04_isoform_exon$strand == -1, ]$start) * 100
pdf("exon_position.pdf")
plot(c(0, 100), c(0, 0.015), type = "n", main = "Distribution of isoform exons along genes", xlab = "gene position", ylab = "density")
lines(density(na.omit(cortex01_GE01_isoform_exon$exon_pos)), col = "red", lty = 1, lwd = 3)
lines(density(na.omit(cortex02_GE02_isoform_exon$exon_pos)), col = "blue", lty = 1, lwd = 3)
lines(density(na.omit(cortex03_GE03_isoform_exon$exon_pos)), col = "green", lty = 1, lwd = 3)
lines(density(na.omit(cortex04_GE04_isoform_exon$exon_pos)), col = "purple", lty = 1, lwd = 3)
legend("bottomright", c("cortex01 vs GE01", "cortex02 vs GE02", "cortex03 vs GE03", "cortex04 vs GE04"), col = c("red", "blue", "green", "purple"), lty = 1, lwd = 5)
dev.off()

# Venn Diagram with average expression level,  No. of exons, and average exon length for MZ twins              
# all isoforms: gene list of different sections of Venn diagram (b: brain01_brain02, c: cortex01_cortex02, g: GE01_GE02)
all_b_only <- geneRPKM[geneRPKM$id %in% setdiff(brain01_brain02_isoform_gene$id, union(cortex01_cortex02_isoform_gene$id, GE01_GE02_isoform_gene$id)), ]
all_c_only <- geneRPKM[geneRPKM$id %in% setdiff(cortex01_cortex02_isoform_gene$id, union(brain01_brain02_isoform_gene$id, GE01_GE02_isoform_gene$id)), ]
all_g_only <- geneRPKM[geneRPKM$id %in% setdiff(GE01_GE02_isoform_gene$id, union(brain01_brain02_isoform_gene$id, cortex01_cortex02_isoform_gene$id)), ]
all_b_c_not_g <- geneRPKM[geneRPKM$id %in% setdiff(intersect(brain01_brain02_isoform_gene$id, cortex01_cortex02_isoform_gene$id), GE01_GE02_isoform_gene$id), ]
all_b_g_not_c <- geneRPKM[geneRPKM$id %in% setdiff(intersect(brain01_brain02_isoform_gene$id, GE01_GE02_isoform_gene$id), cortex01_cortex02_isoform_gene$id), ]
all_c_g_not_b <- geneRPKM[geneRPKM$id %in% setdiff(intersect(cortex01_cortex02_isoform_gene$id, GE01_GE02_isoform_gene$id), brain01_brain02_isoform_gene$id), ]
all_b_c_g <- geneRPKM[geneRPKM$id %in% intersect(intersect(brain01_brain02_isoform_gene$id, cortex01_cortex02_isoform_gene$id), GE01_GE02_isoform_gene$id), ]
all_b_only$Nexon <- Nexon[all_b_only$id, "Nexon"]
all_c_only$Nexon <- Nexon[all_c_only$id, "Nexon"]
all_g_only$Nexon <- Nexon[all_g_only$id, "Nexon"]
all_b_c_not_g$Nexon <- Nexon[all_b_c_not_g$id, "Nexon"]
all_b_g_not_c$Nexon <- Nexon[all_b_g_not_c$id, "Nexon"]
all_c_g_not_b$Nexon <- Nexon[all_c_g_not_b$id, "Nexon"]
all_b_c_g$Nexon <- Nexon[all_b_c_g$id, "Nexon"]
brain01_brain02_isoform_exon$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", brain01_brain02_isoform_exon$V1)))
brain01_brain02_isoform_exon$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", brain01_brain02_isoform_exon$V1)))
brain01_brain02_isoform_exon$exon_length <- abs(brain01_brain02_isoform_exon$exon_end - brain01_brain02_isoform_exon$exon_start)
cortex01_cortex02_isoform_exon$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", cortex01_cortex02_isoform_exon$V1)))
cortex01_cortex02_isoform_exon$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", cortex01_cortex02_isoform_exon$V1)))
cortex01_cortex02_isoform_exon$exon_length <- abs(cortex01_cortex02_isoform_exon$exon_end - cortex01_cortex02_isoform_exon$exon_start)
GE01_GE02_isoform_exon$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", GE01_GE02_isoform_exon$V1)))
GE01_GE02_isoform_exon$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", GE01_GE02_isoform_exon$V1)))
GE01_GE02_isoform_exon$exon_length <- abs(GE01_GE02_isoform_exon$exon_end - GE01_GE02_isoform_exon$exon_start)
all_b_only$exon_length <- sapply(all_b_only$id, function(x) mean(brain01_brain02_isoform_exon[brain01_brain02_isoform_exon$id == x, "exon_length"]))
all_c_only$exon_length <- sapply(all_c_only$id, function(x) mean(cortex01_cortex02_isoform_exon[cortex01_cortex02_isoform_exon$id == x, "exon_length"]))
all_g_only$exon_length <- sapply(all_g_only$id, function(x) mean(GE01_GE02_isoform_exon[GE01_GE02_isoform_exon$id == x, "exon_length"]))
all_b_c_not_g$exon_length <- sapply(all_b_c_not_g$id, function(x) mean(c(brain01_brain02_isoform_exon[brain01_brain02_isoform_exon$id == x, "exon_length"], cortex01_cortex02_isoform_exon[cortex01_cortex02_isoform_exon$id == x, "exon_length"])))
all_b_g_not_c$exon_length <- sapply(all_b_g_not_c$id, function(x) mean(c(brain01_brain02_isoform_exon[brain01_brain02_isoform_exon$id == x, "exon_length"], GE01_GE02_isoform_exon[GE01_GE02_isoform_exon$id == x, "exon_length"])))
all_c_g_not_b$exon_length <- sapply(all_c_g_not_b$id, function(x) mean(c(cortex01_cortex02_isoform_exon[cortex01_cortex02_isoform_exon$id == x, "exon_length"], GE01_GE02_isoform_exon[GE01_GE02_isoform_exon$id == x, "exon_length"])))
all_b_c_g$exon_length <- sapply(all_b_c_g$id, function(x) mean(c(brain01_brain02_isoform_exon[brain01_brain02_isoform_exon$id == x, "exon_length"], cortex01_cortex02_isoform_exon[cortex01_cortex02_isoform_exon$id == x, "exon_length"], GE01_GE02_isoform_exon[GE01_GE02_isoform_exon$id == x, "exon_length"])))
# validated isoforms: gene list of different sections of Venn diagram (b: brain01_brain02, c: cortex01_cortex02, g: GE01_GE02)
valid_b_only <- geneRPKM[geneRPKM$id %in% setdiff(brain01_brain02_isoform_valid_gene$geneID, union(cortex01_cortex02_isoform_valid_gene$geneID, GE01_GE02_isoform_valid_gene$geneID)), ]
valid_c_only <- geneRPKM[geneRPKM$id %in% setdiff(cortex01_cortex02_isoform_valid_gene$geneID, union(brain01_brain02_isoform_valid_gene$geneID, GE01_GE02_isoform_valid_gene$geneID)), ]
valid_g_only <- geneRPKM[geneRPKM$id %in% setdiff(GE01_GE02_isoform_valid_gene$geneID, union(brain01_brain02_isoform_valid_gene$geneID, cortex01_cortex02_isoform_valid_gene$geneID)), ]
valid_b_c_not_g <- geneRPKM[geneRPKM$id %in% setdiff(intersect(brain01_brain02_isoform_valid_gene$geneID, cortex01_cortex02_isoform_valid_gene$geneID), GE01_GE02_isoform_valid_gene$geneID), ]
valid_b_g_not_c <- geneRPKM[geneRPKM$id %in% setdiff(intersect(brain01_brain02_isoform_valid_gene$geneID, GE01_GE02_isoform_valid_gene$geneID), cortex01_cortex02_isoform_valid_gene$geneID), ]
valid_c_g_not_b <- geneRPKM[geneRPKM$id %in% setdiff(intersect(cortex01_cortex02_isoform_valid_gene$geneID, GE01_GE02_isoform_valid_gene$geneID), brain01_brain02_isoform_valid_gene$geneID), ]
valid_b_c_g <- geneRPKM[geneRPKM$id %in% intersect(intersect(brain01_brain02_isoform_valid_gene$geneID, cortex01_cortex02_isoform_valid_gene$geneID), GE01_GE02_isoform_valid_gene$geneID), ]
valid_b_only$Nexon <- Nexon[valid_b_only$id, "Nexon"]
valid_c_only$Nexon <- Nexon[valid_c_only$id, "Nexon"]
valid_g_only$Nexon <- Nexon[valid_g_only$id, "Nexon"]
valid_b_c_not_g$Nexon <- Nexon[valid_b_c_not_g$id, "Nexon"]
valid_b_g_not_c$Nexon <- Nexon[valid_b_g_not_c$id, "Nexon"]
valid_c_g_not_b$Nexon <- Nexon[valid_c_g_not_b$id, "Nexon"]
valid_b_c_g$Nexon <- Nexon[valid_b_c_g$id, "Nexon"]
brain01_brain02_isoform_valid_exon$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", brain01_brain02_isoform_valid_exon$exonID)))
brain01_brain02_isoform_valid_exon$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", brain01_brain02_isoform_valid_exon$exonID)))
brain01_brain02_isoform_valid_exon$exon_length <- abs(brain01_brain02_isoform_valid_exon$exon_end - brain01_brain02_isoform_valid_exon$exon_start)
cortex01_cortex02_isoform_valid_exon$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", cortex01_cortex02_isoform_valid_exon$exonID)))
cortex01_cortex02_isoform_valid_exon$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", cortex01_cortex02_isoform_valid_exon$exonID)))
cortex01_cortex02_isoform_valid_exon$exon_length <- abs(cortex01_cortex02_isoform_valid_exon$exon_end - cortex01_cortex02_isoform_valid_exon$exon_start)
GE01_GE02_isoform_valid_exon$exon_start <- as.numeric(gsub("-[0-9]+<[-_ENSG0-9]+", "", gsub("chr[0-9XY]+:", "", GE01_GE02_isoform_valid_exon$exonID)))
GE01_GE02_isoform_valid_exon$exon_end <- as.numeric(gsub("<[-_ENSG0-9]+", "", gsub("chr[0-9XY:]+-", "", GE01_GE02_isoform_valid_exon$exonID)))
GE01_GE02_isoform_valid_exon$exon_length <- abs(GE01_GE02_isoform_valid_exon$exon_end - GE01_GE02_isoform_valid_exon$exon_start)
brain01_brain02_isoform_valid_exon$geneID <- as.character(brain01_brain02_isoform_valid_exon$geneID)
cortex01_cortex02_isoform_valid_exon$geneID <- as.character(cortex01_cortex02_isoform_valid_exon$geneID)
GE01_GE02_isoform_valid_exon$geneID <- as.character(GE01_GE02_isoform_valid_exon$geneID)
valid_b_only$exon_length <- sapply(valid_b_only$id, function(x) mean(brain01_brain02_isoform_valid_exon[brain01_brain02_isoform_valid_exon$geneID == x, "exon_length"]))
valid_c_only$exon_length <- sapply(valid_c_only$id, function(x) mean(cortex01_cortex02_isoform_valid_exon[cortex01_cortex02_isoform_valid_exon$geneID == x, "exon_length"]))
valid_g_only$exon_length <- sapply(valid_g_only$id, function(x) mean(GE01_GE02_isoform_valid_exon[GE01_GE02_isoform_valid_exon$geneID == x, "exon_length"]))
valid_b_c_not_g$exon_length <- sapply(valid_b_c_not_g$id, function(x) mean(c(brain01_brain02_isoform_valid_exon[brain01_brain02_isoform_valid_exon$geneID == x, "exon_length"], cortex01_cortex02_isoform_valid_exon[cortex01_cortex02_isoform_valid_exon$geneID == x, "exon_length"])))
valid_b_g_not_c$exon_length <- sapply(valid_b_g_not_c$id, function(x) mean(c(brain01_brain02_isoform_valid_exon[brain01_brain02_isoform_valid_exon$geneID == x, "exon_length"], GE01_GE02_isoform_valid_exon[GE01_GE02_isoform_valid_exon$geneID == x, "exon_length"])))
valid_c_g_not_b$exon_length <- sapply(valid_c_g_not_b$id, function(x) mean(c(cortex01_cortex02_isoform_valid_exon[cortex01_cortex02_isoform_valid_exon$geneID == x, "exon_length"], GE01_GE02_isoform_valid_exon[GE01_GE02_isoform_valid_exon$geneID == x, "exon_length"])))
valid_b_c_g$exon_length <- sapply(valid_b_c_g$id, function(x) mean(c(brain01_brain02_isoform_valid_exon[brain01_brain02_isoform_valid_exon$geneID == x, "exon_length"], cortex01_cortex02_isoform_valid_exon[cortex01_cortex02_isoform_valid_exon$geneID == x, "exon_length"], GE01_GE02_isoform_valid_exon[GE01_GE02_isoform_valid_exon$geneID == x, "exon_length"])))
isoform_all <- list(brain01_brain02 = brain01_brain02_isoform_gene$id, cortex01_cortex02 = cortex01_cortex02_isoform_gene$id, GE01_GE02 = GE01_GE02_isoform_gene$id)
isoform_valid <- list(brain01_brain02 = as.character(brain01_brain02_isoform_valid_gene$gene), cortex01_cortex02 = as.character(cortex01_cortex02_isoform_valid_gene$gene), GE01_GE02 = as.character(GE01_GE02_isoform_valid_gene$gene))
isoform_all_N <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5), N = c(nrow(all_b_only), nrow(all_c_only), nrow(all_g_only), nrow(all_b_c_not_g), nrow(all_c_g_not_b), nrow(all_b_g_not_c), nrow(all_b_c_g))) 
isoform_valid_N <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5), N = c(nrow(valid_b_only), nrow(valid_c_only), nrow(valid_g_only), nrow(valid_b_c_not_g), nrow(valid_c_g_not_b), nrow(valid_b_g_not_c), nrow(valid_b_c_g))) 
isoform_all_rpkm <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), RPKM = c(mean(c(all_b_only$brain01, all_b_only$brain02)), mean(c(all_c_only$cortex01, all_c_only$cortex02)), mean(c(all_g_only$ge01, all_g_only$ge02)), mean(c(all_b_c_not_g$brain01, all_b_c_not_g$brain02, all_b_c_not_g$cortex01, all_b_c_not_g$cortex02)), mean(c(all_c_g_not_b$cortex01, all_c_g_not_b$cortex02, all_c_g_not_b$ge01, all_c_g_not_b$ge02)), mean(c(all_b_g_not_c$brain01, all_b_g_not_c$brain02, all_b_g_not_c$ge01, all_b_g_not_c$ge02)), mean(c(all_b_c_g$brain01, all_b_c_g$brain02, all_b_c_g$cortex01, all_b_c_g$cortex02, all_b_c_g$ge01, all_b_c_g$ge02)), with(geneRPKM, mean(c(brain01[brain01 > 0.1], cortex01[cortex01 > 0.1], ge01[ge01 > 0.1], brain02[brain02 > 0.1], cortex02[cortex02 > 0.1], ge02[ge02 > 0.1]))))) 
isoform_valid_rpkm <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), RPKM = c(mean(c(valid_b_only$brain01, valid_b_only$brain02)), mean(c(valid_c_only$cortex01, valid_c_only$cortex02)), mean(c(valid_g_only$ge01, valid_g_only$ge02)), mean(c(valid_b_c_not_g$brain01, valid_b_c_not_g$brain02, valid_b_c_not_g$cortex01, valid_b_c_not_g$cortex02)), mean(c(valid_c_g_not_b$cortex01, valid_c_g_not_b$cortex02, valid_c_g_not_b$ge01, valid_c_g_not_b$ge02)), mean(c(valid_b_g_not_c$brain01, valid_b_g_not_c$brain02, valid_b_g_not_c$ge01, valid_b_g_not_c$ge02)), mean(c(valid_b_c_g$brain01, valid_b_c_g$brain02, valid_b_c_g$cortex01, valid_b_c_g$cortex02, valid_b_c_g$ge01, valid_b_c_g$ge02)), with(geneRPKM, mean(c(brain01[brain01 > 0.1], cortex01[cortex01 > 0.1], ge01[ge01 > 0.1], brain02[brain02 > 0.1], cortex02[cortex02 > 0.1], ge02[ge02 > 0.1]))))) 
isoform_all_Nexon <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), No.exon = c(mean(all_b_only$Nexon), mean(all_c_only$Nexon), mean(all_g_only$Nexon), mean(all_b_c_not_g$Nexon), mean(all_c_g_not_b$Nexon), mean(all_b_g_not_c$Nexon), mean(all_b_c_g$Nexon), mean(Nexon[geneRPKM[(geneRPKM$cortex01 + geneRPKM$cortex02 + geneRPKM$cortex03 + geneRPKM$cortex04 + geneRPKM$ge01 + geneRPKM$ge02 + geneRPKM$ge03 + geneRPKM$ge04) > 0.05,]$id, "Nexon"]))) 
isoform_valid_Nexon <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), No.exon = c(mean(valid_b_only$Nexon), mean(valid_c_only$Nexon), mean(valid_g_only$Nexon), mean(valid_b_c_not_g$Nexon), mean(valid_c_g_not_b$Nexon), mean(valid_b_g_not_c$Nexon), mean(valid_b_c_g$Nexon), mean(Nexon[geneRPKM[(geneRPKM$cortex01 + geneRPKM$cortex02 + geneRPKM$cortex03 + geneRPKM$cortex04 + geneRPKM$ge01 + geneRPKM$ge02 + geneRPKM$ge03 + geneRPKM$ge04) > 0.05,]$id, "Nexon"]))) 
exon_length <- read.delim("~/hg19/hg19v65_exons_for_genes.length", head = F, as.is = T)
exon_length$gene <- gsub("chr[0-9XY:-]+<[-]*1_", "", exon_length$V1)
isoform_all_exon_length <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), exon_length = c(mean(all_b_only$exon_length), mean(all_c_only$exon_length), mean(all_g_only$exon_length), mean(all_b_c_not_g$exon_length), mean(all_c_g_not_b$exon_length), mean(all_b_g_not_c$exon_length), mean(all_b_c_g$exon_length), mean(exon_length[exon_length$gene %in% geneRPKM[(geneRPKM$cortex01 + geneRPKM$cortex02 + geneRPKM$cortex03 + geneRPKM$cortex04 + geneRPKM$ge01 + geneRPKM$ge02 + geneRPKM$ge03 + geneRPKM$ge04) > 0.05,]$id, "V2"]))) 
isoform_valid_exon_length <- data.frame(x= c(0.2, 0.85, 0.5, 0.5, 0.7, 0.3, 0.5, 0.05), y = c(0.7, 0.7, 0.15, 0.8, 0.4, 0.4, 0.5, 0.05), exon_length = c(mean(valid_b_only$exon_length), mean(valid_c_only$exon_length), mean(valid_g_only$exon_length), mean(valid_b_c_not_g$exon_length), mean(valid_c_g_not_b$exon_length), mean(valid_b_g_not_c$exon_length), mean(valid_b_c_g$exon_length), mean(exon_length[exon_length$gene %in% geneRPKM[(geneRPKM$cortex01 + geneRPKM$cortex02 + geneRPKM$cortex03 + geneRPKM$cortex04 + geneRPKM$ge01 + geneRPKM$ge02 + geneRPKM$ge03 + geneRPKM$ge04) > 0.05,]$id, "V2"]))) 

# No. of isoforms
pdf("venn_isoform.pdf")
(bubble_all_N <- ggplot(isoform_all_N) + 
   geom_point(aes(x, y, size = N), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(50, 1500), range = c(5, 20)) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_all <- venn.diagram(isoform_all, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "No. of isoform genes", main.cex = 1.5)
grid.draw(venn_all)
(bubble_valid_N <- ggplot(isoform_valid_N) + 
   geom_point(aes(x, y, size = N), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(50, 1500), range = c(5, 20)) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_valid <- venn.diagram(isoform_valid, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "No. of validated isoform genes", main.cex = 1.5)
grid.draw(venn_valid)
dev.off()
# average expression level
pdf("venn_rpkm.pdf")
(bubble_all_rpkm <- ggplot(isoform_all_rpkm) + 
   geom_point(aes(x, y, size = RPKM), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size_area(max_size = 30, breaks = c(0.5, 2, 15)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_all <- venn.diagram(isoform_all, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average gene RPKM of isoforms", main.cex = 1.5)
grid.draw(venn_all)
(bubble_valid_rpkm <- ggplot(isoform_valid_rpkm) + 
   geom_point(aes(x, y, size = RPKM), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size_area(max_size = 30, breaks = c(0.5, 2, 15)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_valid <- venn.diagram(isoform_valid, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average gene RPKM of validated isoforms", main.cex = 1.5)
grid.draw(venn_valid)
dev.off()
# No. of exons
pdf("venn_Nexons.pdf")
(bubble_all_Nexon <- ggplot(isoform_all_Nexon) + 
   geom_point(aes(x, y, size = No.exon), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(9, 15), range = c(5, 15), breaks = c(11, 12)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_all <- venn.diagram(isoform_all, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average No. of exons of isoforms", main.cex = 1.5)
grid.draw(venn_all)
(bubble_valid_Nexon <- ggplot(isoform_valid_Nexon) + 
   geom_point(aes(x, y, size = No.exon), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(9, 15), range = c(5, 15), breaks = c(11, 12, 13)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_valid <- venn.diagram(isoform_valid, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average No. of exons of validated isoforms", main.cex = 1.5)
grid.draw(venn_valid)
dev.off()
# average exon length 
pdf("venn_exon_length.pdf")
(bubble_all_exon_length <- ggplot(isoform_all_exon_length) + 
   geom_point(aes(x, y, size = exon_length), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(100, 400), range = c(5, 20), breaks = c(200, 300, 400)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_all <- venn.diagram(isoform_all, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average length of isoform exons", main.cex = 1.5)
grid.draw(venn_all)
(bubble_valid_exon_length <- ggplot(isoform_valid_exon_length) + 
   geom_point(aes(x, y, size = exon_length), color = "red") + 
   scale_x_continuous(limits = c(0, 1)) + 
   scale_y_continuous(limits = c(0, 1)) + 
   scale_size(limits = c(100, 400), range = c(5, 20), breaks = c(200, 300, 400)) + 
   geom_text(data = NULL, x = 0.06, y = 0.14, label = "all expressed genes", size = 3.5) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = c(0.9, 0.15), legend.title=element_blank(), legend.text = element_text(size = 15), legend.key = element_rect(fill = 'white')))
venn_valid <- venn.diagram(isoform_valid, filename = NULL, fill = rep("white", 3), cex = 0, alpha = 0, col = c("red", "blue", "green"), main = "Average length of validated isoform exons", main.cex = 1.5)
grid.draw(venn_valid)
dev.off()

rm(ensembl, junction_exon, exon_length)
rm(list = ls(patter = "^all_"))
rm(list = ls(patter = "^valid_"))
save.image(file = "FetalBrain_isoform_valid.Rdata")


