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

cortex_GE_summary <- rbind(cortex01_GE01_summary, cortex02_GE02_summary, cortex03_GE03_summary, cortex04_GE04_summary)
write.table(cortex_GE_summary, file = "cortex_GE_isoform_valid_summary.txt", sep = "\t", quote = F)

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

individual_summary <- rbind(brain01_brain02_summary, cortex01_cortex02_summary, GE01_GE02_summary, cortex03_cortex04_summary, GE03_GE04_summary)
write.table(individual_summary, file = "individual_isoform_valid_summary.txt", sep = "\t", quote = F)
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
venn_cortex_GE <- venn.diagram(isoform_valid_cortex_GE, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of cortex vs GE isoform_valids", main.cex = 2)
pdf("venn_cortex_GE_isoform_valid.pdf")
plot.new()
grid.draw(venn_cortex_GE)
dev.off()

# HuFNSC01 vs HuFNSC02
isoform_valid_HuFNSC01_HuFNSC02 <- list(brain = brain01_brain02_isoform_valid_gene$geneID, cortex = cortex01_cortex02_isoform_valid_gene$geneID, GE = GE01_GE02_isoform_valid_gene$geneID)
venn_HuFNSC01_HuFNSC02 <- venn.diagram(isoform_valid_HuFNSC01_HuFNSC02, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of HuFNSC01 vs HuFNSC02 isoform_valids", main.cex = 2)
pdf("venn_HuFNSC01_HuFNSC02_isoform_valid.pdf")
plot.new()
grid.draw(venn_HuFNSC01_HuFNSC02)
dev.off()

# HuFNSC03 vs HuFNSC04
isoform_valid_HuFNSC03_HuFNSC04 <- list(cortex = cortex03_cortex04_isoform_valid_gene$geneID, GE = GE03_GE04_isoform_valid_gene$geneID)
venn_HuFNSC03_HuFNSC04 <- venn.diagram(isoform_valid_HuFNSC03_HuFNSC04, filename = NULL, fill = c("blue", "green"), main = "Venn diagram of HuFNSC03 vs HuFNSC04 isoform_valids", main.cex = 2)
pdf("venn_HuFNSC03_HuFNSC04_isoform_valid.pdf")
plot.new()
grid.draw(venn_HuFNSC03_HuFNSC04)
dev.off()

save.image("FetalBrain_isoform_valid.Rdata")


