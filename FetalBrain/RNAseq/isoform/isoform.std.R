##################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
# isoforms: DEfine on exons & exonRPKM < cutoff*geneRPKM for one sample & > cutoff2 for the other & expressed in both samples (geneRPKM > RPKMmin) & gene not DE
# cutoff = 0.01; cutoff2 = 0.1; RPKMmin = 0.01
# DEfine: FDR = 0.01; rmin = 0.005; Nmin = 25
##################################################################################
setwd("~/FetalBrain/RNAseq/isoform/")
source("/home/lli/HirstLab/Pipeline/R/isoform.R")

# isoform between cortex and GE
lib1='A03473'; cell1='Cortex'; donor1='HuFNSC01';
lib2='A03474'; cell2='GE'; donor2='HuFNSC01';
cortex01_GE01 <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2)
cortex01_GE01_summary <- cortex01_GE01$summary
cortex01_GE01_isoform_exon <- cortex01_GE01$isoform_exon
cortex01_GE01_isoform_gene <- cortex01_GE01$isoform_gene
rm(cortex01_GE01, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A03475'; cell1='Cortex'; donor1='HuFNSC02';
lib2='A03476'; cell2='GE'; donor2='HuFNSC02';
cortex02_GE02 <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2)
cortex02_GE02_summary <- cortex02_GE02$summary
cortex02_GE02_isoform_exon <- cortex02_GE02$isoform_exon
cortex02_GE02_isoform_gene <- cortex02_GE02$isoform_gene
rm(cortex02_GE02, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A04599'; cell1='Cortex'; donor1='HuFNSC03';
lib2='A15295'; cell2='GE'; donor2='HuFNSC03';
cortex03_GE03 <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2)
cortex03_GE03_summary <- cortex03_GE03$summary
cortex03_GE03_isoform_exon <- cortex03_GE03$isoform_exon
cortex03_GE03_isoform_gene <- cortex03_GE03$isoform_gene
rm(cortex03_GE03, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A15298'; cell1='Cortex'; donor1='HuFNSC04';
lib2='A15299'; cell2='GE'; donor2='HuFNSC04';
cortex04_GE04 <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2)
cortex04_GE04_summary <- cortex04_GE04$summary
cortex04_GE04_isoform_exon <- cortex04_GE04$isoform_exon
cortex04_GE04_isoform_gene <- cortex04_GE04$isoform_gene
rm(cortex04_GE04, lib1, lib2, cell1, cell2, donor1, donor2)

cortex_GE_summary <- rbind(cortex01_GE01_summary, cortex02_GE02_summary, cortex03_GE03_summary, cortex04_GE04_summary)
rm(cortex01_GE01_summary, cortex02_GE02_summary, cortex03_GE03_summary, cortex04_GE04_summary)
write.table(cortex_GE_summary, file = "cortex_GE_isoform_summary.txt", sep = "\t", quote = F)

# isoform between individuals
lib1='A03484'; cell1='Brain'; donor1='HuFNSC01';
lib2='A07825'; cell2='Brain'; donor2='HuFNSC02';
brain01_brain02 <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2)
brain01_brain02_summary <- brain01_brain02$summary
brain01_brain02_isoform_exon <- brain01_brain02$isoform_exon
brain01_brain02_isoform_gene <- brain01_brain02$isoform_gene
rm(brain01_brain02, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A03473'; cell1='Cortex'; donor1='HuFNSC01';
lib2='A03475'; cell2='Cortex'; donor2='HuFNSC02';
cortex01_cortex02 <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2)
cortex01_cortex02_summary <- cortex01_cortex02$summary
cortex01_cortex02_isoform_exon <- cortex01_cortex02$isoform_exon
cortex01_cortex02_isoform_gene <- cortex01_cortex02$isoform_gene
rm(cortex01_cortex02, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A03474'; cell1='GE'; donor1='HuFNSC01';
lib2='A03476'; cell2='GE'; donor2='HuFNSC02';
GE01_GE02 <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2)
GE01_GE02_summary <- GE01_GE02$summary
GE01_GE02_isoform_exon <- GE01_GE02$isoform_exon
GE01_GE02_isoform_gene <- GE01_GE02$isoform_gene
rm(GE01_GE02, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A04599'; cell1='Cortex'; donor1='HuFNSC03';
lib2='A15298'; cell2='Cortex'; donor2='HuFNSC04';
cortex03_cortex04 <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2)
cortex03_cortex04_summary <- cortex03_cortex04$summary
cortex03_cortex04_isoform_exon <- cortex03_cortex04$isoform_exon
cortex03_cortex04_isoform_gene <- cortex03_cortex04$isoform_gene
rm(cortex03_cortex04, lib1, lib2, cell1, cell2, donor1, donor2)

lib1='A15295'; cell1='GE'; donor1='HuFNSC03';
lib2='A15299'; cell2='GE'; donor2='HuFNSC04';
GE03_GE04 <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2)
GE03_GE04_summary <- GE03_GE04$summary
GE03_GE04_isoform_exon <- GE03_GE04$isoform_exon
GE03_GE04_isoform_gene <- GE03_GE04$isoform_gene
rm(GE03_GE04, lib1, lib2, cell1, cell2, donor1, donor2)

individual_summary <- rbind(brain01_brain02_summary, cortex01_cortex02_summary, GE01_GE02_summary, cortex03_cortex04_summary, GE03_GE04_summary)
rm(brain01_brain02_summary, cortex01_cortex02_summary, GE01_GE02_summary, cortex03_cortex04_summary, GE03_GE04_summary)
write.table(individual_summary, file = "individual_isoform_summary.txt", sep = "\t", quote = F)
save.image("FetalBrain_isoform.Rdata")
rm(list = ls())

##################################################################################
setwd("~/快盘/FetalBrain/RNAseq/isoform/")
load("FetalBrain_isoform.Rdata")
source("~/HirstLab/Pipeline/enrich.R")
library(VennDiagram)
library(ggplot2)
load("~/快盘/hg19/hg19v65_genes.Rdata")

# cortex vs GE
isoform_cortex_GE <- list(HuFNSC01 = cortex01_GE01_isoform_gene$id, HuFNSC02 = cortex02_GE02_isoform_gene$id, HuFNSC03 = cortex03_GE03_isoform_gene$id, HuFNSC04 = cortex04_GE04_isoform_gene$id)
venn_cortex_GE <- venn.diagram(isoform_cortex_GE, filename = NULL, fill = c("red", "blue", "green", "purple"), main = "Venn diagram of cortex vs GE isoforms", main.cex = 2)
pdf("venn_cortex_GE_isoform.pdf")
plot.new()
grid.draw(venn_cortex_GE)
dev.off()
isoform_cortex_GE_gene <- c(cortex01_GE01_isoform_gene$id, cortex02_GE02_isoform_gene$id, cortex03_GE03_isoform_gene$id, cortex04_GE04_isoform_gene$id)
isoform_cortex_GE_gene_dup <- ensembl[unique(isoform_cortex_GE_gene[duplicated(isoform_cortex_GE_gene)]), ] 
isoform_cortex_GE_gene <- unique(isoform_cortex_GE_gene)
isoform_cortex_GE_exon <- c(cortex01_GE01_isoform_exon$V1, cortex02_GE02_isoform_exon$V1, cortex03_GE03_isoform_exon$V1, cortex04_GE04_isoform_exon$V1)
isoform_cortex_GE_exon_dup <-isoform_cortex_GE_exon[duplicated(isoform_cortex_GE_exon)] 
isoform_cortex_GE_exon <- unique(isoform_cortex_GE_exon)
write.table(isoform_cortex_GE_gene_dup, file = "./enrich/isoform_cortex_GE_gene_dup.txt", sep = "\t", quote = F, row.names = F, col.names = F)
isoform_cortex_GE_enrich <- enrich(name = "isoform_cortex_GE", fdr = 0.01, p = "FDR", erminej = F, height = 7, width = 9)
(isoform_cortex_GE_enrich$figure)

# HuFNSC01 vs HuFNSC02
isoform_HuFNSC01_HuFNSC02 <- list(brain = brain01_brain02_isoform_gene$id, cortex = cortex01_cortex02_isoform_gene$id, GE = GE01_GE02_isoform_gene$id)
venn_HuFNSC01_HuFNSC02 <- venn.diagram(isoform_HuFNSC01_HuFNSC02, filename = NULL, fill = c("green", "red", "blue"), main = "Venn diagram of HuFNSC01 vs HuFNSC02 isoforms", main.cex = 2)
pdf("venn_HuFNSC01_HuFNSC02_isoform.pdf")
plot.new()
grid.draw(venn_HuFNSC01_HuFNSC02)
dev.off()
isoform_HuFNSC01_HuFNSC02_gene_shared <- intersect(brain01_brain02_isoform_gene$id, intersect(cortex01_cortex02_isoform_gene$id, GE01_GE02_isoform_gene$id))
write.table(ensembl[brain01_brain02_isoform_gene$id, ], file = "./enrich/isoform_brain01_brain02.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(ensembl[cortex01_cortex02_isoform_gene$id, ], file = "./enrich/isoform_cortex01_cortex02.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(ensembl[GE01_GE02_isoform_gene$id, ], file = "./enrich/isoform_GE01_GE02.txt", sep = "\t", quote = F, row.names = F, col.names = F)
isoform_brain01_brain02_enrich <- enrich(name = "isoform_brain01_brain02", fdr = 0.01, p = "FDR", erminej = F, width = 9, height = 6)
(isoform_brain01_brain02_enrich$figure)
isoform_cortex01_cortex02_enrich <- enrich(name = "isoform_cortex01_cortex02", fdr = 0.01, p = "FDR", erminej = F, height = 6, width = 9)
(isoform_cortex01_cortex02_enrich$figure)
isoform_GE01_GE02_enrich <- enrich(name = "isoform_GE01_GE02", fdr = 0.01, p = "FDR", erminej = F, height = 7, width = 9)
(isoform_GE01_GE02_enrich$figure)

# HuFNSC03 vs HuFNSC04
isoform_HuFNSC03_HuFNSC04 <- list(cortex = cortex03_cortex04_isoform_gene$id, GE = GE03_GE04_isoform_gene$id)
venn_HuFNSC03_HuFNSC04 <- venn.diagram(isoform_HuFNSC03_HuFNSC04, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of HuFNSC03 vs HuFNSC04 isoforms", main.cex = 2)
pdf("venn_HuFNSC03_HuFNSC04_isoform.pdf")
plot.new()
grid.draw(venn_HuFNSC03_HuFNSC04)
dev.off()
isoform_HuFNSC03_HuFNSC04_gene_shared <- intersect(cortex03_cortex04_isoform_gene$id, GE03_GE04_isoform_gene$id)
write.table(ensembl[cortex03_cortex04_isoform_gene$id, ], file = "./enrich/isoform_cortex03_cortex04.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(ensembl[GE03_GE04_isoform_gene$id, ], file = "./enrich/isoform_GE03_GE04.txt", sep = "\t", quote = F, row.names = F, col.names = F)
isoform_cortex03_cortex04_enrich <- enrich(name = "isoform_cortex03_cortex04", fdr = 0.01, p = "FDR", erminej = F, width = 9, height = 4)
(isoform_cortex03_cortex04_enrich$figure)
isoform_GE03_GE04_enrich <- enrich(name = "isoform_GE03_GE04", fdr = 0.01, p = "FDR", erminej = F, width = 9, height = 3)
(isoform_GE03_GE04_enrich$figure)

rm(ensembl)
save.image("FetalBrain_isoform.Rdata")


