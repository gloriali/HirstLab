##################################################################################
# run on xhost
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
setwd("~/FetalBrain/RNAseq/isoform/")
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'
libs=c("A03484", "A07825", 
       "A03473", "A03475", "A04599", "A15298",
       "A03474", "A03476", "A15295", "A15299")
colnames <- c("brain01", "brain02", 
              "cortex01", "cortex02", "cortex03", "cortex04", 
              "ge01", "ge02", "ge03", "ge04")

lib <- libs[1]
exon_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.exn.A.rpkm",sep=""))
pc_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.A.rpkm.pc",sep=""),row.names=1)
fetalBrain_gene <- as.data.frame(matrix(data = NA, nrow = nrow(pc_lib), ncol = length(libs), dimnames = list(rownames(pc_lib), colnames)))
fetalBrain_exon <- as.data.frame(matrix(data = NA, nrow = nrow(exon_lib), ncol = length(libs), dimnames = list(paste0(exon_lib$V1, "_", exon_lib$V2), colnames)))

pdf("percent_ave.pdf")
plot(range(0,0.1),range(0,0.3),type='n',main="ECDF percentage of average RPKM",xlab="percentage",ylab="cdf")
i=1
for(lib in libs){
  exon_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.exn.A.rpkm",sep=""))
  pc_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.A.rpkm.pc",sep=""),row.names=1)
  fetalBrain_gene[, i] <- pc_lib$V3
  fetalBrain_exon[, i] <- exon_lib$V4
  exon_lib$ave<-pc_lib[exon_lib$V2,3]
  exon_lib<-exon_lib[exon_lib$ave!=0,]
  exon_lib$percent<-exon_lib$V4/exon_lib$ave
  lines(ecdf(exon_lib$percent[exon_lib$percent!=0]),col=i,lty=1)
  i=i+1
}
legend("bottomright",libs,col=c(1:length(libs)),cex=0.8,lty=1)
dev.off()
save(fetalBrain_gene, fetalBrain_exon, file = "fetalBrain_rpkm.Rdata")

##################################################################################
##################################################################################

# isoforms:DEfine on exons 
# & expressed in both samples(gene RPKM > Rmin)
# & RPKM(exon) <= 1%*RPKM(ave. exons of the gene) for one sample& >= 10% for the other 
# & gene not DE
setwd("~/FetalBrain/RNAseq/isoform/")
load("fetalBrain_rpkm.Rdata")
col <- c("Gene id", "rpkm1", "rpkm2", "p-value", "Multiple testing corrected p-value")
cutoff = 0.01
cutoff2 = 0.1 # sample specific exon: expressed in one sample(>0.1ave.RPKM) & not expressed in the other(<0.01ave.RPKM)
Rmin = 0.005

##################################################################################
## tissue specific isoforms: cortex vs GE
genedir <- "~/FetalBrain/RNAseq/DEfine/gene/cortexge/"
exondir <- "~/FetalBrain/RNAseq/DEfine/exon/"
isoform_cortexge <- matrix(data = NA, nrow = 4, ncol = 6, dimnames = list(c("HuFNSC01", "HuFNSC02", "HuFNSC03", "HuFNSC04"), c("DE_exons", "pc_DE_exons", "gene_expressed", "gene_not_DE", "isoform_exons", "isoform_genes")))

donor1 = "HuFNSC01"; cell1 = "Cortex"; 
donor2 = "HuFNSC01"; cell2 = "GE"; 
cortex_ge_01up <- read.table(paste(exondir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex_ge_01dn <- read.table(paste(exondir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex_ge_01 <- rbind(data.frame(cortex_ge_01up, DE = "up"), data.frame(cortex_ge_01dn, DE = "dn"))
DE_cortex_ge_01up <- read.table(paste(genedir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex_ge_01dn <- read.table(paste(genedir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex_ge_01 <- c(as.character(DE_cortex_ge_01up$V1),as.character(DE_cortex_ge_01dn$V1))
rm(cortex_ge_01up, cortex_ge_01dn, DE_cortex_ge_01up, DE_cortex_ge_01dn)
cortex_ge_01$id <- unlist(strsplit(as.character(cortex_ge_01$V1), '_'))[2*(1:nrow(cortex_ge_01))]
cortex_ge_01$ave2 <- fetalBrain_gene[cortex_ge_01$id, "cortex01"]
cortex_ge_01$ave3 <- fetalBrain_gene[cortex_ge_01$id, "ge01"]
(isoform_cortexge["HuFNSC01", "DE_exons"] <- nrow(cortex_ge_01))
cortex_ge_01 <- na.omit(cortex_ge_01)
(isoform_cortexge["HuFNSC01", "pc_DE_exons"] <- nrow(cortex_ge_01))
cortex_ge_01 <- cortex_ge_01[cortex_ge_01$ave2 > Rmin & cortex_ge_01$ave3 > Rmin, ]
(isoform_cortexge["HuFNSC01", "gene_expressed"] <- nrow(cortex_ge_01))
cortex_ge_01 <- cortex_ge_01[!(cortex_ge_01$id %in% DE_cortex_ge_01), ]
(isoform_cortexge["HuFNSC01", "gene_not_DE"] <- nrow(cortex_ge_01))
cortex_ge_01_isoform <- cortex_ge_01[((cortex_ge_01$V2 <= cutoff * cortex_ge_01$ave2) & (cortex_ge_01$V3 >= cutoff2 * cortex_ge_01$ave3))|
                               ((cortex_ge_01$V2 >= cutoff2 * cortex_ge_01$ave2) & (cortex_ge_01$V3 <= cutoff * cortex_ge_01$ave3)),]
(isoform_cortexge["HuFNSC01", "isoform_exons"] <- nrow(cortex_ge_01_isoform))
cortex_ge_01_isoform_gene <- cortex_ge_01_isoform[!duplicated(cortex_ge_01_isoform$id), ]
(isoform_cortexge["HuFNSC01", "isoform_genes"] <- nrow(cortex_ge_01_isoform_gene))
write.table(cortex_ge_01_isoform, file = "cortex_ge_01_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(cortex_ge_01_isoform_gene, file = "cortex_ge_01_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

donor1 = "HuFNSC02"; cell1 = "Cortex"; 
donor2 = "HuFNSC02"; cell2 = "GE"; 
cortex_ge_02up <- read.table(paste(exondir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex_ge_02dn <- read.table(paste(exondir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex_ge_02 <- rbind(data.frame(cortex_ge_02up, DE = "up"), data.frame(cortex_ge_02dn, DE = "dn"))
DE_cortex_ge_02up <- read.table(paste(genedir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex_ge_02dn <- read.table(paste(genedir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex_ge_02 <- c(as.character(DE_cortex_ge_02up$V1),as.character(DE_cortex_ge_02dn$V1))
rm(cortex_ge_02up, cortex_ge_02dn, DE_cortex_ge_02up, DE_cortex_ge_02dn)
cortex_ge_02$id <- unlist(strsplit(as.character(cortex_ge_02$V1), '_'))[2*(1:nrow(cortex_ge_02))]
cortex_ge_02$ave2 <- fetalBrain_gene[cortex_ge_02$id, "cortex02"]
cortex_ge_02$ave3 <- fetalBrain_gene[cortex_ge_02$id, "ge02"]
(isoform_cortexge["HuFNSC02", "DE_exons"] <- nrow(cortex_ge_02))
cortex_ge_02 <- na.omit(cortex_ge_02)
(isoform_cortexge["HuFNSC02", "pc_DE_exons"] <- nrow(cortex_ge_02))
cortex_ge_02 <- cortex_ge_02[cortex_ge_02$ave2 > Rmin & cortex_ge_02$ave3 > Rmin, ]
(isoform_cortexge["HuFNSC02", "gene_expressed"] <- nrow(cortex_ge_02))
cortex_ge_02 <- cortex_ge_02[!(cortex_ge_02$id %in% DE_cortex_ge_02), ]
(isoform_cortexge["HuFNSC02", "gene_not_DE"] <- nrow(cortex_ge_02))
cortex_ge_02_isoform <- cortex_ge_02[((cortex_ge_02$V2 <= cutoff * cortex_ge_02$ave2) & (cortex_ge_02$V3 >= cutoff2 * cortex_ge_02$ave3))|
                                       ((cortex_ge_02$V2 >= cutoff2 * cortex_ge_02$ave2) & (cortex_ge_02$V3 <= cutoff * cortex_ge_02$ave3)),]
(isoform_cortexge["HuFNSC02", "isoform_exons"] <- nrow(cortex_ge_02_isoform))
cortex_ge_02_isoform_gene <- cortex_ge_02_isoform[!duplicated(cortex_ge_02_isoform$id), ]
(isoform_cortexge["HuFNSC02", "isoform_genes"] <- nrow(cortex_ge_02_isoform_gene))
write.table(cortex_ge_02_isoform, file = "cortex_ge_02_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(cortex_ge_02_isoform_gene, file = "cortex_ge_02_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

donor1 = "HuFNSC03"; cell1 = "Cortex"; 
donor2 = "HuFNSC03"; cell2 = "GE"; 
cortex_ge_03up <- read.table(paste(exondir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex_ge_03dn <- read.table(paste(exondir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex_ge_03 <- rbind(data.frame(cortex_ge_03up, DE = "up"), data.frame(cortex_ge_03dn, DE = "dn"))
DE_cortex_ge_03up <- read.table(paste(genedir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex_ge_03dn <- read.table(paste(genedir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex_ge_03 <- c(as.character(DE_cortex_ge_03up$V1),as.character(DE_cortex_ge_03dn$V1))
rm(cortex_ge_03up, cortex_ge_03dn, DE_cortex_ge_03up, DE_cortex_ge_03dn)
cortex_ge_03$id <- unlist(strsplit(as.character(cortex_ge_03$V1), '_'))[2*(1:nrow(cortex_ge_03))]
cortex_ge_03$ave2 <- fetalBrain_gene[cortex_ge_03$id, "cortex03"]
cortex_ge_03$ave3 <- fetalBrain_gene[cortex_ge_03$id, "ge03"]
(isoform_cortexge["HuFNSC03", "DE_exons"] <- nrow(cortex_ge_03))
cortex_ge_03 <- na.omit(cortex_ge_03)
(isoform_cortexge["HuFNSC03", "pc_DE_exons"] <- nrow(cortex_ge_03))
cortex_ge_03 <- cortex_ge_03[cortex_ge_03$ave2 > Rmin & cortex_ge_03$ave3 > Rmin, ]
(isoform_cortexge["HuFNSC03", "gene_expressed"] <- nrow(cortex_ge_03))
cortex_ge_03 <- cortex_ge_03[!(cortex_ge_03$id %in% DE_cortex_ge_03), ]
(isoform_cortexge["HuFNSC03", "gene_not_DE"] <- nrow(cortex_ge_03))
cortex_ge_03_isoform <- cortex_ge_03[((cortex_ge_03$V2 <= cutoff * cortex_ge_03$ave2) & (cortex_ge_03$V3 >= cutoff2 * cortex_ge_03$ave3))|
                                       ((cortex_ge_03$V2 >= cutoff2 * cortex_ge_03$ave2) & (cortex_ge_03$V3 <= cutoff * cortex_ge_03$ave3)),]
(isoform_cortexge["HuFNSC03", "isoform_exons"] <- nrow(cortex_ge_03_isoform))
cortex_ge_03_isoform_gene <- cortex_ge_03_isoform[!duplicated(cortex_ge_03_isoform$id), ]
(isoform_cortexge["HuFNSC03", "isoform_genes"] <- nrow(cortex_ge_03_isoform_gene))
write.table(cortex_ge_03_isoform, file = "cortex_ge_03_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(cortex_ge_03_isoform_gene, file = "cortex_ge_03_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

donor1 = "HuFNSC04"; cell1 = "Cortex"; 
donor2 = "HuFNSC04"; cell2 = "GE"; 
cortex_ge_04up <- read.table(paste(exondir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex_ge_04dn <- read.table(paste(exondir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex_ge_04 <- rbind(data.frame(cortex_ge_04up, DE = "up"), data.frame(cortex_ge_04dn, DE = "dn"))
DE_cortex_ge_04up <- read.table(paste(genedir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex_ge_04dn <- read.table(paste(genedir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex_ge_04 <- c(as.character(DE_cortex_ge_04up$V1),as.character(DE_cortex_ge_04dn$V1))
rm(cortex_ge_04up, cortex_ge_04dn, DE_cortex_ge_04up, DE_cortex_ge_04dn)
cortex_ge_04$id <- unlist(strsplit(as.character(cortex_ge_04$V1), '_'))[2*(1:nrow(cortex_ge_04))]
cortex_ge_04$ave2 <- fetalBrain_gene[cortex_ge_04$id, "cortex04"]
cortex_ge_04$ave3 <- fetalBrain_gene[cortex_ge_04$id, "ge04"]
(isoform_cortexge["HuFNSC04", "DE_exons"] <- nrow(cortex_ge_04))
cortex_ge_04 <- na.omit(cortex_ge_04)
(isoform_cortexge["HuFNSC04", "pc_DE_exons"] <- nrow(cortex_ge_04))
cortex_ge_04 <- cortex_ge_04[cortex_ge_04$ave2 > Rmin & cortex_ge_04$ave3 > Rmin, ]
(isoform_cortexge["HuFNSC04", "gene_expressed"] <- nrow(cortex_ge_04))
cortex_ge_04 <- cortex_ge_04[!(cortex_ge_04$id %in% DE_cortex_ge_04), ]
(isoform_cortexge["HuFNSC04", "gene_not_DE"] <- nrow(cortex_ge_04))
cortex_ge_04_isoform <- cortex_ge_04[((cortex_ge_04$V2 <= cutoff * cortex_ge_04$ave2) & (cortex_ge_04$V3 >= cutoff2 * cortex_ge_04$ave3))|
                                       ((cortex_ge_04$V2 >= cutoff2 * cortex_ge_04$ave2) & (cortex_ge_04$V3 <= cutoff * cortex_ge_04$ave3)),]
(isoform_cortexge["HuFNSC04", "isoform_exons"] <- nrow(cortex_ge_04_isoform))
cortex_ge_04_isoform_gene <- cortex_ge_04_isoform[!duplicated(cortex_ge_04_isoform$id), ]
(isoform_cortexge["HuFNSC04", "isoform_genes"] <- nrow(cortex_ge_04_isoform_gene))
write.table(cortex_ge_04_isoform, file = "cortex_ge_04_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(cortex_ge_04_isoform_gene, file = "cortex_ge_04_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

cortex_ge_isoform_gene <- rbind(cortex_ge_01_isoform_gene, cortex_ge_02_isoform_gene, cortex_ge_03_isoform_gene, cortex_ge_04_isoform_gene)
cortex_ge_isoform_gene <- cortex_ge_isoform_gene[duplicated(cortex_ge_isoform_gene$id), ]
cortex_ge_isoform_gene <- cortex_ge_isoform_gene[!duplicated(cortex_ge_isoform_gene$id), ]
cortex_ge_isoform <- rbind(cortex_ge_01_isoform, cortex_ge_02_isoform, cortex_ge_03_isoform, cortex_ge_04_isoform)
cortex_ge_isoform <- cortex_ge_isoform[(cortex_ge_isoform$id %in% cortex_ge_isoform_gene$id), ]
(nrow(cortex_ge_isoform))
(nrow(cortex_ge_isoform_gene))
write.table(cortex_ge_isoform, file = "cortex_ge_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(cortex_ge_isoform_gene, file = "cortex_ge_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
cortex_ge_01only_isoform_gene <- cortex_ge_01_isoform_gene[cortex_ge_01_isoform_gene$id %in% setdiff(cortex_ge_01_isoform_gene$id, c(cortex_ge_02_isoform_gene$id, cortex_ge_03_isoform_gene$id, cortex_ge_04_isoform_gene$id)), ]
write.table(cortex_ge_01only_isoform_gene, file = "cortex_ge_01only_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
cortex_ge_02only_isoform_gene <- cortex_ge_02_isoform_gene[cortex_ge_02_isoform_gene$id %in% setdiff(cortex_ge_02_isoform_gene$id, c(cortex_ge_01_isoform_gene$id, cortex_ge_03_isoform_gene$id, cortex_ge_04_isoform_gene$id)), ]
write.table(cortex_ge_02only_isoform_gene, file = "cortex_ge_02only_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
cortex_ge_03only_isoform_gene <- cortex_ge_03_isoform_gene[cortex_ge_03_isoform_gene$id %in% setdiff(cortex_ge_03_isoform_gene$id, c(cortex_ge_01_isoform_gene$id, cortex_ge_02_isoform_gene$id, cortex_ge_04_isoform_gene$id)), ]
write.table(cortex_ge_03only_isoform_gene, file = "cortex_ge_03only_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
cortex_ge_04only_isoform_gene <- cortex_ge_04_isoform_gene[cortex_ge_04_isoform_gene$id %in% setdiff(cortex_ge_04_isoform_gene$id, c(cortex_ge_01_isoform_gene$id, cortex_ge_02_isoform_gene$id, cortex_ge_03_isoform_gene$id)), ]
write.table(cortex_ge_04only_isoform_gene, file = "cortex_ge_04only_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

rm(fetalBrain_gene, fetalBrain_exon, cell1, cell2, donor1, donor2, cortex_ge_01, cortex_ge_02, cortex_ge_03, cortex_ge_04)
save.image(file = "fetalBrain_isoform_tissue.Rdata")

library(VennDiagram)
library(venneuler)
isoform_gene_tissue <- list(HuFNSC01 = cortex_ge_01_isoform_gene$id, HuFNSC02 = cortex_ge_02_isoform_gene$id, HuFNSC03 = cortex_ge_03_isoform_gene$id, HuFNSC04 = cortex_ge_04_isoform_gene$id)
pdf("~/快盘/FetalBrain/RNAseq/isoform/venn_isoform_gene_tissue.pdf")
plot.new()
venn_isoform_gene_tissue <- venn.diagram(isoform_gene_tissue, filename = NULL, fill = c("red", "blue", "green", "yellow"), main = "Venn diagram of cortex vs GE isoforms", main.cex = 2)
grid.draw(venn_isoform_gene_tissue)

venneuler_isoform_gene_tissue <- venneuler(data.frame(elements = c(cortex_ge_01_isoform_gene$id, cortex_ge_02_isoform_gene$id, cortex_ge_03_isoform_gene$id, cortex_ge_04_isoform_gene$id), 
                                                 sets = c(rep("HuFNSC01", nrow(cortex_ge_01_isoform_gene)), rep("HuFNSC02", nrow(cortex_ge_02_isoform_gene)), rep("HuFNSC03", nrow(cortex_ge_03_isoform_gene)), rep("HuFNSC04", nrow(cortex_ge_04_isoform_gene)))))
plot(venneuler_isoform_gene_tissue)
mtext("Venn diagram of cortex vs GE isoforms", side = 3)
dev.off()

##################################################################################
##################################################################################

# isoforms:DEfine on exons 
# & expressed in both samples(gene RPKM > Rmin)
# & RPKM(exon) < 1%*RPKM(ave. exons of the gene) for one sample& > 10% for the other 
# & gene not DE
setwd("~/FetalBrain/RNAseq/isoform/")
load("fetalBrain_rpkm.Rdata")
col <- c("Gene id", "rpkm1", "rpkm2", "p-value", "Multiple testing corrected p-value")
cutoff = 0.01
cutoff2 = 0.1 # sample specific exon: expressed in one sample(>0.1ave.RPKM) & not expressed in the other(<0.01ave.RPKM)
Rmin = 0.005

##################################################################################
## individual specific isoforms: HuFNSC01 vs HuFNSC02 & HuFNSC03 vs HuFNSC04
genedir <- "~/FetalBrain/RNAseq/DEfine/gene/individual/"
exondir <- "~/FetalBrain/RNAseq/DEfine/exon/"
isoform_individual <- matrix(data = NA, nrow = 5, ncol = 6, dimnames = list(c("brain01_02", "cortex01_02", "GE01_02", "cortex03_04", "GE03_04"), c("DE_exons", "pc_DE_exons", "gene_expressed", "gene_not_DE", "isoform_exons", "isoform_genes")))

donor1 = "HuFNSC01"; cell1 = "Brain"; 
donor2 = "HuFNSC02"; cell2 = "Brain"; 
brain01_02up <- read.table(paste(exondir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
brain01_02dn <- read.table(paste(exondir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
brain01_02 <- rbind(data.frame(brain01_02up, DE = "up"), data.frame(brain01_02dn, DE = "dn"))
DE_brain01_02up <- read.table(paste(genedir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_brain01_02dn <- read.table(paste(genedir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_brain01_02 <- c(as.character(DE_brain01_02up$V1),as.character(DE_brain01_02dn$V1))
rm(brain01_02up, brain01_02dn, DE_brain01_02up, DE_brain01_02dn)
brain01_02$id <- unlist(strsplit(as.character(brain01_02$V1), '_'))[2*(1:nrow(brain01_02))]
brain01_02$ave2 <- fetalBrain_gene[brain01_02$id, "brain01"]
brain01_02$ave3 <- fetalBrain_gene[brain01_02$id, "brain02"]
(isoform_individual["brain01_02", "DE_exons"] <- nrow(brain01_02))
brain01_02 <- na.omit(brain01_02)
(isoform_individual["brain01_02", "pc_DE_exons"] <- nrow(brain01_02))
brain01_02 <- brain01_02[brain01_02$ave2 > Rmin & brain01_02$ave3 > Rmin, ]
(isoform_individual["brain01_02", "gene_expressed"] <- nrow(brain01_02))
brain01_02 <- brain01_02[!(brain01_02$id %in% DE_brain01_02), ]
(isoform_individual["brain01_02", "gene_not_DE"] <- nrow(brain01_02))
brain01_02_isoform <- brain01_02[((brain01_02$V2 <= cutoff * brain01_02$ave2) & (brain01_02$V3 >= cutoff2 * brain01_02$ave3))|
                                   ((brain01_02$V2 >= cutoff2 * brain01_02$ave2) & (brain01_02$V3 <= cutoff * brain01_02$ave3)),]
(isoform_individual["brain01_02", "isoform_exons"] <- nrow(brain01_02_isoform))
brain01_02_isoform_gene <- brain01_02_isoform[!duplicated(brain01_02_isoform$id), ]
(isoform_individual["brain01_02", "isoform_genes"] <- nrow(brain01_02_isoform_gene))
write.table(brain01_02_isoform, file = "brain01_02_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(brain01_02_isoform_gene, file = "brain01_02_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

donor1 = "HuFNSC01"; cell1 = "Cortex"; 
donor2 = "HuFNSC02"; cell2 = "Cortex"; 
cortex01_02up <- read.table(paste(exondir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex01_02dn <- read.table(paste(exondir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex01_02 <- rbind(data.frame(cortex01_02up, DE = "up"), data.frame(cortex01_02dn, DE = "dn"))
DE_cortex01_02up <- read.table(paste(genedir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex01_02dn <- read.table(paste(genedir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex01_02 <- c(as.character(DE_cortex01_02up$V1),as.character(DE_cortex01_02dn$V1))
rm(cortex01_02up, cortex01_02dn, DE_cortex01_02up, DE_cortex01_02dn)
cortex01_02$id <- unlist(strsplit(as.character(cortex01_02$V1), '_'))[2*(1:nrow(cortex01_02))]
cortex01_02$ave2 <- fetalBrain_gene[cortex01_02$id, "cortex01"]
cortex01_02$ave3 <- fetalBrain_gene[cortex01_02$id, "cortex02"]
(isoform_individual["cortex01_02", "DE_exons"] <- nrow(cortex01_02))
cortex01_02 <- na.omit(cortex01_02)
(isoform_individual["cortex01_02", "pc_DE_exons"] <- nrow(cortex01_02))
cortex01_02 <- cortex01_02[cortex01_02$ave2 > Rmin & cortex01_02$ave3 > Rmin, ]
(isoform_individual["cortex01_02", "gene_expressed"] <- nrow(cortex01_02))
cortex01_02 <- cortex01_02[!(cortex01_02$id %in% DE_cortex01_02), ]
(isoform_individual["cortex01_02", "gene_not_DE"] <- nrow(cortex01_02))
cortex01_02_isoform <- cortex01_02[((cortex01_02$V2 <= cutoff * cortex01_02$ave2) & (cortex01_02$V3 >= cutoff2 * cortex01_02$ave3))|
                                     ((cortex01_02$V2 >= cutoff2 * cortex01_02$ave2) & (cortex01_02$V3 <= cutoff * cortex01_02$ave3)),]
(isoform_individual["cortex01_02", "isoform_exons"] <- nrow(cortex01_02_isoform))
cortex01_02_isoform_gene <- cortex01_02_isoform[!duplicated(cortex01_02_isoform$id), ]
(isoform_individual["cortex01_02", "isoform_genes"] <- nrow(cortex01_02_isoform_gene))
write.table(cortex01_02_isoform, file = "cortex01_02_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(cortex01_02_isoform_gene, file = "cortex01_02_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

donor1 = "HuFNSC01"; cell1 = "GE"; 
donor2 = "HuFNSC02"; cell2 = "GE"; 
ge01_02up <- read.table(paste(exondir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
ge01_02dn <- read.table(paste(exondir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
ge01_02 <- rbind(data.frame(ge01_02up, DE = "up"), data.frame(ge01_02dn, DE = "dn"))
DE_ge01_02up <- read.table(paste(genedir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_ge01_02dn <- read.table(paste(genedir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_ge01_02 <- c(as.character(DE_ge01_02up$V1),as.character(DE_ge01_02dn$V1))
rm(ge01_02up, ge01_02dn, DE_ge01_02up, DE_ge01_02dn)
ge01_02$id <- unlist(strsplit(as.character(ge01_02$V1), '_'))[2*(1:nrow(ge01_02))]
ge01_02$ave2 <- fetalBrain_gene[ge01_02$id, "ge01"]
ge01_02$ave3 <- fetalBrain_gene[ge01_02$id, "ge02"]
(isoform_individual["GE01_02", "DE_exons"] <- nrow(ge01_02))
ge01_02 <- na.omit(ge01_02)
(isoform_individual["GE01_02", "pc_DE_exons"] <- nrow(ge01_02))
ge01_02 <- ge01_02[ge01_02$ave2 > Rmin & ge01_02$ave3 > Rmin, ]
(isoform_individual["GE01_02", "gene_expressed"] <- nrow(ge01_02))
ge01_02 <- ge01_02[!(ge01_02$id %in% DE_ge01_02), ]
(isoform_individual["GE01_02", "gene_not_DE"] <- nrow(ge01_02))
ge01_02_isoform <- ge01_02[((ge01_02$V2 <= cutoff * ge01_02$ave2) & (ge01_02$V3 >= cutoff2 * ge01_02$ave3))|
                             ((ge01_02$V2 >= cutoff2 * ge01_02$ave2) & (ge01_02$V3 <= cutoff * ge01_02$ave3)),]
(isoform_individual["GE01_02", "isoform_exons"] <- nrow(ge01_02_isoform))
ge01_02_isoform_gene <- ge01_02_isoform[!duplicated(ge01_02_isoform$id), ]
(isoform_individual["GE01_02", "isoform_genes"] <- nrow(ge01_02_isoform_gene))
write.table(ge01_02_isoform, file = "ge01_02_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(ge01_02_isoform_gene, file = "ge01_02_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

donor1 = "HuFNSC03"; cell1 = "Cortex"; 
donor2 = "HuFNSC04"; cell2 = "Cortex"; 
cortex03_04up <- read.table(paste(exondir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex03_04dn <- read.table(paste(exondir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
cortex03_04 <- rbind(data.frame(cortex03_04up, DE = "up"), data.frame(cortex03_04dn, DE = "dn"))
DE_cortex03_04up <- read.table(paste(genedir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex03_04dn <- read.table(paste(genedir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_cortex03_04 <- c(as.character(DE_cortex03_04up$V1),as.character(DE_cortex03_04dn$V1))
rm(cortex03_04up, cortex03_04dn, DE_cortex03_04up, DE_cortex03_04dn)
cortex03_04$id <- unlist(strsplit(as.character(cortex03_04$V1), '_'))[2*(1:nrow(cortex03_04))]
cortex03_04$ave2 <- fetalBrain_gene[cortex03_04$id, "cortex01"]
cortex03_04$ave3 <- fetalBrain_gene[cortex03_04$id, "cortex02"]
(isoform_individual["cortex03_04", "DE_exons"] <- nrow(cortex03_04))
cortex03_04 <- na.omit(cortex03_04)
(isoform_individual["cortex03_04", "pc_DE_exons"] <- nrow(cortex03_04))
cortex03_04 <- cortex03_04[cortex03_04$ave2 > Rmin & cortex03_04$ave3 > Rmin, ]
(isoform_individual["cortex03_04", "gene_expressed"] <- nrow(cortex03_04))
cortex03_04 <- cortex03_04[!(cortex03_04$id %in% DE_cortex03_04), ]
(isoform_individual["cortex03_04", "gene_not_DE"] <- nrow(cortex03_04))
cortex03_04_isoform <- cortex03_04[((cortex03_04$V2 <= cutoff * cortex03_04$ave2) & (cortex03_04$V3 >= cutoff2 * cortex03_04$ave3))|
                                     ((cortex03_04$V2 >= cutoff2 * cortex03_04$ave2) & (cortex03_04$V3 <= cutoff * cortex03_04$ave3)),]
(isoform_individual["cortex03_04", "isoform_exons"] <- nrow(cortex03_04_isoform))
cortex03_04_isoform_gene <- cortex03_04_isoform[!duplicated(cortex03_04_isoform$id), ]
(isoform_individual["cortex03_04", "isoform_genes"] <- nrow(cortex03_04_isoform_gene))
write.table(cortex03_04_isoform, file = "cortex03_04_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(cortex03_04_isoform_gene, file = "cortex03_04_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

donor1 = "HuFNSC03"; cell1 = "GE"; 
donor2 = "HuFNSC04"; cell2 = "GE"; 
ge03_04up <- read.table(paste(exondir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
ge03_04dn <- read.table(paste(exondir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
ge03_04 <- rbind(data.frame(ge03_04up, DE = "up"), data.frame(ge03_04dn, DE = "dn"))
DE_ge03_04up <- read.table(paste(genedir, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_ge03_04dn <- read.table(paste(genedir, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_0.01.rmin_0.005.Nmin_25", sep=""))
DE_ge03_04 <- c(as.character(DE_ge03_04up$V1),as.character(DE_ge03_04dn$V1))
rm(ge03_04up, ge03_04dn, DE_ge03_04up, DE_ge03_04dn)
ge03_04$id <- unlist(strsplit(as.character(ge03_04$V1), '_'))[2*(1:nrow(ge03_04))]
ge03_04$ave2 <- fetalBrain_gene[ge03_04$id, "ge01"]
ge03_04$ave3 <- fetalBrain_gene[ge03_04$id, "ge02"]
(isoform_individual["GE03_04", "DE_exons"] <- nrow(ge03_04))
ge03_04 <- na.omit(ge03_04)
(isoform_individual["GE03_04", "pc_DE_exons"] <- nrow(ge03_04))
ge03_04 <- ge03_04[ge03_04$ave2 > Rmin & ge03_04$ave3 > Rmin, ]
(isoform_individual["GE03_04", "gene_expressed"] <- nrow(ge03_04))
ge03_04 <- ge03_04[!(ge03_04$id %in% DE_ge03_04), ]
(isoform_individual["GE03_04", "gene_not_DE"] <- nrow(ge03_04))
ge03_04_isoform <- ge03_04[((ge03_04$V2 <= cutoff * ge03_04$ave2) & (ge03_04$V3 >= cutoff2 * ge03_04$ave3))|
                             ((ge03_04$V2 >= cutoff2 * ge03_04$ave2) & (ge03_04$V3 <= cutoff * ge03_04$ave3)),]
(isoform_individual["GE03_04", "isoform_exons"] <- nrow(ge03_04_isoform))
ge03_04_isoform_gene <- ge03_04_isoform[!duplicated(ge03_04_isoform$id), ]
(isoform_individual["GE03_04", "isoform_genes"] <- nrow(ge03_04_isoform_gene))
write.table(ge03_04_isoform, file = "ge03_04_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(ge03_04_isoform_gene, file = "ge03_04_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

HuFNSC01_02_isoform_gene <- merge(brain01_02_isoform_gene, merge(cortex01_02_isoform_gene, ge01_02_isoform_gene, by = "id"), by = "id")
HuFNSC01_02_isoform <- rbind(brain01_02_isoform, cortex01_02_isoform, ge01_02_isoform)
HuFNSC01_02_isoform <- HuFNSC01_02_isoform[(HuFNSC01_02_isoform$id %in% HuFNSC01_02_isoform_gene$id), ]
(nrow(HuFNSC01_02_isoform))
(nrow(HuFNSC01_02_isoform_gene))
write.table(HuFNSC01_02_isoform, file = "HuFNSC01_02_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(HuFNSC01_02_isoform_gene, file = "HuFNSC01_02_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
HuFNSC01_02_brain_only_isoform_gene <- brain01_02_isoform_gene[brain01_02_isoform_gene$id %in% setdiff(brain01_02_isoform_gene$id, union(cortex01_02_isoform_gene$id, ge01_02_isoform_gene$id)),]
write.table(HuFNSC01_02_brain_only_isoform_gene, file = "HuFNSC01_02_brain_only_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
HuFNSC01_02_cortex_only_isoform_gene <- cortex01_02_isoform_gene[cortex01_02_isoform_gene$id %in% setdiff(cortex01_02_isoform_gene$id, union(brain01_02_isoform_gene$id, ge01_02_isoform_gene$id)),]
write.table(HuFNSC01_02_cortex_only_isoform_gene, file = "HuFNSC01_02_cortex_only_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
HuFNSC01_02_ge_only_isoform_gene <- ge01_02_isoform_gene[ge01_02_isoform_gene$id %in% setdiff(ge01_02_isoform_gene$id, union(brain01_02_isoform_gene$id, cortex01_02_isoform_gene$id)),]
write.table(HuFNSC01_02_ge_only_isoform_gene, file = "HuFNSC01_02_ge_only_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
HuFNSC01_02_brain_cortex_isoform_gene <- brain01_02_isoform_gene[brain01_02_isoform_gene$id %in% setdiff(intersect(brain01_02_isoform_gene$id, cortex01_02_isoform_gene$id), ge01_02_isoform_gene$id),]
write.table(HuFNSC01_02_brain_cortex_isoform_gene, file = "HuFNSC01_02_brain_cortex_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
HuFNSC01_02_brain_ge_isoform_gene <- brain01_02_isoform_gene[brain01_02_isoform_gene$id %in% setdiff(intersect(brain01_02_isoform_gene$id, ge01_02_isoform_gene$id), cortex01_02_isoform_gene$id),]
write.table(HuFNSC01_02_brain_ge_isoform_gene, file = "HuFNSC01_02_brain_ge_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
HuFNSC01_02_cortex_ge_isoform_gene <- cortex01_02_isoform_gene[cortex01_02_isoform_gene$id %in% setdiff(intersect(cortex01_02_isoform_gene$id, ge01_02_isoform_gene$id), brain01_02_isoform_gene$id),]
write.table(HuFNSC01_02_cortex_ge_isoform_gene, file = "HuFNSC01_02_cortex_ge_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

HuFNSC03_04_isoform_gene <- rbind(cortex03_04_isoform_gene, ge03_04_isoform_gene)
HuFNSC03_04_isoform_gene <- HuFNSC03_04_isoform_gene[duplicated(HuFNSC03_04_isoform_gene$id), ]
HuFNSC03_04_isoform_gene <- HuFNSC03_04_isoform_gene[!duplicated(HuFNSC03_04_isoform_gene$id), ]
HuFNSC03_04_isoform <- rbind(cortex03_04_isoform, ge03_04_isoform)
HuFNSC03_04_isoform <- HuFNSC03_04_isoform[(HuFNSC03_04_isoform$id %in% HuFNSC03_04_isoform_gene$id), ]
(nrow(HuFNSC03_04_isoform))
(nrow(HuFNSC03_04_isoform_gene))
write.table(HuFNSC03_04_isoform, file = "HuFNSC03_04_isoform.txt", sep = "\t", quote=F, row.names=F)
write.table(HuFNSC03_04_isoform_gene, file = "HuFNSC03_04_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
HuFNSC03_04_cortex_only_isoform_gene <- cortex03_04_isoform_gene[cortex03_04_isoform_gene$id %in% setdiff(cortex03_04_isoform_gene$id, ge03_04_isoform_gene$id),]
write.table(HuFNSC03_04_cortex_only_isoform_gene, file = "HuFNSC03_04_cortex_only_isoform_gene.txt", sep = "\t", quote=F, row.names=F)
HuFNSC03_04_ge_only_isoform_gene <- ge03_04_isoform_gene[ge03_04_isoform_gene$id %in% setdiff(ge03_04_isoform_gene$id, cortex03_04_isoform_gene$id),]
write.table(HuFNSC03_04_ge_only_isoform_gene, file = "HuFNSC03_04_ge_only_isoform_gene.txt", sep = "\t", quote=F, row.names=F)

rm(fetalBrain_gene, fetalBrain_exon, cell1, cell2, donor1, donor2, brain01_02, cortex01_02, cortex03_04, ge01_02, ge03_04)
save.image(file = "fetalBrain_isoform_individual.Rdata")

library(VennDiagram)
library(venneuler)
isoform_gene_individual12 <- list(brain = brain01_02_isoform_gene$id, cortex = cortex01_02_isoform_gene$id, GE = ge01_02_isoform_gene$id)
isoform_gene_individual34 <- list(cortex = cortex03_04_isoform_gene$id, GE = ge03_04_isoform_gene$id)
pdf("~/快盘/FetalBrain/RNAseq/isoform/venn_isoform_gene_individual.pdf")
plot.new()
venn_isoform_gene_individual12 <- venn.diagram(isoform_gene_individual12, filename = NULL, fill = c("red", "blue", "green"), main = "Venn diagram of HuFNSC01 vs HuFNSC02 isoforms")
grid.draw(venn_isoform_gene_individual12)

venneuler_isoform_gene_01_02 <- venneuler(data.frame(elements = c(brain01_02_isoform_gene$id, cortex01_02_isoform_gene$id, ge01_02_isoform_gene$id), 
                                                     sets = c(rep("brain", nrow(brain01_02_isoform_gene)), rep("cortex", nrow(cortex01_02_isoform_gene)), rep("GE", nrow(ge01_02_isoform_gene)))))
plot(venneuler_isoform_gene_01_02)
mtext("Venn diagram of HuFNSC01 vs HuFNSC02 isoforms", side = 3)
plot.new()
venn_isoform_gene_individual34 <- venn.diagram(isoform_gene_individual34, filename = NULL, fill = c("red", "blue"), main = "Venn diagram of HuFNSC03 vs HuFNSC04 isoforms")
grid.draw(venn_isoform_gene_individual34)
venneuler_isoform_gene_03_04 <- venneuler(data.frame(elements = c(cortex03_04_isoform_gene$id, ge03_04_isoform_gene$id), 
                                                     sets = c(rep("cortex", nrow(cortex03_04_isoform_gene)), rep("GE", nrow(ge03_04_isoform_gene)))))
plot(venneuler_isoform_gene_03_04)
mtext("Venn diagram of HuFNSC03 vs HuFNSC04 isoforms", side = 3)
dev.off()

save(isoform_cortexge, isoform_individual, file = "fetalBrain_isoform_summary.Rdata")

##################################################################################
##################################################################################







