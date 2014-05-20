# read RNAseq data
setwd("~/FetalBrain/RNAseq/rpkm/")
RNA_cortex01<-read.table(file="A03473.G.A.rpkm.pc",row.names=1,head=F)
RNA_cortex02<-read.table(file="A03475.G.A.rpkm.pc",row.names=1,head=F)
# RNA_cortex03<-read.table(file="A04599.G.A.rpkm.pc",row.names=1,head=F)
# RNA_cortex04<-read.table(file="A15298.G.A.rpkm.pc",row.names=1,head=F)
RNA_ge01<-read.table(file="A03474.G.A.rpkm.pc",row.names=1,head=F)
RNA_ge02<-read.table(file="A03476.G.A.rpkm.pc",row.names=1,head=F)
# RNA_ge03<-read.table(file="A15295.G.A.rpkm.pc",row.names=1,head=F)
# RNA_ge04<-read.table(file="A15299.G.A.rpkm.pc",row.names=1,head=F)
RNA_brain01<-read.table(file="A03484.G.A.rpkm.pc",row.names=1,head=F)
RNA_brain02<-read.table(file="A07825.G.A.rpkm.pc",row.names=1,head=F)

# read all DMRs with methylation level
setwd("~/FetalBrain/MeDIPMRE/CRF/")
diff_brain<-read.csv("diff_brain.csv",as.is=T,fill=T,head=T)
diff_cortex<-read.csv("diff_cortex.csv",as.is=T,fill=T,head=T)
diff_ge<-read.csv("diff_ge.csv",as.is=T,fill=T,head=T)
rownames(diff_brain) <- paste0(diff_brain$chr,":",diff_brain$start,"-",diff_brain$end)
rownames(diff_cortex) <- paste0(diff_cortex$chr,":",diff_cortex$start,"-",diff_cortex$end)
rownames(diff_ge) <- paste0(diff_ge$chr,":",diff_ge$start,"-",diff_ge$end)

# read genebody DMRs
setwd("~/FetalBrain/MeDIPMRE/DMRs/intersect/")
brain<-read.table(file="DMRbrain_genebody_pc.txt",head=F)
cortex<-read.table(file="DMRcortex_genebody_pc.txt",head=F)
ge<-read.table(file="DMRge_genebody_pc.txt",head=F)

#get expression and methylation level 
brain$expression01<-RNA_brain01[brain$V9,]$V3
brain$expression02<-RNA_brain02[brain$V9,]$V3
brain$methylation01<-diff_brain[brain$V4,]$brain01
brain$methylation02<-diff_brain[brain$V4,]$brain02
brain$logfold_methylation <- log2(brain$methylation01/brain$methylation02)
brain$logfold_expression <- log2(brain$expression01/brain$expression02)
cortex$expression01<-RNA_cortex01[cortex$V9,]$V3
cortex$expression02<-RNA_cortex02[cortex$V9,]$V3
cortex$methylation01<-diff_cortex[cortex$V4,]$cortex01
cortex$methylation02<-diff_cortex[cortex$V4,]$cortex02
cortex$logfold_methylation <- log2(cortex$methylation01/cortex$methylation02)
cortex$logfold_expression <- log2(cortex$expression01/cortex$expression02)
ge$expression01<-RNA_ge01[ge$V9,]$V3
ge$expression02<-RNA_ge02[ge$V9,]$V3
ge$methylation01<-diff_ge[ge$V4,]$ge01
ge$methylation02<-diff_ge[ge$V4,]$ge02
ge$logfold_methylation <- log2(ge$methylation01/ge$methylation02)
ge$logfold_expression <- log2(ge$expression01/ge$expression02)
pdf("genebodyDMRs-expression.pdf")
smoothScatter(brain$logfold_methylation,brain$logfold_expression, main="genebodyDMRs-expression brain")
abline(h=0)
smoothScatter(cortex$logfold_methylation,cortex$logfold_expression, main="genebodyDMRs-expression cortex")
abline(h=0)
smoothScatter(ge$logfold_methylation,ge$logfold_expression, main="genebodyDMRs-expression GE")
abline(h=0)
dev.off()

