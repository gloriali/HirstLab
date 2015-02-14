setwd("/Users/epigenomics_lab_02/rotation/CRF")
cortex01<-read.table(file="JOC75JOC75_mCRF.bed",as.is=T,fill=T)
cortex02<-read.table(file="JOC77JOC77_mCRF.bed",as.is=T,fill=T)
brain01<-read.table(file="JOC79JOC79_mCRF.bed",as.is=T,fill=T)
brain02<-read.table(file="JOC80JOC80_mCRF.bed",as.is=T,fill=T)
ge01<-read.table(file="JOC76JOC76_mCRF.bed",as.is=T,fill=T)
ge02<-read.table(file="JOC78JOC78_mCRF.bed",as.is=T,fill=T)
save.image("CRF.Rdata")

############################################################################################################################################

#general density compare
load("CRF.Rdata")
xrange<-range(0,1)
yrange<-range(density(cortex01$methyl,adjust=15)$y,density(cortex02$methyl,adjust=15)$y,density(brain01$methyl,adjust=15)$y,density(brain02$methyl,adjust=15)$y,density(ge01$methyl,adjust=15)$y,density(ge02$methyl,adjust=15)$y)
pdf("densitycompare.pdf")
plot(xrange,yrange,type="n",xlab="Methylation level",ylab="density",main="MethylCRF")
lines(density(brain01$methyl,adjust=15),col=1,pch=1,lty=1)
lines(density(brain02$methyl,adjust=15),col=2,pch=2,lty=1)
lines(density(cortex01$methyl,adjust=15),col=3,pch=3,lty=1)
lines(density(cortex02$methyl,adjust=15),col=4,pch=4,lty=1)
lines(density(ge01$methyl,adjust=15),col=5,pch=5,lty=1)
lines(density(ge02$methyl,adjust=15),col=6,pch=6,lty=1)
legend("topleft",c("brain01","brain02","cortex01","cortex02","ge01","ge02"),col=c(1:6),cex=0.8,lty=c(1,1),pch=c(1:6))
dev.off()

############################################################################################################################################

#bin the genome (1kB bins)
setwd("/Users/epigenomics_lab_02/rotation/")
chrlen<-read.csv("chrlen_hg19.csv",fill=T,as.is=T,header=T,row.names=1)
bin<-data.frame(chr=rep(chrlen$chr,times=(chrlen$length%/%1000+1)))
bin$start<-NA
bin$end<-NA
start<-seq(1,by=1000,length.out=(chrlen$length[1]%/%1000+1))
for(i in 2:24){
  start<-c(start,seq(1,by=1000,length.out=(chrlen$length[i]%/%1000+1)))
}
bin$start<-start
bin$end<-bin$start+999
bin$end[chrlen$X1000e]<-chrlen$length
write.csv(bin,file="bin1000.csv",quote=F)

#methyl-level summary in bins
load("CRF.Rdata")
chrlen<-read.csv("chrlen_hg19.csv",row.names=1)
bin1000<-read.csv("bin1000.csv")

chr<-as.character(brain01$chr)
chr<-gsub("chr","",chr)
startbin<-chrlen[chr,5]  #first bin of the chr
start<-brain01$start
bin<-startbin+start%/%1000    #row number 
summary<-data.frame(id=1:nrow(brain01),bin=bin,chr=chr,start=start,startbin=startbin,count=1,brain01=brain01$methyl,brain02=brain02$methyl,cortex01=cortex01$methyl,cortex02=cortex02$methyl,ge01=ge01$methyl,ge02=ge02$methyl)
summary<-na.omit(summary)
crf<-data.frame(count=1,brain01=summary$brain01,brain02=summary$brain02,cortex01=summary$cortex01,cortex02=summary$cortex02,ge01=summary$ge01,ge02=summary$ge02)
crf<-aggregate(crf,by=list(summary$bin),sum)
colnames(crf)<-c("bin","count","brain01","brain02","cortex01","cortex02","ge01","ge02")
crf$chr<-bin1000$chr[crf$bin]
crf$start<-bin1000$start[crf$bin]
crf$end<-bin1000$end[crf$bin]
crfaverage<-data.frame(chr=crf$chr,start=crf$start,end=crf$end,bin=crf$bin,brain01=(crf$brain01/crf$count),brain02=(crf$brain02/crf$count),cortex01=(crf$cortex01/crf$count),cortex02=(crf$cortex02/crf$count),ge01=(crf$ge01/crf$count),ge02=(crf$ge02/crf$count))
write.csv(summary,file="summary.csv",quote=F)
write.csv(crf,file="crf.csv",quote=F)
write.csv(crfaverage,file="crfaverage.csv",quote=F)

############################################################################################################################################
#################################################should NOT average among different tissues!!!##############################################
crf01<-data.frame(brain=crfaverage$brain01,cortex=crfaverage$cortex01,ge=crfaverage$ge01)
crf02<-data.frame(brain=crfaverage$brain02,cortex=crfaverage$cortex02,ge=crfaverage$ge02)
crfaverage$m01<-rowMeans(crf01)
crfaverage$m02<-rowMeans(crf02)
library(genefilter)
crfaverage$sd01<-rowSds(crf01)
crfaverage$sd02<-rowSds(crf02)
crfaverage$logfold<-log2(crfaverage$m01/crfaverage$m02)
crfaverage$ttest=(crfaverage$m01-crfaverage$m02)/sqrt(crfaverage$sd01^2/3 + crfaverage$sd02^2/3)
crfaverage$pval=2*(1-pt(abs(crfaverage$ttest),4))
diff<-na.omit(crfaverage)
diff<-diff[(abs(diff$logfold)>1)&(diff$pval<0.01),]

pdf("volcano_smooth.pdf")
smoothScatter(crfaverage$logfold,-log10(crfaverage$pval),main="MethylCRF_differential",xlab="log2FoldChange",ylab="-log10(pvalue)")
points(diff$logfold,-log10(diff$pval),col = "red", pch = 20,cex=0.35)
abline(h=2)
abline(v=1)
abline(v=-1)
dev.off()
pdf("volcano.pdf")
plot(crfaverage$logfold,-log10(crfaverage$pval),main="MethylCRF_differential",xlab="log2FoldChange",ylab="-log10(pvalue)",pch = 20,cex=0.15)
points(diff$logfold,-log10(diff$pval),col = "red", pch = 20,cex=0.35)
abline(h=2)
abline(v=1)
abline(v=-1)
dev.off()

write.csv(crfaverage,file="crfaverage.csv",quote=F)
write.csv(diff,file="differential.csv",quote=F)
write.table(diff,file="differential.cod",sep="\t",quote=F)
save.image("differential.Rdata")
#################################################should NOT average among different tissues!!!##############################################
############################################################################################################################################

#methylation level plot check
pdf(file="crf_bin.pdf")
plot(xrange,yrange,type="n",main="CRF density",xlab="CRF",ylab="density")
lines(density(crfaverage$brain01),col=1,pch=1,lty=1)
lines(density(crfaverage$brain02),col=2,pch=2,lty=1)
lines(density(crfaverage$cortex01),col=3,pch=3,lty=1)
lines(density(crfaverage$cortex02),col=4,pch=4,lty=1)
lines(density(crfaverage$ge01),col=5,pch=5,lty=1)
lines(density(crfaverage$ge02),col=6,pch=6,lty=1)
legend("topleft",c("brain01","brain02","cortex01","cortex02","GE01","GE02"),col=c(1:6),pch=c(1:6),lty=1)
dev.off()
pdf(file="methycorrelarion.pdf")
plot(crfaverage$brain01,crfaverage$brain02,pch=20,cex=0.2)
plot(crfaverage$cortex01,crfaverage$cortex02,pch=20,cex=0.2)
plot(crfaverage$ge01,crfaverage$ge02,pch=20,cex=0.2)
dev.off()
pdf(file="methycorrelarion_smooth.pdf")
smoothScatter(crfaverage$brain01,crfaverage$brain02)
smoothScatter(crfaverage$cortex01,crfaverage$cortex02)
smoothScatter(crfaverage$ge01,crfaverage$ge02)
dev.off()
pdf(file="MA_smooth.pdf")
smoothScatter((log(crfaverage$brain01)+log(crfaverage$brain02))/2,log(crfaverage$brain01/crfaverage$brain02))
abline(h=1)
abline(h=-1)
smoothScatter((log(crfaverage$cortex01)+log(crfaverage$cortex02))/2,log(crfaverage$cortex01/crfaverage$cortex02))
abline(h=1)
abline(h=-1)
smoothScatter((log(crfaverage$ge01)+log(crfaverage$ge02))/2,log(crfaverage$ge01/crfaverage$ge02))
abline(h=1)
abline(h=-1)
dev.off()
pdf(file="MA.pdf")
plot((log(crfaverage$brain01)+log(crfaverage$brain02))/2,log(crfaverage$brain01/crfaverage$brain02),pch=20,cex=0.2)
abline(h=1)
abline(h=-1)
plot((log(crfaverage$cortex01)+log(crfaverage$cortex02))/2,log(crfaverage$cortex01/crfaverage$cortex02),pch=20,cex=0.2)
abline(h=1)
abline(h=-1)
plot((log(crfaverage$ge01)+log(crfaverage$ge02))/2,log(crfaverage$ge01/crfaverage$ge02),pch=20,cex=0.2)
abline(h=1)
abline(h=-1)
dev.off()
#differential methylated bins (fold change>2  &  at least one sample CRF>0.5)
crfaverage$logfold_brain<-log2(crfaverage$brain01/crfaverage$brain02)
crfaverage$logfold_cortex<-log2(crfaverage$cortex01/crfaverage$cortex02)
crfaverage$logfold_ge<-log2(crfaverage$ge01/crfaverage$ge02)
write.csv(crfaverage,file="crfaverage.csv",quote=F,row.names=F)
save.image(crfaverage,file="crfaverage.Rdata")

fold_brain<-crfaverage[abs(crfaverage$logfold_brain)>1,]
fold_brain$cortex01<-NULL
fold_brain$cortex02<-NULL
fold_brain$ge01<-NULL
fold_brain$ge02<-NULL
fold_brain$logfold_cortex<-NULL
fold_brain$logfold_ge<-NULL
fold_cortex<-crfaverage[abs(crfaverage$logfold_cortex)>1,]
fold_cortex$brain01<-NULL
fold_cortex$brain02<-NULL
fold_cortex$ge01<-NULL
fold_cortex$ge02<-NULL
fold_cortex$logfold_brain<-NULL
fold_cortex$logfold_ge<-NULL
fold_ge<-crfaverage[abs(crfaverage$logfold_ge)>1,]
fold_ge$brain01<-NULL
fold_ge$brain02<-NULL
fold_ge$cortex01<-NULL
fold_ge$cortex02<-NULL
fold_ge$logfold_brain<-NULL
fold_ge$logfold_cortex<-NULL
save(fold_brain,fold_cortex,fold_ge,file="fold_change.Rdata")

xrange<-range(fold_brain$logfold_brain,fold_cortex$logfold_cortex,fold_ge$logfold_ge)
yrange<-range(density(fold_brain$logfold_brain)$y,density(fold_cortex$logfold_cortex)$y,density(fold_ge$logfold_ge)$y)
pdf(file="Asymmetry.pdf")
plot(xrange,yrange,type="n",main="Global DNA Methylation Asymmetry",xlab="log2(MethylCRF01/MethylCRF02)",ylab="density")
lines(density(fold_brain$logfold_brain),col=1,lty=1,lwd=3)
lines(density(fold_cortex$logfold_cortex),col=2,lty=1,lwd=3)
lines(density(fold_ge$logfold_ge),col=3,lty=1,lwd=3)
legend("topright",c("brain","cortex","GE"),col=c(1:3),lty=1,lwd=5)
dev.off()
pdf("hist.pdf")
hist(fold_brain$logfold_brain,main="brain",breaks=100)
hist(fold_cortex$logfold_cortex,main="cortex",breaks=100)
hist(fold_ge$logfold_ge,main="GE",breaks=100)
dev.off()

#CRF cut off=0.5 (methylated in one sample and unmethylated in the other)
diff_brain<-fold_brain[((fold_brain$brain01>0.5)|(fold_brain$brain02>0.5)),]
diff_cortex<-fold_cortex[((fold_cortex$cortex01>0.5)|(fold_cortex$cortex02>0.5)),]
diff_ge<-fold_ge[((fold_ge$ge01>0.5)|(fold_ge$ge02>0.5)),]
write.csv(diff_brain,file="diff_brain.csv",quote=F,row.names=F)
write.csv(diff_cortex,file="diff_cortex.csv",quote=F,row.names=F)
write.csv(diff_ge,file="diff_ge.csv",quote=F,row.names=F)
DMR_brain<-diff_brain[,2:4]
DMR_brain$ID<-paste(as.character(DMR_brain$chr),":",as.character(DMR_brain$start),"-",as.character(DMR_brain$end),sep="")
DMR_cortex<-diff_cortex[,2:4]
DMR_cortex$ID<-paste(as.character(DMR_cortex$chr),":",as.character(DMR_cortex$start),"-",as.character(DMR_cortex$end),sep="")
DMR_ge<-diff_ge[,2:4]
DMR_ge$ID<-paste(as.character(DMR_ge$chr),":",as.character(DMR_ge$start),"-",as.character(DMR_ge$end),sep="")
write.table(DMR_brain,file="DMR_brain.bed",sep="\t",row.names=F,col.names=F,quote=F)
write.table(DMR_cortex,file="DMR_cortex.bed",sep="\t",row.names=F,col.names=F,quote=F)
write.table(DMR_ge,file="DMR_ge.bed",sep="\t",row.names=F,col.names=F,quote=F)
save(diff_brain,diff_cortex,diff_ge,file="diff_foldchange.Rdata")

setwd("~/FetalBrain/MeDIPMRE/CRF/")
load("diff_foldchange.Rdata")
xrange<-range(diff_brain$logfold_brain,diff_cortex$logfold_cortex,diff_ge$logfold_ge)
yrange<-range(density(diff_brain$logfold_brain)$y,density(diff_cortex$logfold_cortex)$y,density(diff_ge$logfold_ge)$y)
pdf(file="Asymmetry_v2.pdf")
plot(xrange,yrange,type="n",main="Global DNA Methylation Asymmetry",xlab="log2(MethylCRF01/MethylCRF02)",ylab="density")
lines(density(diff_brain$logfold_brain),col=1,lty=1,lwd=3)
lines(density(diff_cortex$logfold_cortex),col=2,lty=1,lwd=3)
lines(density(diff_ge$logfold_ge),col=3,lty=1,lwd=3)
legend("topright",c("brain","cortex","GE"),col=c(1:3),lty=1,lwd=5)
dev.off()

pos_brain<-diff_brain[diff_brain$logfold_brain>0,]
pos_cortex<-diff_cortex[diff_cortex$logfold_cortex>0,]
pos_ge<-diff_ge[diff_ge$logfold_ge>0,]
neg_brain<-diff_brain[diff_brain$logfold_brain<0,]
neg_cortex<-diff_cortex[diff_cortex$logfold_cortex<0,]
neg_ge<-diff_ge[diff_ge$logfold_ge<0,]
DMR_brain_cortex<-merge(DMR_brain,DMR_cortex)
write.table(DMR_brain_cortex,file="DMR_brain_cortex.bed",sep="\t",row.names=F,col.names=F,quote=F)
brain_cortex_pos<-merge(pos_brain,pos_cortex)
DMR_brain_cortex_pos<-brain_cortex_pos[,2:4]
DMR_brain_cortex_pos$ID<-paste(as.character(DMR_brain_cortex_pos$chr),":",as.character(DMR_brain_cortex_pos$start),"-",as.character(DMR_brain_cortex_pos$end),sep="")
write.table(DMR_brain_cortex_pos,file="DMR_brain_cortex_pos.bed",sep="\t",row.names=F,col.names=F,quote=F)
brain_cortex_neg<-merge(neg_brain,neg_cortex)
DMR_brain_cortex_neg<-brain_cortex_neg[,2:4]
DMR_brain_cortex_neg$ID<-paste(as.character(DMR_brain_cortex_neg$chr),":",as.character(DMR_brain_cortex_neg$start),"-",as.character(DMR_brain_cortex_neg$end),sep="")
write.table(DMR_brain_cortex_neg,file="DMR_brain_cortex_neg.bed",sep="\t",row.names=F,col.names=F,quote=F)
save(pos_brain,pos_cortex,pos_ge,neg_brain,neg_cortex,neg_ge,brain_cortex,brain_cortex_pos,brain_cortex_neg,file="overlap.Rdata")
pos<-matrix(c(298,2053,3962,2820552),nrow=2,byrow=T)
fisher.test(pos)
neg<-matrix(c(41,417,1141,2825266),nrow=2,byrow=T)
fisher.test(neg)
############################################################################################################################################

#annotation
setwd("/Users/epigenomics_lab_02/rotation/CRF")
load("diff_foldchange.Rdata")
brain<-data.frame(id=diff_brain$X, chr=diff_brain$chr,start=diff_brain$start,end=diff_brain$end,strand="+")
write.table(brain,file="differential_brain.txt",quote=F,sep="\t",row.names=F,col.names=F)
cortex<-data.frame(id=diff_cortex$X, chr=diff_cortex$chr,start=diff_cortex$start,end=diff_cortex$end,strand="+")
write.table(cortex,file="differential_cortex.txt",quote=F,sep="\t",row.names=F,col.names=F)
ge<-data.frame(id=diff_ge$X, chr=diff_ge$chr,start=diff_ge$start,end=diff_ge$end,strand="+")
write.table(ge,file="differential_ge.txt",quote=F,sep="\t",row.names=F,col.names=F)
#perl annotatePeaks.pl /home/lli/MethylCRF/differential_brain.txt hg19 > /home/lli/MethylCRF/Eannotate_brain.txt -cTSS /home/lli/homer/data/genomes/hg19/ensembl_tss.txt -noann
#perl annotatePeaks.pl /home/lli/MethylCRF/differential_cortex.txt hg19 > /home/lli/MethylCRF/Eannotate_cortex.txt -cTSS /home/lli/homer/data/genomes/hg19/ensembl_tss.txt -noann
#perl annotatePeaks.pl /home/lli/MethylCRF/differential_ge.txt hg19 > /home/lli/MethylCRF/Eannotate_ge.txt -cTSS /home/lli/homer/data/genomes/hg19/ensembl_tss.txt -noann

#expression difference
setwd("/Users/epigenomics_lab_02/rotation/RNAseq/rpkm.pc")
rpkm<-read.csv("rpkm1.csv",row.names=1)
setwd("~/FetalBrain/MeDIPMRE/CRF")
annotate_brain<-read.table("Eannotate_brain.txt",as.is=T,fill=T,head=T,sep="\t")
colnames(annotate_brain)[1]<-"PeakID"
annotate_brain<-annotate_brain[,colSums(is.na(annotate_brain)) == 0]
annotate_brain$Strand<-NULL
annotate_brain$Peak.Score<-NULL
annotate_brain$Ensembl<-gsub("_[0-9,a-z]+","",annotate_brain$Nearest.PromoterID,perl=T,ignore.case=T)
annotate_cortex<-read.table("Eannotate_cortex.txt",as.is=T,fill=T,head=T,sep="\t")
colnames(annotate_cortex)[1]<-"PeakID"
annotate_cortex<-annotate_cortex[,colSums(is.na(annotate_cortex)) == 0]
annotate_cortex$Strand<-NULL
annotate_cortex$Peak.Score<-NULL
annotate_cortex$Ensembl<-gsub("_[0-9,a-z]+","",annotate_cortex$Nearest.PromoterID,perl=T,ignore.case=T)
annotate_ge<-read.table("Eannotate_ge.txt",as.is=T,fill=T,head=T,sep="\t")
colnames(annotate_ge)[1]<-"PeakID"
annotate_ge<-annotate_ge[,colSums(is.na(annotate_ge)) == 0]
annotate_ge$Strand<-NULL
annotate_ge$Peak.Score<-NULL
annotate_ge$Ensembl<-gsub("_[0-9,a-z]+","",annotate_ge$Nearest.PromoterID,perl=T,ignore.case=T)

diff_brain<-read.csv("diff_brain.csv",as.is=T,fill=T,head=T,row.names=1)
diff_cortex<-read.csv("diff_cortex.csv",as.is=T,fill=T,head=T,row.names=1)
diff_ge<-read.csv("diff_ge.csv",as.is=T,fill=T,head=T,row.names=1)

#brain
expression_brain<-annotate_brain[(abs(annotate_brain$Distance.to.TSS)<2000),]
expression_brain$brain01<-rpkm[expression_brain$Ensembl,9]
expression_brain$brain02<-rpkm[expression_brain$Ensembl,10]
expression_brain<-na.omit(expression_brain)
expression_brain$logfold<-log2(expression_brain$brain01/expression_brain$brain02)
expression_brain$z<-abs(expression_brain$brain01-expression_brain$brain02)/sqrt(expression_brain$brain01+(174797872/137879052)*expression_brain$brain02)
expression_brain$CRF01<-diff_brain[as.character(expression_brain$PeakID),5]
expression_brain$CRF02<-diff_brain[as.character(expression_brain$PeakID),6]
expression_brain$CRF_logfold<-diff_brain[as.character(expression_brain$PeakID),7]
expression_brain<-expression_brain[is.finite(expression_brain$logfold),]
#cortex
expression_cortex<-annotate_cortex[(abs(annotate_cortex$Distance.to.TSS)<2000),]
expression_cortex$cortex01<-rpkm[expression_cortex$Ensembl,1]
expression_cortex$cortex02<-rpkm[expression_cortex$Ensembl,2]
expression_cortex<-na.omit(expression_cortex)
expression_cortex$logfold<-log2(expression_cortex$cortex01/expression_cortex$cortex02)
expression_cortex$z<-abs(expression_cortex$cortex01-expression_cortex$cortex02)/sqrt(expression_cortex$cortex01+(154389698/182310872)*expression_cortex$cortex02)
expression_cortex$CRF01<-diff_cortex[as.character(expression_cortex$PeakID),5]
expression_cortex$CRF02<-diff_cortex[as.character(expression_cortex$PeakID),6]
expression_cortex$CRF_logfold<-diff_cortex[as.character(expression_cortex$PeakID),7]
expression_cortex<-expression_cortex[is.finite(expression_cortex$logfold),]
#ge
expression_ge<-annotate_ge[(abs(annotate_ge$Distance.to.TSS)<2000),]
expression_ge$ge01<-rpkm[expression_ge$Ensembl,1]
expression_ge$ge02<-rpkm[expression_ge$Ensembl,2]
expression_ge<-na.omit(expression_ge)
expression_ge$logfold<-log2(expression_ge$ge01/expression_ge$ge02)
expression_ge$z<-abs(expression_ge$ge01-expression_ge$ge02)/sqrt(expression_ge$ge01+(191064090/219123224)*expression_ge$ge02)
expression_ge$CRF01<-diff_ge[as.character(expression_ge$PeakID),5]
expression_ge$CRF02<-diff_ge[as.character(expression_ge$PeakID),6]
expression_ge$CRF_logfold<-diff_ge[as.character(expression_ge$PeakID),7]
expression_ge<-expression_ge[is.finite(expression_ge$logfold),]
write.csv(expression_brain,file="expression_brain2000.csv",row.names=F)
write.csv(expression_cortex,file="expression_cortex2000.csv",row.names=F)
write.csv(expression_ge,file="expression_ge2000.csv",row.names=F)
#differential expressed
diff_ex_brain<-expression_brain[((abs(expression_brain$logfold)>1)&(expression_brain$z>0.75)),]
diff_ex_cortex<-expression_cortex[((abs(expression_cortex$logfold)>1)&(expression_cortex$z>0.75)),]
diff_ex_ge<-expression_ge[((abs(expression_ge$logfold)>1)&(expression_ge$z>0.75)),]
diff_ex_brain<-diff_ex_brain[((diff_ex_brain$logfold*diff_ex_brain$CRF_logfold)<0),]
diff_ex_cortex<-diff_ex_cortex[((diff_ex_cortex$logfold*diff_ex_cortex$CRF_logfold)<0),]
diff_ex_ge<-diff_ex_ge[((diff_ex_ge$logfold*diff_ex_ge$CRF_logfold)<0),]
#get gene information
setwd("/Users/epigenomics_lab_02/FetalBrain/Annotation_hg19/")
gene<-read.table(file="hg19v65_genes",as.is=T,head=F)
colnames(gene)<-c("Ensembl","chr","start","end","strand","type","name","description")
gene<-gene[!duplicated(gene$Ensembl),]
rownames(gene)<-gene$Ensembl
diff_ex_brain$Name<-gene[diff_ex_brain$Ensembl,7]
diff_ex_brain$Description<-gene[diff_ex_brain$Ensembl,8]
#diff_ex_brain$Gchr<-gene[diff_ex_brain$Ensembl,2]
diff_ex_brain$geneStart<-gene[diff_ex_brain$Ensembl,3]
diff_ex_brain$geneEnd<-gene[diff_ex_brain$Ensembl,4]
diff_ex_brain$ucsc<-paste(diff_ex_brain$Chr,":",diff_ex_brain$Start,"-",diff_ex_brain$End,sep="")
diff_ex_cortex$Name<-gene[diff_ex_cortex$Ensembl,7]
diff_ex_cortex$Description<-gene[diff_ex_cortex$Ensembl,8]
#diff_ex_cortex$Gchr<-gene[diff_ex_cortex$Ensembl,2]
diff_ex_cortex$geneStart<-gene[diff_ex_cortex$Ensembl,3]
diff_ex_cortex$geneEnd<-gene[diff_ex_cortex$Ensembl,4]
diff_ex_cortex$ucsc<-paste(diff_ex_cortex$Chr,":",diff_ex_cortex$Start,"-",diff_ex_cortex$End,sep="")
diff_ex_ge$Name<-gene[diff_ex_ge$Ensembl,7]
diff_ex_ge$Description<-gene[diff_ex_ge$Ensembl,8]
#diff_ex_ge$Gchr<-gene[diff_ex_ge$Ensembl,2]
diff_ex_ge$geneStart<-gene[diff_ex_ge$Ensembl,3]
diff_ex_ge$geneEnd<-gene[diff_ex_ge$Ensembl,4]
diff_ex_ge$ucsc<-paste(diff_ex_ge$Chr,":",diff_ex_ge$Start,"-",diff_ex_ge$End,sep="")

diff_ex_brain<-diff_ex_brain[order(diff_ex_brain$CRF_logfold,decreasing=T),]
diff_ex_cortex<-diff_ex_cortex[order(diff_ex_cortex$CRF_logfold,decreasing=T),]
diff_ex_ge<-diff_ex_ge[order(diff_ex_ge$CRF_logfold,decreasing=T),]
setwd("/Users/epigenomics_lab_02/rotation/CRF")
write.csv(diff_ex_brain,file="diff_ex_brain2000.csv",row.names=F)
write.csv(diff_ex_cortex,file="diff_ex_cortex2000.csv",row.names=F)
write.csv(diff_ex_ge,file="diff_ex_ge2000.csv",row.names=F)
save.image("annotation2000.Rdata")

# get promoter DMRs: 
setwd("~/FetalBrain/MeDIPMRE/CRF")
annotate_brain<-read.table("Eannotate_brain.txt",as.is=T,fill=T,head=T,sep="\t")
colnames(annotate_brain)[1]<-"PeakID"
annotate_brain<-annotate_brain[,colSums(is.na(annotate_brain)) == 0]
annotate_brain$Strand<-NULL
annotate_brain$Peak.Score<-NULL
annotate_brain$Ensembl<-gsub("_[0-9,a-z]+","",annotate_brain$Nearest.PromoterID,perl=T,ignore.case=T)
annotate_cortex<-read.table("Eannotate_cortex.txt",as.is=T,fill=T,head=T,sep="\t")
colnames(annotate_cortex)[1]<-"PeakID"
annotate_cortex<-annotate_cortex[,colSums(is.na(annotate_cortex)) == 0]
annotate_cortex$Strand<-NULL
annotate_cortex$Peak.Score<-NULL
annotate_cortex$Ensembl<-gsub("_[0-9,a-z]+","",annotate_cortex$Nearest.PromoterID,perl=T,ignore.case=T)
annotate_ge<-read.table("Eannotate_ge.txt",as.is=T,fill=T,head=T,sep="\t")
colnames(annotate_ge)[1]<-"PeakID"
annotate_ge<-annotate_ge[,colSums(is.na(annotate_ge)) == 0]
annotate_ge$Strand<-NULL
annotate_ge$Peak.Score<-NULL
annotate_ge$Ensembl<-gsub("_[0-9,a-z]+","",annotate_ge$Nearest.PromoterID,perl=T,ignore.case=T)
promoter_brain<-annotate_brain[(abs(annotate_brain$Distance.to.TSS)<2000),]
promoter_cortex<-annotate_cortex[(abs(annotate_cortex$Distance.to.TSS)<2000),]
promoter_ge<-annotate_ge[(abs(annotate_ge$Distance.to.TSS)<2000),]
promoter_brain$category<-gsub("ENSG[0-9]+_","",promoter_brain$Nearest.PromoterID,perl=T,ignore.case=F)
promoter_cortex$category<-gsub("ENSG[0-9]+_","",promoter_cortex$Nearest.PromoterID,perl=T,ignore.case=F)
promoter_ge$category<-gsub("ENSG[0-9]+_","",promoter_ge$Nearest.PromoterID,perl=T,ignore.case=F)
write.table(data.frame(chr=promoter_brain$Chr,start=promoter_brain$Start,end=promoter_brain$End,id=promoter_brain$Nearest.PromoterID),file="DMR_promoter_brain.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(data.frame(chr=promoter_cortex$Chr,start=promoter_cortex$Start,end=promoter_cortex$End,id=promoter_cortex$Nearest.PromoterID),file="DMR_promoter_cortex.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(data.frame(chr=promoter_ge$Chr,start=promoter_ge$Start,end=promoter_ge$End,id=promoter_ge$Nearest.PromoterID),file="DMR_promoter_ge.bed",sep="\t",col.names=F,row.names=F,quote=F)
promoter_pc_brain<-promoter_brain[promoter_brain$category=="protein_coding",]
write.table(data.frame(chr=promoter_pc_brain$Chr,start=promoter_pc_brain$Start,end=promoter_pc_brain$End,id=promoter_pc_brain$Nearest.PromoterID),file="DMR_promoter_pc_brain.bed",sep="\t",col.names=F,row.names=F,quote=F)
promoter_pc_cortex<-promoter_cortex[promoter_cortex$category=="protein_coding",]
write.table(data.frame(chr=promoter_pc_cortex$Chr,start=promoter_pc_cortex$Start,end=promoter_pc_cortex$End,id=promoter_pc_cortex$Nearest.PromoterID),file="DMR_promoter_pc_cortex.bed",sep="\t",col.names=F,row.names=F,quote=F)
promoter_pc_ge<-promoter_ge[promoter_ge$category=="protein_coding",]
write.table(data.frame(chr=promoter_pc_ge$Chr,start=promoter_pc_ge$Start,end=promoter_pc_ge$End,id=promoter_pc_ge$Nearest.PromoterID),file="DMR_promoter_pc_ge.bed",sep="\t",col.names=F,row.names=F,quote=F)
save.image("promoterDMR.Rdata")

############################################################################################################################################
############################################################################################################################################

#bin around TSS (+/- 1.5kB)
load("CRF.Rdata")
CpG<-data.frame(id=rownames(brain01),chr=brain01$chr,start=brain01$start,end=brain01$end,strand="+")
write.table(CpG,file="CpG.txt",quote=F,sep="\t",row.names=F,col.names=F)
save(CpG,file="CpG.Rdata")
CpG1<-CpG[1:6000000,]
CpG2<-CpG[6000001:12000000,]
CpG3<-CpG[12000001:18000000,]
CpG4<-CpG[18000001:24000000,]
CpG5<-CpG[24000001:nrow(CpG),]
write.table(CpG1,file="CpG1.txt",quote=F,sep="\t",row.names=F,col.names=F)
write.table(CpG2,file="CpG2.txt",quote=F,sep="\t",row.names=F,col.names=F)
write.table(CpG3,file="CpG3.txt",quote=F,sep="\t",row.names=F,col.names=F)
write.table(CpG4,file="CpG4.txt",quote=F,sep="\t",row.names=F,col.names=F)
write.table(CpG5,file="CpG5.txt",quote=F,sep="\t",row.names=F,col.names=F)

#perl annotatePeaks.pl /home/lli/MethylCRF/binTSS/CpG1.txt hg19> /home/lli/MethylCRF/binTSS/Eannotate_CpG1.txt -cTSS /home/lli/homer/data/genomes/hg19/ensembl_tss.txt -noann
#perl annotatePeaks.pl /home/lli/MethylCRF/binTSS/CpG2.txt hg19> /home/lli/MethylCRF/binTSS/Eannotate_CpG2.txt -cTSS /home/lli/homer/data/genomes/hg19/ensembl_tss.txt -noann
#perl annotatePeaks.pl /home/lli/MethylCRF/binTSS/CpG3.txt hg19> /home/lli/MethylCRF/binTSS/Eannotate_CpG3.txt -cTSS /home/lli/homer/data/genomes/hg19/ensembl_tss.txt -noann
#perl annotatePeaks.pl /home/lli/MethylCRF/binTSS/CpG4.txt hg19> /home/lli/MethylCRF/binTSS/Eannotate_CpG4.txt -cTSS /home/lli/homer/data/genomes/hg19/ensembl_tss.txt -noann
#perl annotatePeaks.pl /home/lli/MethylCRF/binTSS/CpG5.txt hg19> /home/lli/MethylCRF/binTSS/Eannotate_CpG5.txt -cTSS /home/lli/homer/data/genomes/hg19/ensembl_tss.txt -noann

setwd("/home/lli/MethylCRF/binTSS")
#perl CpG.pl (convert annotation file to csv file)
annotate1<-read.csv("Eannotate_CpG1.csv",as.is=T,head=T)
annotate2<-read.csv("Eannotate_CpG2.csv",as.is=T,head=T)
annotate3<-read.csv("Eannotate_CpG3.csv",as.is=T,head=T)
annotate4<-read.csv("Eannotate_CpG4.csv",as.is=T,head=T)
annotate5<-read.csv("Eannotate_CpG5.csv",as.is=T,head=T)
annotate1<-na.omit(annotate1)
annotate2<-na.omit(annotate2)
annotate3<-na.omit(annotate3)
annotate4<-na.omit(annotate4)
annotate5<-na.omit(annotate5)
data<-list(annotate1,annotate2,annotate3,annotate4,annotate5)
annotate<-do.call("rbind", data)
#annotate<-annotate[order(annotate$promoterID),]
write.csv(annotate,"annotate_CpG.csv",row.names=F,quote=T)
save(annotate,file="annotate.Rdata")

annotate<-annotate[(abs(annotate$distance_to_TSS)<1500),]
annotate$Ensembl<-gsub("_[0-9a-z]+","",annotate$promoterID,perl=T,ignore.case=T)
annotate<-na.omit(annotate)
annotate<-annotate[grepl("E",annotate$Ensembl),]
annotate<-annotate[order(annotate$Ensembl),]
write.csv(annotate,"annotate_CpG.csv",row.names=F)
annotate<-read.csv("annotate_CpG.csv",as.is=T,head=T)
setwd("/home/lli/MethylCRF")
load("CRF.Rdata")
setwd("/home/lli/MethylCRF/binTSS")
crf<-data.frame(count=1,brain01=brain01[annotate$ID,5],brain02=brain02[annotate$ID,5],cortex01=cortex01[annotate$ID,5],cortex02=cortex02[annotate$ID,5],ge01=ge01[annotate$ID,5],ge02=ge02[annotate$ID,5])
crf<-aggregate(crf,by=list(annotate$Ensembl),sum)
colnames(crf)<-c("Ensembl","count","brain01","brain02","cortex01","cortex02","ge01","ge02")
write.csv(crf,file="crf.csv",quote=F)
crfaverage<-data.frame(Ensembl=crf$Ensembl,brain01=(crf$brain01/crf$count),brain02=(crf$brain02/crf$count),cortex01=(crf$cortex01/crf$count),cortex02=(crf$cortex02/crf$count),ge01=(crf$ge01/crf$count),ge02=(crf$ge02/crf$count))
crfaverage$logfold_brain<-log2(crfaverage$brain01/crfaverage$brain02)
crfaverage$logfold_cortex<-log2(crfaverage$cortex01/crfaverage$cortex02)
crfaverage$logfold_ge<-log2(crfaverage$ge01/crfaverage$ge02)
write.csv(crfaverage,file="crfaverage.csv",row.names=F)
save(crfaverage,file="crfaverage.Rdata")
#differential methylated (fold change >2)
setwd("/Users/epigenomics_lab_02/rotation/CRF/binTSS")
load("crfaverage.Rdata")
diff_brain<-crfaverage[abs(crfaverage$logfold_brain)>1,]
diff_brain$cortex01<-NULL
diff_brain$cortex02<-NULL
diff_brain$ge01<-NULL
diff_brain$ge02<-NULL
diff_brain$logfold_cortex<-NULL
diff_brain$logfold_ge<-NULL
diff_cortex<-crfaverage[abs(crfaverage$logfold_cortex)>1,]
diff_cortex$brain01<-NULL
diff_cortex$brain02<-NULL
diff_cortex$ge01<-NULL
diff_cortex$ge02<-NULL
diff_cortex$logfold_brain<-NULL
diff_cortex$logfold_ge<-NULL
diff_ge<-crfaverage[abs(crfaverage$logfold_ge)>1,]
diff_ge$brain01<-NULL
diff_ge$brain02<-NULL
diff_ge$cortex01<-NULL
diff_ge$cortex02<-NULL
diff_ge$logfold_brain<-NULL
diff_ge$logfold_cortex<-NULL
#CRF cut off=0.5 (methylated in one sample and unmethylated in the other)
diff_brain<-diff_brain[((diff_brain$brain01>0.5)|(diff_brain$brain02>0.5)),]
diff_cortex<-diff_cortex[((diff_cortex$cortex01>0.5)|(diff_cortex$cortex02>0.5)),]
diff_ge<-diff_ge[((diff_ge$ge01>0.5)|(diff_ge$ge02>0.5)),]
write.csv(diff_brain,file="diff_brain.csv",row.names=F)
write.csv(diff_cortex,file="diff_cortex.csv",row.names=F)
write.csv(diff_ge,file="diff_ge.csv",row.names=F)
save(diff_brain,diff_cortex,diff_ge,file="diff_foldchange.Rdata")
#expression analysis
setwd("/home/lli/RNAseq/rpkm.pc")
rpkm<-read.csv("rpkm1.csv",row.names=1)
setwd("/home/lli/MethylCRF/binTSS")
expression_brain<-diff_brain
expression_brain$brain01ex<-rpkm[expression_brain$Ensembl,9]
expression_brain$brain02ex<-rpkm[expression_brain$Ensembl,10]
expression_brain<-na.omit(expression_brain)
expression_brain$logfoldex<-log2(expression_brain$brain01ex/expression_brain$brain02ex)
expression_brain$z<-abs(expression_brain$brain01ex-expression_brain$brain02ex)/sqrt(expression_brain$brain01ex+(174797872/137879052)*expression_brain$brain02ex)
expression_brain<-expression_brain[is.finite(expression_brain$logfoldex),]
expression_cortex<-diff_cortex
expression_cortex$cortex01ex<-rpkm[expression_cortex$Ensembl,9]
expression_cortex$cortex02ex<-rpkm[expression_cortex$Ensembl,10]
expression_cortex<-na.omit(expression_cortex)
expression_cortex$logfoldex<-log2(expression_cortex$cortex01ex/expression_cortex$cortex02ex)
expression_cortex$z<-abs(expression_cortex$cortex01ex-expression_cortex$cortex02ex)/sqrt(expression_cortex$cortex01ex+(154389698/182310872)*expression_cortex$cortex02ex)
expression_cortex<-expression_cortex[is.finite(expression_cortex$logfoldex),]
expression_ge<-diff_ge
expression_ge$ge01ex<-rpkm[expression_ge$Ensembl,9]
expression_ge$ge02ex<-rpkm[expression_ge$Ensembl,10]
expression_ge<-na.omit(expression_ge)
expression_ge$logfoldex<-log2(expression_ge$ge01ex/expression_ge$ge02ex)
expression_ge$z<-abs(expression_ge$ge01ex-expression_ge$ge02ex)/sqrt(expression_ge$ge01ex+(191064090/219123224)*expression_ge$ge02ex)
expression_ge<-expression_ge[is.finite(expression_ge$logfoldex),]
write.csv(expression_brain,file="expression_brain.csv",row.names=F)
write.csv(expression_cortex,file="expression_cortex.csv",row.names=F)
write.csv(expression_ge,file="expression_ge.csv",row.names=F)
#differential expressed
diff_ex_brain<-expression_brain[((abs(expression_brain$logfoldex)>1)&(expression_brain$z>0.5)),]
diff_ex_cortex<-expression_cortex[((abs(expression_cortex$logfoldex)>1)&(expression_cortex$z>0.5)),]
diff_ex_ge<-expression_ge[((abs(expression_ge$logfoldex)>1)&(expression_ge$z>0.5)),]
gene<-read.csv(file="hg19v65_genes.csv",as.is=T,head=F)
colnames(gene)<-c("Ensembl","chr","start","end","strand","type","name","description")
gene<-gene[!duplicated(gene$Ensembl),]
rownames(gene)<-gene$Ensembl
diff_ex_brain$Name<-gene[diff_ex_brain$Ensembl,7]
diff_ex_brain$Description<-gene[diff_ex_brain$Ensembl,8]
diff_ex_brain$geneStart<-gene[diff_ex_brain$Ensembl,3]
diff_ex_brain$geneEnd<-gene[diff_ex_brain$Ensembl,4]
diff_ex_brain$ucsc<-paste(diff_ex_brain$Chr,":",diff_ex_brain$Start,"-",diff_ex_brain$End,sep="")
diff_ex_cortex$Name<-gene[diff_ex_cortex$Ensembl,7]
diff_ex_cortex$Description<-gene[diff_ex_cortex$Ensembl,8]
diff_ex_cortex$geneStart<-gene[diff_ex_cortex$Ensembl,3]
diff_ex_cortex$geneEnd<-gene[diff_ex_cortex$Ensembl,4]
diff_ex_cortex$ucsc<-paste(diff_ex_cortex$Chr,":",diff_ex_cortex$Start,"-",diff_ex_cortex$End,sep="")
diff_ex_ge$Name<-gene[diff_ex_ge$Ensembl,7]
diff_ex_ge$Description<-gene[diff_ex_ge$Ensembl,8]
diff_ex_ge$geneStart<-gene[diff_ex_ge$Ensembl,3]
diff_ex_ge$geneEnd<-gene[diff_ex_ge$Ensembl,4]
diff_ex_ge$ucsc<-paste(diff_ex_ge$Chr,":",diff_ex_ge$Start,"-",diff_ex_ge$End,sep="")
write.csv(diff_ex_brain,file="diff_ex_brain.csv",row.names=F)
write.csv(diff_ex_cortex,file="diff_ex_cortex.csv",row.names=F)
write.csv(diff_ex_ge,file="diff_ex_ge.csv",row.names=F)
save.image("diff_ex.Rdata")

############################################################################################################################################
############################################################################################################################################

#CpG by CpG approach
load("CRF.Rdata")
brain<-data.frame(chr=brain01$chr,start=brain01$start,end=brain01$end,id=brain01$id,brain01=brain01$methyl,brain02=brain02$methyl)
brain$logfold_brain<-log2(brain$brain01/brain$brain02)
cortex<-data.frame(chr=cortex01$chr,start=cortex01$start,end=cortex01$end,id=cortex01$id,cortex01=cortex01$methyl,cortex02=cortex02$methyl)
cortex$logfold_cortex<-log2(cortex$cortex01/cortex$cortex02)
ge<-data.frame(chr=ge01$chr,start=ge01$start,end=ge01$end,id=ge01$id,ge01=ge01$methyl,ge02=ge02$methyl)
ge$logfold_ge<-log2(ge$ge01/ge$ge02)
write.csv(brain,file="brain.csv",quote=F,row.names=F)
write.csv(cortex,file="cortex.csv",quote=F,row.names=F)
write.csv(ge,file="ge.csv",quote=F,row.names=F)
save.image("CpG.Rdata")
fold_brain<-brain[(abs(brain$logfold_brain)>1),]
fold_cortex<-cortex[(abs(cortex$logfold_cortex)>1),]
fold_ge<-ge[(abs(ge$logfold_ge)>1),]
write.csv(fold_brain,file="fold_brain.csv",quote=F,row.names=F)
write.csv(fold_cortex,file="fold_cortex.csv",quote=F,row.names=F)
write.csv(fold_ge,file="fold_ge.csv",quote=F,row.names=F)
save(fold_brain,fold_cortex,fold_ge,file="fold.Rdata")

pos_brain<-fold_brain[fold_brain$logfold_brain>0,]
pos_cortex<-fold_cortex[fold_cortex$logfold_cortex>0,]
pos_ge<-fold_ge[fold_ge$logfold_ge>0,]
neg_brain<-fold_brain[fold_brain$logfold_brain<0,]
neg_cortex<-fold_cortex[fold_cortex$logfold_cortex<0,]
neg_ge<-fold_ge[fold_ge$logfold_ge<0,]
brain_cortex<-merge(pos_brain,pos_cortex)
overlap<-merge(brain_cortex,neg_ge)
save(pos_brain,pos_cortex,pos_ge,neg_brain,neg_cortex,neg_ge,brain_cortex,overlap,file="overlap.Rdata")

xrange<-range(density(brain$logfold)$x,density(cortex$logfold)$x,density(ge$logfold)$x)
yrange<-range(density(brain$logfold)$y,density(cortex$logfold)$y,density(ge$logfold)$y)
#yrange=range(0,0.5)
pdf(file="foldchange_all.pdf")
plot(xrange,yrange,type="n",main="log fold change density",xlab="log2 fold change",ylab="density")
lines(density(brain$logfold),col=1,pch=1,lty=1)
lines(density(cortex$logfold),col=2,pch=2,lty=1)
lines(density(ge$logfold),col=3,pch=3,lty=1)
legend("topleft",c("brain","cortex","GE"),col=c(1:3),pch=c(1:3),lty=1)
dev.off()
xrange<-range(density(fold_brain$logfold)$x,density(fold_cortex$logfold)$x,density(fold_ge$logfold)$x)
yrange<-range(density(fold_brain$logfold)$y,density(fold_cortex$logfold)$y,density(fold_ge$logfold)$y)
pdf(file="foldchange.pdf")
plot(xrange,yrange,type="n",main="log fold change density",xlab="log2 fold change",ylab="density")
lines(density(fold_brain$logfold),col=1,pch=1,lty=1)
lines(density(fold_cortex$logfold),col=2,pch=2,lty=1)
lines(density(fold_ge$logfold),col=3,pch=3,lty=1)
legend("topleft",c("brain","cortex","GE"),col=c(1:3),pch=c(1:3),lty=1)
dev.off()
pdf("cdf.pdf")
plot(ecdf(brain$logfold),verticals=T,col=1,main="CDF log fold change",xlab="log2 fold change",ylab="CDF")
lines(ecdf(cortex$logfold),verticals=T,col=2)
lines(ecdf(ge$logfold),verticals=T,col=3)
abline(h=0.5)
legend("topleft",c("brain","cortex","GE"),col=c(1:3),lty=1)
dev.off()

diff_brain<-brain[abs(brain$logfold)>2,]
diff_cortex<-cortex[abs(cortex$logfold)>2,]
diff_ge<-ge[abs(ge$logfold)>2,]
diff_brain<-diff_brain[((diff_brain$brain01>0.5)|(diff_brain$brain02>0.5)),]
diff_cortex<-diff_cortex[((diff_cortex$cortex01>0.5)|(diff_cortex$cortex02>0.5)),]
diff_ge<-diff_ge[((diff_ge$ge01>0.5)|(diff_ge$ge02>0.5)),]
write.csv(diff_brain,file="diff_brain.csv",row.names=F)
write.csv(diff_cortex,file="diff_cortex.csv",row.names=F)
write.csv(diff_ge,file="diff_ge.csv",row.names=F)
save.image("diff_foldchange.Rdata")

############################################################################################################################################
############################################################################################################################################

#sliding window approach
load("CRF.Rdata")
chrlen<-read.csv("chrlen_hg19.csv",row.names=1)
bin200<-read.csv("bin200.csv")
#summary to 200bp bins
chr<-as.character(brain01$chr)
chr<-gsub("chr","",chr)
methyl<-brain01$methyl
startbin<-chrlen[chr,5]  #first bin of the chr
start<-brain01$start
bin<-startbin+start%/%200    #row number 
summary<-data.frame(id=1:nrow(brain01),bin=bin,chr=chr,start=start,startbin=startbin,count=1,brain01=brain01$methyl,brain02=brain02$methyl,cortex01=cortex01$methyl,cortex02=cortex02$methyl,ge01=ge01$methyl,ge02=ge02$methyl)
summary<-na.omit(summary)
crf<-data.frame(count=1,brain01=summary$brain01,brain02=summary$brain02,cortex01=summary$cortex01,cortex02=summary$cortex02,ge01=summary$ge01,ge02=summary$ge02)
crf<-aggregate(crf,by=list(summary$bin),sum)
colnames(crf)<-c("bin","count","brain01","brain02","cortex01","cortex02","ge01","ge02")
rownames(crf)<-as.character(crf$bin)
#crf$chr<-bin200$chr[crf$bin]
#crf$start<-bin200$start[crf$bin]
#crf$end<-bin200$end[crf$bin]
bin200$count<-1
bin200$brain01<-0
bin200$brain02<-0
bin200$cortex01<-0
bin200$cortex02<-0
bin200$ge01<-0
bin200$ge02<-0
bin200$count[crf$bin]<-crf[as.character(crf$bin),2]
bin200$brain01[crf$bin]<-crf[as.character(crf$bin),3]
bin200$brain02[crf$bin]<-crf[as.character(crf$bin),4]
bin200$cortex01[crf$bin]<-crf[as.character(crf$bin),5]
bin200$cortex02[crf$bin]<-crf[as.character(crf$bin),6]
bin200$ge01[crf$bin]<-crf[as.character(crf$bin),7]
bin200$ge02[crf$bin]<-crf[as.character(crf$bin),8]
crf<-bin200
crfaverage<-data.frame(chr=crf$chr,start=crf$start,end=crf$end,bin=crf$X,brain01=(crf$brain01/crf$count),brain02=(crf$brain02/crf$count),cortex01=(crf$cortex01/crf$count),cortex02=(crf$cortex02/crf$count),ge01=(crf$ge01/crf$count),ge02=(crf$ge02/crf$count))
save(summary,file="summary.Rdata")
save(crf,file="crf.Rdata")
save(crfaverage,file="crfaverage.Rdata")

#sliding window
install.packages("zoo",lib="/home/lli/R")
require(zoo,lib="/home/lli/R")
crfaverage$chr<-gsub("chr","",crfaverage$chr)
crfaverage$chr<-gsub("X","23",crfaverage$chr)
crfaverage$chr<-gsub("Y","24",crfaverage$chr)
crfaverage$chr<-as.numeric(crfaverage$chr)
slide<-rollmean(crfaverage,5)
slide<-as.data.frame(slide)
slide<-slide[(slide$chr%%1==0),]
slide$start<-slide$start-400
slide$end<-slide$end+400
slide$chr<-as.character(slide$chr)
slide$chr<-gsub("23","X",slide$chr)
slide$chr<-gsub("24","Y",slide$chr)
slide$chr<-paste("chr",slide$chr,sep="")
crfaverage<-slide
save(crfaverage,file="crfaverage.Rdata")
