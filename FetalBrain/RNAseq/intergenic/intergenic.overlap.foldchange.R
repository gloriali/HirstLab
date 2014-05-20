setwd("/Users/epigenomics_lab_02/FetalBrain/RNAseq/intergenic/")
# brain12up=read.table("UP.Brain-HuFNSC01_Brain-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_0")
# brain12dn=read.table("DN.Brain-HuFNSC01_Brain-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_0")
# brainDefine=rbind(brain12up,brain12dn)

pdf("plot2D.BetweenTwins.pdf")
lib1='A03484'
lib2='A07825'
posReg=read.table(paste(lib1,'_',lib2,'.pos.intersect.bed',sep=""),as.is=T,fill=T)
posReg$ID=paste(posReg$V1,":",posReg$V2,"-",posReg$V3,"<1",sep="")
posRPKM=read.table(paste(lib1,'_',lib2,'.pos.bed',sep=""),as.is=T,fill=T)
pos=data.frame(ID=posReg$ID,chr=posReg$V1,start=posReg$V2,end=posReg$V3,strand="+",RPKM1=posRPKM$V4,RPKM2=posRPKM$V8)
negReg=read.table(paste(lib1,'_',lib2,'.neg.intersect.bed',sep=""),as.is=T,fill=T)
negReg$ID=paste(negReg$V1,":",negReg$V2,"-",negReg$V3,"<1",sep="")
negRPKM=read.table(paste(lib1,'_',lib2,'.neg.bed',sep=""),as.is=T,fill=T)
neg=data.frame(ID=negReg$ID,chr=negReg$V1,start=negReg$V2,end=negReg$V3,strand="-",RPKM1=negRPKM$V4,RPKM2=negRPKM$V8)
brain12=rbind(pos,neg)
brain12$LFC=log2(brain12$RPKM1/brain12$RPKM2)
smoothScatter(log10(brain12$RPKM1),log10(brain12$RPKM2),main=paste(lib1,'_',lib2,sep=""))
# q5=quantile(brain12$LFC,probs=0.05)
# q95=quantile(brain12$LFC,probs=0.95)
# brain12=brain12[brain12$LFC<q5|brain12$LFC>q95,]
# intersectq5=intersect(brainDefine$V1,brain12$ID)
q1=quantile(brain12$LFC,probs=0.01)
q99=quantile(brain12$LFC,probs=0.99)
brain12=brain12[brain12$LFC<q1|brain12$LFC>q99,]
points(log10(brain12$RPKM1),log10(brain12$RPKM2),col="red",pch=20,cex=0.5)
# rownames(brain12)=brain12$ID
# intersectq1=brain12[intersect(brainDefine$V1,brain12$ID),]

lib1='A03473'
lib2='A03475'
posReg=read.table(paste(lib1,'_',lib2,'.pos.intersect.bed',sep=""),as.is=T,fill=T)
posReg$ID=paste(posReg$V1,":",posReg$V2,"-",posReg$V3,"<1",sep="")
posRPKM=read.table(paste(lib1,'_',lib2,'.pos.bed',sep=""),as.is=T,fill=T)
pos=data.frame(ID=posReg$ID,chr=posReg$V1,start=posReg$V2,end=posReg$V3,strand="+",RPKM1=posRPKM$V4,RPKM2=posRPKM$V8)
negReg=read.table(paste(lib1,'_',lib2,'.neg.intersect.bed',sep=""),as.is=T,fill=T)
negReg$ID=paste(negReg$V1,":",negReg$V2,"-",negReg$V3,"<1",sep="")
negRPKM=read.table(paste(lib1,'_',lib2,'.neg.bed',sep=""),as.is=T,fill=T)
neg=data.frame(ID=negReg$ID,chr=negReg$V1,start=negReg$V2,end=negReg$V3,strand="-",RPKM1=negRPKM$V4,RPKM2=negRPKM$V8)
cortex12=rbind(pos,neg)
cortex12$LFC=log2(cortex12$RPKM1/cortex12$RPKM2)
smoothScatter(log10(cortex12$RPKM1),log10(cortex12$RPKM2),main=paste(lib1,'_',lib2,sep=""))
q1=quantile(cortex12$LFC,probs=0.01)
q99=quantile(cortex12$LFC,probs=0.99)
cortex12=cortex12[cortex12$LFC<q1|cortex12$LFC>q99,]
points(log10(cortex12$RPKM1),log10(cortex12$RPKM2),col="red",pch=20,cex=0.5)

lib1='A03474'
lib2='A03476'
posReg=read.table(paste(lib1,'_',lib2,'.pos.intersect.bed',sep=""),as.is=T,fill=T)
posReg$ID=paste(posReg$V1,":",posReg$V2,"-",posReg$V3,"<1",sep="")
posRPKM=read.table(paste(lib1,'_',lib2,'.pos.bed',sep=""),as.is=T,fill=T)
pos=data.frame(ID=posReg$ID,chr=posReg$V1,start=posReg$V2,end=posReg$V3,strand="+",RPKM1=posRPKM$V4,RPKM2=posRPKM$V8)
negReg=read.table(paste(lib1,'_',lib2,'.neg.intersect.bed',sep=""),as.is=T,fill=T)
negReg$ID=paste(negReg$V1,":",negReg$V2,"-",negReg$V3,"<1",sep="")
negRPKM=read.table(paste(lib1,'_',lib2,'.neg.bed',sep=""),as.is=T,fill=T)
neg=data.frame(ID=negReg$ID,chr=negReg$V1,start=negReg$V2,end=negReg$V3,strand="-",RPKM1=negRPKM$V4,RPKM2=negRPKM$V8)
ge12=rbind(pos,neg)
ge12$LFC=log2(ge12$RPKM1/ge12$RPKM2)
smoothScatter(log10(ge12$RPKM1),log10(ge12$RPKM2),main=paste(lib1,'_',lib2,sep=""))
q1=quantile(ge12$LFC,probs=0.01)
q99=quantile(ge12$LFC,probs=0.99)
ge12=ge12[ge12$LFC<q1|ge12$LFC>q99,]
points(log10(ge12$RPKM1),log10(ge12$RPKM2),col="red",pch=20,cex=0.5)

lib1='A04599'
lib2='A15298'
posReg=read.table(paste(lib1,'_',lib2,'.pos.intersect.bed',sep=""),as.is=T,fill=T)
posReg$ID=paste(posReg$V1,":",posReg$V2,"-",posReg$V3,"<1",sep="")
posRPKM=read.table(paste(lib1,'_',lib2,'.pos.bed',sep=""),as.is=T,fill=T)
pos=data.frame(ID=posReg$ID,chr=posReg$V1,start=posReg$V2,end=posReg$V3,strand="+",RPKM1=posRPKM$V4,RPKM2=posRPKM$V8)
negReg=read.table(paste(lib1,'_',lib2,'.neg.intersect.bed',sep=""),as.is=T,fill=T)
negReg$ID=paste(negReg$V1,":",negReg$V2,"-",negReg$V3,"<1",sep="")
negRPKM=read.table(paste(lib1,'_',lib2,'.neg.bed',sep=""),as.is=T,fill=T)
neg=data.frame(ID=negReg$ID,chr=negReg$V1,start=negReg$V2,end=negReg$V3,strand="-",RPKM1=negRPKM$V4,RPKM2=negRPKM$V8)
cortex34=rbind(pos,neg)
cortex34$LFC=log2(cortex34$RPKM1/cortex34$RPKM2)
smoothScatter(log10(cortex34$RPKM1),log10(cortex34$RPKM2),main=paste(lib1,'_',lib2,sep=""))
q1=quantile(cortex34$LFC,probs=0.01)
q99=quantile(cortex34$LFC,probs=0.99)
cortex34=cortex34[cortex34$LFC<q1|cortex34$LFC>q99,]
points(log10(cortex34$RPKM1),log10(cortex34$RPKM2),col="red",pch=20,cex=0.5)

lib1='A15295'
lib2='A15299'
posReg=read.table(paste(lib1,'_',lib2,'.pos.intersect.bed',sep=""),as.is=T,fill=T)
posReg$ID=paste(posReg$V1,":",posReg$V2,"-",posReg$V3,"<1",sep="")
posRPKM=read.table(paste(lib1,'_',lib2,'.pos.bed',sep=""),as.is=T,fill=T)
pos=data.frame(ID=posReg$ID,chr=posReg$V1,start=posReg$V2,end=posReg$V3,strand="+",RPKM1=posRPKM$V4,RPKM2=posRPKM$V8)
negReg=read.table(paste(lib1,'_',lib2,'.neg.intersect.bed',sep=""),as.is=T,fill=T)
negReg$ID=paste(negReg$V1,":",negReg$V2,"-",negReg$V3,"<1",sep="")
negRPKM=read.table(paste(lib1,'_',lib2,'.neg.bed',sep=""),as.is=T,fill=T)
neg=data.frame(ID=negReg$ID,chr=negReg$V1,start=negReg$V2,end=negReg$V3,strand="-",RPKM1=negRPKM$V4,RPKM2=negRPKM$V8)
ge34=rbind(pos,neg)
ge34$LFC=log2(ge34$RPKM1/ge34$RPKM2)
smoothScatter(log10(ge34$RPKM1),log10(ge34$RPKM2),main=paste(lib1,'_',lib2,sep=""))
q1=quantile(ge34$LFC,probs=0.01)
q99=quantile(ge34$LFC,probs=0.99)
ge34=ge34[ge34$LFC<q1|ge34$LFC>q99,]
points(log10(ge34$RPKM1),log10(ge34$RPKM2),col="red",pch=20,cex=0.5)
dev.off()

pdf("CDF_log2fold.pdf") # plot before filter with 1% cut-off
plot(range(-3,3),range(0,1),type="n",main="Ecdf(log2 fold change)",xlab="LFC",ylab="cdf")
lines(ecdf(brain12$LFC),col=1)
lines(ecdf(cortex12$LFC),col=2)
lines(ecdf(ge12$LFC),col=3)
lines(ecdf(cortex34$LFC),col=4)
lines(ecdf(ge34$LFC),col=5)
abline(h=0.5)
abline(v=0)
legend("bottomright",c("brain01_brain02","cortex01_cortex02","GE01_GE02","cortex03_cortex04","GE03_GE04"),col=(1:5),lty=1)
dev.off()

write.csv(brain12,file="brain12.csv",quote=F,row.names=F)
write.csv(cortex12,file="cortex12.csv",quote=F,row.names=F)
write.csv(ge12,file="ge12.csv",quote=F,row.names=F)
write.csv(cortex34,file="cortex34.csv",quote=F,row.names=F)
write.csv(ge34,file="ge34.csv",quote=F,row.names=F)


