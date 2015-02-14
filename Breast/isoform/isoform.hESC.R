##################################################################################
# getting the cutoff for expressing and non-expressing
# run on xhost
# /gsc/software/linux-x86_64-centos5/R-2.15.2/bin/R
setwd("~/REMC/hESC/")
dirIn='/projects/epigenomics/ep50/external/jqc.1.7.6/'
libs=c("CD184",
       "hESC_Derived_CD56plus_Ectoderm_Cultured_Cells",
       "hESC_Derived_CD56plus_Mesoderm_Cultured_Cells"
       )

pdf("percent_ave_including0_zoom.pdf")
plot(range(0,0.3),range(0,0.5),type='n',main="ECDF percentage of average RPKM",xlab="exon RPKM/average RPKM",ylab="cdf")
i=1
for(lib in libs){
  exon_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.exn.A.rpkm",sep=""))
  pc_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.A.rpkm.pc",sep=""),row.names=1)
  exon_lib$ave<-pc_lib[exon_lib$V2,3]
  exon_lib<-exon_lib[exon_lib$ave!=0,]
  exon_lib$percent<-exon_lib$V4/exon_lib$ave
  #lines(ecdf(exon_lib$percent[exon_lib$percent!=0]),col=i,lty=1)
  lines(ecdf(exon_lib$percent),col=i,lty=1)
  i=i+1
}
abline(v=0.01)
abline(v=0.1)
legend("bottomright",c("Endoderm","Ectoderm","Mesoderm"),col=c(1:length(libs)),cex=0.8,lty=1)
dev.off()
##################################################################################


##################################################################################
##isoforms:DEfine on exons & RPKM(exon) < 1%*RPKM(ave. exons of the gene) for one sample& > 10% for the other & expressed in both samples(ave.RPKM>Rmin)
##################################################################################

##################################################################################
## breast: tissue specific (vHMEC; myo; lum; stem). 
##################################################################################
# run on xhost
# /gsc/software/linux-x86_64-centos5/R-2.15.2/bin/R
setwd("~/REMC/hESC/")
col<-c("Gene id", "rpkm1", "rpkm2", "p-value", "Multiple testing corrected p-value")
dirIn='/projects/epigenomics/ep50/external/jqc.1.7.6/'
cutoff=0.01
cutoff2=0.1 # sample specific exon: expressed in one sample(>0.1ave.RPKM) & not expressed in the other(<0.01ave.RPKM)
Rmin=0.005

lib1='hESC_Derived_CD56plus_Ectoderm_Cultured_Cells'; cell1='ecto'; donor1='hESC';
lib2='CD184'; cell2='endo'; donor2='hESC';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
ecto_endo_up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
ecto_endo_up$id<-unlist(strsplit(as.character(ecto_endo_up$V1),'_'))[2*(1:nrow(ecto_endo_up))]
ecto_endo_up$ave2<-pc_lib1[ecto_endo_up$id,3]
ecto_endo_up$ave3<-pc_lib2[ecto_endo_up$id,3]
ecto_endo_up<-na.omit(ecto_endo_up)
ecto_endo_up<-ecto_endo_up[(ecto_endo_up$V2<=cutoff*ecto_endo_up$ave2&ecto_endo_up$V3>=cutoff2*ecto_endo_up$ave3)|(ecto_endo_up$V2>=cutoff2*ecto_endo_up$ave2&ecto_endo_up$V3<=cutoff*ecto_endo_up$ave3),]
ecto_endo_up_gene<-unique(ecto_endo_up$id)
ecto_endo_up_isoform<-ecto_endo_up[ecto_endo_up$ave2>Rmin&ecto_endo_up$ave3>Rmin,]
write.csv(ecto_endo_up,file="ecto_endo_up.csv",quote=F,row.names=F)
write.csv(ecto_endo_up_gene,file="ecto_endo_up_gene.csv",quote=F,row.names=F)
write.csv(ecto_endo_up_isoform,file="ecto_endo_up_isoform.csv",quote=F,row.names=F)
ecto_endo_dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
ecto_endo_dn$id<-unlist(strsplit(as.character(ecto_endo_dn$V1),'_'))[2*(1:nrow(ecto_endo_dn))]
ecto_endo_dn$ave2<-pc_lib1[ecto_endo_dn$id,3]
ecto_endo_dn$ave3<-pc_lib2[ecto_endo_dn$id,3]
ecto_endo_dn<-na.omit(ecto_endo_dn)
ecto_endo_dn<-ecto_endo_dn[(ecto_endo_dn$V2<=cutoff*ecto_endo_dn$ave2&ecto_endo_dn$V3>=cutoff2*ecto_endo_dn$ave3)|(ecto_endo_dn$V2>=cutoff2*ecto_endo_dn$ave2&ecto_endo_dn$V3<=cutoff*ecto_endo_dn$ave3),]
ecto_endo_dn_gene<-unique(ecto_endo_dn$id)
ecto_endo_dn_isoform<-ecto_endo_dn[ecto_endo_dn$ave2>Rmin&ecto_endo_dn$ave3>Rmin,]
write.csv(ecto_endo_dn,file="ecto_endo_dn.csv",quote=F,row.names=F)
write.csv(ecto_endo_dn_gene,file="ecto_endo_dn_gene.csv",quote=F,row.names=F)
write.csv(ecto_endo_dn_isoform,file="ecto_endo_dn_isoform.csv",quote=F,row.names=F)
ecto_endo_isoform<-unique(c(ecto_endo_up_isoform$id,ecto_endo_dn_isoform$id))
write.csv(ecto_endo_isoform,file="ecto_endo_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/hESC/gene/")
DE_ecto_endo_up<-read.table("UP.Endo_Ecto.FDR_0.01.rmin_0.005.Nmin_25")
DE_ecto_endo_dn<-read.table("DN.Endo_Ecto.FDR_0.01.rmin_0.005.Nmin_25")
DE_ecto_endo<-c(as.character(DE_ecto_endo_up$V1),as.character(DE_ecto_endo_dn$V1))
setwd("~/REMC/hESC/")
ecto_endo_isoform_only<-ecto_endo_isoform[!is.element(ecto_endo_isoform,DE_ecto_endo)]
write.csv(ecto_endo_isoform_only,file="ecto_endo_isoform_only.csv",quote=F,row.names=F)

lib1='hESC_Derived_CD56plus_Ectoderm_Cultured_Cells'; cell1='ecto'; donor1='hESC';
lib2='hESC_Derived_CD56plus_Mesoderm_Cultured_Cells'; cell2='meso'; donor2='hESC';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
ecto_meso_up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
ecto_meso_up$id<-unlist(strsplit(as.character(ecto_meso_up$V1),'_'))[2*(1:nrow(ecto_meso_up))]
ecto_meso_up$ave2<-pc_lib1[ecto_meso_up$id,3]
ecto_meso_up$ave3<-pc_lib2[ecto_meso_up$id,3]
ecto_meso_up<-na.omit(ecto_meso_up)
ecto_meso_up<-ecto_meso_up[(ecto_meso_up$V2<=cutoff*ecto_meso_up$ave2&ecto_meso_up$V3>=cutoff2*ecto_meso_up$ave3)|(ecto_meso_up$V2>=cutoff2*ecto_meso_up$ave2&ecto_meso_up$V3<=cutoff*ecto_meso_up$ave3),]
ecto_meso_up_gene<-unique(ecto_meso_up$id)
ecto_meso_up_isoform<-ecto_meso_up[ecto_meso_up$ave2>Rmin&ecto_meso_up$ave3>Rmin,]
write.csv(ecto_meso_up,file="ecto_meso_up.csv",quote=F,row.names=F)
write.csv(ecto_meso_up_gene,file="ecto_meso_up_gene.csv",quote=F,row.names=F)
write.csv(ecto_meso_up_isoform,file="ecto_meso_up_isoform.csv",quote=F,row.names=F)
ecto_meso_dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
ecto_meso_dn$id<-unlist(strsplit(as.character(ecto_meso_dn$V1),'_'))[2*(1:nrow(ecto_meso_dn))]
ecto_meso_dn$ave2<-pc_lib1[ecto_meso_dn$id,3]
ecto_meso_dn$ave3<-pc_lib2[ecto_meso_dn$id,3]
ecto_meso_dn<-na.omit(ecto_meso_dn)
ecto_meso_dn<-ecto_meso_dn[(ecto_meso_dn$V2<=cutoff*ecto_meso_dn$ave2&ecto_meso_dn$V3>=cutoff2*ecto_meso_dn$ave3)|(ecto_meso_dn$V2>=cutoff2*ecto_meso_dn$ave2&ecto_meso_dn$V3<=cutoff*ecto_meso_dn$ave3),]
ecto_meso_dn_gene<-unique(ecto_meso_dn$id)
ecto_meso_dn_isoform<-ecto_meso_dn[ecto_meso_dn$ave2>Rmin&ecto_meso_dn$ave3>Rmin,]
write.csv(ecto_meso_dn,file="ecto_meso_dn.csv",quote=F,row.names=F)
write.csv(ecto_meso_dn_gene,file="ecto_meso_dn_gene.csv",quote=F,row.names=F)
write.csv(ecto_meso_dn_isoform,file="ecto_meso_dn_isoform.csv",quote=F,row.names=F)
ecto_meso_isoform<-unique(c(ecto_meso_up_isoform$id,ecto_meso_dn_isoform$id))
write.csv(ecto_meso_isoform,file="ecto_meso_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/hESC/gene/")
DE_ecto_meso_up<-read.table("UP.Ecto_Meso.FDR_0.01.rmin_0.005.Nmin_25")
DE_ecto_meso_dn<-read.table("DN.Ecto_Meso.FDR_0.01.rmin_0.005.Nmin_25")
DE_ecto_meso<-c(as.character(DE_ecto_meso_up$V1),as.character(DE_ecto_meso_dn$V1))
setwd("~/REMC/hESC/")
ecto_meso_isoform_only<-ecto_meso_isoform[!is.element(ecto_meso_isoform,DE_ecto_meso)]
write.csv(ecto_meso_isoform_only,file="ecto_meso_isoform_only.csv",quote=F,row.names=F)

lib1='CD184'; cell1='endo'; donor1='hESC';
lib2='hESC_Derived_CD56plus_Mesoderm_Cultured_Cells'; cell2='meso'; donor2='hESC';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
endo_meso_up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
endo_meso_up$id<-unlist(strsplit(as.character(endo_meso_up$V1),'_'))[2*(1:nrow(endo_meso_up))]
endo_meso_up$ave2<-pc_lib1[endo_meso_up$id,3]
endo_meso_up$ave3<-pc_lib2[endo_meso_up$id,3]
endo_meso_up<-na.omit(endo_meso_up)
endo_meso_up<-endo_meso_up[(endo_meso_up$V2<=cutoff*endo_meso_up$ave2&endo_meso_up$V3>=cutoff2*endo_meso_up$ave3)|(endo_meso_up$V2>=cutoff2*endo_meso_up$ave2&endo_meso_up$V3<=cutoff*endo_meso_up$ave3),]
endo_meso_up_gene<-unique(endo_meso_up$id)
endo_meso_up_isoform<-endo_meso_up[endo_meso_up$ave2>Rmin&endo_meso_up$ave3>Rmin,]
write.csv(endo_meso_up,file="endo_meso_up.csv",quote=F,row.names=F)
write.csv(endo_meso_up_gene,file="endo_meso_up_gene.csv",quote=F,row.names=F)
write.csv(endo_meso_up_isoform,file="endo_meso_up_isoform.csv",quote=F,row.names=F)
endo_meso_dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
endo_meso_dn$id<-unlist(strsplit(as.character(endo_meso_dn$V1),'_'))[2*(1:nrow(endo_meso_dn))]
endo_meso_dn$ave2<-pc_lib1[endo_meso_dn$id,3]
endo_meso_dn$ave3<-pc_lib2[endo_meso_dn$id,3]
endo_meso_dn<-na.omit(endo_meso_dn)
endo_meso_dn<-endo_meso_dn[(endo_meso_dn$V2<=cutoff*endo_meso_dn$ave2&endo_meso_dn$V3>=cutoff2*endo_meso_dn$ave3)|(endo_meso_dn$V2>=cutoff2*endo_meso_dn$ave2&endo_meso_dn$V3<=cutoff*endo_meso_dn$ave3),]
endo_meso_dn_gene<-unique(endo_meso_dn$id)
endo_meso_dn_isoform<-endo_meso_dn[endo_meso_dn$ave2>Rmin&endo_meso_dn$ave3>Rmin,]
write.csv(endo_meso_dn,file="endo_meso_dn.csv",quote=F,row.names=F)
write.csv(endo_meso_dn_gene,file="endo_meso_dn_gene.csv",quote=F,row.names=F)
write.csv(endo_meso_dn_isoform,file="endo_meso_dn_isoform.csv",quote=F,row.names=F)
endo_meso_isoform<-unique(c(endo_meso_up_isoform$id,endo_meso_dn_isoform$id))
write.csv(endo_meso_isoform,file="endo_meso_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/hESC/gene/")
DE_endo_meso_up<-read.table("UP.Endo_Meso.FDR_0.01.rmin_0.005.Nmin_25")
DE_endo_meso_dn<-read.table("DN.Endo_Meso.FDR_0.01.rmin_0.005.Nmin_25")
DE_endo_meso<-c(as.character(DE_endo_meso_up$V1),as.character(DE_endo_meso_dn$V1))
setwd("~/REMC/hESC/")
endo_meso_isoform_only<-endo_meso_isoform[!is.element(endo_meso_isoform,DE_endo_meso)]
write.csv(endo_meso_isoform_only,file="endo_meso_isoform_only.csv",quote=F,row.names=F)

rm(lib1,cell1,donor1,lib2,cell2,donor2,pc_lib1,pc_lib2)
save.image(file="REMC.isoform.hESC.RData")

#####################################################################################################
# get gene information
setwd("~/hg19/")
gene<-read.table(file="hg19v65_genes",as.is=T,head=F)
colnames(gene)<-c("Ensembl","chr","start","end","strand","type","name","description")
gene<-gene[!duplicated(gene$Ensembl),]
rownames(gene)<-gene$Ensembl

setwd("~/REMC/hESC/")
load("REMC.isoform.hESC.RData")
ecto_endo_isoform_only<-data.frame(ID=ecto_endo_isoform_only,Name=gene[ecto_endo_isoform_only,7],Description=gene[ecto_endo_isoform_only,8])
write.csv(ecto_endo_isoform_only,file="ecto_endo_isoform_only.csv",quote=F,row.names=F)
ecto_meso_isoform_only<-data.frame(ID=ecto_meso_isoform_only,Name=gene[ecto_meso_isoform_only,7],Description=gene[ecto_meso_isoform_only,8])
write.csv(ecto_meso_isoform_only,file="ecto_meso_isoform_only.csv",quote=F,row.names=F)
endo_meso_isoform_only<-data.frame(ID=endo_meso_isoform_only,Name=gene[endo_meso_isoform_only,7],Description=gene[endo_meso_isoform_only,8])
write.csv(endo_meso_isoform_only,file="endo_meso_isoform_only.csv",quote=F,row.names=F)

save.image(file="REMC.isoform.hESC.RData")



