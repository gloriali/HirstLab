##################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
##isoforms:DEfine on exons & RPKM(exon) < 1%*RPKM(ave. exons of the gene) for one sample& > 10% for the other & expressed in both samples(ave.RPKM>Rmin)
# isoform between myo RM084 and H1
##################################################################################

setwd("~/REMC/breast/H1/")
col<-c("Gene id", "rpkm1", "rpkm2", "p-value", "Multiple testing corrected p-value")
dirIn1='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
dirIn2='/projects/epigenomics/ep50/external/jqc.1.7.6/'; 
cutoff=0.01
cutoff2=0.1 # sample specific exon: expressed in one sample(>0.1ave.RPKM) & not expressed in the other(<0.01ave.RPKM)
Rmin=0.005

lib1='A17919'; cell1='myo'; donor1='RM084';
lib2='H1_r1a'; cell2='H1'; donor2='H1_r1a';
pc_lib1<-read.table(paste(dirIn1,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn2,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
myo_H1_r1a_84up<-read.table(paste("./exon/UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_H1_r1a_84up$id<-unlist(strsplit(as.character(myo_H1_r1a_84up$V1),'_'))[2*(1:nrow(myo_H1_r1a_84up))]
myo_H1_r1a_84up$ave2<-pc_lib1[myo_H1_r1a_84up$id,3]
myo_H1_r1a_84up$ave3<-pc_lib2[myo_H1_r1a_84up$id,3]
myo_H1_r1a_84up<-na.omit(myo_H1_r1a_84up)
myo_H1_r1a_84up<-myo_H1_r1a_84up[(myo_H1_r1a_84up$V2<=cutoff*myo_H1_r1a_84up$ave2&myo_H1_r1a_84up$V3>=cutoff2*myo_H1_r1a_84up$ave3)|(myo_H1_r1a_84up$V2>=cutoff2*myo_H1_r1a_84up$ave2&myo_H1_r1a_84up$V3<=cutoff*myo_H1_r1a_84up$ave3),]
myo_H1_r1a_84up_gene<-unique(myo_H1_r1a_84up$id)
myo_H1_r1a_84up_isoform<-myo_H1_r1a_84up[myo_H1_r1a_84up$ave2>Rmin&myo_H1_r1a_84up$ave3>Rmin,]
myo_H1_r1a_84dn<-read.table(paste("./exon/DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_H1_r1a_84dn$id<-unlist(strsplit(as.character(myo_H1_r1a_84dn$V1),'_'))[2*(1:nrow(myo_H1_r1a_84dn))]
myo_H1_r1a_84dn$ave2<-pc_lib1[myo_H1_r1a_84dn$id,3]
myo_H1_r1a_84dn$ave3<-pc_lib2[myo_H1_r1a_84dn$id,3]
myo_H1_r1a_84dn<-na.omit(myo_H1_r1a_84dn)
myo_H1_r1a_84dn<-myo_H1_r1a_84dn[(myo_H1_r1a_84dn$V2<=cutoff*myo_H1_r1a_84dn$ave2&myo_H1_r1a_84dn$V3>=cutoff2*myo_H1_r1a_84dn$ave3)|(myo_H1_r1a_84dn$V2>=cutoff2*myo_H1_r1a_84dn$ave2&myo_H1_r1a_84dn$V3<=cutoff*myo_H1_r1a_84dn$ave3),]
myo_H1_r1a_84dn_gene<-unique(myo_H1_r1a_84dn$id)
myo_H1_r1a_84dn_isoform<-myo_H1_r1a_84dn[myo_H1_r1a_84dn$ave2>Rmin&myo_H1_r1a_84dn$ave3>Rmin,]
myo_H1_r1a_84_isoform<-rbind(myo_H1_r1a_84up_isoform,myo_H1_r1a_84dn_isoform)
# exclude gene level DE
DE_myo_H1_r1a_84up<-read.table("./gene/UP.H1_r1aRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_H1_r1a_84dn<-read.table("./gene/DN.H1_r1aRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_H1_r1a_84<-c(as.character(DE_myo_H1_r1a_84up$V1),as.character(DE_myo_H1_r1a_84dn$V1))
myo_H1_r1a_84_isoform_only<-myo_H1_r1a_84_isoform[!(myo_H1_r1a_84_isoform$id %in% DE_myo_H1_r1a_84), ]
