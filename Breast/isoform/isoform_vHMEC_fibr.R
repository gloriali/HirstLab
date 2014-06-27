##################################################################################
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
##isoforms:DEfine on exons & RPKM(exon) < 1%*RPKM(ave. exons of the gene) for one sample& > 10% for the other & expressed in both samples(ave.RPKM>Rmin)
# isoform between fibr RM070 and vHMEC
##################################################################################

setwd("~/REMC/breast/vHMEC_fibr/")
col<-c("Gene id", "rpkm1", "rpkm2", "p-value", "Multiple testing corrected p-value")
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'; 
cutoff=0.01
cutoff2=0.1 # sample specific exon: expressed in one sample(>0.1ave.RPKM) & not expressed in the other(<0.01ave.RPKM)
Rmin=0.005

lib1='A18760'; cell1='fibr'; donor1='RM070';
lib2='HS2263'; cell2='vHMEC'; donor2='RM035';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
fibr070_vHMEC035up<-read.table(paste("./exon/UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
fibr070_vHMEC035up$id<-unlist(strsplit(as.character(fibr070_vHMEC035up$V1),'_'))[2*(1:nrow(fibr070_vHMEC035up))]
fibr070_vHMEC035up$ave2<-pc_lib1[fibr070_vHMEC035up$id,3]
fibr070_vHMEC035up$ave3<-pc_lib2[fibr070_vHMEC035up$id,3]
fibr070_vHMEC035up<-na.omit(fibr070_vHMEC035up)
fibr070_vHMEC035up<-fibr070_vHMEC035up[(fibr070_vHMEC035up$V2<=cutoff*fibr070_vHMEC035up$ave2&fibr070_vHMEC035up$V3>=cutoff2*fibr070_vHMEC035up$ave3)|(fibr070_vHMEC035up$V2>=cutoff2*fibr070_vHMEC035up$ave2&fibr070_vHMEC035up$V3<=cutoff*fibr070_vHMEC035up$ave3),]
fibr070_vHMEC035up_isoform<-fibr070_vHMEC035up[fibr070_vHMEC035up$ave2>Rmin&fibr070_vHMEC035up$ave3>Rmin,]
fibr070_vHMEC035dn<-read.table(paste("./exon/DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
fibr070_vHMEC035dn$id<-unlist(strsplit(as.character(fibr070_vHMEC035dn$V1),'_'))[2*(1:nrow(fibr070_vHMEC035dn))]
fibr070_vHMEC035dn$ave2<-pc_lib1[fibr070_vHMEC035dn$id,3]
fibr070_vHMEC035dn$ave3<-pc_lib2[fibr070_vHMEC035dn$id,3]
fibr070_vHMEC035dn<-na.omit(fibr070_vHMEC035dn)
fibr070_vHMEC035dn<-fibr070_vHMEC035dn[(fibr070_vHMEC035dn$V2<=cutoff*fibr070_vHMEC035dn$ave2&fibr070_vHMEC035dn$V3>=cutoff2*fibr070_vHMEC035dn$ave3)|(fibr070_vHMEC035dn$V2>=cutoff2*fibr070_vHMEC035dn$ave2&fibr070_vHMEC035dn$V3<=cutoff*fibr070_vHMEC035dn$ave3),]
fibr070_vHMEC035dn_isoform<-fibr070_vHMEC035dn[fibr070_vHMEC035dn$ave2>Rmin&fibr070_vHMEC035dn$ave3>Rmin,]
fibr070_vHMEC035_isoform<-rbind(fibr070_vHMEC035up_isoform,fibr070_vHMEC035dn_isoform)
nrow(fibr070_vHMEC035_isoform)
# exclude gene level DE
DE_fibr070_vHMEC035up<-read.table("./gene/UP.fibrRM070_vHMECRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_fibr070_vHMEC035dn<-read.table("./gene/DN.fibrRM070_vHMECRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_fibr070_vHMEC035<-c(as.character(DE_fibr070_vHMEC035up$V1),as.character(DE_fibr070_vHMEC035dn$V1))
length(DE_fibr070_vHMEC035)
fibr070_vHMEC035_isoform<-fibr070_vHMEC035_isoform[!(fibr070_vHMEC035_isoform$id %in% DE_fibr070_vHMEC035), ]
nrow(fibr070_vHMEC035_isoform)
fibr070_vHMEC035_isoform_gene <- fibr070_vHMEC035_isoform[!duplicated(fibr070_vHMEC035_isoform$id), ]
nrow(fibr070_vHMEC035_isoform_gene)
write.table(fibr070_vHMEC035_isoform, file = "fibr070_vHMEC035_isoform.txt", sep = "\t", quote = F, row.names = F)
write.table(fibr070_vHMEC035_isoform_gene, file = "fibr070_vHMEC035_isoform_gene.txt", sep = "\t", quote = F, row.names = F)

lib1='A18760'; cell1='fibr'; donor1='RM070';
lib2='A18761'; cell2='vHMEC'; donor2='RM071';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
fibr070_vHMEC071up<-read.table(paste("./exon/UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
fibr070_vHMEC071up$id<-unlist(strsplit(as.character(fibr070_vHMEC071up$V1),'_'))[2*(1:nrow(fibr070_vHMEC071up))]
fibr070_vHMEC071up$ave2<-pc_lib1[fibr070_vHMEC071up$id,3]
fibr070_vHMEC071up$ave3<-pc_lib2[fibr070_vHMEC071up$id,3]
fibr070_vHMEC071up<-na.omit(fibr070_vHMEC071up)
fibr070_vHMEC071up<-fibr070_vHMEC071up[(fibr070_vHMEC071up$V2<=cutoff*fibr070_vHMEC071up$ave2&fibr070_vHMEC071up$V3>=cutoff2*fibr070_vHMEC071up$ave3)|(fibr070_vHMEC071up$V2>=cutoff2*fibr070_vHMEC071up$ave2&fibr070_vHMEC071up$V3<=cutoff*fibr070_vHMEC071up$ave3),]
fibr070_vHMEC071up_isoform<-fibr070_vHMEC071up[fibr070_vHMEC071up$ave2>Rmin&fibr070_vHMEC071up$ave3>Rmin,]
fibr070_vHMEC071dn<-read.table(paste("./exon/DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
fibr070_vHMEC071dn$id<-unlist(strsplit(as.character(fibr070_vHMEC071dn$V1),'_'))[2*(1:nrow(fibr070_vHMEC071dn))]
fibr070_vHMEC071dn$ave2<-pc_lib1[fibr070_vHMEC071dn$id,3]
fibr070_vHMEC071dn$ave3<-pc_lib2[fibr070_vHMEC071dn$id,3]
fibr070_vHMEC071dn<-na.omit(fibr070_vHMEC071dn)
fibr070_vHMEC071dn<-fibr070_vHMEC071dn[(fibr070_vHMEC071dn$V2<=cutoff*fibr070_vHMEC071dn$ave2&fibr070_vHMEC071dn$V3>=cutoff2*fibr070_vHMEC071dn$ave3)|(fibr070_vHMEC071dn$V2>=cutoff2*fibr070_vHMEC071dn$ave2&fibr070_vHMEC071dn$V3<=cutoff*fibr070_vHMEC071dn$ave3),]
fibr070_vHMEC071dn_isoform<-fibr070_vHMEC071dn[fibr070_vHMEC071dn$ave2>Rmin&fibr070_vHMEC071dn$ave3>Rmin,]
fibr070_vHMEC071_isoform<-rbind(fibr070_vHMEC071up_isoform,fibr070_vHMEC071dn_isoform)
nrow(fibr070_vHMEC071_isoform)
# exclude gene level DE
DE_fibr070_vHMEC071up<-read.table("./gene/UP.fibrRM070_vHMECRM071.FDR_0.015.rmin_0.005.Nmin_25")
DE_fibr070_vHMEC071dn<-read.table("./gene/DN.fibrRM070_vHMECRM071.FDR_0.015.rmin_0.005.Nmin_25")
DE_fibr070_vHMEC071<-c(as.character(DE_fibr070_vHMEC071up$V1),as.character(DE_fibr070_vHMEC071dn$V1))
length(DE_fibr070_vHMEC071)
fibr070_vHMEC071_isoform<-fibr070_vHMEC071_isoform[!(fibr070_vHMEC071_isoform$id %in% DE_fibr070_vHMEC071), ]
nrow(fibr070_vHMEC071_isoform)
fibr070_vHMEC071_isoform_gene <- fibr070_vHMEC071_isoform[!duplicated(fibr070_vHMEC071_isoform$id), ]
nrow(fibr070_vHMEC071_isoform_gene)
write.table(fibr070_vHMEC071_isoform, file = "fibr070_vHMEC071_isoform.txt", sep = "\t", quote = F, row.names = F)
write.table(fibr070_vHMEC071_isoform_gene, file = "fibr070_vHMEC071_isoform_gene.txt", sep = "\t", quote = F, row.names = F)

lib1='A18472'; cell1='fibr'; donor1='RM071';
lib2='HS2263'; cell2='vHMEC'; donor2='RM035';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
fibr071_vHMEC035up<-read.table(paste("./exon/UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
fibr071_vHMEC035up$id<-unlist(strsplit(as.character(fibr071_vHMEC035up$V1),'_'))[2*(1:nrow(fibr071_vHMEC035up))]
fibr071_vHMEC035up$ave2<-pc_lib1[fibr071_vHMEC035up$id,3]
fibr071_vHMEC035up$ave3<-pc_lib2[fibr071_vHMEC035up$id,3]
fibr071_vHMEC035up<-na.omit(fibr071_vHMEC035up)
fibr071_vHMEC035up<-fibr071_vHMEC035up[(fibr071_vHMEC035up$V2<=cutoff*fibr071_vHMEC035up$ave2&fibr071_vHMEC035up$V3>=cutoff2*fibr071_vHMEC035up$ave3)|(fibr071_vHMEC035up$V2>=cutoff2*fibr071_vHMEC035up$ave2&fibr071_vHMEC035up$V3<=cutoff*fibr071_vHMEC035up$ave3),]
fibr071_vHMEC035up_isoform<-fibr071_vHMEC035up[fibr071_vHMEC035up$ave2>Rmin&fibr071_vHMEC035up$ave3>Rmin,]
fibr071_vHMEC035dn<-read.table(paste("./exon/DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
fibr071_vHMEC035dn$id<-unlist(strsplit(as.character(fibr071_vHMEC035dn$V1),'_'))[2*(1:nrow(fibr071_vHMEC035dn))]
fibr071_vHMEC035dn$ave2<-pc_lib1[fibr071_vHMEC035dn$id,3]
fibr071_vHMEC035dn$ave3<-pc_lib2[fibr071_vHMEC035dn$id,3]
fibr071_vHMEC035dn<-na.omit(fibr071_vHMEC035dn)
fibr071_vHMEC035dn<-fibr071_vHMEC035dn[(fibr071_vHMEC035dn$V2<=cutoff*fibr071_vHMEC035dn$ave2&fibr071_vHMEC035dn$V3>=cutoff2*fibr071_vHMEC035dn$ave3)|(fibr071_vHMEC035dn$V2>=cutoff2*fibr071_vHMEC035dn$ave2&fibr071_vHMEC035dn$V3<=cutoff*fibr071_vHMEC035dn$ave3),]
fibr071_vHMEC035dn_isoform<-fibr071_vHMEC035dn[fibr071_vHMEC035dn$ave2>Rmin&fibr071_vHMEC035dn$ave3>Rmin,]
fibr071_vHMEC035_isoform<-rbind(fibr071_vHMEC035up_isoform,fibr071_vHMEC035dn_isoform)
nrow(fibr071_vHMEC035_isoform)
# exclude gene level DE
DE_fibr071_vHMEC035up<-read.table("./gene/UP.fibrRM071_vHMECRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_fibr071_vHMEC035dn<-read.table("./gene/DN.fibrRM071_vHMECRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_fibr071_vHMEC035<-c(as.character(DE_fibr071_vHMEC035up$V1),as.character(DE_fibr071_vHMEC035dn$V1))
length(DE_fibr071_vHMEC035)
fibr071_vHMEC035_isoform<-fibr071_vHMEC035_isoform[!(fibr071_vHMEC035_isoform$id %in% DE_fibr071_vHMEC035), ]
nrow(fibr071_vHMEC035_isoform)
fibr071_vHMEC035_isoform_gene <- fibr071_vHMEC035_isoform[!duplicated(fibr071_vHMEC035_isoform$id), ]
nrow(fibr071_vHMEC035_isoform_gene)
write.table(fibr071_vHMEC035_isoform, file = "fibr071_vHMEC035_isoform.txt", sep = "\t", quote = F, row.names = F)
write.table(fibr071_vHMEC035_isoform_gene, file = "fibr071_vHMEC035_isoform_gene.txt", sep = "\t", quote = F, row.names = F)

