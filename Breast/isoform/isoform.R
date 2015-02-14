##################################################################################
# clutering on exon RPKM
# run on xhost
# /gsc/software/linux-x86_64-centos5/R-2.15.2/bin/R
# Function to color branches
colbranches <- function(n, col)
{
  a <- attributes(n) # Find the attributes of current node
  # Color edges with requested color
  attr(n, "edgePar") <- c(a$edgePar, list(col=col, lwd=3))
  n # Don't forget to return the node!
}

coltissue<-data.frame(tissue=c("fibr","lum","myo","stem-like","vHMEC"),
                      col=c("darkorchid2","red","springgreen","royalblue1","greenyellow"))

setwd("~/REMC/breast/")
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'
libs=c("A01029",
       "A01030",
       "A01031",
       "A17918",
       "A17919",
       "A17920",
       "A18760",
       "HS1187",
       "HS1188",
       "HS2263",
       "A18472",
       "A18761")
ID<-c("lumRM080","myoRM080","stem-likeRM080","lumRM084","myoRM084","stem-likeRM084","fibrRM070","lumRM035","myoRM035","vHMECRM035","fibrRM071","vHMECRM071")
lib=libs[1]
# exon_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.exn.A.rpkm",sep=""),as.is=T)
pc_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.A.rpkm.nc",sep=""),as.is=T)
# exon_lib<-exon_lib[is.element(exon_lib$V2,pc_lib$V1),]

# exon<-matrix(data=NA,ncol=length(libs),nrow=nrow(exon_lib))
gene<-matrix(data=NA,ncol=length(libs),nrow=nrow(pc_lib))

# No.of expressed exons 
# nexpress<-c(rep(0,times=length(libs)))
# colnames(exon)<-ID

# # myo lum expressed exons
# myolum<-c(8,1,4,9,2,5)
# libs<-libs[myolum]
# ID<-ID[myolum]
# express<-data.frame(lum35=0,lum80=0,lum84=0,myo35=0,myo80=0,myo84=0,ID=exon_lib$V1,gene=exon_lib$V2)

i=1
for(lib in libs){
#   exon_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.exn.A.rpkm",sep=""),as.is=T)
  pc_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.A.rpkm.nc",sep=""),as.is=T)
#   exon_lib<-exon_lib[is.element(exon_lib$V2,pc_lib$V1),]
#   exon[,i]<-exon_lib$V4
  gene[,i]<-exon_lib$V3
  
# # myo lum expressed exons
#   pc_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.A.rpkm.pc",sep=""),row.names=1)
#   exon_lib$ave<-pc_lib[exon_lib$V2,3]
#   exon_lib$ave[is.na(exon_lib$ave)]<-0
#   express[(exon_lib$V4>=0.1*exon_lib$ave&exon_lib$ave>=0.05),i]<-1
  
# # No.of expressed exons 
#   exon_lib<-na.omit(exon_lib)
#   exon_lib<-exon_lib[(exon_lib$V4>=0.1*exon_lib$ave&exon_lib$ave>=0.05),]
#   nexpress[i]<-nrow(exon_lib)

  i=i+1
}
# colnames(exon)<-ID
# write.table(exon,file="exons_all_libraries.txt",quote=F,sep="\t",col.names=F,row.names=F)
colnames(gene)<-ID
write.table(gene,file="nc_all_libraries.txt",quote=F,sep="\t",col.names=F,row.names=F)

# # No.of expressed exons 
# nexpress<-data.frame(Library=libs,Nexons=nexpress)
# nexpress$ID<-c("lumRM080","myoRM080","stem-likeRM080","lumRM084","myoRM084","stem-likeRM084","fibrRM070","lumRM035","myoRM035","vHMECRM035","fibrRM071","vHMECRM071")
# write.csv(nexpress,file="N.expressed.exons.csv",quote=F,row.names=F)
# 
# # myo lum expressed exons
# express$lum<-express$lum35&express$lum80&express$lum84
# express$myo<-express$myo35&express$myo80&express$myo84
# save(express,file="~/REMC/breast/bismark/express.Rdata")

# c<-cor(exon, method="spearman")
# write.table(c,file="exons_spearman.txt",quote=F,sep="\t",col.names=F,row.names=F)
c<-cor(gene, method="spearman")
write.table(c,file="nc_spearman.txt",quote=F,sep="\t",col.names=F,row.names=F)
d<-as.dist(1-c)
hc<-as.dendrogram(hclust(d,method="complete"))

# hc[[1]][[1]] = dendrapply(hc[[1]][[1]], colbranches, "greenyellow")
# hc[[1]][[2]] = dendrapply(hc[[1]][[2]], colbranches, "red")
# hc[[2]][[1]] = dendrapply(hc[[2]][[1]], colbranches, "darkorchid2")
# hc[[2]][[2]][[1]] = dendrapply(hc[[2]][[2]][[1]], colbranches, "springgreen")
# hc[[2]][[2]][[2]] = dendrapply(hc[[2]][[2]][[2]], colbranches, "royalblue1")
# hc[[2]][[2]][[2]][[2]][[1]] = dendrapply(hc[[2]][[2]][[2]][[2]][[1]], colbranches, "springgreen")

# pdf("cluster.exon.nc.pdf")
# par(mar=c(5,5,2,8))
# plot(hc,main="Exon RPKM clustering",horiz=TRUE)
# dev.off()
pdf("cluster.gene.nc.pdf")
par(mar=c(5,5,2,8))
plot(hc,main="nc-gene RPKM clustering",horiz=TRUE)
dev.off()
##################################################################################


##################################################################################
# getting the cutoff for expressing and non-expressing
# run on xhost
# /gsc/software/linux-x86_64-centos5/R-2.15.2/bin/R
setwd("~/REMC/breast/")
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'
libs=c("A01029",
       "A01030",
       "A01031",
       "A17918",
       "A17919",
       "A17920",
       "A18760",
       "HS1187",
       "HS1188",
       "HS2263",
       "A18472",
       "A18761")

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
legend("bottomright",libs,col=c(1:length(libs)),cex=0.8,lty=1)
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
setwd("~/REMC/breast/tissue/")
col<-c("Gene id", "rpkm1", "rpkm2", "p-value", "Multiple testing corrected p-value")
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'
cutoff=0.01
cutoff2=0.1 # sample specific exon: expressed in one sample(>0.1ave.RPKM) & not expressed in the other(<0.01ave.RPKM)
Rmin=0.005

lib1='HS2263'; cell1='vHMEC'; donor1='RM035';
lib2='HS1188'; cell2='myo'; donor2='RM035';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
vHMEC_myo_35up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
vHMEC_myo_35up$id<-unlist(strsplit(as.character(vHMEC_myo_35up$V1),'_'))[2*(1:nrow(vHMEC_myo_35up))]
vHMEC_myo_35up$ave2<-pc_lib1[vHMEC_myo_35up$id,3]
vHMEC_myo_35up$ave3<-pc_lib2[vHMEC_myo_35up$id,3]
vHMEC_myo_35up<-na.omit(vHMEC_myo_35up)
vHMEC_myo_35up<-vHMEC_myo_35up[(vHMEC_myo_35up$V2<=cutoff*vHMEC_myo_35up$ave2&vHMEC_myo_35up$V3>=cutoff2*vHMEC_myo_35up$ave3)|(vHMEC_myo_35up$V2>=cutoff2*vHMEC_myo_35up$ave2&vHMEC_myo_35up$V3<=cutoff*vHMEC_myo_35up$ave3),]
vHMEC_myo_35up_gene<-unique(vHMEC_myo_35up$id)
vHMEC_myo_35up_isoform<-vHMEC_myo_35up[vHMEC_myo_35up$ave2>Rmin&vHMEC_myo_35up$ave3>Rmin,]
write.csv(vHMEC_myo_35up,file="vHMEC_myo_35up.csv",quote=F,row.names=F)
write.csv(vHMEC_myo_35up_gene,file="vHMEC_myo_35up_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_myo_35up_isoform,file="vHMEC_myo_35up_isoform.csv",quote=F,row.names=F)
vHMEC_myo_35dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
vHMEC_myo_35dn$id<-unlist(strsplit(as.character(vHMEC_myo_35dn$V1),'_'))[2*(1:nrow(vHMEC_myo_35dn))]
vHMEC_myo_35dn$ave2<-pc_lib1[vHMEC_myo_35dn$id,3]
vHMEC_myo_35dn$ave3<-pc_lib2[vHMEC_myo_35dn$id,3]
vHMEC_myo_35dn<-na.omit(vHMEC_myo_35dn)
vHMEC_myo_35dn<-vHMEC_myo_35dn[(vHMEC_myo_35dn$V2<=cutoff*vHMEC_myo_35dn$ave2&vHMEC_myo_35dn$V3>=cutoff2*vHMEC_myo_35dn$ave3)|(vHMEC_myo_35dn$V2>=cutoff2*vHMEC_myo_35dn$ave2&vHMEC_myo_35dn$V3<=cutoff*vHMEC_myo_35dn$ave3),]
vHMEC_myo_35dn_gene<-unique(vHMEC_myo_35dn$id)
vHMEC_myo_35dn_isoform<-vHMEC_myo_35dn[vHMEC_myo_35dn$ave2>Rmin&vHMEC_myo_35dn$ave3>Rmin,]
write.csv(vHMEC_myo_35dn,file="vHMEC_myo_35dn.csv",quote=F,row.names=F)
write.csv(vHMEC_myo_35dn_gene,file="vHMEC_myo_35dn_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_myo_35dn_isoform,file="vHMEC_myo_35dn_isoform.csv",quote=F,row.names=F)
vHMEC_myo_35_isoform<-unique(c(vHMEC_myo_35up_isoform$id,vHMEC_myo_35dn_isoform$id))
write.csv(vHMEC_myo_35_isoform,file="vHMEC_myo_35_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_vHMEC_myo_35up<-read.table("UP.myoRM035_vHMECRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_vHMEC_myo_35dn<-read.table("DN.myoRM035_vHMECRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_vHMEC_myo_35<-c(as.character(DE_vHMEC_myo_35up$V1),as.character(DE_vHMEC_myo_35dn$V1))
setwd("~/REMC/breast/tissue/")
vHMEC_myo_35_isoform_only<-vHMEC_myo_35_isoform[!is.element(vHMEC_myo_35_isoform,DE_vHMEC_myo_35)]
write.csv(vHMEC_myo_35_isoform_only,file="vHMEC_myo_35_isoform_only.csv",quote=F,row.names=F)

lib1='HS2263'; cell1='vHMEC'; donor1='RM035';
lib2='HS1187'; cell2='lum'; donor2='RM035';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
vHMEC_lum_35up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
vHMEC_lum_35up$id<-unlist(strsplit(as.character(vHMEC_lum_35up$V1),'_'))[2*(1:nrow(vHMEC_lum_35up))]
vHMEC_lum_35up$ave2<-pc_lib1[vHMEC_lum_35up$id,3]
vHMEC_lum_35up$ave3<-pc_lib2[vHMEC_lum_35up$id,3]
vHMEC_lum_35up<-na.omit(vHMEC_lum_35up)
vHMEC_lum_35up<-vHMEC_lum_35up[(vHMEC_lum_35up$V2<=cutoff*vHMEC_lum_35up$ave2&vHMEC_lum_35up$V3>=cutoff2*vHMEC_lum_35up$ave3)|(vHMEC_lum_35up$V2>=cutoff2*vHMEC_lum_35up$ave2&vHMEC_lum_35up$V3<=cutoff*vHMEC_lum_35up$ave3),]
vHMEC_lum_35up_gene<-unique(vHMEC_lum_35up$id)
vHMEC_lum_35up_isoform<-vHMEC_lum_35up[vHMEC_lum_35up$ave2>Rmin&vHMEC_lum_35up$ave3>Rmin,]
write.csv(vHMEC_lum_35up,file="vHMEC_lum_35up.csv",quote=F,row.names=F)
write.csv(vHMEC_lum_35up_gene,file="vHMEC_lum_35up_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_lum_35up_isoform,file="vHMEC_lum_35up_isoform.csv",quote=F,row.names=F)
vHMEC_lum_35dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
vHMEC_lum_35dn$id<-unlist(strsplit(as.character(vHMEC_lum_35dn$V1),'_'))[2*(1:nrow(vHMEC_lum_35dn))]
vHMEC_lum_35dn$ave2<-pc_lib1[vHMEC_lum_35dn$id,3]
vHMEC_lum_35dn$ave3<-pc_lib2[vHMEC_lum_35dn$id,3]
vHMEC_lum_35dn<-na.omit(vHMEC_lum_35dn)
vHMEC_lum_35dn<-vHMEC_lum_35dn[(vHMEC_lum_35dn$V2<=cutoff*vHMEC_lum_35dn$ave2&vHMEC_lum_35dn$V3>=cutoff2*vHMEC_lum_35dn$ave3)|(vHMEC_lum_35dn$V2>=cutoff2*vHMEC_lum_35dn$ave2&vHMEC_lum_35dn$V3<=cutoff*vHMEC_lum_35dn$ave3),]
vHMEC_lum_35dn_gene<-unique(vHMEC_lum_35dn$id)
vHMEC_lum_35dn_isoform<-vHMEC_lum_35dn[vHMEC_lum_35dn$ave2>Rmin&vHMEC_lum_35dn$ave3>Rmin,]
write.csv(vHMEC_lum_35dn,file="vHMEC_lum_35dn.csv",quote=F,row.names=F)
write.csv(vHMEC_lum_35dn_gene,file="vHMEC_lum_35dn_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_lum_35dn_isoform,file="vHMEC_lum_35dn_isoform.csv",quote=F,row.names=F)
vHMEC_lum_35_isoform<-unique(c(vHMEC_lum_35up_isoform$id,vHMEC_lum_35dn_isoform$id))
write.csv(vHMEC_lum_35_isoform,file="vHMEC_lum_35_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_vHMEC_lum_35up<-read.table("UP.lumRM035_vHMECRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_vHMEC_lum_35dn<-read.table("DN.lumRM035_vHMECRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_vHMEC_lum_35<-c(as.character(DE_vHMEC_lum_35up$V1),as.character(DE_vHMEC_lum_35dn$V1))
setwd("~/REMC/breast/tissue/")
vHMEC_lum_35_isoform_only<-vHMEC_lum_35_isoform[!is.element(vHMEC_lum_35_isoform,DE_vHMEC_lum_35)]
write.csv(vHMEC_lum_35_isoform_only,file="vHMEC_lum_35_isoform_only.csv",quote=F,row.names=F)

lib1='HS1188'; cell1='myo'; donor1='RM035';
lib2='HS1187'; cell2='lum'; donor2='RM035';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
myo_lum_35up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_lum_35up$id<-unlist(strsplit(as.character(myo_lum_35up$V1),'_'))[2*(1:nrow(myo_lum_35up))]
myo_lum_35up$ave2<-pc_lib1[myo_lum_35up$id,3]
myo_lum_35up$ave3<-pc_lib2[myo_lum_35up$id,3]
myo_lum_35up<-na.omit(myo_lum_35up)
myo_lum_35up<-myo_lum_35up[(myo_lum_35up$V2<=cutoff*myo_lum_35up$ave2&myo_lum_35up$V3>=cutoff2*myo_lum_35up$ave3)|(myo_lum_35up$V2>=cutoff2*myo_lum_35up$ave2&myo_lum_35up$V3<=cutoff*myo_lum_35up$ave3),]
myo_lum_35up_gene<-unique(myo_lum_35up$id)
myo_lum_35up_isoform<-myo_lum_35up[myo_lum_35up$ave2>Rmin&myo_lum_35up$ave3>Rmin,]
write.csv(myo_lum_35up,file="myo_lum_35up.csv",quote=F,row.names=F)
write.csv(myo_lum_35up_gene,file="myo_lum_35up_gene.csv",quote=F,row.names=F)
write.csv(myo_lum_35up_isoform,file="myo_lum_35up_isoform.csv",quote=F,row.names=F)
myo_lum_35dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_lum_35dn$id<-unlist(strsplit(as.character(myo_lum_35dn$V1),'_'))[2*(1:nrow(myo_lum_35dn))]
myo_lum_35dn$ave2<-pc_lib1[myo_lum_35dn$id,3]
myo_lum_35dn$ave3<-pc_lib2[myo_lum_35dn$id,3]
myo_lum_35dn<-na.omit(myo_lum_35dn)
myo_lum_35dn<-myo_lum_35dn[(myo_lum_35dn$V2<=cutoff*myo_lum_35dn$ave2&myo_lum_35dn$V3>=cutoff2*myo_lum_35dn$ave3)|(myo_lum_35dn$V2>=cutoff2*myo_lum_35dn$ave2&myo_lum_35dn$V3<=cutoff*myo_lum_35dn$ave3),]
myo_lum_35dn_gene<-unique(myo_lum_35dn$id)
myo_lum_35dn_isoform<-myo_lum_35dn[myo_lum_35dn$ave2>Rmin&myo_lum_35dn$ave3>Rmin,]
write.csv(myo_lum_35dn,file="myo_lum_35dn.csv",quote=F,row.names=F)
write.csv(myo_lum_35dn_gene,file="myo_lum_35dn_gene.csv",quote=F,row.names=F)
write.csv(myo_lum_35dn_isoform,file="myo_lum_35dn_isoform.csv",quote=F,row.names=F)
myo_lum_35_isoform<-unique(c(myo_lum_35up_isoform$id,myo_lum_35dn_isoform$id))
write.csv(myo_lum_35_isoform,file="myo_lum_35_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_myo_lum_35up<-read.table("UP.lumRM035_myoRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_lum_35dn<-read.table("DN.lumRM035_myoRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_lum_35<-c(as.character(DE_myo_lum_35up$V1),as.character(DE_myo_lum_35dn$V1))
setwd("~/REMC/breast/tissue/")
myo_lum_35_isoform_only<-myo_lum_35_isoform[!is.element(myo_lum_35_isoform,DE_myo_lum_35)]
write.csv(myo_lum_35_isoform_only,file="myo_lum_35_isoform_only.csv",quote=F,row.names=F)

lib1='A01030'; cell1='myo'; donor1='RM080';
lib2='A01029'; cell2='lum'; donor2='RM080';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
myo_lum_80up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_lum_80up$id<-unlist(strsplit(as.character(myo_lum_80up$V1),'_'))[2*(1:nrow(myo_lum_80up))]
myo_lum_80up$ave2<-pc_lib1[myo_lum_80up$id,3]
myo_lum_80up$ave3<-pc_lib2[myo_lum_80up$id,3]
myo_lum_80up<-na.omit(myo_lum_80up)
myo_lum_80up<-myo_lum_80up[(myo_lum_80up$V2<=cutoff*myo_lum_80up$ave2&myo_lum_80up$V3>=cutoff2*myo_lum_80up$ave3)|(myo_lum_80up$V2>=cutoff2*myo_lum_80up$ave2&myo_lum_80up$V3<=cutoff*myo_lum_80up$ave3),]
myo_lum_80up_gene<-unique(myo_lum_80up$id)
myo_lum_80up_isoform<-myo_lum_80up[myo_lum_80up$ave2>Rmin&myo_lum_80up$ave3>Rmin,]
write.csv(myo_lum_80up,file="myo_lum_80up.csv",quote=F,row.names=F)
write.csv(myo_lum_80up_gene,file="myo_lum_80up_gene.csv",quote=F,row.names=F)
write.csv(myo_lum_80up_isoform,file="myo_lum_80up_isoform.csv",quote=F,row.names=F)
myo_lum_80dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_lum_80dn$id<-unlist(strsplit(as.character(myo_lum_80dn$V1),'_'))[2*(1:nrow(myo_lum_80dn))]
myo_lum_80dn$ave2<-pc_lib1[myo_lum_80dn$id,3]
myo_lum_80dn$ave3<-pc_lib2[myo_lum_80dn$id,3]
myo_lum_80dn<-na.omit(myo_lum_80dn)
myo_lum_80dn<-myo_lum_80dn[(myo_lum_80dn$V2<=cutoff*myo_lum_80dn$ave2&myo_lum_80dn$V3>=cutoff2*myo_lum_80dn$ave3)|(myo_lum_80dn$V2>=cutoff2*myo_lum_80dn$ave2&myo_lum_80dn$V3<=cutoff*myo_lum_80dn$ave3),]
myo_lum_80dn_gene<-unique(myo_lum_80dn$id)
myo_lum_80dn_isoform<-myo_lum_80dn[myo_lum_80dn$ave2>Rmin&myo_lum_80dn$ave3>Rmin,]
write.csv(myo_lum_80dn,file="myo_lum_80dn.csv",quote=F,row.names=F)
write.csv(myo_lum_80dn_gene,file="myo_lum_80dn_gene.csv",quote=F,row.names=F)
write.csv(myo_lum_80dn_isoform,file="myo_lum_80dn_isoform.csv",quote=F,row.names=F)
myo_lum_80_isoform<-unique(c(myo_lum_80up_isoform$id,myo_lum_80dn_isoform$id))
write.csv(myo_lum_80_isoform,file="myo_lum_80_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_myo_lum_80up<-read.table("UP.lumRM080_myoRM080.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_lum_80dn<-read.table("DN.lumRM080_myoRM080.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_lum_80<-c(as.character(DE_myo_lum_80up$V1),as.character(DE_myo_lum_80dn$V1))
setwd("~/REMC/breast/tissue/")
myo_lum_80_isoform_only<-myo_lum_80_isoform[!is.element(myo_lum_80_isoform,DE_myo_lum_80)]
write.csv(myo_lum_80_isoform_only,file="myo_lum_80_isoform_only.csv",quote=F,row.names=F)

lib1='A17919'; cell1='myo'; donor1='RM084';
lib2='A17918'; cell2='lum'; donor2='RM084';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
myo_lum_84up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_lum_84up$id<-unlist(strsplit(as.character(myo_lum_84up$V1),'_'))[2*(1:nrow(myo_lum_84up))]
myo_lum_84up$ave2<-pc_lib1[myo_lum_84up$id,3]
myo_lum_84up$ave3<-pc_lib2[myo_lum_84up$id,3]
myo_lum_84up<-na.omit(myo_lum_84up)
myo_lum_84up<-myo_lum_84up[(myo_lum_84up$V2<=cutoff*myo_lum_84up$ave2&myo_lum_84up$V3>=cutoff2*myo_lum_84up$ave3)|(myo_lum_84up$V2>=cutoff2*myo_lum_84up$ave2&myo_lum_84up$V3<=cutoff*myo_lum_84up$ave3),]
myo_lum_84up_gene<-unique(myo_lum_84up$id)
myo_lum_84up_isoform<-myo_lum_84up[myo_lum_84up$ave2>Rmin&myo_lum_84up$ave3>Rmin,]
write.csv(myo_lum_84up,file="myo_lum_84up.csv",quote=F,row.names=F)
write.csv(myo_lum_84up_gene,file="myo_lum_84up_gene.csv",quote=F,row.names=F)
write.csv(myo_lum_84up_isoform,file="myo_lum_84up_isoform.csv",quote=F,row.names=F)
myo_lum_84dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_lum_84dn$id<-unlist(strsplit(as.character(myo_lum_84dn$V1),'_'))[2*(1:nrow(myo_lum_84dn))]
myo_lum_84dn$ave2<-pc_lib1[myo_lum_84dn$id,3]
myo_lum_84dn$ave3<-pc_lib2[myo_lum_84dn$id,3]
myo_lum_84dn<-na.omit(myo_lum_84dn)
myo_lum_84dn<-myo_lum_84dn[(myo_lum_84dn$V2<=cutoff*myo_lum_84dn$ave2&myo_lum_84dn$V3>=cutoff2*myo_lum_84dn$ave3)|(myo_lum_84dn$V2>=cutoff2*myo_lum_84dn$ave2&myo_lum_84dn$V3<=cutoff*myo_lum_84dn$ave3),]
myo_lum_84dn_gene<-unique(myo_lum_84dn$id)
myo_lum_84dn_isoform<-myo_lum_84dn[myo_lum_84dn$ave2>Rmin&myo_lum_84dn$ave3>Rmin,]
write.csv(myo_lum_84dn,file="myo_lum_84dn.csv",quote=F,row.names=F)
write.csv(myo_lum_84dn_gene,file="myo_lum_84dn_gene.csv",quote=F,row.names=F)
write.csv(myo_lum_84dn_isoform,file="myo_lum_84dn_isoform.csv",quote=F,row.names=F)
myo_lum_84_isoform<-unique(c(myo_lum_84up_isoform$id,myo_lum_84dn_isoform$id))
write.csv(myo_lum_84_isoform,file="myo_lum_84_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_myo_lum_84up<-read.table("UP.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_lum_84dn<-read.table("DN.lumRM084_myoRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_lum_84<-c(as.character(DE_myo_lum_84up$V1),as.character(DE_myo_lum_84dn$V1))
setwd("~/REMC/breast/tissue/")
myo_lum_84_isoform_only<-myo_lum_84_isoform[!is.element(myo_lum_84_isoform,DE_myo_lum_84)]
write.csv(myo_lum_84_isoform_only,file="myo_lum_84_isoform_only.csv",quote=F,row.names=F)

lib1='A01030'; cell1='myo'; donor1='RM080';
lib2='A01031'; cell2='stem'; donor2='RM080';
myo_stem_80up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_stem_80up$id<-unlist(strsplit(as.character(myo_stem_80up$V1),'_'))[2*(1:nrow(myo_stem_80up))]
myo_stem_80up$ave2<-pc_lib1[myo_stem_80up$id,3]
myo_stem_80up$ave3<-pc_lib2[myo_stem_80up$id,3]
myo_stem_80up<-na.omit(myo_stem_80up)
myo_stem_80up<-myo_stem_80up[(myo_stem_80up$V2<=cutoff*myo_stem_80up$ave2&myo_stem_80up$V3>=cutoff2*myo_stem_80up$ave3)|(myo_stem_80up$V2>=cutoff2*myo_stem_80up$ave2&myo_stem_80up$V3<=cutoff*myo_stem_80up$ave3),]
myo_stem_80up_gene<-unique(myo_stem_80up$id)
myo_stem_80up_isoform<-myo_stem_80up[myo_stem_80up$ave2>Rmin&myo_stem_80up$ave3>Rmin,]
write.csv(myo_stem_80up,file="myo_stem_80up.csv",quote=F,row.names=F)
write.csv(myo_stem_80up_gene,file="myo_stem_80up_gene.csv",quote=F,row.names=F)
write.csv(myo_stem_80up_isoform,file="myo_stem_80up_isoform.csv",quote=F,row.names=F)
myo_stem_80dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_stem_80dn$id<-unlist(strsplit(as.character(myo_stem_80dn$V1),'_'))[2*(1:nrow(myo_stem_80dn))]
myo_stem_80dn$ave2<-pc_lib1[myo_stem_80dn$id,3]
myo_stem_80dn$ave3<-pc_lib2[myo_stem_80dn$id,3]
myo_stem_80dn<-na.omit(myo_stem_80dn)
myo_stem_80dn<-myo_stem_80dn[(myo_stem_80dn$V2<=cutoff*myo_stem_80dn$ave2&myo_stem_80dn$V3>=cutoff2*myo_stem_80dn$ave3)|(myo_stem_80dn$V2>=cutoff2*myo_stem_80dn$ave2&myo_stem_80dn$V3<=cutoff*myo_stem_80dn$ave3),]
myo_stem_80dn_gene<-unique(myo_stem_80dn$id)
myo_stem_80dn_isoform<-myo_stem_80dn[myo_stem_80dn$ave2>Rmin&myo_stem_80dn$ave3>Rmin,]
write.csv(myo_stem_80dn,file="myo_stem_80dn.csv",quote=F,row.names=F)
write.csv(myo_stem_80dn_gene,file="myo_stem_80dn_gene.csv",quote=F,row.names=F)
write.csv(myo_stem_80dn_isoform,file="myo_stem_80dn_isoform.csv",quote=F,row.names=F)
myo_stem_80_isoform<-unique(c(myo_stem_80up_isoform$id,myo_stem_80dn_isoform$id))
write.csv(myo_stem_80_isoform,file="myo_stem_80_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_myo_stem_80up<-read.table("UP.myoRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_stem_80dn<-read.table("DN.myoRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_stem_80<-c(as.character(DE_myo_stem_80up$V1),as.character(DE_myo_stem_80dn$V1))
setwd("~/REMC/breast/tissue/")
myo_stem_80_isoform_only<-myo_stem_80_isoform[!is.element(myo_stem_80_isoform,DE_myo_stem_80)]
write.csv(myo_stem_80_isoform_only,file="myo_stem_80_isoform_only.csv",quote=F,row.names=F)

lib1='A17919'; cell1='myo'; donor1='RM084';
lib2='A17920'; cell2='stem'; donor2='RM084';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
myo_stem_84up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_stem_84up$id<-unlist(strsplit(as.character(myo_stem_84up$V1),'_'))[2*(1:nrow(myo_stem_84up))]
myo_stem_84up$ave2<-pc_lib1[myo_stem_84up$id,3]
myo_stem_84up$ave3<-pc_lib2[myo_stem_84up$id,3]
myo_stem_84up<-na.omit(myo_stem_84up)
myo_stem_84up<-myo_stem_84up[(myo_stem_84up$V2<=cutoff*myo_stem_84up$ave2&myo_stem_84up$V3>=cutoff2*myo_stem_84up$ave3)|(myo_stem_84up$V2>=cutoff2*myo_stem_84up$ave2&myo_stem_84up$V3<=cutoff*myo_stem_84up$ave3),]
myo_stem_84up_gene<-unique(myo_stem_84up$id)
myo_stem_84up_isoform<-myo_stem_84up[myo_stem_84up$ave2>Rmin&myo_stem_84up$ave3>Rmin,]
write.csv(myo_stem_84up,file="myo_stem_84up.csv",quote=F,row.names=F)
write.csv(myo_stem_84up_gene,file="myo_stem_84up_gene.csv",quote=F,row.names=F)
write.csv(myo_stem_84up_isoform,file="myo_stem_84up_isoform.csv",quote=F,row.names=F)
myo_stem_84dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_stem_84dn$id<-unlist(strsplit(as.character(myo_stem_84dn$V1),'_'))[2*(1:nrow(myo_stem_84dn))]
myo_stem_84dn$ave2<-pc_lib1[myo_stem_84dn$id,3]
myo_stem_84dn$ave3<-pc_lib2[myo_stem_84dn$id,3]
myo_stem_84dn<-na.omit(myo_stem_84dn)
myo_stem_84dn<-myo_stem_84dn[(myo_stem_84dn$V2<=cutoff*myo_stem_84dn$ave2&myo_stem_84dn$V3>=cutoff2*myo_stem_84dn$ave3)|(myo_stem_84dn$V2>=cutoff2*myo_stem_84dn$ave2&myo_stem_84dn$V3<=cutoff*myo_stem_84dn$ave3),]
myo_stem_84dn_gene<-unique(myo_stem_84dn$id)
myo_stem_84dn_isoform<-myo_stem_84dn[myo_stem_84dn$ave2>Rmin&myo_stem_84dn$ave3>Rmin,]
write.csv(myo_stem_84dn,file="myo_stem_84dn.csv",quote=F,row.names=F)
write.csv(myo_stem_84dn_gene,file="myo_stem_84dn_gene.csv",quote=F,row.names=F)
write.csv(myo_stem_84dn_isoform,file="myo_stem_84dn_isoform.csv",quote=F,row.names=F)
myo_stem_84_isoform<-unique(c(myo_stem_84up_isoform$id,myo_stem_84dn_isoform$id))
write.csv(myo_stem_84_isoform,file="myo_stem_84_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_myo_stem_84up<-read.table("UP.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_stem_84dn<-read.table("DN.myoRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_stem_84<-c(as.character(DE_myo_stem_84up$V1),as.character(DE_myo_stem_84dn$V1))
setwd("~/REMC/breast/tissue/")
myo_stem_84_isoform_only<-myo_stem_84_isoform[!is.element(myo_stem_84_isoform,DE_myo_stem_84)]
write.csv(myo_stem_84_isoform_only,file="myo_stem_84_isoform_only.csv",quote=F,row.names=F)

lib1='A01029'; cell1='lum'; donor1='RM080';
lib2='A01031'; cell2='stem'; donor2='RM080';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
lum_stem_80up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
lum_stem_80up$id<-unlist(strsplit(as.character(lum_stem_80up$V1),'_'))[2*(1:nrow(lum_stem_80up))]
lum_stem_80up$ave2<-pc_lib1[lum_stem_80up$id,3]
lum_stem_80up$ave3<-pc_lib2[lum_stem_80up$id,3]
lum_stem_80up<-na.omit(lum_stem_80up)
lum_stem_80up<-lum_stem_80up[(lum_stem_80up$V2<=cutoff*lum_stem_80up$ave2&lum_stem_80up$V3>=cutoff2*lum_stem_80up$ave3)|(lum_stem_80up$V2>=cutoff2*lum_stem_80up$ave2&lum_stem_80up$V3<=cutoff*lum_stem_80up$ave3),]
lum_stem_80up_gene<-unique(lum_stem_80up$id)
lum_stem_80up_isoform<-lum_stem_80up[lum_stem_80up$ave2>Rmin&lum_stem_80up$ave3>Rmin,]
write.csv(lum_stem_80up,file="lum_stem_80up.csv",quote=F,row.names=F)
write.csv(lum_stem_80up_gene,file="lum_stem_80up_gene.csv",quote=F,row.names=F)
write.csv(lum_stem_80up_isoform,file="lum_stem_80up_isoform.csv",quote=F,row.names=F)
lum_stem_80dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
lum_stem_80dn$id<-unlist(strsplit(as.character(lum_stem_80dn$V1),'_'))[2*(1:nrow(lum_stem_80dn))]
lum_stem_80dn$ave2<-pc_lib1[lum_stem_80dn$id,3]
lum_stem_80dn$ave3<-pc_lib2[lum_stem_80dn$id,3]
lum_stem_80dn<-na.omit(lum_stem_80dn)
lum_stem_80dn<-lum_stem_80dn[(lum_stem_80dn$V2<=cutoff*lum_stem_80dn$ave2&lum_stem_80dn$V3>=cutoff2*lum_stem_80dn$ave3)|(lum_stem_80dn$V2>=cutoff2*lum_stem_80dn$ave2&lum_stem_80dn$V3<=cutoff*lum_stem_80dn$ave3),]
lum_stem_80dn_gene<-unique(lum_stem_80dn$id)
lum_stem_80dn_isoform<-lum_stem_80dn[lum_stem_80dn$ave2>Rmin&lum_stem_80dn$ave3>Rmin,]
write.csv(lum_stem_80dn,file="lum_stem_80dn.csv",quote=F,row.names=F)
write.csv(lum_stem_80dn_gene,file="lum_stem_80dn_gene.csv",quote=F,row.names=F)
write.csv(lum_stem_80dn_isoform,file="lum_stem_80dn_isoform.csv",quote=F,row.names=F)
lum_stem_80_isoform<-unique(c(lum_stem_80up_isoform$id,lum_stem_80dn_isoform$id))
write.csv(lum_stem_80_isoform,file="lum_stem_80_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_lum_stem_80up<-read.table("UP.lumRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25")
DE_lum_stem_80dn<-read.table("DN.lumRM080_stemRM080.FDR_0.015.rmin_0.005.Nmin_25")
DE_lum_stem_80<-c(as.character(DE_lum_stem_80up$V1),as.character(DE_lum_stem_80dn$V1))
setwd("~/REMC/breast/tissue/")
lum_stem_80_isoform_only<-lum_stem_80_isoform[!is.element(lum_stem_80_isoform,DE_lum_stem_80)]
write.csv(lum_stem_80_isoform_only,file="lum_stem_80_isoform_only.csv",quote=F,row.names=F)

lib1='A17918'; cell1='lum'; donor1='RM084';
lib2='A17920'; cell2='stem'; donor2='RM084';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
lum_stem_84up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
lum_stem_84up$id<-unlist(strsplit(as.character(lum_stem_84up$V1),'_'))[2*(1:nrow(lum_stem_84up))]
lum_stem_84up$ave2<-pc_lib1[lum_stem_84up$id,3]
lum_stem_84up$ave3<-pc_lib2[lum_stem_84up$id,3]
lum_stem_84up<-na.omit(lum_stem_84up)
lum_stem_84up<-lum_stem_84up[(lum_stem_84up$V2<=cutoff*lum_stem_84up$ave2&lum_stem_84up$V3>=cutoff2*lum_stem_84up$ave3)|(lum_stem_84up$V2>=cutoff2*lum_stem_84up$ave2&lum_stem_84up$V3<=cutoff*lum_stem_84up$ave3),]
lum_stem_84up_gene<-unique(lum_stem_84up$id)
lum_stem_84up_isoform<-lum_stem_84up[lum_stem_84up$ave2>Rmin&lum_stem_84up$ave3>Rmin,]
write.csv(lum_stem_84up,file="lum_stem_84up.csv",quote=F,row.names=F)
write.csv(lum_stem_84up_gene,file="lum_stem_84up_gene.csv",quote=F,row.names=F)
write.csv(lum_stem_84up_isoform,file="lum_stem_84up_isoform.csv",quote=F,row.names=F)
lum_stem_84dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
lum_stem_84dn$id<-unlist(strsplit(as.character(lum_stem_84dn$V1),'_'))[2*(1:nrow(lum_stem_84dn))]
lum_stem_84dn$ave2<-pc_lib1[lum_stem_84dn$id,3]
lum_stem_84dn$ave3<-pc_lib2[lum_stem_84dn$id,3]
lum_stem_84dn<-na.omit(lum_stem_84dn)
lum_stem_84dn<-lum_stem_84dn[(lum_stem_84dn$V2<=cutoff*lum_stem_84dn$ave2&lum_stem_84dn$V3>=cutoff2*lum_stem_84dn$ave3)|(lum_stem_84dn$V2>=cutoff2*lum_stem_84dn$ave2&lum_stem_84dn$V3<=cutoff*lum_stem_84dn$ave3),]
lum_stem_84dn_gene<-unique(lum_stem_84dn$id)
lum_stem_84dn_isoform<-lum_stem_84dn[lum_stem_84dn$ave2>Rmin&lum_stem_84dn$ave3>Rmin,]
write.csv(lum_stem_84dn,file="lum_stem_84dn.csv",quote=F,row.names=F)
write.csv(lum_stem_84dn_gene,file="lum_stem_84dn_gene.csv",quote=F,row.names=F)
write.csv(lum_stem_84dn_isoform,file="lum_stem_84dn_isoform.csv",quote=F,row.names=F)
lum_stem_84_isoform<-unique(c(lum_stem_84up_isoform$id,lum_stem_84dn_isoform$id))
write.csv(lum_stem_84_isoform,file="lum_stem_84_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_lum_stem_84up<-read.table("UP.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_lum_stem_84dn<-read.table("DN.lumRM084_stemRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_lum_stem_84<-c(as.character(DE_lum_stem_84up$V1),as.character(DE_lum_stem_84dn$V1))
setwd("~/REMC/breast/tissue/")
lum_stem_84_isoform_only<-lum_stem_84_isoform[!is.element(lum_stem_84_isoform,DE_lum_stem_84)]
write.csv(lum_stem_84_isoform_only,file="lum_stem_84_isoform_only.csv",quote=F,row.names=F)

lib1='A18761'; cell1='vHMEC'; donor1='RM071';
lib2='A18472'; cell2='fibro'; donor2='RM071';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
vHMEC_fibro_71up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
vHMEC_fibro_71up$id<-unlist(strsplit(as.character(vHMEC_fibro_71up$V1),'_'))[2*(1:nrow(vHMEC_fibro_71up))]
vHMEC_fibro_71up$ave2<-pc_lib1[vHMEC_fibro_71up$id,3]
vHMEC_fibro_71up$ave3<-pc_lib2[vHMEC_fibro_71up$id,3]
vHMEC_fibro_71up<-na.omit(vHMEC_fibro_71up)
vHMEC_fibro_71up<-vHMEC_fibro_71up[(vHMEC_fibro_71up$V2<=cutoff*vHMEC_fibro_71up$ave2&vHMEC_fibro_71up$V3>=cutoff2*vHMEC_fibro_71up$ave3)|(vHMEC_fibro_71up$V2>=cutoff2*vHMEC_fibro_71up$ave2&vHMEC_fibro_71up$V3<=cutoff*vHMEC_fibro_71up$ave3),]
vHMEC_fibro_71up_gene<-unique(vHMEC_fibro_71up$id)
vHMEC_fibro_71up_isoform<-vHMEC_fibro_71up[vHMEC_fibro_71up$ave2>Rmin&vHMEC_fibro_71up$ave3>Rmin,]
write.csv(vHMEC_fibro_71up,file="vHMEC_fibro_71up.csv",quote=F,row.names=F)
write.csv(vHMEC_fibro_71up_gene,file="vHMEC_fibro_71up_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_fibro_71up_isoform,file="vHMEC_fibro_71up_isoform.csv",quote=F,row.names=F)
vHMEC_fibro_71dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
vHMEC_fibro_71dn$id<-unlist(strsplit(as.character(vHMEC_fibro_71dn$V1),'_'))[2*(1:nrow(vHMEC_fibro_71dn))]
vHMEC_fibro_71dn$ave2<-pc_lib1[vHMEC_fibro_71dn$id,3]
vHMEC_fibro_71dn$ave3<-pc_lib2[vHMEC_fibro_71dn$id,3]
vHMEC_fibro_71dn<-na.omit(vHMEC_fibro_71dn)
vHMEC_fibro_71dn<-vHMEC_fibro_71dn[(vHMEC_fibro_71dn$V2<=cutoff*vHMEC_fibro_71dn$ave2&vHMEC_fibro_71dn$V3>=cutoff2*vHMEC_fibro_71dn$ave3)|(vHMEC_fibro_71dn$V2>=cutoff2*vHMEC_fibro_71dn$ave2&vHMEC_fibro_71dn$V3<=cutoff*vHMEC_fibro_71dn$ave3),]
vHMEC_fibro_71dn_gene<-unique(vHMEC_fibro_71dn$id)
vHMEC_fibro_71dn_isoform<-vHMEC_fibro_71dn[vHMEC_fibro_71dn$ave2>Rmin&vHMEC_fibro_71dn$ave3>Rmin,]
write.csv(vHMEC_fibro_71dn,file="vHMEC_fibro_71dn.csv",quote=F,row.names=F)
write.csv(vHMEC_fibro_71dn_gene,file="vHMEC_fibro_71dn_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_fibro_71dn_isoform,file="vHMEC_fibro_71dn_isoform.csv",quote=F,row.names=F)
vHMEC_fibro_71_isoform<-unique(c(vHMEC_fibro_71up_isoform$id,vHMEC_fibro_71dn_isoform$id))
write.csv(vHMEC_fibro_71_isoform,file="vHMEC_fibro_71_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_vHMEC_fibro_71up<-read.table("UP.fibrRM071_vHMECRM071.FDR_0.015.rmin_0.005.Nmin_25")
DE_vHMEC_fibro_71dn<-read.table("DN.fibrRM071_vHMECRM071.FDR_0.015.rmin_0.005.Nmin_25")
DE_vHMEC_fibro_71<-c(as.character(DE_vHMEC_fibro_71up$V1),as.character(DE_vHMEC_fibro_71dn$V1))
setwd("~/REMC/breast/tissue/")
vHMEC_fibro_71_isoform_only<-vHMEC_fibro_71_isoform[!is.element(vHMEC_fibro_71_isoform,DE_vHMEC_fibro_71)]
write.csv(vHMEC_fibro_71_isoform_only,file="vHMEC_fibro_71_isoform_only.csv",quote=F,row.names=F)

rm(lib1,cell1,donor1,lib2,cell2,donor2,pc_lib1,pc_lib2)
save.image(file="REMC.isoform.breast.RData")

##################################################################################

setwd("~/REMC/breast/tissue/")
load("REMC.isoform.breast.RData")

# side by side barplot comparing isoforms vs DE genes
# library(ggplot2)
# library(reshape)
# df = melt(data.frame(Isoform=c(3204,2963,1829,3240,3038.5,3254), 
#                      DE_gene=c(994,1474,1845,1211,1339.5,1522), 
#                      tissue=c("vHMEC_myo", "vHMEC_lum","vHMEC_fibr","myo_lum","myo_stem-like","lum_stem-like")))
# ggplot(df, aes(tissue, value, fill=variable)) + geom_bar(position="dodge")
summary<-data.frame(tissue=c("vHMEC_myo", "vHMEC_lum","vHMEC_fibr","myo_lum","myo_stem-like","lum_stem-like"),
                    Isoform=c(3204,2963,1829,2349,2408.5,2387), 
                    DE_gene=c(994,1474,1845,1211,1339.5,1522),
                    Isoform_max=c(3204,2963,1829,2381,2427,2429),
                    Isoform_min=c(3204,2963,1829,2325,2390,2345),
                    DE_max=c(994,1474,1845,1255,1401,1569),
                    DE_min=c(994,1474,1845,1171,1278,1475),
                    Nexons1=c(101826.5,101826.5,101826.5,113100.6667,113100.6667,112374.3333),
                    Nexons2=c(113100.6667,112374.3333,102248.5,112374.3333,112586.5,112586.5))
error.bar <- function(x, y, upper, lower, length=0.1,...){
  arrows(x,upper,x,lower,angle=90,code=3,length=length,lwd=1.8,...)
}
yy<-matrix(c(summary$Isoform,summary$DE_gene),2,6,byrow=T)
ee_up<-matrix(c(summary$Isoform_max,summary$DE_max),2,6,byrow=T)
ee_low<-matrix(c(summary$Isoform_min,summary$DE_min),2,6,byrow=T)
eNexons1_up<-c(106163,106163,106163,117809,117809,117201)
eNexons1_low<-c(97490,97490,97490,109462,109462,108666)
eNexons2_up<-c(117809,117201,104419,117201,114900,114900)
eNexons2_low<-c(109462,108666,100078,108666,110273,110273)

coltissue<-data.frame(tissue=c("fibr","lum","myo","stem-like","vHMEC"),
                      col=c(rgb(110,190,30, maxColorValue = 255), rgb(200,50,0, maxColorValue = 255),rgb(50,200,50, maxColorValue = 255),rgb(0,100,200, maxColorValue = 255),rgb(150,70,50, maxColorValue = 255)))
colNexons1<-c(as.character(coltissue$col[5]),as.character(coltissue$col[5]),as.character(coltissue$col[5]),as.character(coltissue$col[3]),as.character(coltissue$col[3]),as.character(coltissue$col[2]))
colNexons2<-c(as.character(coltissue$col[3]),as.character(coltissue$col[2]),as.character(coltissue$col[1]),as.character(coltissue$col[2]),as.character(coltissue$col[4]),as.character(coltissue$col[4]))

pdf("~/快盘/REMC/figures/isoform_DE.pdf",width=13,height=8)
par(mar = c(5,5,2,15), xpd=TRUE,cex.axis=1.2,cex.lab=1.5)
barx<-barplot(yy, beside=TRUE,col=c("blue","red"), ylim=c(0,5000),cex.names=1.2,names.arg=summary$tissue, axis.lty=1, xlab="Cell type", ylab="No. of genes")  
error.bar(barx,yy,ee_up,ee_low)
par(new = T,cex.axis=1.2,cex.lab=1.5)
plot(c(0.2,1.5,2.6,3.9,5.1,6.35),summary$Nexons1,col=colNexons1,axes=F,xlab=NA,ylab=NA,type="p",pch=16,cex=1.5,ylim=c(0,120000),xlim=c(0,7))
error.bar(c(0.2,1.5,2.6,3.9,5.1,6.35),summary$Nexon1,eNexons1_up,eNexons1_low,col=colNexons1,length=0.05)
points(c(0.55,1.85,2.95,4.25,5.45,6.75),summary$Nexons2,col=colNexons2,xlab=NA,ylab=NA,type="p",pch=16,cex=1.5)
error.bar(c(0.55,1.85,2.95,4.25,5.45,6.75),summary$Nexons2,eNexons2_up,eNexons2_low,col=colNexons2,length=0.05)
axis(side = 4,cex=1.2)
mtext(side = 4, line = 3, "No. of expressed exons",cex=1.5)
legend("topright",c("DE isoforms","DE genes"),inset=c(-0.33,0),col=c("blue","red"),lty=c(1,1),lwd=c(8,8),cex=1.1,pt.cex=1.2)
legend("right",as.character(coltissue$tissue),inset=c(-0.33,0),col=as.character(coltissue$col),pch=19,cex=1.2,title="No. of expressed exons",pt.cex=1.2)
dev.off()

# hiearchal clustering
isoform<-union(union(union(union(vHMEC_myo_35_isoform_only,vHMEC_lum_35_isoform_only),union(vHMEC_fibro_71_isoform_only,myo_lum_35_isoform_only)),union(union(myo_lum_80_isoform_only,myo_lum_84_isoform_only),union(myo_stem_80_isoform_only,myo_stem_84_isoform_only))),union(lum_stem_80_isoform_only,lum_stem_84_isoform_only))
isoform<-data.frame(ID=isoform,vHMEC_myo_35=0,vHMEC_lum_35=0,vHMEC_fibro_71=0,myo_lum_35=0,myo_lum_80=0,myo_lum_84=0,myo_stem_80=0,myo_stem_84=0,lum_stem_80=0,lum_stem_84=0)
rownames(isoform)<-isoform$ID
isoform[vHMEC_myo_35_isoform_only,2]<-1
isoform[vHMEC_lum_35_isoform_only,3]<-1
isoform[vHMEC_fibro_71_isoform_only,4]<-1
isoform[myo_lum_35_isoform_only,5]<-1
isoform[myo_lum_80_isoform_only,6]<-1
isoform[myo_lum_84_isoform_only,7]<-1
isoform[myo_stem_80_isoform_only,8]<-1
isoform[myo_stem_84_isoform_only,9]<-1
isoform[lum_stem_80_isoform_only,10]<-1
isoform[lum_stem_84_isoform_only,11]<-1
c<-cor(as.matrix(isoform[,2:11]), method="spearman")
d<-as.dist(1-c)
pdf("dendrogram.pdf")
plot(hclust(d,method="ward"),main="Isoform clustering",xlab=NA)
dev.off()

# summarize by tissue type
vHMEC_myo_up=as.character(vHMEC_myo_35up$V1)
vHMEC_myo_up_gene=vHMEC_myo_35up_gene
vHMEC_myo_dn=as.character(vHMEC_myo_35dn$V1)
vHMEC_myo_dn_gene=vHMEC_myo_35dn_gene
vHMEC_myo_isoform=vHMEC_myo_35_isoform
vHMEC_myo_isoform_only=vHMEC_myo_35_isoform_only
write.csv(vHMEC_myo_isoform_only,file="vHMEC_myo_isoform_only.csv",quote=F,row.names=F)

vHMEC_lum_up=as.character(vHMEC_lum_35up$V1)
vHMEC_lum_up_gene=vHMEC_lum_35up_gene
vHMEC_lum_dn=as.character(vHMEC_lum_35dn$V1)
vHMEC_lum_dn_gene=vHMEC_lum_35dn_gene
vHMEC_lum_isoform=vHMEC_lum_35_isoform
vHMEC_lum_isoform_only=vHMEC_lum_35_isoform_only
write.csv(vHMEC_lum_isoform_only,file="vHMEC_lum_isoform_only.csv",quote=F,row.names=F)

vHMEC_fibro_up=as.character(vHMEC_fibro_71up$V1)
vHMEC_fibro_up_gene=vHMEC_fibro_71up_gene
vHMEC_fibro_dn=as.character(vHMEC_fibro_71dn$V1)
vHMEC_fibro_dn_gene=vHMEC_fibro_71dn_gene
vHMEC_fibro_isoform=vHMEC_fibro_71_isoform
vHMEC_fibro_isoform_only=vHMEC_fibro_71_isoform_only
write.csv(vHMEC_fibro_isoform_only,file="vHMEC_fibro_isoform_only.csv",quote=F,row.names=F)

myo_lum_up=intersect(as.character(myo_lum_35up$V1),intersect(as.character(myo_lum_80up$V1),as.character(myo_lum_84up$V1)))
myo_lum_up_gene=intersect(myo_lum_35up_gene,intersect(myo_lum_80up_gene,myo_lum_84up_gene))
myo_lum_dn=intersect(as.character(myo_lum_35dn$V1),intersect(as.character(myo_lum_80dn$V1),as.character(myo_lum_84dn$V1)))
myo_lum_dn_gene=intersect(myo_lum_35dn_gene,intersect(myo_lum_80dn_gene,myo_lum_84dn_gene))
myo_lum_isoform=intersect(myo_lum_35_isoform,intersect(myo_lum_80_isoform,myo_lum_84_isoform))
myo_lum_isoform_only=intersect(myo_lum_35_isoform_only,intersect(myo_lum_80_isoform_only,myo_lum_84_isoform_only))
write.csv(myo_lum_isoform_only,file="myo_lum_isoform_only.csv",quote=F,row.names=F)

myo_stem_up=intersect(as.character(myo_stem_80up$V1),as.character(myo_stem_84up$V1))
myo_stem_up_gene=intersect(myo_stem_80up_gene,myo_stem_84up_gene)
myo_stem_dn=intersect(as.character(myo_stem_80dn$V1),as.character(myo_stem_84dn$V1))
myo_stem_dn_gene=intersect(myo_stem_80dn_gene,myo_stem_84dn_gene)
myo_stem_isoform=intersect(myo_stem_80_isoform,myo_stem_84_isoform)
myo_stem_isoform_only=intersect(myo_stem_80_isoform_only,myo_stem_84_isoform_only)
write.csv(myo_stem_isoform_only,file="myo_stem_isoform_only.csv",quote=F,row.names=F)

lum_stem_up=intersect(as.character(lum_stem_80up$V1),as.character(lum_stem_84up$V1))
lum_stem_up_gene=intersect(lum_stem_80up_gene,lum_stem_84up_gene)
lum_stem_dn=intersect(as.character(lum_stem_80dn$V1),as.character(lum_stem_84dn$V1))
lum_stem_dn_gene=intersect(lum_stem_80dn_gene,lum_stem_84dn_gene)
lum_stem_isoform=intersect(lum_stem_80_isoform,lum_stem_84_isoform)
lum_stem_isoform_only=intersect(lum_stem_80_isoform_only,lum_stem_84_isoform_only)
write.csv(lum_stem_isoform_only,file="lum_stem_isoform_only.csv",quote=F,row.names=F)

myo_isoform<-intersect(intersect(vHMEC_myo_isoform,myo_lum_isoform),myo_stem_isoform)
myo_isoform_only<-intersect(intersect(vHMEC_myo_isoform_only,myo_lum_isoform_only),myo_stem_isoform_only)
lum_isoform<-intersect(intersect(vHMEC_lum_isoform,myo_lum_isoform),lum_stem_isoform)
lum_isoform_only<-intersect(intersect(vHMEC_lum_isoform_only,myo_lum_isoform_only),lum_stem_isoform_only)

write.csv(myo_isoform_only,file="myo_isoform_only.csv",quote=F,row.names=F)
write.csv(lum_isoform_only,file="lum_isoform_only.csv",quote=F,row.names=F)
write.csv(vHMEC_myo_up,file="vHMEC_myo_up.csv",quote=F,row.names=F)
write.csv(vHMEC_myo_up_gene,file="vHMEC_myo_up_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_myo_dn,file="vHMEC_myo_dn.csv",quote=F,row.names=F)
write.csv(vHMEC_myo_dn_gene,file="vHMEC_myo_dn_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_lum_up,file="vHMEC_lum_up.csv",quote=F,row.names=F)
write.csv(vHMEC_lum_up_gene,file="vHMEC_lum_up_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_lum_dn,file="vHMEC_lum_dn.csv",quote=F,row.names=F)
write.csv(vHMEC_lum_dn_gene,file="vHMEC_lum_dn_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_fibro_up,file="vHMEC_fibro_up.csv",quote=F,row.names=F)
write.csv(vHMEC_fibro_up_gene,file="vHMEC_fibro_up_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_fibro_dn,file="vHMEC_fibro_dn.csv",quote=F,row.names=F)
write.csv(vHMEC_fibro_dn_gene,file="vHMEC_fibro_dn_gene.csv",quote=F,row.names=F)
write.csv(myo_lum_up,file="myo_lum_up.csv",quote=F,row.names=F)
write.csv(myo_lum_up_gene,file="myo_lum_up_gene.csv",quote=F,row.names=F)
write.csv(myo_lum_dn,file="myo_lum_dn.csv",quote=F,row.names=F)
write.csv(myo_lum_dn_gene,file="myo_lum_dn_gene.csv",quote=F,row.names=F)
write.csv(myo_stem_up,file="myo_stem_up.csv",quote=F,row.names=F)
write.csv(myo_stem_up_gene,file="myo_stem_up_gene.csv",quote=F,row.names=F)
write.csv(myo_stem_dn,file="myo_stem_dn.csv",quote=F,row.names=F)
write.csv(myo_stem_dn_gene,file="myo_stem_dn_gene.csv",quote=F,row.names=F)
write.csv(lum_stem_up,file="lum_stem_up.csv",quote=F,row.names=F)
write.csv(lum_stem_up_gene,file="lum_stem_up_gene.csv",quote=F,row.names=F)
write.csv(lum_stem_dn,file="lum_stem_dn.csv",quote=F,row.names=F)
write.csv(lum_stem_dn_gene,file="lum_stem_dn_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_myo_isoform,file="vHMEC_myo_isoform.csv",quote=F,row.names=F)
write.csv(vHMEC_lum_isoform,file="vHMEC_lum_isoform.csv",quote=F,row.names=F)
write.csv(vHMEC_fibro_isoform,file="vHMEC_fibro_isoform.csv",quote=F,row.names=F)
write.csv(myo_lum_isoform,file="myo_lum_isoform.csv",quote=F,row.names=F)
write.csv(myo_stem_isoform,file="myo_stem_isoform.csv",quote=F,row.names=F)
write.csv(lum_stem_isoform,file="lum_stem_isoform.csv",quote=F,row.names=F)
write.csv(myo_isoform,file="myo_isoform.csv",quote=F,row.names=F)
write.csv(lum_isoform,file="lum_isoform.csv",quote=F,row.names=F)

save.image(file="REMC.isoform.breast.RData")

##################################################################################

setwd("~/REMC/breast/tissue/")
load("REMC.isoform.breast.RData")

#get gene information
setwd("~/hg19/")
gene<-read.table(file="hg19v65_genes",as.is=T,head=F)
colnames(gene)<-c("Ensembl","chr","start","end","strand","type","name","description")
gene<-gene[!duplicated(gene$Ensembl),]
rownames(gene)<-gene$Ensembl
setwd("~/REMC/breast/tissue/")

vHMEC_myo_isoform<-data.frame(ID=vHMEC_myo_isoform,Name=gene[vHMEC_myo_isoform,7],Description=gene[vHMEC_myo_isoform,8])
vHMEC_myo_isoform_only<-data.frame(ID=vHMEC_myo_isoform_only,Name=gene[vHMEC_myo_isoform_only,7],Description=gene[vHMEC_myo_isoform_only,8])
write.csv(vHMEC_myo_isoform,file="vHMEC_myo_isoform.csv",quote=F,row.names=F)
write.csv(vHMEC_myo_isoform_only,file="vHMEC_myo_isoform_only.csv",quote=F,row.names=F)

vHMEC_lum_isoform<-data.frame(ID=vHMEC_lum_isoform,Name=gene[vHMEC_lum_isoform,7],Description=gene[vHMEC_lum_isoform,8])
vHMEC_lum_isoform_only<-data.frame(ID=vHMEC_lum_isoform_only,Name=gene[vHMEC_lum_isoform_only,7],Description=gene[vHMEC_lum_isoform_only,8])
write.csv(vHMEC_lum_isoform,file="vHMEC_lum_isoform.csv",quote=F,row.names=F)
write.csv(vHMEC_lum_isoform_only,file="vHMEC_lum_isoform_only.csv",quote=F,row.names=F)

vHMEC_fibro_isoform<-data.frame(ID=vHMEC_fibro_isoform,Name=gene[vHMEC_fibro_isoform,7],Description=gene[vHMEC_fibro_isoform,8])
vHMEC_fibro_isoform_only<-data.frame(ID=vHMEC_fibro_isoform_only,Name=gene[vHMEC_fibro_isoform_only,7],Description=gene[vHMEC_fibro_isoform_only,8])
write.csv(vHMEC_fibro_isoform,file="vHMEC_fibro_isoform.csv",quote=F,row.names=F)
write.csv(vHMEC_fibro_isoform_only,file="vHMEC_fibro_isoform_only.csv",quote=F,row.names=F)

myo_lum_isoform<-data.frame(ID=myo_lum_isoform,Name=gene[myo_lum_isoform,7],Description=gene[myo_lum_isoform,8])
myo_lum_isoform_only<-data.frame(ID=myo_lum_isoform_only,Name=gene[myo_lum_isoform_only,7],Description=gene[myo_lum_isoform_only,8])
write.csv(myo_lum_isoform,file="myo_lum_isoform.csv",quote=F,row.names=F)
write.csv(myo_lum_isoform_only,file="myo_lum_isoform_only.csv",quote=F,row.names=F)

myo_stem_isoform<-data.frame(ID=myo_stem_isoform,Name=gene[myo_stem_isoform,7],Description=gene[myo_stem_isoform,8])
myo_stem_isoform_only<-data.frame(ID=myo_stem_isoform_only,Name=gene[myo_stem_isoform_only,7],Description=gene[myo_stem_isoform_only,8])
write.csv(myo_stem_isoform,file="myo_stem_isoform.csv",quote=F,row.names=F)
write.csv(myo_stem_isoform_only,file="myo_stem_isoform_only.csv",quote=F,row.names=F)

lum_stem_isoform<-data.frame(ID=lum_stem_isoform,Name=gene[lum_stem_isoform,7],Description=gene[lum_stem_isoform,8])
lum_stem_isoform_only<-data.frame(ID=lum_stem_isoform_only,Name=gene[lum_stem_isoform_only,7],Description=gene[lum_stem_isoform_only,8])
write.csv(lum_stem_isoform,file="lum_stem_isoform.csv",quote=F,row.names=F)
write.csv(lum_stem_isoform_only,file="lum_stem_isoform_only.csv",quote=F,row.names=F)

save.image(file="REMC.isoform.breast.RData")

##################################################################################
# annotate exons
setwd("~/REMC/breast/tissue/")
load("REMC.isoform.breast.RData")

# vHMEC vs myo 
exonannot<-function(gene){
  exonsup<-as.character(vHMEC_myo_35up_isoform[vHMEC_myo_35up_isoform$id==gene,]$V1)
  exonsup<-paste(unlist(strsplit(exonsup,'<'))[2*(1:length(exonsup))-1],sep="",collapse=";")
  exonsdn<-as.character(vHMEC_myo_35dn_isoform[vHMEC_myo_35dn_isoform$id==gene,]$V1)
  exonsdn<-paste(unlist(strsplit(exonsdn,'<'))[2*(1:length(exonsdn))-1],sep="",collapse=";")
  return(c(exonsup,exonsdn))
}
vexonannot<-Vectorize(exonannot,vectorize.args="gene",SIMPLIFY=T)
y<-vexonannot(vHMEC_myo_isoform_only$ID)
vHMEC_myo_isoform_only$up<-y[1,]
vHMEC_myo_isoform_only$dn<-y[2,]
write.table(vHMEC_myo_isoform_only,file="vHMEC_myo_isoform_only.txt",sep="\t",row.names=F,col.names=T)

# vHMEC vs lum 
exonannot<-function(gene){
  exonsup<-as.character(vHMEC_lum_35up_isoform[vHMEC_lum_35up_isoform$id==gene,]$V1)
  exonsup<-paste(unlist(strsplit(exonsup,'<'))[2*(1:length(exonsup))-1],sep="",collapse=";")
  exonsdn<-as.character(vHMEC_lum_35dn_isoform[vHMEC_lum_35dn_isoform$id==gene,]$V1)
  exonsdn<-paste(unlist(strsplit(exonsdn,'<'))[2*(1:length(exonsdn))-1],sep="",collapse=";")
  return(c(exonsup,exonsdn))
}
vexonannot<-Vectorize(exonannot,vectorize.args="gene",SIMPLIFY=T)
y<-vexonannot(vHMEC_lum_isoform_only$ID)
vHMEC_lum_isoform_only$up<-y[1,]
vHMEC_lum_isoform_only$dn<-y[2,]
write.table(vHMEC_lum_isoform_only,file="vHMEC_lum_isoform_only.txt",sep="\t",row.names=F,col.names=T)

# vHMEC vs fibro 
exonannot<-function(gene){
  exonsup<-as.character(vHMEC_fibro_71up_isoform[vHMEC_fibro_71up_isoform$id==gene,]$V1)
  exonsup<-paste(unlist(strsplit(exonsup,'<'))[2*(1:length(exonsup))-1],sep="",collapse=";")
  exonsdn<-as.character(vHMEC_fibro_71dn_isoform[vHMEC_fibro_71dn_isoform$id==gene,]$V1)
  exonsdn<-paste(unlist(strsplit(exonsdn,'<'))[2*(1:length(exonsdn))-1],sep="",collapse=";")
  return(c(exonsup,exonsdn))
}
vexonannot<-Vectorize(exonannot,vectorize.args="gene",SIMPLIFY=T)
y<-vexonannot(vHMEC_fibro_isoform_only$ID)
vHMEC_fibro_isoform_only$up<-y[1,]
vHMEC_fibro_isoform_only$dn<-y[2,]
write.table(vHMEC_fibro_isoform_only,file="vHMEC_fibro_isoform_only.txt",sep="\t",row.names=F,col.names=T)

# myo vs stem
exonannot<-function(gene){
  exonsup80<-as.character(myo_stem_80up_isoform[myo_stem_80up_isoform$id==gene,]$V1)
  exonsdn80<-as.character(myo_stem_80dn_isoform[myo_stem_80dn_isoform$id==gene,]$V1)
  exonsup84<-as.character(myo_stem_84up_isoform[myo_stem_84up_isoform$id==gene,]$V1)
  exonsdn84<-as.character(myo_stem_84dn_isoform[myo_stem_84dn_isoform$id==gene,]$V1)
  exonsup<-unique(c(exonsup80,exonsup84))
  exonsdn<-unique(c(exonsdn80,exonsdn84))
  exonsup<-paste(unlist(strsplit(exonsup,'<'))[2*(1:length(exonsup))-1],sep="",collapse=";")
  exonsdn<-paste(unlist(strsplit(exonsdn,'<'))[2*(1:length(exonsdn))-1],sep="",collapse=";")
  return(c(exonsup,exonsdn))
}
vexonannot<-Vectorize(exonannot,vectorize.args="gene",SIMPLIFY=T)
y<-vexonannot(myo_stem_isoform_only$ID)
myo_stem_isoform_only$up<-y[1,]
myo_stem_isoform_only$dn<-y[2,]
write.table(myo_stem_isoform_only,file="myo_stem_isoform_only.txt",sep="\t",row.names=F,col.names=T)

# lum vs stem
exonannot<-function(gene){
  exonsup80<-as.character(lum_stem_80up_isoform[lum_stem_80up_isoform$id==gene,]$V1)
  exonsdn80<-as.character(lum_stem_80dn_isoform[lum_stem_80dn_isoform$id==gene,]$V1)
  exonsup84<-as.character(lum_stem_84up_isoform[lum_stem_84up_isoform$id==gene,]$V1)
  exonsdn84<-as.character(lum_stem_84dn_isoform[lum_stem_84dn_isoform$id==gene,]$V1)
  exonsup<-unique(c(exonsup80,exonsup84))
  exonsdn<-unique(c(exonsdn80,exonsdn84))
  exonsup<-paste(unlist(strsplit(exonsup,'<'))[2*(1:length(exonsup))-1],sep="",collapse=";")
  exonsdn<-paste(unlist(strsplit(exonsdn,'<'))[2*(1:length(exonsdn))-1],sep="",collapse=";")
  return(c(exonsup,exonsdn))
}
vexonannot<-Vectorize(exonannot,vectorize.args="gene",SIMPLIFY=T)
y<-vexonannot(lum_stem_isoform_only$ID)
lum_stem_isoform_only$up<-y[1,]
lum_stem_isoform_only$dn<-y[2,]
write.table(lum_stem_isoform_only,file="lum_stem_isoform_only.txt",sep="\t",row.names=F,col.names=T)

# myo vs lum
exonannot<-function(gene){
  exonsup35<-as.character(myo_lum_35up_isoform[myo_lum_35up_isoform$id==gene,]$V1)
  exonsdn35<-as.character(myo_lum_35dn_isoform[myo_lum_35dn_isoform$id==gene,]$V1)
  exonsup80<-as.character(myo_lum_80up_isoform[myo_lum_80up_isoform$id==gene,]$V1)
  exonsdn80<-as.character(myo_lum_80dn_isoform[myo_lum_80dn_isoform$id==gene,]$V1)
  exonsup84<-as.character(myo_lum_84up_isoform[myo_lum_84up_isoform$id==gene,]$V1)
  exonsdn84<-as.character(myo_lum_84dn_isoform[myo_lum_84dn_isoform$id==gene,]$V1)
  exonsup<-unique(c(exonsup35,exonsup80,exonsup84))
  exonsdn<-unique(c(exonsdn35,exonsdn80,exonsdn84))
  exonsup<-paste(unlist(strsplit(exonsup,'<'))[2*(1:length(exonsup))-1],sep="",collapse=";")
  exonsdn<-paste(unlist(strsplit(exonsdn,'<'))[2*(1:length(exonsdn))-1],sep="",collapse=";")
  return(c(exonsup,exonsdn))
}
vexonannot<-Vectorize(exonannot,vectorize.args="gene",SIMPLIFY=T)
y<-vexonannot(myo_lum_isoform_only$ID)
myo_lum_isoform_only$up<-y[1,]
myo_lum_isoform_only$dn<-y[2,]
write.table(myo_lum_isoform_only,file="myo_lum_isoform_only.txt",sep="\t",row.names=F,col.names=T)

save.image(file="REMC.isoform.breast.RData")

##################################################################################
## breast: individual specific (RM035,070,071,080,084). 
##################################################################################
# run on xhost
# /gsc/software/linux-x86_64-centos5/R-2.15.2/bin/R
setwd("~/REMC/breast/individual/")
col<-c("Gene id", "rpkm1", "rpkm2", "p-value", "Multiple testing corrected p-value")
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'
cutoff=0.01
cutoff2=0.1
Rmin=0.005

lib1='HS2263'; cell1='vHMEC'; donor1='RM035';
lib2='A18761'; cell2='vHMEC'; donor2='RM071';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
vHMEC_35_71up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
vHMEC_35_71up$id<-unlist(strsplit(as.character(vHMEC_35_71up$V1),'_'))[2*(1:nrow(vHMEC_35_71up))]
vHMEC_35_71up$ave2<-pc_lib1[vHMEC_35_71up$id,3]
vHMEC_35_71up$ave3<-pc_lib2[vHMEC_35_71up$id,3]
vHMEC_35_71up<-na.omit(vHMEC_35_71up)
vHMEC_35_71up<-vHMEC_35_71up[(vHMEC_35_71up$V2<=cutoff*vHMEC_35_71up$ave2&vHMEC_35_71up$V3>=cutoff2*vHMEC_35_71up$ave3)|(vHMEC_35_71up$V2>=cutoff2*vHMEC_35_71up$ave2&vHMEC_35_71up$V3<=cutoff*vHMEC_35_71up$ave3),]
vHMEC_35_71up_gene<-unique(vHMEC_35_71up$id)
vHMEC_35_71up_isoform<-vHMEC_35_71up[vHMEC_35_71up$ave2>Rmin&vHMEC_35_71up$ave3>Rmin,]
write.csv(vHMEC_35_71up,file="vHMEC_35_71up.csv",quote=F,row.names=F)
write.csv(vHMEC_35_71up_gene,file="vHMEC_35_71up_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_35_71up_isoform,file="vHMEC_35_71up_isoform.csv",quote=F,row.names=F)
vHMEC_35_71dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
vHMEC_35_71dn$id<-unlist(strsplit(as.character(vHMEC_35_71dn$V1),'_'))[2*(1:nrow(vHMEC_35_71dn))]
vHMEC_35_71dn$ave2<-pc_lib1[vHMEC_35_71dn$id,3]
vHMEC_35_71dn$ave3<-pc_lib2[vHMEC_35_71dn$id,3]
vHMEC_35_71dn<-na.omit(vHMEC_35_71dn)
vHMEC_35_71dn<-vHMEC_35_71dn[(vHMEC_35_71dn$V2<=cutoff*vHMEC_35_71dn$ave2&vHMEC_35_71dn$V3>=cutoff2*vHMEC_35_71dn$ave3)|(vHMEC_35_71dn$V2>=cutoff2*vHMEC_35_71dn$ave2&vHMEC_35_71dn$V3<=cutoff*vHMEC_35_71dn$ave3),]
vHMEC_35_71dn_gene<-unique(vHMEC_35_71dn$id)
vHMEC_35_71dn_isoform<-vHMEC_35_71dn[vHMEC_35_71dn$ave2>Rmin&vHMEC_35_71dn$ave3>Rmin,]
write.csv(vHMEC_35_71dn,file="vHMEC_35_71dn.csv",quote=F,row.names=F)
write.csv(vHMEC_35_71dn_gene,file="vHMEC_35_71dn_gene.csv",quote=F,row.names=F)
write.csv(vHMEC_35_71dn_isoform,file="vHMEC_35_71dn_isoform.csv",quote=F,row.names=F)
vHMEC_35_71_isoform<-unique(c(vHMEC_35_71up_isoform$id,vHMEC_35_71dn_isoform$id))
write.csv(vHMEC_35_71_isoform,file="vHMEC_35_71_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_vHMEC_35_71up<-read.table("UP.vHMECRM035_vHMECRM071.FDR_0.015.rmin_0.005.Nmin_25")
DE_vHMEC_35_71dn<-read.table("DN.vHMECRM035_vHMECRM071.FDR_0.015.rmin_0.005.Nmin_25")
DE_vHMEC_35_71<-c(as.character(DE_vHMEC_35_71up$V1),as.character(DE_vHMEC_35_71dn$V1))
setwd("~/REMC/breast/individual/")
vHMEC_35_71_isoform_only<-vHMEC_35_71_isoform[!is.element(vHMEC_35_71_isoform,DE_vHMEC_35_71)]
write.csv(vHMEC_35_71_isoform_only,file="vHMEC_35_71_isoform_only.csv",quote=F,row.names=F)

lib1='HS1187'; cell1='lum'; donor1='RM035';
lib2='A01029'; cell2='lum'; donor2='RM080';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
lum_35_80up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
lum_35_80up$id<-unlist(strsplit(as.character(lum_35_80up$V1),'_'))[2*(1:nrow(lum_35_80up))]
lum_35_80up$ave2<-pc_lib1[lum_35_80up$id,3]
lum_35_80up$ave3<-pc_lib2[lum_35_80up$id,3]
lum_35_80up<-na.omit(lum_35_80up)
lum_35_80up<-lum_35_80up[(lum_35_80up$V2<=cutoff*lum_35_80up$ave2&lum_35_80up$V3>=cutoff2*lum_35_80up$ave3)|(lum_35_80up$V2>=cutoff2*lum_35_80up$ave2&lum_35_80up$V3<=cutoff*lum_35_80up$ave3),]
lum_35_80up_gene<-unique(lum_35_80up$id)
lum_35_80up_isoform<-lum_35_80up[lum_35_80up$ave2>Rmin&lum_35_80up$ave3>Rmin,]
write.csv(lum_35_80up,file="lum_35_80up.csv",quote=F,row.names=F)
write.csv(lum_35_80up_gene,file="lum_35_80up_gene.csv",quote=F,row.names=F)
write.csv(lum_35_80up_isoform,file="lum_35_80up_isoform.csv",quote=F,row.names=F)
lum_35_80dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
lum_35_80dn$id<-unlist(strsplit(as.character(lum_35_80dn$V1),'_'))[2*(1:nrow(lum_35_80dn))]
lum_35_80dn$ave2<-pc_lib1[lum_35_80dn$id,3]
lum_35_80dn$ave3<-pc_lib2[lum_35_80dn$id,3]
lum_35_80dn<-na.omit(lum_35_80dn)
lum_35_80dn<-lum_35_80dn[(lum_35_80dn$V2<=cutoff*lum_35_80dn$ave2&lum_35_80dn$V3>=cutoff2*lum_35_80dn$ave3)|(lum_35_80dn$V2>=cutoff2*lum_35_80dn$ave2&lum_35_80dn$V3<=cutoff*lum_35_80dn$ave3),]
lum_35_80dn_gene<-unique(lum_35_80dn$id)
lum_35_80dn_isoform<-lum_35_80dn[lum_35_80dn$ave2>Rmin&lum_35_80dn$ave3>Rmin,]
write.csv(lum_35_80dn,file="lum_35_80dn.csv",quote=F,row.names=F)
write.csv(lum_35_80dn_gene,file="lum_35_80dn_gene.csv",quote=F,row.names=F)
write.csv(lum_35_80dn_isoform,file="lum_35_80dn_isoform.csv",quote=F,row.names=F)
lum_35_80_isoform<-unique(c(lum_35_80up_isoform$id,lum_35_80dn_isoform$id))
write.csv(lum_35_80_isoform,file="lum_35_80_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_lum_35_80up<-read.table("UP.lumRM080_lumRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_lum_35_80dn<-read.table("DN.lumRM080_lumRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_lum_35_80<-c(as.character(DE_lum_35_80up$V1),as.character(DE_lum_35_80dn$V1))
setwd("~/REMC/breast/individual/")
lum_35_80_isoform_only<-lum_35_80_isoform[!is.element(lum_35_80_isoform,DE_lum_35_80)]
write.csv(lum_35_80_isoform_only,file="lum_35_80_isoform_only.csv",quote=F,row.names=F)

lib1='HS1188'; cell1='myo'; donor1='RM035';
lib2='A01030'; cell2='myo'; donor2='RM080';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
myo_35_80up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_35_80up$id<-unlist(strsplit(as.character(myo_35_80up$V1),'_'))[2*(1:nrow(myo_35_80up))]
myo_35_80up$ave2<-pc_lib1[myo_35_80up$id,3]
myo_35_80up$ave3<-pc_lib2[myo_35_80up$id,3]
myo_35_80up<-na.omit(myo_35_80up)
myo_35_80up<-myo_35_80up[(myo_35_80up$V2<=cutoff*myo_35_80up$ave2&myo_35_80up$V3>=cutoff2*myo_35_80up$ave3)|(myo_35_80up$V2>=cutoff2*myo_35_80up$ave2&myo_35_80up$V3<=cutoff*myo_35_80up$ave3),]
myo_35_80up_gene<-unique(myo_35_80up$id)
myo_35_80up_isoform<-myo_35_80up[myo_35_80up$ave2>Rmin&myo_35_80up$ave3>Rmin,]
write.csv(myo_35_80up,file="myo_35_80up.csv",quote=F,row.names=F)
write.csv(myo_35_80up_gene,file="myo_35_80up_gene.csv",quote=F,row.names=F)
write.csv(myo_35_80up_isoform,file="myo_35_80up_isoform.csv",quote=F,row.names=F)
myo_35_80dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_35_80dn$id<-unlist(strsplit(as.character(myo_35_80dn$V1),'_'))[2*(1:nrow(myo_35_80dn))]
myo_35_80dn$ave2<-pc_lib1[myo_35_80dn$id,3]
myo_35_80dn$ave3<-pc_lib2[myo_35_80dn$id,3]
myo_35_80dn<-na.omit(myo_35_80dn)
myo_35_80dn<-myo_35_80dn[(myo_35_80dn$V2<=cutoff*myo_35_80dn$ave2&myo_35_80dn$V3>=cutoff2*myo_35_80dn$ave3)|(myo_35_80dn$V2>=cutoff2*myo_35_80dn$ave2&myo_35_80dn$V3<=cutoff*myo_35_80dn$ave3),]
myo_35_80dn_gene<-unique(myo_35_80dn$id)
myo_35_80dn_isoform<-myo_35_80dn[myo_35_80dn$ave2>Rmin&myo_35_80dn$ave3>Rmin,]
write.csv(myo_35_80dn,file="myo_35_80dn.csv",quote=F,row.names=F)
write.csv(myo_35_80dn_gene,file="myo_35_80dn_gene.csv",quote=F,row.names=F)
write.csv(myo_35_80dn_isoform,file="myo_35_80dn_isoform.csv",quote=F,row.names=F)
myo_35_80_isoform<-unique(c(myo_35_80up_isoform$id,myo_35_80dn_isoform$id))
write.csv(myo_35_80_isoform,file="myo_35_80_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_myo_35_80up<-read.table("UP.myoRM080_myoRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_35_80dn<-read.table("DN.myoRM080_myoRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_35_80<-c(as.character(DE_myo_35_80up$V1),as.character(DE_myo_35_80dn$V1))
setwd("~/REMC/breast/individual/")
myo_35_80_isoform_only<-myo_35_80_isoform[!is.element(myo_35_80_isoform,DE_myo_35_80)]
write.csv(myo_35_80_isoform_only,file="myo_35_80_isoform_only.csv",quote=F,row.names=F)

lib1='HS1187'; cell1='lum'; donor1='RM035';
lib2='A17918'; cell2='lum'; donor2='RM084';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
lum_35_84up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
lum_35_84up$id<-unlist(strsplit(as.character(lum_35_84up$V1),'_'))[2*(1:nrow(lum_35_84up))]
lum_35_84up$ave2<-pc_lib1[lum_35_84up$id,3]
lum_35_84up$ave3<-pc_lib2[lum_35_84up$id,3]
lum_35_84up<-na.omit(lum_35_84up)
lum_35_84up<-lum_35_84up[(lum_35_84up$V2<=cutoff*lum_35_84up$ave2&lum_35_84up$V3>=cutoff2*lum_35_84up$ave3)|(lum_35_84up$V2>=cutoff2*lum_35_84up$ave2&lum_35_84up$V3<=cutoff*lum_35_84up$ave3),]
lum_35_84up_gene<-unique(lum_35_84up$id)
lum_35_84up_isoform<-lum_35_84up[lum_35_84up$ave2>Rmin&lum_35_84up$ave3>Rmin,]
write.csv(lum_35_84up,file="lum_35_84up.csv",quote=F,row.names=F)
write.csv(lum_35_84up_gene,file="lum_35_84up_gene.csv",quote=F,row.names=F)
write.csv(lum_35_84up_isoform,file="lum_35_84up_isoform.csv",quote=F,row.names=F)
lum_35_84dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
lum_35_84dn$id<-unlist(strsplit(as.character(lum_35_84dn$V1),'_'))[2*(1:nrow(lum_35_84dn))]
lum_35_84dn$ave2<-pc_lib1[lum_35_84dn$id,3]
lum_35_84dn$ave3<-pc_lib2[lum_35_84dn$id,3]
lum_35_84dn<-na.omit(lum_35_84dn)
lum_35_84dn<-lum_35_84dn[(lum_35_84dn$V2<=cutoff*lum_35_84dn$ave2&lum_35_84dn$V3>=cutoff2*lum_35_84dn$ave3)|(lum_35_84dn$V2>=cutoff2*lum_35_84dn$ave2&lum_35_84dn$V3<=cutoff*lum_35_84dn$ave3),]
lum_35_84dn_gene<-unique(lum_35_84dn$id)
lum_35_84dn_isoform<-lum_35_84dn[lum_35_84dn$ave2>Rmin&lum_35_84dn$ave3>Rmin,]
write.csv(lum_35_84dn,file="lum_35_84dn.csv",quote=F,row.names=F)
write.csv(lum_35_84dn_gene,file="lum_35_84dn_gene.csv",quote=F,row.names=F)
write.csv(lum_35_84dn_isoform,file="lum_35_84dn_isoform.csv",quote=F,row.names=F)
lum_35_84_isoform<-unique(c(lum_35_84up_isoform$id,lum_35_84dn_isoform$id))
write.csv(lum_35_84_isoform,file="lum_35_84_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_lum_35_84up<-read.table("UP.lumRM084_lumRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_lum_35_84dn<-read.table("DN.lumRM084_lumRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_lum_35_84<-c(as.character(DE_lum_35_84up$V1),as.character(DE_lum_35_84dn$V1))
setwd("~/REMC/breast/individual/")
lum_35_84_isoform_only<-lum_35_84_isoform[!is.element(lum_35_84_isoform,DE_lum_35_84)]
write.csv(lum_35_84_isoform_only,file="lum_35_84_isoform_only.csv",quote=F,row.names=F)

lib1='HS1188'; cell1='myo'; donor1='RM035';
lib2='A17919'; cell2='myo'; donor2='RM084';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
myo_35_84up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_35_84up$id<-unlist(strsplit(as.character(myo_35_84up$V1),'_'))[2*(1:nrow(myo_35_84up))]
myo_35_84up$ave2<-pc_lib1[myo_35_84up$id,3]
myo_35_84up$ave3<-pc_lib2[myo_35_84up$id,3]
myo_35_84up<-na.omit(myo_35_84up)
myo_35_84up<-myo_35_84up[(myo_35_84up$V2<=cutoff*myo_35_84up$ave2&myo_35_84up$V3>=cutoff2*myo_35_84up$ave3)|(myo_35_84up$V2>=cutoff2*myo_35_84up$ave2&myo_35_84up$V3<=cutoff*myo_35_84up$ave3),]
myo_35_84up_gene<-unique(myo_35_84up$id)
myo_35_84up_isoform<-myo_35_84up[myo_35_84up$ave2>Rmin&myo_35_84up$ave3>Rmin,]
write.csv(myo_35_84up,file="myo_35_84up.csv",quote=F,row.names=F)
write.csv(myo_35_84up_gene,file="myo_35_84up_gene.csv",quote=F,row.names=F)
write.csv(myo_35_84up_isoform,file="myo_35_84up_isoform.csv",quote=F,row.names=F)
myo_35_84dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_35_84dn$id<-unlist(strsplit(as.character(myo_35_84dn$V1),'_'))[2*(1:nrow(myo_35_84dn))]
myo_35_84dn$ave2<-pc_lib1[myo_35_84dn$id,3]
myo_35_84dn$ave3<-pc_lib2[myo_35_84dn$id,3]
myo_35_84dn<-na.omit(myo_35_84dn)
myo_35_84dn<-myo_35_84dn[(myo_35_84dn$V2<=cutoff*myo_35_84dn$ave2&myo_35_84dn$V3>=cutoff2*myo_35_84dn$ave3)|(myo_35_84dn$V2>=cutoff2*myo_35_84dn$ave2&myo_35_84dn$V3<=cutoff*myo_35_84dn$ave3),]
myo_35_84dn_gene<-unique(myo_35_84dn$id)
myo_35_84dn_isoform<-myo_35_84dn[myo_35_84dn$ave2>Rmin&myo_35_84dn$ave3>Rmin,]
write.csv(myo_35_84dn,file="myo_35_84dn.csv",quote=F,row.names=F)
write.csv(myo_35_84dn_gene,file="myo_35_84dn_gene.csv",quote=F,row.names=F)
write.csv(myo_35_84dn_isoform,file="myo_35_84dn_isoform.csv",quote=F,row.names=F)
myo_35_84_isoform<-unique(c(myo_35_84up_isoform$id,myo_35_84dn_isoform$id))
write.csv(myo_35_84_isoform,file="myo_35_84_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_myo_35_84up<-read.table("UP.myoRM084_myoRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_35_84dn<-read.table("DN.myoRM084_myoRM035.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_35_84<-c(as.character(DE_myo_35_84up$V1),as.character(DE_myo_35_84dn$V1))
setwd("~/REMC/breast/individual/")
myo_35_84_isoform_only<-myo_35_84_isoform[!is.element(myo_35_84_isoform,DE_myo_35_84)]
write.csv(myo_35_84_isoform_only,file="myo_35_84_isoform_only.csv",quote=F,row.names=F)

lib1='A18760'; cell1='fibro'; donor1='RM070';
lib2='A18472'; cell2='fibro'; donor2='RM071';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
fibro_70_71up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
fibro_70_71up$id<-unlist(strsplit(as.character(fibro_70_71up$V1),'_'))[2*(1:nrow(fibro_70_71up))]
fibro_70_71up$ave2<-pc_lib1[fibro_70_71up$id,3]
fibro_70_71up$ave3<-pc_lib2[fibro_70_71up$id,3]
fibro_70_71up<-na.omit(fibro_70_71up)
fibro_70_71up<-fibro_70_71up[(fibro_70_71up$V2<=cutoff*fibro_70_71up$ave2&fibro_70_71up$V3>=cutoff2*fibro_70_71up$ave3)|(fibro_70_71up$V2>=cutoff2*fibro_70_71up$ave2&fibro_70_71up$V3<=cutoff*fibro_70_71up$ave3),]
fibro_70_71up_gene<-unique(fibro_70_71up$id)
fibro_70_71up_isoform<-fibro_70_71up[fibro_70_71up$ave2>Rmin&fibro_70_71up$ave3>Rmin,]
write.csv(fibro_70_71up,file="fibro_70_71up.csv",quote=F,row.names=F)
write.csv(fibro_70_71up_gene,file="fibro_70_71up_gene.csv",quote=F,row.names=F)
write.csv(fibro_70_71up_isoform,file="fibro_70_71up_isoform.csv",quote=F,row.names=F)
fibro_70_71dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
fibro_70_71dn$id<-unlist(strsplit(as.character(fibro_70_71dn$V1),'_'))[2*(1:nrow(fibro_70_71dn))]
fibro_70_71dn$ave2<-pc_lib1[fibro_70_71dn$id,3]
fibro_70_71dn$ave3<-pc_lib2[fibro_70_71dn$id,3]
fibro_70_71dn<-na.omit(fibro_70_71dn)
fibro_70_71dn<-fibro_70_71dn[(fibro_70_71dn$V2<=cutoff*fibro_70_71dn$ave2&fibro_70_71dn$V3>=cutoff2*fibro_70_71dn$ave3)|(fibro_70_71dn$V2>=cutoff2*fibro_70_71dn$ave2&fibro_70_71dn$V3<=cutoff*fibro_70_71dn$ave3),]
fibro_70_71dn_gene<-unique(fibro_70_71dn$id)
fibro_70_71dn_isoform<-fibro_70_71dn[fibro_70_71dn$ave2>Rmin&fibro_70_71dn$ave3>Rmin,]
write.csv(fibro_70_71dn,file="fibro_70_71dn.csv",quote=F,row.names=F)
write.csv(fibro_70_71dn_gene,file="fibro_70_71dn_gene.csv",quote=F,row.names=F)
write.csv(fibro_70_71dn_isoform,file="fibro_70_71dn_isoform.csv",quote=F,row.names=F)
fibro_70_71_isoform<-unique(c(fibro_70_71up_isoform$id,fibro_70_71dn_isoform$id))
write.csv(fibro_70_71_isoform,file="fibro_70_71_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_fibro_70_71up<-read.table("UP.fibrRM070_fibrRM071.FDR_0.015.rmin_0.005.Nmin_25")
DE_fibro_70_71dn<-read.table("DN.fibrRM070_fibrRM071.FDR_0.015.rmin_0.005.Nmin_25")
DE_fibro_70_71<-c(as.character(DE_fibro_70_71up$V1),as.character(DE_fibro_70_71dn$V1))
setwd("~/REMC/breast/individual/")
fibro_70_71_isoform_only<-fibro_70_71_isoform[!is.element(fibro_70_71_isoform,DE_fibro_70_71)]
write.csv(fibro_70_71_isoform_only,file="fibro_70_71_isoform_only.csv",quote=F,row.names=F)

lib1='A01029'; cell1='lum'; donor1='RM080';
lib2='A17918'; cell2='lum'; donor2='RM084';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
lum_80_84up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
lum_80_84up$id<-unlist(strsplit(as.character(lum_80_84up$V1),'_'))[2*(1:nrow(lum_80_84up))]
lum_80_84up$ave2<-pc_lib1[lum_80_84up$id,3]
lum_80_84up$ave3<-pc_lib2[lum_80_84up$id,3]
lum_80_84up<-na.omit(lum_80_84up)
lum_80_84up<-lum_80_84up[(lum_80_84up$V2<=cutoff*lum_80_84up$ave2&lum_80_84up$V3>=cutoff2*lum_80_84up$ave3)|(lum_80_84up$V2>=cutoff2*lum_80_84up$ave2&lum_80_84up$V3<=cutoff*lum_80_84up$ave3),]
lum_80_84up_gene<-unique(lum_80_84up$id)
lum_80_84up_isoform<-lum_80_84up[lum_80_84up$ave2>Rmin&lum_80_84up$ave3>Rmin,]
write.csv(lum_80_84up,file="lum_80_84up.csv",quote=F,row.names=F)
write.csv(lum_80_84up_gene,file="lum_80_84up_gene.csv",quote=F,row.names=F)
write.csv(lum_80_84up_isoform,file="lum_80_84up_isoform.csv",quote=F,row.names=F)
lum_80_84dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
lum_80_84dn$id<-unlist(strsplit(as.character(lum_80_84dn$V1),'_'))[2*(1:nrow(lum_80_84dn))]
lum_80_84dn$ave2<-pc_lib1[lum_80_84dn$id,3]
lum_80_84dn$ave3<-pc_lib2[lum_80_84dn$id,3]
lum_80_84dn<-na.omit(lum_80_84dn)
lum_80_84dn<-lum_80_84dn[(lum_80_84dn$V2<=cutoff*lum_80_84dn$ave2&lum_80_84dn$V3>=cutoff2*lum_80_84dn$ave3)|(lum_80_84dn$V2>=cutoff2*lum_80_84dn$ave2&lum_80_84dn$V3<=cutoff*lum_80_84dn$ave3),]
lum_80_84dn_gene<-unique(lum_80_84dn$id)
lum_80_84dn_isoform<-lum_80_84dn[lum_80_84dn$ave2>Rmin&lum_80_84dn$ave3>Rmin,]
write.csv(lum_80_84dn,file="lum_80_84dn.csv",quote=F,row.names=F)
write.csv(lum_80_84dn_gene,file="lum_80_84dn_gene.csv",quote=F,row.names=F)
write.csv(lum_80_84dn_isoform,file="lum_80_84dn_isoform.csv",quote=F,row.names=F)
lum_80_84_isoform<-unique(c(lum_80_84up_isoform$id,lum_80_84dn_isoform$id))
write.csv(lum_80_84_isoform,file="lum_80_84_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_lum_80_84up<-read.table("UP.lumRM080_lumRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_lum_80_84dn<-read.table("DN.lumRM080_lumRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_lum_80_84<-c(as.character(DE_lum_80_84up$V1),as.character(DE_lum_80_84dn$V1))
setwd("~/REMC/breast/individual/")
lum_80_84_isoform_only<-lum_80_84_isoform[!is.element(lum_80_84_isoform,DE_lum_80_84)]
write.csv(lum_80_84_isoform_only,file="lum_80_84_isoform_only.csv",quote=F,row.names=F)

lib1='A01030'; cell1='myo'; donor1='RM080';
lib2='A17919'; cell2='myo'; donor2='RM084';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
myo_80_84up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_80_84up$id<-unlist(strsplit(as.character(myo_80_84up$V1),'_'))[2*(1:nrow(myo_80_84up))]
myo_80_84up$ave2<-pc_lib1[myo_80_84up$id,3]
myo_80_84up$ave3<-pc_lib2[myo_80_84up$id,3]
myo_80_84up<-na.omit(myo_80_84up)
myo_80_84up<-myo_80_84up[(myo_80_84up$V2<=cutoff*myo_80_84up$ave2&myo_80_84up$V3>=cutoff2*myo_80_84up$ave3)|(myo_80_84up$V2>=cutoff2*myo_80_84up$ave2&myo_80_84up$V3<=cutoff*myo_80_84up$ave3),]
myo_80_84up_gene<-unique(myo_80_84up$id)
myo_80_84up_isoform<-myo_80_84up[myo_80_84up$ave2>Rmin&myo_80_84up$ave3>Rmin,]
write.csv(myo_80_84up,file="myo_80_84up.csv",quote=F,row.names=F)
write.csv(myo_80_84up_gene,file="myo_80_84up_gene.csv",quote=F,row.names=F)
write.csv(myo_80_84up_isoform,file="myo_80_84up_isoform.csv",quote=F,row.names=F)
myo_80_84dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
myo_80_84dn$id<-unlist(strsplit(as.character(myo_80_84dn$V1),'_'))[2*(1:nrow(myo_80_84dn))]
myo_80_84dn$ave2<-pc_lib1[myo_80_84dn$id,3]
myo_80_84dn$ave3<-pc_lib2[myo_80_84dn$id,3]
myo_80_84dn<-na.omit(myo_80_84dn)
myo_80_84dn<-myo_80_84dn[(myo_80_84dn$V2<=cutoff*myo_80_84dn$ave2&myo_80_84dn$V3>=cutoff2*myo_80_84dn$ave3)|(myo_80_84dn$V2>=cutoff2*myo_80_84dn$ave2&myo_80_84dn$V3<=cutoff*myo_80_84dn$ave3),]
myo_80_84dn_gene<-unique(myo_80_84dn$id)
myo_80_84dn_isoform<-myo_80_84dn[myo_80_84dn$ave2>Rmin&myo_80_84dn$ave3>Rmin,]
write.csv(myo_80_84dn,file="myo_80_84dn.csv",quote=F,row.names=F)
write.csv(myo_80_84dn_gene,file="myo_80_84dn_gene.csv",quote=F,row.names=F)
write.csv(myo_80_84dn_isoform,file="myo_80_84dn_isoform.csv",quote=F,row.names=F)
myo_80_84_isoform<-unique(c(myo_80_84up_isoform$id,myo_80_84dn_isoform$id))
write.csv(myo_80_84_isoform,file="myo_80_84_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_myo_80_84up<-read.table("UP.myoRM080_myoRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_80_84dn<-read.table("DN.myoRM080_myoRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_myo_80_84<-c(as.character(DE_myo_80_84up$V1),as.character(DE_myo_80_84dn$V1))
setwd("~/REMC/breast/individual/")
myo_80_84_isoform_only<-myo_80_84_isoform[!is.element(myo_80_84_isoform,DE_myo_80_84)]
write.csv(myo_80_84_isoform_only,file="myo_80_84_isoform_only.csv",quote=F,row.names=F)

lib1='A01031'; cell1='stem'; donor1='RM080';
lib2='A17920'; cell2='stem'; donor2='RM084';
pc_lib1<-read.table(paste(dirIn,lib1,"/coverage/",lib1,".G.A.rpkm.pc",sep=""),row.names=1)
pc_lib2<-read.table(paste(dirIn,lib2,"/coverage/",lib2,".G.A.rpkm.pc",sep=""),row.names=1)
stem_80_84up<-read.table(paste("UP.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
stem_80_84up$id<-unlist(strsplit(as.character(stem_80_84up$V1),'_'))[2*(1:nrow(stem_80_84up))]
stem_80_84up$ave2<-pc_lib1[stem_80_84up$id,3]
stem_80_84up$ave3<-pc_lib2[stem_80_84up$id,3]
stem_80_84up<-na.omit(stem_80_84up)
stem_80_84up<-stem_80_84up[(stem_80_84up$V2<=cutoff*stem_80_84up$ave2&stem_80_84up$V3>=cutoff2*stem_80_84up$ave3)|(stem_80_84up$V2>=cutoff2*stem_80_84up$ave2&stem_80_84up$V3<=cutoff*stem_80_84up$ave3),]
stem_80_84up_gene<-unique(stem_80_84up$id)
stem_80_84up_isoform<-stem_80_84up[stem_80_84up$ave2>Rmin&stem_80_84up$ave3>Rmin,]
write.csv(stem_80_84up,file="stem_80_84up.csv",quote=F,row.names=F)
write.csv(stem_80_84up_gene,file="stem_80_84up_gene.csv",quote=F,row.names=F)
write.csv(stem_80_84up_isoform,file="stem_80_84up_isoform.csv",quote=F,row.names=F)
stem_80_84dn<-read.table(paste("DN.",cell1,"-",donor1,"_",cell2,"-",donor2,".FDR_0.015.rmin_0.005.Nmin_25",sep=""))
stem_80_84dn$id<-unlist(strsplit(as.character(stem_80_84dn$V1),'_'))[2*(1:nrow(stem_80_84dn))]
stem_80_84dn$ave2<-pc_lib1[stem_80_84dn$id,3]
stem_80_84dn$ave3<-pc_lib2[stem_80_84dn$id,3]
stem_80_84dn<-na.omit(stem_80_84dn)
stem_80_84dn<-stem_80_84dn[(stem_80_84dn$V2<=cutoff*stem_80_84dn$ave2&stem_80_84dn$V3>=cutoff2*stem_80_84dn$ave3)|(stem_80_84dn$V2>=cutoff2*stem_80_84dn$ave2&stem_80_84dn$V3<=cutoff*stem_80_84dn$ave3),]
stem_80_84dn_gene<-unique(stem_80_84dn$id)
stem_80_84dn_isoform<-stem_80_84dn[stem_80_84dn$ave2>Rmin&stem_80_84dn$ave3>Rmin,]
write.csv(stem_80_84dn,file="stem_80_84dn.csv",quote=F,row.names=F)
write.csv(stem_80_84dn_gene,file="stem_80_84dn_gene.csv",quote=F,row.names=F)
write.csv(stem_80_84dn_isoform,file="stem_80_84dn_isoform.csv",quote=F,row.names=F)
stem_80_84_isoform<-unique(c(stem_80_84up_isoform$id,stem_80_84dn_isoform$id))
write.csv(stem_80_84_isoform,file="stem_80_84_isoform.csv",quote=F,row.names=F)
# exclude gene level DE
setwd("~/REMC/breast/gene/")
DE_stem_80_84up<-read.table("UP.stemRM080_stemRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_stem_80_84dn<-read.table("DN.stemRM080_stemRM084.FDR_0.015.rmin_0.005.Nmin_25")
DE_stem_80_84<-c(as.character(DE_stem_80_84up$V1),as.character(DE_stem_80_84dn$V1))
setwd("~/REMC/breast/individual/")
stem_80_84_isoform_only<-stem_80_84_isoform[!is.element(stem_80_84_isoform,DE_stem_80_84)]
write.csv(stem_80_84_isoform_only,file="stem_80_84_isoform_only.csv",quote=F,row.names=F)

rm(lib1,cell1,donor1,lib2,cell2,donor2,pc_lib1,pc_lib2)
save.image(file="REMC.isoform.breast.RData")

##################################################################################

setwd("~/REMC/breast/individual/")
load("REMC.isoform.breast.RData")

# hiearchal clustering
isoform<-union(union(union(union(vHMEC_35_71_isoform_only,lum_35_80_isoform_only),
                           union(myo_35_80_isoform_only,lum_35_84_isoform_only)),
                     union(union(myo_35_84_isoform_only,fibro_70_71_isoform_only),
                           union(lum_80_84_isoform_only,myo_80_84_isoform_only))),
               stem_80_84_isoform_only)
isoform<-data.frame(ID=isoform,vHMEC_35_71=0,lum_35_80=0,
                    myo_35_80=0,lum_35_84=0,
                    myo_35_84=0,fibro_70_71=0,
                    lum_80_84=0,myo_80_84=0,
                    stem_80_84=0)
rownames(isoform)<-isoform$ID
isoform[vHMEC_35_71_isoform_only,2]<-1
isoform[lum_35_80_isoform_only,3]<-1
isoform[myo_35_80_isoform_only,4]<-1
isoform[lum_35_84_isoform_only,5]<-1
isoform[myo_35_84_isoform_only,6]<-1
isoform[fibro_70_71_isoform_only,7]<-1
isoform[lum_80_84_isoform_only,8]<-1
isoform[myo_80_84_isoform_only,9]<-1
isoform[stem_80_84_isoform_only,10]<-1
c<-cor(as.matrix(isoform[,2:10]), method="spearman")
d<-as.dist(1-c)
pdf("dendrogram.pdf")
plot(hclust(d,method="ward"),main="Isoform clustering",xlab=NA)
dev.off()

# individual summary
RM035_71up<-as.character(vHMEC_35_71up$V1)
RM035_71up_gene<-vHMEC_35_71up_gene
RM035_71dn<-as.character(vHMEC_35_71dn$V1)
RM035_71dn_gene<-vHMEC_35_71dn_gene
RM035_71_isoform<-vHMEC_35_71_isoform
RM035_71_isoform_only<-vHMEC_35_71_isoform_only
write.csv(RM035_71_isoform_only,file="RM035_71_isoform_only.csv",quote=F,row.names=F)

RM035_80up<-intersect(as.character(lum_35_80up$V1),as.character(myo_35_80up$V1))
RM035_80up_gene<-intersect(lum_35_80up_gene,myo_35_80up_gene)
RM035_80dn<-intersect(as.character(lum_35_80dn$V1),as.character(myo_35_80dn$V1))
RM035_80dn_gene<-intersect(lum_35_80dn_gene,myo_35_80dn_gene)
RM035_80_isoform<-intersect(lum_35_80_isoform,myo_35_80_isoform)
RM035_80_isoform_only<-intersect(lum_35_80_isoform_only,myo_35_80_isoform_only)
write.csv(RM035_80_isoform_only,file="RM035_80_isoform_only.csv",quote=F,row.names=F)

RM035_84up<-intersect(as.character(lum_35_84up$V1),as.character(myo_35_84up$V1))
RM035_84up_gene<-intersect(lum_35_84up_gene,myo_35_84up_gene)
RM035_84dn<-intersect(as.character(lum_35_84dn$V1),as.character(myo_35_84dn$V1))
RM035_84dn_gene<-intersect(lum_35_84dn_gene,myo_35_84dn_gene)
RM035_84_isoform<-intersect(lum_35_84_isoform,myo_35_84_isoform)
RM035_84_isoform_only<-intersect(lum_35_84_isoform_only,myo_35_84_isoform_only)
write.csv(RM035_84_isoform_only,file="RM035_84_isoform_only.csv",quote=F,row.names=F)

RM070_71up<-as.character(fibro_70_71up$V1)
RM070_71up_gene<-fibro_70_71up_gene
RM070_71dn<-as.character(fibro_70_71dn$V1)
RM070_71dn_gene<-fibro_70_71dn_gene
RM070_71_isoform<-fibro_70_71_isoform
RM070_71_isoform_only<-fibro_70_71_isoform_only
write.csv(RM070_71_isoform_only,file="RM070_71_isoform_only.csv",quote=F,row.names=F)

RM080_84up<-intersect(intersect(as.character(lum_80_84up$V1),as.character(myo_80_84up$V1)),as.character(stem_80_84up$V1))
RM080_84up_gene<-intersect(intersect(lum_80_84up_gene,myo_80_84up_gene),stem_80_84up_gene)
RM080_84dn<-intersect(intersect(as.character(lum_80_84dn$V1),as.character(myo_80_84dn$V1)),as.character(stem_80_84dn$V1))
RM080_84dn_gene<-intersect(intersect(lum_80_84dn_gene,myo_80_84dn_gene),stem_80_84dn_gene)
RM080_84_isoform<-intersect(intersect(lum_80_84_isoform,myo_80_84_isoform),stem_80_84_isoform)
RM080_84_isoform_only<-intersect(intersect(lum_80_84_isoform_only,myo_80_84_isoform_only),stem_80_84_isoform_only)
write.csv(RM080_84_isoform_only,file="RM080_84_isoform_only.csv",quote=F,row.names=F)

write.csv(RM035_71up,file="RM035_71up.csv",quote=F,row.names=F)
write.csv(RM035_71up_gene,file="RM035_71up_gene.csv",quote=F,row.names=F)
write.csv(RM035_71dn,file="RM035_71dn.csv",quote=F,row.names=F)
write.csv(RM035_71dn_gene,file="RM035_71dn_gene.csv",quote=F,row.names=F)
write.csv(RM035_80up,file="RM035_80up.csv",quote=F,row.names=F)
write.csv(RM035_80up_gene,file="RM035_80up_gene.csv",quote=F,row.names=F)
write.csv(RM035_80dn,file="RM035_80dn.csv",quote=F,row.names=F)
write.csv(RM035_80dn_gene,file="RM035_80dn_gene.csv",quote=F,row.names=F)
write.csv(RM035_84up,file="RM035_84up.csv",quote=F,row.names=F)
write.csv(RM035_84up_gene,file="RM035_84up_gene.csv",quote=F,row.names=F)
write.csv(RM035_84dn,file="RM035_84dn.csv",quote=F,row.names=F)
write.csv(RM035_84dn_gene,file="RM035_84dn_gene.csv",quote=F,row.names=F)
write.csv(RM070_71up,file="RM070_71up.csv",quote=F,row.names=F)
write.csv(RM070_71up_gene,file="RM070_71up_gene.csv",quote=F,row.names=F)
write.csv(RM070_71dn,file="RM070_71dn.csv",quote=F,row.names=F)
write.csv(RM070_71dn_gene,file="RM070_71dn_gene.csv",quote=F,row.names=F)
write.csv(RM080_84up,file="RM080_84up.csv",quote=F,row.names=F)
write.csv(RM080_84up_gene,file="RM080_84up_gene.csv",quote=F,row.names=F)
write.csv(RM080_84dn,file="RM080_84dn.csv",quote=F,row.names=F)
write.csv(RM080_84dn_gene,file="RM080_84dn_gene.csv",quote=F,row.names=F)
write.csv(RM035_71_isoform,file="RM035_71_isoform.csv",quote=F,row.names=F)
write.csv(RM035_80_isoform,file="RM035_80_isoform.csv",quote=F,row.names=F)
write.csv(RM035_84_isoform,file="RM035_84_isoform.csv",quote=F,row.names=F)
write.csv(RM070_71_isoform,file="RM070_71_isoform.csv",quote=F,row.names=F)
write.csv(RM080_84_isoform,file="RM080_84_isoform.csv",quote=F,row.names=F)

# RM080 vs (RM084+RM035)
RM080vs84lum_RM080vs35lum<-intersect(lum_80_84_isoform_only,lum_35_80_isoform_only)
RM080vs84myo_RM080vs35myo<-intersect(myo_80_84_isoform_only,myo_35_80_isoform_only)
RM080vs84_35<-intersect(intersect(RM080vs84lum_RM080vs35lum,RM080vs84myo_RM080vs35myo),stem_80_84_isoform_only)
write.csv(RM080vs84lum_RM080vs35lum,file="RM080vs84lum_RM080vs35lum.csv",quote=F,row.names=F)
write.csv(RM080vs84myo_RM080vs35myo,file="RM080vs84myo_RM080vs35myo.csv",quote=F,row.names=F)
write.csv(RM080vs84_35,file="RM080vs84_35.csv",quote=F,row.names=F)

#get gene information
setwd("~/hg19/")
gene<-read.table(file="hg19v65_genes",as.is=T,head=F)
colnames(gene)<-c("Ensembl","chr","start","end","strand","type","name","description")
gene<-gene[!duplicated(gene$Ensembl),]
rownames(gene)<-gene$Ensembl
setwd("~/REMC/breast/individual/")

RM035_84_isoform_only<-data.frame(ID=RM035_84_isoform_only,Name=gene[RM035_84_isoform_only,7],Description=gene[RM035_84_isoform_only,8])
write.csv(RM035_84_isoform_only,file="RM035_84_isoform_only.csv",quote=F,row.names=F)
RM080vs84_35<-data.frame(ID=RM080vs84_35,Name=gene[RM080vs84_35,7],Description=gene[RM080vs84_35,8])
write.csv(RM080vs84_35,file="RM080vs84_35.csv",quote=F,row.names=F)

save.image(file="REMC.isoform.breast.RData")
