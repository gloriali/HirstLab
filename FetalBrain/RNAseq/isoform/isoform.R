##################################################################################
# getting the cutoff for expressing and non-expressing
# run on xhost
# /gsc/software/linux-x86_64-centos5/R-3.0.2/bin/R CMD BATCH 
setwd("~/FetalBrain/RNAseq/DEfine/exon/")
dirIn='/projects/epigenomics/ep50/internal/jqc.1.7.6/'
libs=c("A03484",
  "A07825",
  "A03473",
  "A03475",
  "A04599",
  "A15298",
  "A03474",
  "A03476",
  "A15295",
  "A15299")

pdf("percent_ave.pdf")
plot(range(0,0.1),range(0,0.3),type='n',main="ECDF percentage of average RPKM",xlab="percentage",ylab="cdf")
i=1
for(lib in libs){
  exon_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.exn.A.rpkm",sep=""))
  pc_lib<-read.table(paste(dirIn,lib,"/coverage/",lib,".G.A.rpkm.pc",sep=""),row.names=1)
  exon_lib$ave<-pc_lib[exon_lib$V2,3]
  exon_lib<-exon_lib[exon_lib$ave!=0,]
  exon_lib$percent<-exon_lib$V4/exon_lib$ave
  lines(ecdf(exon_lib$percent[exon_lib$percent!=0]),col=i,lty=1)
  i=i+1
}
legend("bottomright",libs,col=c(1:length(libs)),cex=0.8,lty=1)
dev.off()
##################################################################################
##################################################################################
##isoforms:DEfine on exons & RPKM(exon) < 1/10*RPKM(ave. exons of the gene)
setwd("~/FetalBrain/RNAseq/rpkm/")
brain1<-read.table("A03484.G.A.rpkm.pc",row.names=1)
brain2<-read.table("A07825.G.A.rpkm.pc",row.names=1)
cortex1<-read.table("A03473.G.A.rpkm.pc",row.names=1)
cortex2<-read.table("A03475.G.A.rpkm.pc",row.names=1)
cortex3<-read.table("A04599.G.A.rpkm.pc",row.names=1)
cortex4<-read.table("A15298.G.A.rpkm.pc",row.names=1)
ge1<-read.table("A03474.G.A.rpkm.pc",row.names=1)
ge2<-read.table("A03476.G.A.rpkm.pc",row.names=1)
ge3<-read.table("A15295.G.A.rpkm.pc",row.names=1)
ge4<-read.table("A15299.G.A.rpkm.pc",row.names=1)
# brain1<-read.table("A03484.G.A.rpkm.nc",row.names=1)
# brain2<-read.table("A07825.G.A.rpkm.nc",row.names=1)
# cortex1<-read.table("A03473.G.A.rpkm.nc",row.names=1)
# cortex2<-read.table("A03475.G.A.rpkm.nc",row.names=1)
# cortex3<-read.table("A04599.G.A.rpkm.nc",row.names=1)
# cortex4<-read.table("A15298.G.A.rpkm.nc",row.names=1)
# ge1<-read.table("A03474.G.A.rpkm.nc",row.names=1)
# ge2<-read.table("A03476.G.A.rpkm.nc",row.names=1)
# ge3<-read.table("A15295.G.A.rpkm.nc",row.names=1)
# ge4<-read.table("A15299.G.A.rpkm.nc",row.names=1)
# brain1<-read.table("A03484.G.exn.A.rpkm")
# brain2<-read.table("A07825.G.exn.A.rpkm")
# cortex1<-read.table("A03473.G.exn.A.rpkm")
# cortex2<-read.table("A03475.G.exn.A.rpkm")
# cortex3<-read.table("A04599.G.exn.A.rpkm")
# cortex4<-read.table("A15298.G.exn.A.rpkm")
# ge1<-read.table("A03474.G.exn.A.rpkm")
# ge2<-read.table("A03476.G.exn.A.rpkm")
# ge3<-read.table("A15295.G.exn.A.rpkm")
# ge4<-read.table("A15299.G.exn.A.rpkm")

setwd("/Users/epigenomics_lab_02/FetalBrain/RNAseq/DEfine/exon/")
col<-c("Gene id", "rpkm1", "rpkm2", "p-value", "Multiple testing corrected p-value")
cutoff=0.01
Rmin=0.005

#between tissue: cortex Vs GE
cortexge1up<-read.table("UP.Cortex-HuFNSC01_GE-HuFNSC01.FDR_0.01.rmin_0.005.Nmin_25")
cortexge1up$id<-unlist(strsplit(as.character(cortexge1up$V1),'_'))[2*(1:nrow(cortexge1up))]
cortexge1up$ave2<-cortex1[cortexge1up$id,3]
cortexge1up$ave3<-ge1[cortexge1up$id,3]
cortexge1up<-na.omit(cortexge1up)
cortexge1up<-cortexge1up[(cortexge1up$V2<=cutoff*cortexge1up$ave2&cortexge1up$V3>=cutoff*cortexge1up$ave3)|(cortexge1up$V2>=cutoff*cortexge1up$ave2&cortexge1up$V3<=cutoff*cortexge1up$ave3),]
cortexge1up_gene<-unique(cortexge1up$id)
cortexge1up_isoform<-cortexge1up[cortexge1up$ave2>Rmin&cortexge1up$ave3>Rmin,]
write.csv(cortexge1up,file="cortexge1up.csv",quote=F,row.names=F)
write.csv(cortexge1up_gene,file="cortexge1up_gene.csv",quote=F,row.names=F)
write.csv(cortexge1up_isoform,file="cortexge1up_isoform.csv",quote=F,row.names=F)
cortexge1dn<-read.table("DN.Cortex-HuFNSC01_GE-HuFNSC01.FDR_0.01.rmin_0.005.Nmin_25")
cortexge1dn$id<-unlist(strsplit(as.character(cortexge1dn$V1),'_'))[2*(1:nrow(cortexge1dn))]
cortexge1dn$ave2<-cortex1[cortexge1dn$id,3]
cortexge1dn$ave3<-ge1[cortexge1dn$id,3]
cortexge1dn<-na.omit(cortexge1dn)
cortexge1dn<-cortexge1dn[(cortexge1dn$V2<=cutoff*cortexge1dn$ave2&cortexge1dn$V3>=cutoff*cortexge1dn$ave3)|(cortexge1dn$V2>=cutoff*cortexge1dn$ave2&cortexge1dn$V3<=cutoff*cortexge1dn$ave3),]
cortexge1dn_gene<-unique(cortexge1dn$id)
cortexge1dn_isoform<-cortexge1dn[cortexge1dn$ave2>Rmin&cortexge1dn$ave3>Rmin,]
write.csv(cortexge1dn,file="cortexge1dn.csv",quote=F,row.names=F)
write.csv(cortexge1dn_gene,file="cortexge1dn_gene.csv",quote=F,row.names=F)
write.csv(cortexge1dn_isoform,file="cortexge1dn_isoform.csv",quote=F,row.names=F)
cortexge1_isoform<-unique(c(cortexge1up_isoform$id,cortexge1dn_isoform$id))
write.csv(cortexge1_isoform,file="cortexge1_isoform.csv",quote=F,row.names=F)

cortexge2up<-read.table("UP.Cortex-HuFNSC02_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
cortexge2up$id<-unlist(strsplit(as.character(cortexge2up$V1),'_'))[2*(1:nrow(cortexge2up))]
cortexge2up$ave2<-cortex1[cortexge2up$id,3]
cortexge2up$ave3<-ge1[cortexge2up$id,3]
cortexge2up<-na.omit(cortexge2up)
cortexge2up<-cortexge2up[(cortexge2up$V2<=cutoff*cortexge2up$ave2&cortexge2up$V3>=cutoff*cortexge2up$ave3)|(cortexge2up$V2>=cutoff*cortexge2up$ave2&cortexge2up$V3<=cutoff*cortexge2up$ave3),]
cortexge2up_gene<-unique(cortexge2up$id)
cortexge2up_isoform<-cortexge2up[cortexge2up$ave2>Rmin&cortexge2up$ave3>Rmin,]
write.csv(cortexge2up,file="cortexge2up.csv",quote=F,row.names=F)
write.csv(cortexge2up_gene,file="cortexge2up_gene.csv",quote=F,row.names=F)
write.csv(cortexge2up_isoform,file="cortexge2up_isoform.csv",quote=F,row.names=F)
cortexge2dn<-read.table("DN.Cortex-HuFNSC02_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
cortexge2dn$id<-unlist(strsplit(as.character(cortexge2dn$V1),'_'))[2*(1:nrow(cortexge2dn))]
cortexge2dn$ave2<-cortex1[cortexge2dn$id,3]
cortexge2dn$ave3<-ge1[cortexge2dn$id,3]
cortexge2dn<-na.omit(cortexge2dn)
cortexge2dn<-cortexge2dn[(cortexge2dn$V2<=cutoff*cortexge2dn$ave2&cortexge2dn$V3>=cutoff*cortexge2dn$ave3)|(cortexge2dn$V2>=cutoff*cortexge2dn$ave2&cortexge2dn$V3<=cutoff*cortexge2dn$ave3),]
cortexge2dn_gene<-unique(cortexge2dn$id)
cortexge2dn_isoform<-cortexge2dn[cortexge2dn$ave2>Rmin&cortexge2dn$ave3>Rmin,]
write.csv(cortexge2dn,file="cortexge2dn.csv",quote=F,row.names=F)
write.csv(cortexge2dn_gene,file="cortexge2dn_gene.csv",quote=F,row.names=F)
write.csv(cortexge2dn_isoform,file="cortexge2dn_isoform.csv",quote=F,row.names=F)
cortexge2_isoform<-unique(c(cortexge2up_isoform$id,cortexge2dn_isoform$id))
write.csv(cortexge2_isoform,file="cortexge2_isoform.csv",quote=F,row.names=F)

cortexge3up<-read.table("UP.Cortex-HuFNSC03_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25")
cortexge3up$id<-unlist(strsplit(as.character(cortexge3up$V1),'_'))[2*(1:nrow(cortexge3up))]
cortexge3up$ave2<-cortex1[cortexge3up$id,3]
cortexge3up$ave3<-ge1[cortexge3up$id,3]
cortexge3up<-na.omit(cortexge3up)
cortexge3up<-cortexge3up[(cortexge3up$V2<=cutoff*cortexge3up$ave2&cortexge3up$V3>=cutoff*cortexge3up$ave3)|(cortexge3up$V2>=cutoff*cortexge3up$ave2&cortexge3up$V3<=cutoff*cortexge3up$ave3),]
cortexge3up_gene<-unique(cortexge3up$id)
cortexge3up_isoform<-cortexge3up[cortexge3up$ave2>Rmin&cortexge3up$ave3>Rmin,]
write.csv(cortexge3up,file="cortexge3up.csv",quote=F,row.names=F)
write.csv(cortexge3up_gene,file="cortexge3up_gene.csv",quote=F,row.names=F)
write.csv(cortexge3up_isoform,file="cortexge3up_isoform.csv",quote=F,row.names=F)
cortexge3dn<-read.table("DN.Cortex-HuFNSC03_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25")
cortexge3dn$id<-unlist(strsplit(as.character(cortexge3dn$V1),'_'))[2*(1:nrow(cortexge3dn))]
cortexge3dn$ave2<-cortex1[cortexge3dn$id,3]
cortexge3dn$ave3<-ge1[cortexge3dn$id,3]
cortexge3dn<-na.omit(cortexge3dn)
cortexge3dn<-cortexge3dn[(cortexge3dn$V2<=cutoff*cortexge3dn$ave2&cortexge3dn$V3>=cutoff*cortexge3dn$ave3)|(cortexge3dn$V2>=cutoff*cortexge3dn$ave2&cortexge3dn$V3<=cutoff*cortexge3dn$ave3),]
cortexge3dn_gene<-unique(cortexge3dn$id)
cortexge3dn_isoform<-cortexge3dn[cortexge3dn$ave2>Rmin&cortexge3dn$ave3>Rmin,]
write.csv(cortexge3dn,file="cortexge3dn.csv",quote=F,row.names=F)
write.csv(cortexge3dn_gene,file="cortexge3dn_gene.csv",quote=F,row.names=F)
write.csv(cortexge3dn_isoform,file="cortexge3dn_isoform.csv",quote=F,row.names=F)
cortexge3_isoform<-unique(c(cortexge3up_isoform$id,cortexge3dn_isoform$id))
write.csv(cortexge3_isoform,file="cortexge3_isoform.csv",quote=F,row.names=F)

cortexge4up<-read.table("UP.Cortex-HuFNSC04_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
cortexge4up$id<-unlist(strsplit(as.character(cortexge4up$V1),'_'))[2*(1:nrow(cortexge4up))]
cortexge4up$ave2<-cortex1[cortexge4up$id,3]
cortexge4up$ave3<-ge1[cortexge4up$id,3]
cortexge4up<-na.omit(cortexge4up)
cortexge4up<-cortexge4up[(cortexge4up$V2<=cutoff*cortexge4up$ave2&cortexge4up$V3>=cutoff*cortexge4up$ave3)|(cortexge4up$V2>=cutoff*cortexge4up$ave2&cortexge4up$V3<=cutoff*cortexge4up$ave3),]
cortexge4up_gene<-unique(cortexge4up$id)
cortexge4up_isoform<-cortexge4up[cortexge4up$ave2>Rmin&cortexge4up$ave3>Rmin,]
write.csv(cortexge4up,file="cortexge4up.csv",quote=F,row.names=F)
write.csv(cortexge4up_gene,file="cortexge4up_gene.csv",quote=F,row.names=F)
write.csv(cortexge4up_isoform,file="cortexge4up_isoform.csv",quote=F,row.names=F)
cortexge4dn<-read.table("DN.Cortex-HuFNSC04_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
cortexge4dn$id<-unlist(strsplit(as.character(cortexge4dn$V1),'_'))[2*(1:nrow(cortexge4dn))]
cortexge4dn$ave2<-cortex1[cortexge4dn$id,3]
cortexge4dn$ave3<-ge1[cortexge4dn$id,3]
cortexge4dn<-na.omit(cortexge4dn)
cortexge4dn<-cortexge4dn[(cortexge4dn$V2<=cutoff*cortexge4dn$ave2&cortexge4dn$V3>=cutoff*cortexge4dn$ave3)|(cortexge4dn$V2>=cutoff*cortexge4dn$ave2&cortexge4dn$V3<=cutoff*cortexge4dn$ave3),]
cortexge4dn_gene<-unique(cortexge4dn$id)
cortexge4dn_isoform<-cortexge4dn[cortexge4dn$ave2>Rmin&cortexge4dn$ave3>Rmin,]
write.csv(cortexge4dn,file="cortexge4dn.csv",quote=F,row.names=F)
write.csv(cortexge4dn_gene,file="cortexge4dn_gene.csv",quote=F,row.names=F)
write.csv(cortexge4dn_isoform,file="cortexge4dn_isoform.csv",quote=F,row.names=F)
cortexge4_isoform<-unique(c(cortexge4up_isoform$id,cortexge4dn_isoform$id))
write.csv(cortexge4_isoform,file="cortexge4_isoform.csv",quote=F,row.names=F)

# setwd("/Users/epigenomics_lab_02/FetalBrain/RNAseq/DEfine/exon/0.1ave/Cortex_GE/pc")
# cortexge1up_gene<-c(read.csv("cortexge1up_gene.csv",as.is=T)$x)
# cortexge2up_gene<-c(read.csv("cortexge2up_gene.csv",as.is=T)$x)
# cortexge3up_gene<-c(read.csv("cortexge3up_gene.csv",as.is=T)$x)
# cortexge4up_gene<-c(read.csv("cortexge4up_gene.csv",as.is=T)$x)
# cortexge1dn_gene<-c(read.csv("cortexge1dn_gene.csv",as.is=T)$x)
# cortexge2dn_gene<-c(read.csv("cortexge2dn_gene.csv",as.is=T)$x)
# cortexge3dn_gene<-c(read.csv("cortexge3dn_gene.csv",as.is=T)$x)
# cortexge4dn_gene<-c(read.csv("cortexge4dn_gene.csv",as.is=T)$x)
# 
# cortexge_up_gene<-intersect(intersect(cortexge1up_gene,cortexge2up_gene),intersect(cortexge3up_gene,cortexge4up_gene))
# write.csv(cortexge_up_gene,file="cortexge_up_gene.csv",quote=F,row.names=F)
# cortexge_dn_gene<-intersect(intersect(cortexge1dn_gene,cortexge2dn_gene),intersect(cortexge3dn_gene,cortexge4dn_gene))
# write.csv(cortexge_dn_gene,file="cortexge_dn_gene.csv",quote=F,row.names=F)
# cortexge_gene<-intersect(cortexge_up_gene,cortexge_dn_gene)
# write.csv(cortexge_gene,file="cortexge_gene.csv",quote=F,row.names=F)
# 
# cortexge12_up_gene<-intersect(cortexge1up_gene,cortexge2up_gene)
# cortexge34_up_gene<-intersect(cortexge3up_gene,cortexge4up_gene)
# cortexge12_dn_gene<-intersect(cortexge1dn_gene,cortexge2dn_gene)
# cortexge34_dn_gene<-intersect(cortexge3dn_gene,cortexge4dn_gene)

cortexge_up_gene_all<-c(cortexge1up_gene,cortexge2up_gene,cortexge3up_gene,cortexge4up_gene)
cortexge_up_gene_duplicate<-unique(cortexge_up_gene_all[duplicated(cortexge_up_gene_all)])
cortexge_dn_gene_all<-c(cortexge1dn_gene,cortexge2dn_gene,cortexge3dn_gene,cortexge4dn_gene)
cortexge_dn_gene_duplicate<-unique(cortexge_dn_gene_all[duplicated(cortexge_dn_gene_all)])
cortexge_gene_duplicate<-intersect(cortexge_up_gene_duplicate,cortexge_dn_gene_duplicate)
cortex_specific_gene<-setdiff(cortexge_up_gene_duplicate,cortexge_gene_duplicate)
ge_specific_gene<-setdiff(cortexge_dn_gene_duplicate,cortexge_gene_duplicate)
cortexge_isoform_all<-c(cortexge1_isoform,cortexge2_isoform,cortexge3_isoform,cortexge4_isoform)
cortexge_isoform_duplicate<-unique(cortexge_isoform_all[duplicated(cortexge_isoform_all)])

write.csv(cortexge_up_gene_duplicate,file="cortexge_up_gene_duplicate.csv",quote=F,row.names=F)
write.csv(cortexge_dn_gene_duplicate,file="cortexge_dn_gene_duplicate.csv",quote=F,row.names=F)
write.csv(cortexge_gene_duplicate,file="cortexge_gene_duplicate.csv",quote=F,row.names=F)
write.csv(cortex_specific_gene,file="cortex_specific_gene.csv",quote=F,row.names=F)
write.csv(ge_specific_gene,file="ge_specific_gene.csv",quote=F,row.names=F)
write.csv(cortexge_isoform_duplicate,file="cortexge_isoform_duplicate.csv",quote=F,row.names=F)

save.image(file="FetalBrain.isoform.individual.RData")

#between twins
brain12up<-read.table("UP.Brain-HuFNSC01_Brain-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
brain12dn<-read.table("DN.Brain-HuFNSC01_Brain-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
brain12up$id<-unlist(strsplit(as.character(brain12up$V1),'_'))[2*(1:nrow(brain12up))]
brain12up$ave2<-cortex1[brain12up$id,3]
brain12up$ave3<-ge1[brain12up$id,3]
brain12up<-na.omit(brain12up)
brain12up<-brain12up[(brain12up$V2<=cutoff*brain12up$ave2&brain12up$V3>=cutoff*brain12up$ave3)|(brain12up$V2>=cutoff*brain12up$ave2&brain12up$V3<=cutoff*brain12up$ave3),]
brain12up_gene<-unique(brain12up$id)
brain12up_isoform<-brain12up[brain12up$ave2>Rmin&brain12up$ave3>Rmin,]
write.csv(brain12up,file="brain12up.csv",quote=F,row.names=F)
write.csv(brain12up_gene,file="brain12up_gene.csv",quote=F,row.names=F)
write.csv(brain12up_isoform,file="brain12up_isoform.csv",quote=F,row.names=F)
brain12dn$id<-unlist(strsplit(as.character(brain12dn$V1),'_'))[2*(1:nrow(brain12dn))]
brain12dn$ave2<-cortex1[brain12dn$id,3]
brain12dn$ave3<-ge1[brain12dn$id,3]
brain12dn<-na.omit(brain12dn)
brain12dn<-brain12dn[(brain12dn$V2<=cutoff*brain12dn$ave2&brain12dn$V3>=cutoff*brain12dn$ave3)|(brain12dn$V2>=cutoff*brain12dn$ave2&brain12dn$V3<=cutoff*brain12dn$ave3),]
brain12dn_gene<-unique(brain12dn$id)
brain12dn_isoform<-brain12dn[brain12dn$ave2>Rmin&brain12dn$ave3>Rmin,]
write.csv(brain12dn,file="brain12dn.csv",quote=F,row.names=F)
write.csv(brain12dn_gene,file="brain12dn_gene.csv",quote=F,row.names=F)
write.csv(brain12dn_isoform,file="brain12dn_isoform.csv",quote=F,row.names=F)
brain12_isoform<-unique(c(brain12up_isoform$id,brain12dn_isoform$id))
write.csv(brain12_isoform,file="brain12_isoform.csv",quote=F,row.names=F)

cortex12up<-read.table("UP.Cortex-HuFNSC01_Cortex-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
cortex12dn<-read.table("DN.Cortex-HuFNSC01_Cortex-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
cortex12up$id<-unlist(strsplit(as.character(cortex12up$V1),'_'))[2*(1:nrow(cortex12up))]
cortex12up$ave2<-cortex1[cortex12up$id,3]
cortex12up$ave3<-ge1[cortex12up$id,3]
cortex12up<-na.omit(cortex12up)
cortex12up<-cortex12up[(cortex12up$V2<=cutoff*cortex12up$ave2&cortex12up$V3>=cutoff*cortex12up$ave3)|(cortex12up$V2>=cutoff*cortex12up$ave2&cortex12up$V3<=cutoff*cortex12up$ave3),]
cortex12up_gene<-unique(cortex12up$id)
cortex12up_isoform<-cortex12up[cortex12up$ave2>Rmin&cortex12up$ave3>Rmin,]
write.csv(cortex12up,file="cortex12up.csv",quote=F,row.names=F)
write.csv(cortex12up_gene,file="cortex12up_gene.csv",quote=F,row.names=F)
write.csv(cortex12up_isoform,file="cortex12up_isoform.csv",quote=F,row.names=F)
cortex12dn$id<-unlist(strsplit(as.character(cortex12dn$V1),'_'))[2*(1:nrow(cortex12dn))]
cortex12dn$ave2<-cortex1[cortex12dn$id,3]
cortex12dn$ave3<-ge1[cortex12dn$id,3]
cortex12dn<-na.omit(cortex12dn)
cortex12dn<-cortex12dn[(cortex12dn$V2<=cutoff*cortex12dn$ave2&cortex12dn$V3>=cutoff*cortex12dn$ave3)|(cortex12dn$V2>=cutoff*cortex12dn$ave2&cortex12dn$V3<=cutoff*cortex12dn$ave3),]
cortex12dn_gene<-unique(cortex12dn$id)
cortex12dn_isoform<-cortex12dn[cortex12dn$ave2>Rmin&cortex12dn$ave3>Rmin,]
write.csv(cortex12dn,file="cortex12dn.csv",quote=F,row.names=F)
write.csv(cortex12dn_gene,file="cortex12dn_gene.csv",quote=F,row.names=F)
write.csv(cortex12dn_isoform,file="cortex12dn_isoform.csv",quote=F,row.names=F)
cortex12_isoform<-unique(c(cortex12up_isoform$id,cortex12dn_isoform$id))
write.csv(cortex12_isoform,file="cortex12_isoform.csv",quote=F,row.names=F)

ge12up<-read.table("UP.GE-HuFNSC01_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
ge12dn<-read.table("DN.GE-HuFNSC01_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
ge12up$id<-unlist(strsplit(as.character(ge12up$V1),'_'))[2*(1:nrow(ge12up))]
ge12up$ave2<-cortex1[ge12up$id,3]
ge12up$ave3<-ge1[ge12up$id,3]
ge12up<-na.omit(ge12up)
ge12up<-ge12up[(ge12up$V2<=cutoff*ge12up$ave2&ge12up$V3>=cutoff*ge12up$ave3)|(ge12up$V2>=cutoff*ge12up$ave2&ge12up$V3<=cutoff*ge12up$ave3),]
ge12up_gene<-unique(ge12up$id)
ge12up_isoform<-ge12up[ge12up$ave2>Rmin&ge12up$ave3>Rmin,]
write.csv(ge12up,file="ge12up.csv",quote=F,row.names=F)
write.csv(ge12up_gene,file="ge12up_gene.csv",quote=F,row.names=F)
write.csv(ge12up_isoform,file="ge12up_isoform.csv",quote=F,row.names=F)
ge12dn$id<-unlist(strsplit(as.character(ge12dn$V1),'_'))[2*(1:nrow(ge12dn))]
ge12dn$ave2<-cortex1[ge12dn$id,3]
ge12dn$ave3<-ge1[ge12dn$id,3]
ge12dn<-na.omit(ge12dn)
ge12dn<-ge12dn[(ge12dn$V2<=cutoff*ge12dn$ave2&ge12dn$V3>=cutoff*ge12dn$ave3)|(ge12dn$V2>=cutoff*ge12dn$ave2&ge12dn$V3<=cutoff*ge12dn$ave3),]
ge12dn_gene<-unique(ge12dn$id)
ge12dn_isoform<-ge12dn[ge12dn$ave2>Rmin&ge12dn$ave3>Rmin,]
write.csv(ge12dn,file="ge12dn.csv",quote=F,row.names=F)
write.csv(ge12dn_gene,file="ge12dn_gene.csv",quote=F,row.names=F)
write.csv(ge12dn_isoform,file="ge12dn_isoform.csv",quote=F,row.names=F)
ge12_isoform<-unique(c(ge12up_isoform$id,ge12dn_isoform$id))
write.csv(ge12_isoform,file="ge12_isoform.csv",quote=F,row.names=F)

cortex34up<-read.table("UP.Cortex-HuFNSC03_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
cortex34dn<-read.table("DN.Cortex-HuFNSC03_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
cortex34up$id<-unlist(strsplit(as.character(cortex34up$V1),'_'))[2*(1:nrow(cortex34up))]
cortex34up$ave2<-cortex1[cortex34up$id,3]
cortex34up$ave3<-ge1[cortex34up$id,3]
cortex34up<-na.omit(cortex34up)
cortex34up<-cortex34up[(cortex34up$V2<=cutoff*cortex34up$ave2&cortex34up$V3>=cutoff*cortex34up$ave3)|(cortex34up$V2>=cutoff*cortex34up$ave2&cortex34up$V3<=cutoff*cortex34up$ave3),]
cortex34up_gene<-unique(cortex34up$id)
cortex34up_isoform<-cortex34up[cortex34up$ave2>Rmin&cortex34up$ave3>Rmin,]
write.csv(cortex34up,file="cortex34up.csv",quote=F,row.names=F)
write.csv(cortex34up_gene,file="cortex34up_gene.csv",quote=F,row.names=F)
write.csv(cortex34up_isoform,file="cortex34up_isoform.csv",quote=F,row.names=F)
cortex34dn$id<-unlist(strsplit(as.character(cortex34dn$V1),'_'))[2*(1:nrow(cortex34dn))]
cortex34dn$ave2<-cortex1[cortex34dn$id,3]
cortex34dn$ave3<-ge1[cortex34dn$id,3]
cortex34dn<-na.omit(cortex34dn)
cortex34dn<-cortex34dn[(cortex34dn$V2<=cutoff*cortex34dn$ave2&cortex34dn$V3>=cutoff*cortex34dn$ave3)|(cortex34dn$V2>=cutoff*cortex34dn$ave2&cortex34dn$V3<=cutoff*cortex34dn$ave3),]
cortex34dn_gene<-unique(cortex34dn$id)
cortex34dn_isoform<-cortex34dn[cortex34dn$ave2>Rmin&cortex34dn$ave3>Rmin,]
write.csv(cortex34dn,file="cortex34dn.csv",quote=F,row.names=F)
write.csv(cortex34dn_gene,file="cortex34dn_gene.csv",quote=F,row.names=F)
write.csv(cortex34dn_isoform,file="cortex34dn_isoform.csv",quote=F,row.names=F)
cortex34_isoform<-unique(c(cortex34up_isoform$id,cortex34dn_isoform$id))
write.csv(cortex34_isoform,file="cortex34_isoform.csv",quote=F,row.names=F)

ge34up<-read.table("UP.GE-HuFNSC03_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
ge34dn<-read.table("DN.GE-HuFNSC03_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
ge34up$id<-unlist(strsplit(as.character(ge34up$V1),'_'))[2*(1:nrow(ge34up))]
ge34up$ave2<-cortex1[ge34up$id,3]
ge34up$ave3<-ge1[ge34up$id,3]
ge34up<-na.omit(ge34up)
ge34up<-ge34up[(ge34up$V2<=cutoff*ge34up$ave2&ge34up$V3>=cutoff*ge34up$ave3)|(ge34up$V2>=cutoff*ge34up$ave2&ge34up$V3<=cutoff*ge34up$ave3),]
ge34up_gene<-unique(ge34up$id)
ge34up_isoform<-ge34up[ge34up$ave2>Rmin&ge34up$ave3>Rmin,]
write.csv(ge34up,file="ge34up.csv",quote=F,row.names=F)
write.csv(ge34up_gene,file="ge34up_gene.csv",quote=F,row.names=F)
write.csv(ge34up_isoform,file="ge34up_isoform.csv",quote=F,row.names=F)
ge34dn$id<-unlist(strsplit(as.character(ge34dn$V1),'_'))[2*(1:nrow(ge34dn))]
ge34dn$ave2<-cortex1[ge34dn$id,3]
ge34dn$ave3<-ge1[ge34dn$id,3]
ge34dn<-na.omit(ge34dn)
ge34dn<-ge34dn[(ge34dn$V2<=cutoff*ge34dn$ave2&ge34dn$V3>=cutoff*ge34dn$ave3)|(ge34dn$V2>=cutoff*ge34dn$ave2&ge34dn$V3<=cutoff*ge34dn$ave3),]
ge34dn_gene<-unique(ge34dn$id)
ge34dn_isoform<-ge34dn[ge34dn$ave2>Rmin&ge34dn$ave3>Rmin,]
write.csv(ge34dn,file="ge34dn.csv",quote=F,row.names=F)
write.csv(ge34dn_gene,file="ge34dn_gene.csv",quote=F,row.names=F)
write.csv(ge34dn_isoform,file="ge34dn_isoform.csv",quote=F,row.names=F)
ge34_isoform<-unique(c(ge34up_isoform$id,ge34dn_isoform$id))
write.csv(ge34_isoform,file="ge34_isoform.csv",quote=F,row.names=F)

# brain12<-rbind(brain12up,brain12dn)
# brain12_gene<-unique(c(brain12up_gene,brain12dn_gene))
# write.csv(brain12_gene,file="brain12_gene.csv",quote=F,row.names=F)
# cortex12<-rbind(cortex12up,cortex12dn)
# cortex12_gene<-unique(c(cortex12up_gene,cortex12dn_gene))
# write.csv(cortex12_gene,file="cortex12_gene.csv",quote=F,row.names=F)
# cortex34<-rbind(cortex34up,cortex34dn)
# cortex34_gene<-unique(c(cortex34up_gene,cortex34dn_gene))
# write.csv(cortex34_gene,file="cortex34_gene.csv",quote=F,row.names=F)
# ge12<-rbind(ge12up,ge12dn)
# ge12_gene<-unique(c(ge12up_gene,ge12dn_gene))
# write.csv(ge12_gene,file="ge12_gene.csv",quote=F,row.names=F)
# ge34<-rbind(ge34up,ge34dn)
# ge34_gene<-unique(c(ge34up_gene,ge34dn_gene))
# write.csv(ge34_gene,file="ge34_gene.csv",quote=F,row.names=F)

# setwd("/Users/epigenomics_lab_02/FetalBrain/RNAseq/DEfine/exon/0.1ave/BetweenTwins/pc")
# cortex12up_gene<-c(read.csv("cortex12up_gene.csv",as.is=T)$x)
# cortex12dn_gene<-c(read.csv("cortex12dn_gene.csv",as.is=T)$x)
# cortex34up_gene<-c(read.csv("cortex34up_gene.csv",as.is=T)$x)
# cortex34dn_gene<-c(read.csv("cortex34dn_gene.csv",as.is=T)$x)
# ge12up_gene<-c(read.csv("ge12up_gene.csv",as.is=T)$x)
# ge12dn_gene<-c(read.csv("ge12dn_gene.csv",as.is=T)$x)
# ge34up_gene<-c(read.csv("ge34up_gene.csv",as.is=T)$x)
# ge34dn_gene<-c(read.csv("ge34dn_gene.csv",as.is=T)$x)

# cortex12_ge12_isoform<-intersect(cortex12_isoform,ge12_isoform)
# cortex34_ge34_isoform<-intersect(cortex34_isoform,ge34_isoform)
# cortex12_ge12_cortex34_ge34_isoform<-intersect(cortex12_ge12_isoform,cortex34_ge34_isoform)

# unique_01<-intersect(cortex12up_gene,ge12up_gene)
# unique_02<-intersect(cortex12dn_gene,ge12dn_gene)
# unique_03<-intersect(cortex34up_gene,ge34up_gene)
# unique_04<-intersect(cortex34dn_gene,ge34dn_gene)
# as_12<-intersect(unique_01,unique_02)
# as_34<-intersect(unique_03,unique_04)
# as_12_34<-intersect(as_12,as_34)
# 
# cortex12_cortex34<-intersect(cortex12_gene,cortex34_gene)
# ge12_ge34<-intersect(ge12_gene,ge34_gene)
# brain12_cortex12<-intersect(brain12_gene,cortex12_gene)
# brain12_ge12<-intersect(brain12_gene,ge12_gene)
# cortex12_ge12<-intersect(cortex12_gene,ge12_gene)
# brain12_cortex12_ge12<-intersect(brain12_gene,cortex12_ge12)
# cortex34_ge34<-intersect(cortex34_gene,ge34_gene)
# brain12_cortex12_ge12_cortex34_ge34<-intersect(brain12_cortex12_ge12,cortex34_ge34)

# write.csv(cortex12_ge12_isoform,file="cortex12_ge12_isoform.csv",quote=F,row.names=F)
# write.csv(cortex34_ge34_isoform,file="cortex34_ge34_isoform.csv",quote=F,row.names=F)
# write.csv(cortex12_ge12_cortex34_ge34_isoform,file="cortex12_ge12_cortex34_ge34_isoform.csv",quote=F,row.names=F)
# write.csv(cortex12_cortex34,file="cortex12_cortex34.csv",quote=F,row.names=F)
# write.csv(ge12_ge34,file="ge12_ge34.csv",quote=F,row.names=F)
# write.csv(cortex12_ge12,file="cortex12_ge12.csv",quote=F,row.names=F)
# write.csv(brain12_cortex12,file="brain12_cortex12.csv",quote=F,row.names=F)
# write.csv(brain12_ge12,file="brain12_ge12.csv",quote=F,row.names=F)
# write.csv(brain12_cortex12_ge12,file="brain12_cortex12_ge12.csv",quote=F,row.names=F)
# write.csv(cortex34_ge34,file="cortex34_ge34.csv",quote=F,row.names=F)
# write.csv(brain12_cortex12_ge12_cortex34_ge34,file="brain12_cortex12_ge12_cortex34_ge34.csv",quote=F,row.names=F)
# write.csv(unique_01,file="unique_01.csv",quote=F,row.names=F)
# write.csv(unique_02,file="unique_02.csv",quote=F,row.names=F)
# write.csv(unique_03,file="unique_03.csv",quote=F,row.names=F)
# write.csv(unique_04,file="unique_04.csv",quote=F,row.names=F)
# write.csv(as_12,file="as_12.csv",quote=F,row.names=F)
# write.csv(as_34,file="as_34.csv",quote=F,row.names=F)
# write.csv(as_12_34,file="as_12_34.csv",quote=F,row.names=F)

twin01_02up<-intersect(as.character(brain12up$V1),intersect(as.character(cortex12up$V1),as.character(ge12up$V1)))
twin01_02up_gene<-intersect(brain12up_gene,intersect(cortex12up_gene,ge12up_gene))
twin01_02dn<-intersect(as.character(brain12dn$V1),intersect(as.character(cortex12dn$V1),as.character(ge12dn$V1)))
twin01_02dn_gene<-intersect(brain12dn_gene,intersect(cortex12dn_gene,ge12dn_gene))
twin01_02_isoform<-intersect(brain12_isoform,intersect(cortex12_isoform,ge12_isoform))

twin03_04up<-intersect(as.character(cortex34up$V1),as.character(ge34up$V1))
twin03_04up_gene<-intersect(cortex34up_gene,ge34up_gene)
twin03_04dn<-intersect(as.character(cortex34dn$V1),as.character(ge34dn$V1))
twin03_04dn_gene<-intersect(cortex34dn_gene,ge34dn_gene)
twin03_04_isoform<-intersect(cortex34_isoform,ge34_isoform)

twin12_twin34_isoform<-intersect(twin01_02_isoform,twin03_04_isoform)

write.csv(twin01_02up,file="twin01_02up.csv",quote=F,row.names=F)
write.csv(twin01_02up_gene,file="twin01_02up_gene.csv",quote=F,row.names=F)
write.csv(twin01_02dn,file="twin01_02dn.csv",quote=F,row.names=F)
write.csv(twin01_02dn_gene,file="twin01_02dn_gene.csv",quote=F,row.names=F)
write.csv(twin03_04up,file="twin03_04up.csv",quote=F,row.names=F)
write.csv(twin03_04up_gene,file="twin03_04up_gene.csv",quote=F,row.names=F)
write.csv(twin03_04dn,file="twin03_04dn.csv",quote=F,row.names=F)
write.csv(twin03_04dn_gene,file="twin03_04dn_gene.csv",quote=F,row.names=F)
write.csv(twin01_02_isoform,file="twin01_02_isoform.csv",quote=F,row.names=F)
write.csv(twin03_04_isoform,file="twin03_04_isoform.csv",quote=F,row.names=F)
write.csv(twin12_twin34_isoform,file="twin12_twin34_isoform.csv",quote=F,row.names=F)

save.image(file="FetalBrain.isoform.individual.RData")

##isoforms:DEfine on exons & plot cutoff
# cutoff: 
pdf("ave.distribution.pdf")
plot(range(0,40),range(0,0.06),xlab="RPKM",ylab="density",main="Distribution of ave. RPKM",type='n')
lines(density(cortex1[cortex1$V5>0,3]),col=1,lty=1)
lines(density(cortex2[cortex2$V5>0,3]),col=2,lty=1)
lines(density(cortex3[cortex3$V5>0,3]),col=3,lty=1)
lines(density(cortex4[cortex4$V5>0,3]),col=4,lty=1)
lines(density(ge1[ge1$V5>0,3]),col=1,lty=2)
lines(density(ge2[ge2$V5>0,3]),col=2,lty=2)
lines(density(ge3[ge3$V5>0,3]),col=3,lty=2)
lines(density(ge4[ge4$V5>0,3]),col=4,lty=2)
legend("topright",c("HuFNSC01","HuFNSC02","HuFNSC03","HuFNSC04"),col=c(1:4),cex=0.8,lty=c(1,1))
legend(34,0.05,c("Cortex","GE"),col=c(1,1),cex=0.8,lty=c(1,2))
dev.off()

cortex1$ave<-cortex1pc[cortex1$V2,4]
cortex1<-na.omit(cortex1)
cortex2$ave<-cortex2pc[cortex2$V2,4]
cortex2<-na.omit(cortex2)
cortex3$ave<-cortex3pc[cortex3$V2,4]
cortex3<-na.omit(cortex3)
cortex4$ave<-cortex4pc[cortex4$V2,4]
cortex4<-na.omit(cortex4)
ge1$ave<-ge1pc[ge1$V2,4]
ge1<-na.omit(ge1)
ge2$ave<-ge2pc[ge2$V2,4]
ge2<-na.omit(ge2)
ge3$ave<-ge3pc[ge3$V2,4]
ge3<-na.omit(ge3)
ge4$ave<-ge4pc[ge4$V2,4]
ge4<-na.omit(ge4)

pdf("exon.ditribution.pdf")
plot(range(0,80),range(0,0.2),xlab="RPKM",ylab="density",main="Distribution of exon RPKM",type='n')
lines(density(cortex1[cortex1$V4>0,4]),col=1,lty=1)
lines(density(cortex2[cortex2$V4>0,4]),col=2,lty=1)
lines(density(cortex3[cortex3$V4>0,4]),col=3,lty=1)
lines(density(cortex4[cortex4$V4>0,4]),col=4,lty=1)
lines(density(ge1[ge1$V4>0,4]),col=1,lty=2)
lines(density(ge2[ge2$V4>0,4]),col=2,lty=2)
lines(density(ge3[ge3$V4>0,4]),col=3,lty=2)
lines(density(ge4[ge4$V4>0,4]),col=4,lty=2)
legend("topright",c("HuFNSC01","HuFNSC02","HuFNSC03","HuFNSC04"),col=c(1:4),cex=0.8,lty=c(1,1))
legend(69,0.17,c("Cortex","GE"),col=c(1,1),cex=0.8,lty=c(1,2))
dev.off()

setwd("~/FetalBrain/RNAseq/DEfine/exon/0.1ave/Cortex_GE/pc")
hearing=read.csv("hearing.csv",head=F)
hearing=c(hearing[1,])
hearing=unlist(unique(hearing))
hearing2cortex=cortexge2up[is.element(cortexge2up$id,hearing),]
hearing2ge=cortexge2dn[is.element(cortexge2dn$id,hearing),]
hearing2=rbind(hearing2cortex,hearing2ge)
hearing2=hearing2[order(hearing2$id),]
write.csv(hearing2,file="hearing_HuFNSC02.csv",quote=F,row.names=F)


##################################################################################
setwd("~/FetalBrain/RNAseq/DEfine/")
# individual specific DE vs isoform summary plot
summary<-data.frame(tissue=c("Cortex_01vs02", "GE_01vs02","Cortex_03vs04","GE_03vs04"),
                    Isoform=c(2913,2545,2989,2087), 
                    DE_gene=c(596,173,642,545))
yy<-matrix(c(summary$Isoform,summary$DE_gene),2,4,byrow=T)
pdf("isoform_DE_individual.pdf",width=9,height=6)
par(mar = c(5,5,5,8), xpd=TRUE)
barplot(yy, beside=TRUE,col=c("blue","red"),ylim=c(0,3500),names.arg=summary$tissue, axis.lty=1, xlab=NULL, ylab="#genes",main="Individual specific isoforms and DE genes")  
legend("topright",c("Isoforms","DE genes"),inset=c(-0.1,0),col=c("blue","red"),lty=c(1,1),lwd=c(8,8),cex=0.8)
dev.off()
# tissue specific DE vs isoform summary plot
summary<-data.frame(tissue=c("01_CortexVsGE", "02_CortexVsGE","03_CortexVsGE","04_CortexVsGE"),
                    Isoform=c(2853,2828,3234,1834), 
                    DE_gene=c(911,1228,674,627))
yy<-matrix(c(summary$Isoform,summary$DE_gene),2,4,byrow=T)
pdf("isoform_DE_tissue.pdf",width=9,height=6)
par(mar = c(5,5,5,8), xpd=TRUE)
barplot(yy, beside=TRUE,col=c("blue","red"),ylim=c(0,3500),names.arg=summary$tissue, axis.lty=1, xlab=NULL, ylab="#genes",main="Tissue specific isoforms and DE genes")  
legend("topright",c("Isoforms","DE genes"),inset=c(-0.1,0),col=c("blue","red"),lty=c(1,1),lwd=c(8,8),cex=0.8)
dev.off()


