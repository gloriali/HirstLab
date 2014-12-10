setwd("~/FetalBrain/RNAseq/rpkm/")
library(dplyr)
# gene expression matrix 
cortex01<-read.table(file="A03473.G.A.rpkm.pc",row.names=1,head=F)
cortex02<-read.table(file="A03475.G.A.rpkm.pc",row.names=1,head=F)
cortex03<-read.table(file="A04599.G.A.rpkm.pc",row.names=1,head=F)
cortex04<-read.table(file="A15298.G.A.rpkm.pc",row.names=1,head=F)
ge01<-read.table(file="A03474.G.A.rpkm.pc",row.names=1,head=F)
ge02<-read.table(file="A03476.G.A.rpkm.pc",row.names=1,head=F)
ge03<-read.table(file="A15295.G.A.rpkm.pc",row.names=1,head=F)
ge04<-read.table(file="A15299.G.A.rpkm.pc",row.names=1,head=F)
brain01<-read.table(file="A03484.G.A.rpkm.pc",row.names=1,head=F)
brain02<-read.table(file="A07825.G.A.rpkm.pc",row.names=1,head=F)
IDs<-intersect(rownames(cortex04),rownames(cortex03))
rpkm_pc<-data.frame(Ensembl = IDs, 
                    Cortex.HuFNSC01=cortex01[IDs,2],
                    Cortex.HuFNSC02=cortex02[IDs,2],
                    Cortex.HuFNSC03=cortex03[IDs,2],
                    Cortex.HuFNSC04=cortex04[IDs,2],
                    GE.HuFNSC01=ge01[IDs,2],
                    GE.HuFNSC02=ge02[IDs,2],
                    GE.HuFNSC03=ge03[IDs,2],
                    GE.HuFNSC04=ge04[IDs,2],
                    Brain.HuFNSC01=brain01[IDs,2],
                    Brain.HuFNSC02=brain02[IDs,2])
rownames(rpkm_pc)<-IDs
write.table(rpkm_pc, file="rpkm_pc.txt", sep = "\t", row.names = F, col.names = T, quote=F)
# exon expression matrix
cortex01<-read.table(file="A03473.G.exn.A.rpkm",head=F) %>% distinct(ID = paste0(V1, "_", V2))
cortex02<-read.table(file="A03475.G.exn.A.rpkm",head=F) %>% distinct(ID = paste0(V1, "_", V2))
cortex03<-read.table(file="A04599.G.exn.A.rpkm",head=F) %>% distinct(ID = paste0(V1, "_", V2))
cortex04<-read.table(file="A15298.G.exn.A.rpkm",head=F) %>% distinct(ID = paste0(V1, "_", V2))
ge01<-read.table(file="A03474.G.exn.A.rpkm",head=F) %>% distinct(ID = paste0(V1, "_", V2))
ge02<-read.table(file="A03476.G.exn.A.rpkm",head=F) %>% distinct(ID = paste0(V1, "_", V2))
ge03<-read.table(file="A15295.G.exn.A.rpkm",head=F) %>% distinct(ID = paste0(V1, "_", V2))
ge04<-read.table(file="A15299.G.exn.A.rpkm",head=F) %>% distinct(ID = paste0(V1, "_", V2))
brain01<-read.table(file="A03484.G.exn.A.rpkm",head=F) %>% distinct(ID = paste0(V1, "_", V2))
brain02<-read.table(file="A07825.G.exn.A.rpkm",head=F) %>% distinct(ID = paste0(V1, "_", V2))
rpkm_exon<-data.frame(ID = cortex01$ID, 
                    Cortex.HuFNSC01=cortex01$V4,
                    Cortex.HuFNSC02=cortex02$V4,
                    Cortex.HuFNSC03=cortex03$V4,
                    Cortex.HuFNSC04=cortex04$V4,
                    GE.HuFNSC01=ge01$V4,
                    GE.HuFNSC02=ge02$V4,
                    GE.HuFNSC03=ge03$V4,
                    GE.HuFNSC04=ge04$V4,
                    Brain.HuFNSC01=brain01$V4,
                    Brain.HuFNSC02=brain02$V4)
rownames(rpkm_exon)<-rpkm_exon$ID
write.table(rpkm_exon, file="rpkm_exon.txt", sep = "\t", row.names = F, col.names = T, quote=F)

# RPKM for 5mC regulators 
regulator_5mC <- read.table("~/快盘/hg19/regulators.5mC", head = F, as.is = T)
regulator_5mC <- data.frame(gene = regulator_5mC$V1, rpkm_pc[regulator_5mC$V2,])
write.table(regulator_5mC, file = "~/快盘/FetalBrain/RNAseq/regulator_5mC.rpkm", sep = "\t", row.names = F, col.names = T, quote = F)

##################################################################################################
#expression correlation
pdf(file="ExpressionCorrelation_RPKM1.pdf")
smoothScatter((logrpkm1$cortex01+logrpkm1$cortex02)/2,logrpkm1$cortex01-logrpkm1$cortex02,xlab="A",ylab="M",main="NeuroshpereCortex-BetweenTwins-01&02")
abline(h=1)
abline(h=-1)
smoothScatter((logrpkm1$cortex03+logrpkm1$cortex04)/2,logrpkm1$cortex03-logrpkm1$cortex04,xlab="A",ylab="M",main="NeuroshpereCortex-BetweenTwins-03&04")
abline(h=1)
abline(h=-1)
smoothScatter((logrpkm1$cortex01+logrpkm1$cortex03)/2,logrpkm1$cortex01-logrpkm1$cortex03,xlab="A",ylab="M",main="NeuroshpereCortex-AcrossTwins-01&03")
abline(h=1)
abline(h=-1)
smoothScatter((logrpkm1$cortex01+logrpkm1$cortex04)/2,logrpkm1$cortex01-logrpkm1$cortex04,xlab="A",ylab="M",main="NeuroshpereCortex-AcrossTwins-01&04")
abline(h=1)
abline(h=-1)
smoothScatter((logrpkm1$cortex02+logrpkm1$cortex03)/2,logrpkm1$cortex02-logrpkm1$cortex03,xlab="A",ylab="M",main="NeuroshpereCortex-AcrossTwins-02&03")
abline(h=1)
abline(h=-1)
smoothScatter((logrpkm1$cortex02+logrpkm1$cortex04)/2,logrpkm1$cortex02-logrpkm1$cortex04,xlab="A",ylab="M",main="NeuroshpereCortex-AcrossTwins-02&04")
abline(h=1)
abline(h=-1)

smoothScatter((logrpkm1$ge01+logrpkm1$ge02)/2,logrpkm1$ge01-logrpkm1$ge02,xlab="A",ylab="M",main="NeuroshpereGE-BetweenTwins-01&02")
abline(h=1)
abline(h=-1)
smoothScatter((logrpkm1$ge03+logrpkm1$ge04)/2,logrpkm1$ge03-logrpkm1$ge04,xlab="A",ylab="M",main="NeuroshpereGE-BetweenTwins-03&04")
abline(h=1)
abline(h=-1)
smoothScatter((logrpkm1$ge01+logrpkm1$ge03)/2,logrpkm1$ge01-logrpkm1$ge03,xlab="A",ylab="M",main="NeuroshpereGE-AcrossTwins-01&03")
abline(h=1)
abline(h=-1)
smoothScatter((logrpkm1$ge01+logrpkm1$ge04)/2,logrpkm1$ge01-logrpkm1$ge04,xlab="A",ylab="M",main="NeuroshpereGE-AcrossTwins-01&04")
abline(h=1)
abline(h=-1)
smoothScatter((logrpkm1$ge02+logrpkm1$ge03)/2,logrpkm1$ge02-logrpkm1$ge03,xlab="A",ylab="M",main="NeuroshpereGE-AcrossTwins-02&03")
abline(h=1)
abline(h=-1)
smoothScatter((logrpkm1$ge02+logrpkm1$ge04)/2,logrpkm1$ge02-logrpkm1$ge04,xlab="A",ylab="M",main="NeuroshpereGE-AcrossTwins-02&04")
abline(h=1)
abline(h=-1)

smoothScatter((logrpkm1$brain01+logrpkm1$brain02)/2,logrpkm1$brain01-logrpkm1$brain02,xlab="A",ylab="M",main="Brain-BetweenTwins-01&02")
abline(h=1)
abline(h=-1)
dev.off()

##########################################################################################################
##differential expression:fold change & z-score
setwd("~/FetalBrain/RNAseq/rpkm/")
rpkm1<-read.csv("rpkm1.csv",as.is=T,head=T,row.names=1)
setwd("~/FetalBrain/RNAseq/DE_fold/")
#brain
expression_brain<-data.frame(brain01=rpkm1$brain01,brain02=rpkm1$brain02)
rownames(expression_brain)<-rownames(rpkm1)
expression_brain$logfold<-log2(expression_brain$brain01/expression_brain$brain02)
expression_brain$z<-abs(expression_brain$brain01-expression_brain$brain02)/sqrt(expression_brain$brain01+(174797872/137879052)*expression_brain$brain02)
expression_brain<-expression_brain[is.finite(expression_brain$logfold),]
diff_ex_brain<-expression_brain[((abs(expression_brain$logfold)>1)&(expression_brain$z>0.75)),]
write.csv(diff_ex_brain,file="diff_ex_brain.csv",quote=F)

#cortex
expression_cortex12<-data.frame(cortex01=rpkm1$cortex01,cortex02=rpkm1$cortex02)
rownames(expression_cortex12)<-rownames(rpkm1)
expression_cortex12$logfold<-log2(expression_cortex12$cortex01/expression_cortex12$cortex02)
expression_cortex12$z<-abs(expression_cortex12$cortex01-expression_cortex12$cortex02)/sqrt(expression_cortex12$cortex01+(154389698/182310872)*expression_cortex12$cortex02)
expression_cortex12<-expression_cortex12[is.finite(expression_cortex12$logfold),]
diff_ex_cortex12<-expression_cortex12[((abs(expression_cortex12$logfold)>1)&(expression_cortex12$z>0.75)),]
write.csv(diff_ex_cortex12,file="diff_ex_cortex12.csv",quote=F)

expression_cortex34<-data.frame(cortex03=rpkm1$cortex03,cortex04=rpkm1$cortex04)
rownames(expression_cortex34)<-rownames(rpkm1)
expression_cortex34$logfold<-log2(expression_cortex34$cortex03/expression_cortex34$cortex04)
expression_cortex34$z<-abs(expression_cortex34$cortex03-expression_cortex34$cortex04)/sqrt(expression_cortex34$cortex03+(192334650/450442130)*expression_cortex34$cortex04)
expression_cortex34<-expression_cortex34[is.finite(expression_cortex34$logfold),]
diff_ex_cortex34<-expression_cortex34[((abs(expression_cortex34$logfold)>1)&(expression_cortex34$z>0.75)),]
write.csv(diff_ex_cortex34,file="diff_ex_cortex34.csv",quote=F)

diff_ex_cortex<-intersect(rownames(diff_ex_cortex12),rownames(diff_ex_cortex34))
write.csv(diff_ex_cortex,file="diff_ex_cortex.csv",quote=F)

#ge
expression_ge12<-data.frame(ge01=rpkm1$ge01,ge02=rpkm1$ge02)
rownames(expression_ge12)<-rownames(rpkm1)
expression_ge12$logfold<-log2(expression_ge12$ge01/expression_ge12$ge02)
expression_ge12$z<-abs(expression_ge12$ge01-expression_ge12$ge02)/sqrt(expression_ge12$ge01+(191064090/219123224)*expression_ge12$ge02)
expression_ge12<-expression_ge12[is.finite(expression_ge12$logfold),]
diff_ex_ge12<-expression_ge12[((abs(expression_ge12$logfold)>1)&(expression_ge12$z>0.75)),]
write.csv(diff_ex_ge12,file="diff_ex_ge12.csv",quote=F)

expression_ge34<-data.frame(ge03=rpkm1$ge03,ge04=rpkm1$ge04)
rownames(expression_ge34)<-rownames(rpkm1)
expression_ge34$logfold<-log2(expression_ge34$ge03/expression_ge34$ge04)
expression_ge34$z<-abs(expression_ge34$ge03-expression_ge34$ge04)/sqrt(expression_ge34$ge03+(442028153/386943893)*expression_ge34$ge04)
expression_ge34<-expression_ge34[is.finite(expression_ge34$logfold),]
diff_ex_ge34<-expression_ge34[((abs(expression_ge34$logfold)>1)&(expression_ge34$z>0.75)),]
write.csv(diff_ex_ge34,file="diff_ex_ge34.csv",quote=F)

diff_ex_ge<-intersect(rownames(diff_ex_ge12),rownames(diff_ex_ge34))
write.csv(diff_ex_ge,file="diff_ex_ge.csv",quote=F)

diff_brain_cortex12<-intersect(rownames(diff_ex_cortex12),rownames(diff_ex_brain))
write.csv(diff_brain_cortex12,file="diff_brain_cortex12.csv",quote=F)
diff_brain_ge12<-intersect(rownames(diff_ex_ge12),rownames(diff_ex_brain))
write.csv(diff_brain_ge12,file="diff_brain_ge12.csv",quote=F)
diff_ge12_cortex12<-intersect(rownames(diff_ex_cortex12),rownames(diff_ex_ge12))
write.csv(diff_ge12_cortex12,file="diff_ge12_cortex12.csv",quote=F)
diff_ge34_cortex34<-intersect(rownames(diff_ex_cortex34),rownames(diff_ex_ge34))
write.csv(diff_ge34_cortex34,file="diff_ge34_cortex34.csv",quote=F)

fisher_cortex<-matrix(c(76,767,126,18826),nrow=2)
fisher.test(fisher_cortex)
fisher_ge<-matrix(c(3,331,22,19439),nrow=2)
fisher.test(fisher_ge)
fisher_braincortex<-matrix(c(27,175,797,18796),nrow=2)
fisher.test(fisher_braincortex)
fisher_brainge<-matrix(c(4,21,820,18950),nrow=2)
fisher.test(fisher_brainge)
fisher_cortexge12<-matrix(c(3,22,199,19571),nrow=2)
fisher.test(fisher_cortexge12)
fisher_cortexge34<-matrix(c(137,197,706,18755),nrow=2)
fisher.test(fisher_cortexge34)

save.image("FetalBrain.DE.fold.Rdata")

##################################################################################
##differential expression:DEfine
# between twins
setwd("~/FetalBrain/RNAseq/DEfine/gene/individual/")
col<-c("Gene id", "rpkm1", "rpkm2", "p-value", "Multiple testing corrected p-value")
brain12up<-read.table("UP.Brain-HuFNSC01_Brain-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
rownames(brain12up)<-brain12up$V1
brain12dn<-read.table("DN.Brain-HuFNSC01_Brain-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
rownames(brain12dn)<-brain12dn$V1
cortex12up<-read.table("UP.Cortex-HuFNSC01_Cortex-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
rownames(cortex12up)<-cortex12up$V1
cortex12dn<-read.table("DN.Cortex-HuFNSC01_Cortex-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
rownames(cortex12dn)<-cortex12dn$V1
cortex34up<-read.table("UP.Cortex-HuFNSC03_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
rownames(cortex34up)<-cortex34up$V1
cortex34dn<-read.table("DN.Cortex-HuFNSC03_Cortex-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
rownames(cortex34dn)<-cortex34dn$V1
ge12up<-read.table("UP.GE-HuFNSC01_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
rownames(ge12up)<-ge12up$V1
ge12dn<-read.table("DN.GE-HuFNSC01_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
rownames(ge12dn)<-ge12dn$V1
ge34up<-read.table("UP.GE-HuFNSC03_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
rownames(ge34up)<-ge34up$V1
ge34dn<-read.table("DN.GE-HuFNSC03_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
rownames(ge34dn)<-ge34dn$V1
brain12<-rbind(brain12up,brain12dn)
rownames(brain12)<-brain12$V1
cortex12<-rbind(cortex12up,cortex12dn)
rownames(cortex12)<-cortex12$V1
cortex34<-rbind(cortex34up,cortex34dn)
rownames(cortex34)<-cortex34$V1
ge12<-rbind(ge12up,ge12dn)
rownames(ge12)<-ge12$V1
ge34<-rbind(ge34up,ge34dn)
rownames(ge34)<-ge34$V1

diff_brain_cortex12<-intersect(rownames(cortex12),rownames(brain12))
write.csv(diff_brain_cortex12,file="diff_brain_cortex12_0.05.csv",quote=F)
diff_brain_ge12<-intersect(rownames(ge12),rownames(brain12))
write.csv(diff_brain_ge12,file="diff_brain_ge12_0.05.csv",quote=F)
diff_ge12_cortex12<-intersect(rownames(cortex12),rownames(ge12))
write.csv(diff_ge12_cortex12,file="diff_ge12_cortex12_0.05.csv",quote=F)
diff_ge34_cortex34<-intersect(rownames(cortex34),rownames(ge34))
write.csv(diff_ge34_cortex34,file="diff_ge34_cortex34_0.05.csv",quote=F)
diff_cortex12_cortex34<-intersect(rownames(cortex34),rownames(cortex12))
write.csv(diff_cortex12_cortex34,file="diff_cortex12_cortex34_0.05.csv",quote=F)
diff_ge12_ge34<-intersect(rownames(ge34),rownames(ge12))
write.csv(diff_ge12_ge34,file="diff_ge12_ge34_0.05.csv",quote=F)
diff_all12<-intersect(diff_brain_cortex12,rownames(ge12))
write.csv(diff_all12,file="diff_all12_0.05.csv",quote=F)
diff_all<-intersect(diff_all12,diff_ge34_cortex34)
write.csv(diff_all,file="diff_all_0.05.csv",quote=F)

cortex_GEup_12<-cbind(cortex12up[intersect(cortex12up$V1,ge12up$V1),],ge12up[intersect(cortex12up$V1,ge12up$V1),])
cortex_GEdn_12<-cbind(cortex12dn[intersect(cortex12dn$V1,ge12dn$V1),],ge12dn[intersect(cortex12dn$V1,ge12dn$V1),])
cortex_GEup_34<-cbind(cortex34up[intersect(cortex34up$V1,ge34up$V1),],ge34up[intersect(cortex34up$V1,ge34up$V1),])
cortex_GEdn_34<-cbind(cortex34dn[intersect(cortex34dn$V1,ge34dn$V1),],ge34dn[intersect(cortex34dn$V1,ge34dn$V1),])

#get gene information
setwd("~/hg19/")
gene<-read.table(file="hg19v65_genes",as.is=T,head=F)
colnames(gene)<-c("Ensembl","chr","start","end","strand","type","name","description")
gene<-gene[!duplicated(gene$Ensembl),]
rownames(gene)<-gene$Ensembl
setwd("~/FetalBrain/RNAseq/DEfine/gene/individual/")
cortex_GEup_12$Name=gene[cortex_GEup_12$V1,7];cortex_GEup_12$Description=gene[cortex_GEup_12$V1,8]
cortex_GEup_34$Name=gene[cortex_GEup_34$V1,7];cortex_GEup_34$Description=gene[cortex_GEup_34$V1,8]
cortex_GEdn_12$Name=gene[cortex_GEdn_12$V1,7];cortex_GEdn_12$Description=gene[cortex_GEdn_12$V1,8]
cortex_GEdn_34$Name=gene[cortex_GEdn_34$V1,7];cortex_GEdn_34$Description=gene[cortex_GEdn_34$V1,8]

cortex_GE_12<-rbind(cortex_GEup_12,cortex_GEdn_12)
cortex_GE_34<-rbind(cortex_GEup_34,cortex_GEdn_34)
write.csv(cortex_GE_12,file="cortex_GE_12.csv",quote=F)
write.csv(cortex_GE_34,file="cortex_GE_34.csv",quote=F)

save.image("FetalBrain.DEfine.gene.individual.Rdata")

# between tissues: Cortex VS GE
setwd("~/FetalBrain/RNAseq/DEfine/gene/cortexge/")
col<-c("Gene id", "rpkm1", "rpkm2", "p-value", "Multiple testing corrected p-value")
cortexge1up<-read.table("UP.Cortex-HuFNSC01_GE-HuFNSC01.FDR_0.01.rmin_0.005.Nmin_25")
cortexge1dn<-read.table("DN.Cortex-HuFNSC01_GE-HuFNSC01.FDR_0.01.rmin_0.005.Nmin_25")
cortexge2up<-read.table("UP.Cortex-HuFNSC02_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
cortexge2dn<-read.table("DN.Cortex-HuFNSC02_GE-HuFNSC02.FDR_0.01.rmin_0.005.Nmin_25")
cortexge3up<-read.table("UP.Cortex-HuFNSC03_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25")
cortexge3dn<-read.table("DN.Cortex-HuFNSC03_GE-HuFNSC03.FDR_0.01.rmin_0.005.Nmin_25")
cortexge4up<-read.table("UP.Cortex-HuFNSC04_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")
cortexge4dn<-read.table("DN.Cortex-HuFNSC04_GE-HuFNSC04.FDR_0.01.rmin_0.005.Nmin_25")

cortexge_up<-intersect(intersect(cortexge1up$V1,cortexge2up$V1),intersect(cortexge3up$V1,cortexge4up$V1))
write.csv(cortexge_up,file="cortexge_up.csv",quote=F,row.names=F)
cortexge_dn<-intersect(intersect(cortexge1dn$V1,cortexge2dn$V1),intersect(cortexge3dn$V1,cortexge4dn$V1))
write.csv(cortexge_dn,file="cortexge_dn.csv",quote=F,row.names=F)
#cortexge<-intersect(cortexge_up,cortexge_dn)
cortexge_up12<-intersect(cortexge1up$V1,cortexge2up$V1)
cortexge_up34<-intersect(cortexge3up$V1,cortexge4up$V1)
cortexge_dn12<-intersect(cortexge1dn$V1,cortexge2dn$V1)
cortexge_dn34<-intersect(cortexge3dn$V1,cortexge4dn$V1)
write.csv(cortexge_up12,file="cortexge_up12.csv",quote=F,row.names=F)
write.csv(cortexge_up34,file="cortexge_up34.csv",quote=F,row.names=F)
write.csv(cortexge_dn12,file="cortexge_dn12.csv",quote=F,row.names=F)
write.csv(cortexge_dn34,file="cortexge_dn34.csv",quote=F,row.names=F)
cortexge_up_all<-c(as.character(cortexge1up$V1),as.character(cortexge2up$V1),as.character(cortexge3up$V1),as.character(cortexge4up$V1))
cortexge_up_duplicate<-unique(cortexge_up_all[duplicated(cortexge_up_all)])
cortexge_dn_all<-c(as.character(cortexge1dn$V1),as.character(cortexge2dn$V1),as.character(cortexge3dn$V1),as.character(cortexge4dn$V1))
cortexge_dn_duplicate<-unique(cortexge_dn_all[duplicated(cortexge_dn_all)])
#get gene information
setwd("~/hg19/")
gene<-read.table(file="hg19v65_genes",as.is=T,head=F)
colnames(gene)<-c("Ensembl","chr","start","end","strand","type","name","description")
gene<-gene[!duplicated(gene$Ensembl),]
rownames(gene)<-gene$Ensembl
setwd("~/FetalBrain/RNAseq/DEfine/gene/cortexge/")
cortexge_up_duplicate<-data.frame(ID=cortexge_up_duplicate,Name=gene[cortexge_up_duplicate,7],Description=gene[cortexge_up_duplicate,8])
cortexge_dn_duplicate<-data.frame(ID=cortexge_dn_duplicate,Name=gene[cortexge_dn_duplicate,7],Description=gene[cortexge_dn_duplicate,8])

write.csv(cortexge_up_duplicate,file="cortexge_up_duplicate.csv",quote=F,row.names=F)
write.csv(cortexge_dn_duplicate,file="cortexge_dn_duplicate.csv",quote=F,row.names=F)

save.image("FetalBrain.DEfine.gene.tissue.Rdata")
