setwd("/Users/epigenomics_lab_02/rotation")
dataset<-read.csv("FetalBrain.csv",as.is=T,fill=T)
dataset$X<-NULL
dataset$inx<-NA
inx<-read.csv("inx.csv",as.is=T,fill=T)
rownames(inx)<-as.character(inx$Library)

for(i in 1:nrow(dataset)){
  lib<-dataset$Library[i]
  dataset$inx[i]<-as.character(inx[lib,2])
}
write.csv(dataset,file="FetalBrain_inx.csv")
