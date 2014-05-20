setwd("/Users/epigenomics_lab_02/FetalBrain")
hismod<-read.csv("HistoneMods.csv", as.is=T, fill=T)

for (i in 1: nrow(hismod)){
  id<-hismod$Library[i]
  target<-hismod$Target[i]
  twin<-hismod$Original.Source[i]
  tissue<-hismod$Tissue[i]
  bw<-hismod$bigWig[i]
  if (twin=='HuFNSC01'){
    cl<-'200,50,0'
  }
  else{
    cl<-'50,200,50'
  }
  cat("track",id,"\n")
  cat("shortLabel",id,"\n")
  cat("longLabel",target,twin,tissue,"\n")
  cat("type bigWig\nvisibility full\nmaxHeightPixels 70:70:32\nconfigurable on\nautoScale on\nalwaysZero on\npriority 0.1\n")
  cat("bigDataUrl",bw,"\n")
  cat("color",cl,"\n")
  cat("\n")
}

