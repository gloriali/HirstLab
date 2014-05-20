setwd("/Users/epigenomics_lab_02/rotation/")
chrlen<-read.csv("chrlen_hg19.csv",row.names=1)
bin1000<-read.csv("bin1000.csv")
setwd("/Users/epigenomics_lab_02/rotation/chipseq")
brain01_summit<-read.table("H3K4me3_brain01_summits.bed")
brain02_summit<-read.table("H3K4me3_brain02_summits.bed")
colnames(brain01_summit)<-c("chr","start","end","id","height")
colnames(brain02_summit)<-c("chr","start","end","id","height")
brain01_peak<-read.table("H3K4me3_brain01_peaks.cod",header=T)
brain02_peak<-read.table("H3K4me3_brain02_peaks.cod",header=T)

H3K4me3<-bin1000
H3K4me3$brain_summit01<-0
H3K4me3$brain_summit02<-0
H3K4me3$brain_peak01<-0
H3K4me3$brain_peak02<-0

#summit height
#brain01_summit
chr<-as.character(brain01_summit$chr)
height<-brain01_summit$height
startbin<-chrlen[chr,5]  #first bin of the chr
start<-brain01_summit$start
bin<-startbin+start%/%1000    #row number 
temp<-data.frame(id=1:nrow(brain01_summit),bin=bin,chr=chr,height=height,start=start,startbin=startbin)
temp<-na.omit(temp)
if((temp$start<=H3K4me3$start[temp$bin])||(temp$start>=H3K4me3$end[temp$bin])||(paste("chr",temp$chr,sep="")!=as.character(H3K4me3$chr[temp$bin]))){
  error01<-temp$id
  temp$startbin<-NA
}
temp<-na.omit(temp)
H3K4me3$brain_summit01[temp$bin]<-H3K4me3$brain_summit01[temp$bin]+temp$height
#brain02_summit
chr<-as.character(brain02_summit$chr)
height<-brain02_summit$height
startbin<-chrlen[chr,5]  #first bin of the chr
start<-brain02_summit$start
bin<-startbin+start%/%1000    #row number 
temp<-data.frame(id=1:nrow(brain02_summit),bin=bin,chr=chr,height=height,start=start,startbin=startbin)
temp<-na.omit(temp)
if((temp$start<=H3K4me3$start[temp$bin])||(temp$start>=H3K4me3$end[temp$bin])||(paste("chr",temp$chr,sep="")!=as.character(H3K4me3$chr[temp$bin]))){
  error02<-temp$id
  temp$startbin<-NA
}
temp<-na.omit(temp)
H3K4me3$brain_summit02[temp$bin]<-H3K4me3$brain_summit02[temp$bin]+temp$height

#peak fold enrichment
#brain01_peak
chr<-as.character(brain01_peak$chr)
fold_enrichment<-brain01_peak$fold_enrichment
startbin<-chrlen[chr,5]  #first bin of the chr
start<-brain01_peak$start
end<-brain01_peak$end
binstart<-startbin+start%/%1000    #bin row number of the start point
binend<-startbin+end%/%1000    #bin row number of the end point
temp<-data.frame(id=1:nrow(brain01_peak),binstart=binstart,binend=binend,chr=chr,fold_enrichment=fold_enrichment,start=start,end=end,startbin=startbin)
temp<-na.omit(temp)
if((temp$start<=H3K4me3$start[temp$binstart])||(temp$start>=H3K4me3$end[temp$binstart])||(temp$end<=H3K4me3$start[temp$binend])||(temp$end>=H3K4me3$end[temp$binend])||(paste("chr",temp$chr,sep="")!=as.character(H3K4me3$chr[temp$binstart]))||(paste("chr",temp$chr,sep="")!=as.character(H3K4me3$chr[temp$binend]))){
  error03<-temp$id
  temp$startbin<-NA
}
temp<-na.omit(temp)
for(i in 1:nrow(temp)){
  print(i)
  binstart<-temp$binstart[i]
  binend<-temp$binend[i]
  fold_enrichment<-temp$fold_enrichment[i]
  H3K4me3$brain_peak01[binstart:binend]<-H3K4me3$brain_peak01[binstart:binend]+fold_enrichment
}
#brain02_peak
chr<-as.character(brain02_peak$chr)
fold_enrichment<-brain02_peak$fold_enrichment
startbin<-chrlen[chr,5]  #first bin of the chr
start<-brain02_peak$start
end<-brain02_peak$end
binstart<-startbin+start%/%1000    #bin row number of the start point
binend<-startbin+end%/%1000    #bin row number of the end point
temp<-data.frame(id=1:nrow(brain02_peak),binstart=binstart,binend=binend,chr=chr,fold_enrichment=fold_enrichment,start=start,end=end,startbin=startbin)
temp<-na.omit(temp)
if((temp$start<=H3K4me3$start[temp$binstart])||(temp$start>=H3K4me3$end[temp$binstart])||(temp$end<=H3K4me3$start[temp$binend])||(temp$end>=H3K4me3$end[temp$binend])||(paste("chr",temp$chr,sep="")!=as.character(H3K4me3$chr[temp$binstart]))||(paste("chr",temp$chr,sep="")!=as.character(H3K4me3$chr[temp$binend]))){
  error04<-temp$id
  temp$startbin<-NA
}
temp<-na.omit(temp)
for(i in 1:nrow(temp)){
  print(i)
  binstart<-temp$binstart[i]
  binend<-temp$binend[i]
  fold_enrichment<-temp$fold_enrichment[i]
  H3K4me3$brain_peak02[binstart:binend]<-H3K4me3$brain_peak02[binstart:binend]+fold_enrichment
}
save.image("H3K4me3.Rdata")
write.csv(H3K4me3,"H3K4me3.csv")

H3K4me3$logfold<-log2(H3K4me3$brain_peak01/H3K4me3$brain_peak02)
H3K4me3$logheight<-log2(H3K4me3$brain_summit01/H3K4me3$brain_summit02)
logfold<-na.omit(H3K4me3$logfold)
logheight<-na.omit(H3K4me3$logheight)
