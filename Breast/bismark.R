##################################################################################
# read coverage and methylation call distribution check
# run on xhost
# /gsc/software/linux-x86_64-centos5/R-2.15.2/bin/R
# chr1rc<-read.table(paste(dirIn,"rc/chr1.rc",sep=""))
# chr1fm<-read.table(paste(dirIn,"fm/chr1.fm",sep=""))
# chr1<-data.frame(chr="chr1",position=chr1rc$V1,rc=chr1rc$V2,fm=chr1fm$V2)
# chr2rc<-read.table(paste(dirIn,"rc/chr2.rc",sep=""))
# chr2fm<-read.table(paste(dirIn,"fm/chr2.fm",sep=""))
# chr2<-data.frame(chr="chr2",position=chr2rc$V1,rc=chr2rc$V2,fm=chr2fm$V2)
# chr3rc<-read.table(paste(dirIn,"rc/chr3.rc",sep=""))
# chr3fm<-read.table(paste(dirIn,"fm/chr3.fm",sep=""))
# chr3<-data.frame(chr="chr3",position=chr3rc$V1,rc=chr3rc$V2,fm=chr3fm$V2)
# chr4rc<-read.table(paste(dirIn,"rc/chr4.rc",sep=""))
# chr4fm<-read.table(paste(dirIn,"fm/chr4.fm",sep=""))
# chr4<-data.frame(chr="chr4",position=chr4rc$V1,rc=chr4rc$V2,fm=chr4fm$V2)
# chr5rc<-read.table(paste(dirIn,"rc/chr5.rc",sep=""))
# chr5fm<-read.table(paste(dirIn,"fm/chr5.fm",sep=""))
# chr5<-data.frame(chr="chr5",position=chr5rc$V1,rc=chr5rc$V2,fm=chr5fm$V2)
# chr6rc<-read.table(paste(dirIn,"rc/chr6.rc",sep=""))
# chr6fm<-read.table(paste(dirIn,"fm/chr6.fm",sep=""))
# chr6<-data.frame(chr="chr6",position=chr6rc$V1,rc=chr6rc$V2,fm=chr6fm$V2)
# chr7rc<-read.table(paste(dirIn,"rc/chr7.rc",sep=""))
# chr7fm<-read.table(paste(dirIn,"fm/chr7.fm",sep=""))
# chr7<-data.frame(chr="chr7",position=chr7rc$V1,rc=chr7rc$V2,fm=chr7fm$V2)
# chr8rc<-read.table(paste(dirIn,"rc/chr8.rc",sep=""))
# chr8fm<-read.table(paste(dirIn,"fm/chr8.fm",sep=""))
# chr8<-data.frame(chr="chr8",position=chr8rc$V1,rc=chr8rc$V2,fm=chr8fm$V2)
# chr9rc<-read.table(paste(dirIn,"rc/chr9.rc",sep=""))
# chr9fm<-read.table(paste(dirIn,"fm/chr9.fm",sep=""))
# chr9<-data.frame(chr="chr9",position=chr9rc$V1,rc=chr9rc$V2,fm=chr9fm$V2)
# chr10rc<-read.table(paste(dirIn,"rc/chr10.rc",sep=""))
# chr10fm<-read.table(paste(dirIn,"fm/chr10.fm",sep=""))
# chr10<-data.frame(chr="chr10",position=chr10rc$V1,rc=chr10rc$V2,fm=chr10fm$V2)
# chr11rc<-read.table(paste(dirIn,"rc/chr11.rc",sep=""))
# chr11fm<-read.table(paste(dirIn,"fm/chr11.fm",sep=""))
# chr11<-data.frame(chr="chr11",position=chr11rc$V1,rc=chr11rc$V2,fm=chr11fm$V2)
# chr12rc<-read.table(paste(dirIn,"rc/chr12.rc",sep=""))
# chr12fm<-read.table(paste(dirIn,"fm/chr12.fm",sep=""))
# chr12<-data.frame(chr="chr12",position=chr12rc$V1,rc=chr12rc$V2,fm=chr12fm$V2)
# chr13rc<-read.table(paste(dirIn,"rc/chr13.rc",sep=""))
# chr13fm<-read.table(paste(dirIn,"fm/chr13.fm",sep=""))
# chr13<-data.frame(chr="chr13",position=chr13rc$V1,rc=chr13rc$V2,fm=chr13fm$V2)
# chr14rc<-read.table(paste(dirIn,"rc/chr14.rc",sep=""))
# chr14fm<-read.table(paste(dirIn,"fm/chr14.fm",sep=""))
# chr14<-data.frame(chr="chr14",position=chr14rc$V1,rc=chr14rc$V2,fm=chr14fm$V2)
# chr15rc<-read.table(paste(dirIn,"rc/chr15.rc",sep=""))
# chr15fm<-read.table(paste(dirIn,"fm/chr15.fm",sep=""))
# chr15<-data.frame(chr="chr15",position=chr15rc$V1,rc=chr15rc$V2,fm=chr15fm$V2)
# chr16rc<-read.table(paste(dirIn,"rc/chr16.rc",sep=""))
# chr16fm<-read.table(paste(dirIn,"fm/chr16.fm",sep=""))
# chr16<-data.frame(chr="chr16",position=chr16rc$V1,rc=chr16rc$V2,fm=chr16fm$V2)
# chr17rc<-read.table(paste(dirIn,"rc/chr17.rc",sep=""))
# chr17fm<-read.table(paste(dirIn,"fm/chr17.fm",sep=""))
# chr17<-data.frame(chr="chr17",position=chr17rc$V1,rc=chr17rc$V2,fm=chr17fm$V2)
# chr18rc<-read.table(paste(dirIn,"rc/chr18.rc",sep=""))
# chr18fm<-read.table(paste(dirIn,"fm/chr18.fm",sep=""))
# chr18<-data.frame(chr="chr18",position=chr18rc$V1,rc=chr18rc$V2,fm=chr18fm$V2)
# chr19rc<-read.table(paste(dirIn,"rc/chr19.rc",sep=""))
# chr19fm<-read.table(paste(dirIn,"fm/chr19.fm",sep=""))
# chr19<-data.frame(chr="chr19",position=chr19rc$V1,rc=chr19rc$V2,fm=chr19fm$V2)
# chr20rc<-read.table(paste(dirIn,"rc/chr20.rc",sep=""))
# chr20fm<-read.table(paste(dirIn,"fm/chr20.fm",sep=""))
# chr20<-data.frame(chr="chr20",position=chr20rc$V1,rc=chr20rc$V2,fm=chr20fm$V2)
# chr21rc<-read.table(paste(dirIn,"rc/chr21.rc",sep=""))
# chr21fm<-read.table(paste(dirIn,"fm/chr21.fm",sep=""))
# chr21<-data.frame(chr="chr21",position=chr21rc$V1,rc=chr21rc$V2,fm=chr21fm$V2)
# chr22rc<-read.table(paste(dirIn,"rc/chr22.rc",sep=""))
# chr22fm<-read.table(paste(dirIn,"fm/chr22.fm",sep=""))
# chr22<-data.frame(chr="chr22",position=chr22rc$V1,rc=chr22rc$V2,fm=chr22fm$V2)
# chrXrc<-read.table(paste(dirIn,"rc/chrX.rc",sep=""))
# chrXfm<-read.table(paste(dirIn,"fm/chrX.fm",sep=""))
# chrX<-data.frame(chr="chrX",position=chrXrc$V1,rc=chrXrc$V2,fm=chrXfm$V2)
# chrYrc<-read.table(paste(dirIn,"rc/chrY.rc",sep=""))
# chrYfm<-read.table(paste(dirIn,"fm/chrY.fm",sep=""))
# chrY<-data.frame(chr="chrY",position=chrYrc$V1,rc=chrYrc$V2,fm=chrYfm$V2)

# setwd("~/REMC/breast/bismark/myoRM045/")
# dirIn="/projects/mbilenky/REMC/breast/WGBS/bismark/myoRM045/"
setwd("~/REMC/breast/bismark/lumRM066/")
dirIn="/projects/mbilenky/REMC/breast/WGBS/bismark/lumRM066/"
chrs<-c(1:22)
chrs<-c(as.character(chrs),"X","Y")

# Ecdf plot for read coverage
pdf("coverage.ecdf.pdf")
plot(range(0,50),range(0,1),type="n",main="Ecdf coverage",xlab="coverage",ylab="Ecdf")
i=1
for(chr in chrs){
  chr_rc<-read.table(paste(dirIn,"rc/chr",chr,".rc",sep=""))
  chr_fm<-read.table(paste(dirIn,"fm/chr",chr,".fm",sep=""))
  chrom<-data.frame(chrom=paste("chr",chr,sep=""),position=chr_rc$V1,rc=chr_rc$V2,fm=chr_fm$V2)
  lines(ecdf(chrom$rc),col=i)
  i=i+1
}
abline(v=5)
legend("bottomright",chrs,col=c(1:length(chrs)),cex=0.8,lty=1)
dev.off()

# Ecdf plot for methylation level (read coverage>=5)
pdf("methylation.ecdf.pdf")
plot(range(0,1),range(0,1),type="n",main="Ecdf methylation level(coverage>=5)",xlab="methylation level",ylab="Ecdf")
i=1
for(chr in chrs){
  chr_rc<-read.table(paste(dirIn,"rc/chr",chr,".rc",sep=""))
  chr_fm<-read.table(paste(dirIn,"fm/chr",chr,".fm",sep=""))
  chrom<-data.frame(chrom=chr,position=chr_rc$V1,rc=chr_rc$V2,fm=chr_fm$V2)
  write.table(chrom,file=paste("chr",chr,sep=""),sep="\t",quote=F,row.names=F)
  chrom<-chrom[chrom$rc>=5,]
  lines(ecdf(chrom$fm),col=i)
  i=i+1
}
legend("bottomright",chrs,col=c(1:length(chrs)),cex=0.5,lty=1)
dev.off()

##################################################################################
# get read coverage and methylation call (rc>=5)
# run on xhost
# /gsc/software/linux-x86_64-centos5/R-2.15.2/bin/R
# order exon file by chr
# setwd("~/REMC/breast/bismark/")
# # exon<-read.table("/projects/03/genereg/projects/SOLEXA/resources/Ensembl/hg19v65/hg19v65_exons.3prime_200",as.is=T)
# exon<-read.table("/projects/03/genereg/projects/SOLEXA/resources/Ensembl/hg19v65/hg19v65_exons.5prime_200",as.is=T)
# colnames(exon)<-c("chr","start","end","strand","ID")
# index<-with(exon, order(chr, start))
# exon<-exon[index,]
# rm(index)
# # save.image("hg19v65_exons.3prime_200.Rdata")
# save.image("hg19v65_exons.5prime_200.Rdata")
 
extract<-function(start,end,methyl){
  methylextract<-methyl[methyl$pos>=start&methyl$pos<=end&methyl$rc>=5,]
  count<-nrow(methylextract)
  aveM<-mean(methylextract$m/methylextract$rc)
  y<-c(count,aveM)
  return(y)
}
vextract<-Vectorize(extract,vectorize.args=c("start","end"),SIMPLIFY=T)

extractchr<-function(chr,exon){
  print(chr)
  methyl<-read.table(paste(dirIn,"chr",chr,".m",sep=""),as.is=T,col.names=c("pos","m","rc"))
  y<-vextract(exon[exon$chr==chr,]$start,exon[exon$chr==chr,]$end,methyl)
  return(y)
}
vextractchr<-Vectorize(extractchr,vectorize.args="chr",SIMPLIFY=F)

# setwd("~/REMC/breast/bismark/lumRM066/")
# dirIn="/projects/mbilenky/REMC/breast/WGBS/bismark/lumRM066/m/"
setwd("~/REMC/breast/bismark/myoRM045/")
dirIn="/projects/mbilenky/REMC/breast/WGBS/bismark/myoRM045/m/"

load("~/REMC/breast/bismark/hg19v65_exons.3prime_200.Rdata")
# load("~/REMC/breast/bismark/hg19v65_exons.5prime_200.Rdata")

chrs<-c(as.character(c(1:22)),"X","Y")
y<-vextractchr(chrs,exon)
# y<-vextractchr("1",exon)
y<-do.call("cbind",y)
exon<-exon[is.element(exon$chr,chrs),]
exon$count<-y[1,]
exon$aveM<-y[2,]

save(exon,file="hg19v65_exons.3prime_200.Rdata")
# save(exon,file="hg19v65_exons.5prime_200.Rdata")

# save(exon,file="hg19v65_exons.3prime_200chr1.Rdata")
# save(exon,file="hg19v65_exons.5prime_200chr1.Rdata")

# chr<-as.character(exon$chr[1])
# methyl<-read.table(paste(dirIn,"chr",chr,".m",sep=""),as.is=T)
# for(i in 1:nrow(exon)){
#   if(as.character(exon$chr[i])!=chr){
#     chr<-as.character(exon$chr[i])
#     print(chr)
#     methyl<-read.table(paste(dirIn,"chr",chr,".m",sep=""),as.is=T)
#   }
#   y<-extract(as.character(exon[i,]),methyl)
#   exon$count[i]<-y[1]
#   exon$aveM[i]<-y[2]
# }

# write.table(exon,file="hg19v65_exons.3prime_2001.m",sep="\t",quote=F,row.names=F)
# write.table(exon,file="hg19v65_exons.5prime_2001.m",sep="\t",quote=F,row.names=F)

#################################################################################################################
# identify expressed exons
expressed<-function(gene,start,end){
  return(rownames(exon[exon$gene==gene&((exon$start>=start&exon$start<=end)|(exon$end>=start&exon$end<=end)),]))
}
vexpressed<-Vectorize(expressed)

setwd("~/REMC/breast/bismark/")
# load("hg19v65_exons.3prime_200.Rdata")
load("hg19v65_exons.5prime_200.Rdata")
exon$gene<-sapply(exon$ID,function(x) unlist(strsplit(as.character(x),'_'))[1])
exon$lumExpress<-0
exon$myoExpress<-0

load("~/REMC/breast/bismark/express.Rdata")
lum<-express[express$lum,7:8]
lum$chr<-unlist(strsplit(as.character(lum$ID),':'))[2*(1:nrow(lum))-1]
id<-unlist(strsplit(as.character(lum$ID),':'))[2*(1:nrow(lum))]
lum$strand<-as.numeric(unlist(strsplit(as.character(id),'<'))[2*(1:nrow(lum))])
id<-unlist(strsplit(as.character(id),'<'))[2*(1:nrow(lum))-1]
lum$start<-as.numeric(unlist(strsplit(as.character(id),'-'))[2*(1:nrow(lum))-1])
lum$end<-as.numeric(unlist(strsplit(as.character(id),'-'))[2*(1:nrow(lum))])
myo<-express[express$myo,7:8]
myo$chr<-unlist(strsplit(as.character(myo$ID),':'))[2*(1:nrow(myo))-1]
id<-unlist(strsplit(as.character(myo$ID),':'))[2*(1:nrow(myo))]
myo$strand<-as.numeric(unlist(strsplit(as.character(id),'<'))[2*(1:nrow(myo))])
id<-unlist(strsplit(as.character(id),'<'))[2*(1:nrow(myo))-1]
myo$start<-as.numeric(unlist(strsplit(as.character(id),'-'))[2*(1:nrow(myo))-1])
myo$end<-as.numeric(unlist(strsplit(as.character(id),'-'))[2*(1:nrow(myo))])

exon[unique(unlist(vexpressed(lum$gene,lum$start,lum$end))),]$lumExpress<-1
exon[unique(unlist(vexpressed(myo$gene,myo$start,myo$end))),]$myoExpress<-1

# save(exon,file="myolum_exons.3prime_200.Rdata")
save(exon,file="myolum_exons.5prime_200.Rdata")

#################################################################################################################
# identify cell type specific exons between myo vs lum
setwd("~/REMC/breast/tissue/")
load("REMC.isoform.breast.RData")
myolumgene<-myo_lum_isoform_only   # Ensembl ID of isoforms between myo vs lum
myoexon35<-myo_lum_35up_isoform[is.element(myo_lum_35up_isoform$id,myolumgene$ID),]
myoexon80<-myo_lum_80up_isoform[is.element(myo_lum_80up_isoform$id,myolumgene$ID),]
myoexon84<-myo_lum_84up_isoform[is.element(myo_lum_84up_isoform$id,myolumgene$ID),]
myoexon<-rbind(myoexon35,myoexon80,myoexon84)
myoexon<-myoexon[duplicated(myoexon$V1),]
myoexon$chr<-unlist(strsplit(as.character(myoexon$V1),':'))[2*(1:nrow(myoexon))-1]
id<-unlist(strsplit(as.character(myoexon$V1),':'))[2*(1:nrow(myoexon))]
id<-unlist(strsplit(as.character(id),'_'))[2*(1:nrow(myoexon))-1]
myoexon$strand<-as.numeric(unlist(strsplit(as.character(id),'<'))[2*(1:nrow(myoexon))])
id<-unlist(strsplit(as.character(id),'<'))[2*(1:nrow(myoexon))-1]
myoexon$start<-as.numeric(unlist(strsplit(as.character(id),'-'))[2*(1:nrow(myoexon))-1])
myoexon$end<-as.numeric(unlist(strsplit(as.character(id),'-'))[2*(1:nrow(myoexon))])
lumexon35<-myo_lum_35dn_isoform[is.element(myo_lum_35dn_isoform$id,myolumgene$ID),]
lumexon80<-myo_lum_80dn_isoform[is.element(myo_lum_80dn_isoform$id,myolumgene$ID),]
lumexon84<-myo_lum_84dn_isoform[is.element(myo_lum_84dn_isoform$id,myolumgene$ID),]
lumexon<-rbind(lumexon35,lumexon80,lumexon84)
lumexon<-lumexon[duplicated(lumexon$V1),]
lumexon$chr<-unlist(strsplit(as.character(lumexon$V1),':'))[2*(1:nrow(lumexon))-1]
id<-unlist(strsplit(as.character(lumexon$V1),':'))[2*(1:nrow(lumexon))]
id<-unlist(strsplit(as.character(id),'_'))[2*(1:nrow(lumexon))-1]
lumexon$strand<-as.numeric(unlist(strsplit(as.character(id),'<'))[2*(1:nrow(lumexon))])
id<-unlist(strsplit(as.character(id),'<'))[2*(1:nrow(lumexon))-1]
lumexon$start<-as.numeric(unlist(strsplit(as.character(id),'-'))[2*(1:nrow(lumexon))-1])
lumexon$end<-as.numeric(unlist(strsplit(as.character(id),'-'))[2*(1:nrow(lumexon))])
save(myoexon,lumexon,file="~/REMC/breast/bismark/myo_lum_exons.Rdata")

expressed<-function(gene,start,end){
  return(rownames(exon[exon$gene==gene&((exon$start>=start&exon$start<=end)|(exon$end>=start&exon$end<=end)),]))
}
vexpressed<-Vectorize(expressed)

setwd("~/REMC/breast/bismark/")
# load("myolum_exons.3prime_200.Rdata")
load("myolum_exons.5prime_200.Rdata")
load("myo_lum_exons.Rdata")
exon$myoExon<-0
exon$lumExon<-0

exon[unique(unlist(vexpressed(myoexon$id,myoexon$start,myoexon$end))),]$myoExon<-1
exon[unique(unlist(vexpressed(lumexon$id,lumexon$start,lumexon$end))),]$lumExon<-1

# save(exon,file="myolum_exons.3prime_200.Rdata")
save(exon,file="myolum_exons.5prime_200.Rdata")

#################################################################################################################
# compare methylation between myo vs lum
setwd("~/REMC/breast/bismark/")
load("myolum_exons.3prime_200.Rdata")
exon3prime<-exon
load("myolum_exons.5prime_200.Rdata")
exon5prime<-exon
exon3prime$expression<-"not expressed"
exon3prime[(exon3prime$lumExpress+exon3prime$myoExpress)>0,]$expression<-"non-specific"
exon3prime[exon3prime$lumExon==1,]$expression<-"lum-specific"
exon3prime[exon3prime$myoExon==1,]$expression<-"myo-specific"
exon3prime$expression<-as.factor(exon3prime$expression)
exon5prime$expression<-"not expressed"
exon5prime[(exon5prime$lumExpress+exon5prime$myoExpress)>0,]$expression<-"non-specific"
exon5prime[exon5prime$lumExon==1,]$expression<-"lum-specific"
exon5prime[exon5prime$myoExon==1,]$expression<-"myo-specific"
exon5prime$expression<-as.factor(exon5prime$expression)
save(exon3prime,exon5prime,file="myolum_exons_200.Rdata")

setwd("~/REMC/breast/bismark/")
load("myolum_exons_200.Rdata")
# exon3prime<-exon3prime[exon3prime$chr=="1",]
# exon5prime<-exon3prime[exon5prime$chr=="1",]
chrs<-c(as.character(c(1:22)),"X","Y")
exon3prime<-exon3prime[is.element(exon3prime$chr,chrs),]
exon5prime<-exon5prime[is.element(exon5prime$chr,chrs),]
load("~/REMC/breast/bismark/myoRM045/hg19v65_exons.3prime_200.Rdata")
exon3prime$myo_count<-exon$count
exon3prime$myo_aveM<-exon$aveM
load("~/REMC/breast/bismark/myoRM045/hg19v65_exons.5prime_200.Rdata")
exon5prime$myo_count<-exon$count
exon5prime$myo_aveM<-exon$aveM
load("~/REMC/breast/bismark/lumRM066/hg19v65_exons.3prime_200.Rdata")
exon3prime$lum_count<-exon$count
exon3prime$lum_aveM<-exon$aveM
load("~/REMC/breast/bismark/lumRM066/hg19v65_exons.5prime_200.Rdata")
exon5prime$lum_count<-exon$count
exon5prime$lum_aveM<-exon$aveM
exon3prime$myo_lum<-exon3prime$myo_ave-exon3prime$lum_ave
exon5prime$myo_lum<-exon5prime$myo_ave-exon5prime$lum_ave
save(exon3prime,exon5prime,file="exon_methylation.Rdata")

setwd("~/REMC/breast/bismark/")
load("exon_methylation.Rdata")
exon3prime<-exon3prime[exon3prime$myo_count>=5&exon3prime$lum_count>=5,]
exon5prime<-exon5prime[exon5prime$myo_count>=5&exon5prime$lum_count>=5,]

exon3Specific<-exon3prime[exon3prime$lumExon|exon3prime$myoExon,]
exon5Specific<-exon5prime[exon5prime$lumExon|exon5prime$myoExon,]

pdf("Exons_methylation.pdf")
plot(exon5prime$expression,exon5prime$myo_lum,cex=0.5,xlab=NULL,ylab="methylation:myo-lum",
     main="Exon 5prime methylation")
abline(h=0)
plot(exon3prime$expression,exon3prime$myo_lum,cex=0.5,xlab=NULL,ylab="methylation:myo-lum",
     main="Exon 3prime methylation")
abline(h=0)

# plot(exon5Specific$lumExon,exon5Specific$myo_lum,xlab="0:myo-specific;1:lum-specific",ylab="methylation:myo-lum",
#      main="Cell type specific exon 5prime methylation")
# abline(h=0)
# plot(exon3Specific$lumExon,exon3Specific$myo_lum,xlab="0:myo-specific;1:lum-specific",ylab="methylation:myo-lum",
#      main="Cell type specific exon 3prime methylation")
# abline(h=0)

smoothScatter(exon5prime[exon5prime$expression=="non-specific",]$lum_aveM,exon5prime[exon5prime$expression=="non-specific",]$myo_aveM,
              xlab="lum",ylab="myo",main="5prime methylation level")
points(exon5Specific$lum_aveM,exon5Specific$myo_aveM,col=exon5Specific$expression,pch=19,cex=0.5)
abline(a=0,b=1)
legend("bottomright",c("lum-specific","myo-specific","non-specific"),col=c("black","red","blue"),pch=19)
smoothScatter(exon3prime[exon3prime$expression=="non-specific",]$lum_aveM,exon3prime[exon3prime$expression=="non-specific",]$myo_aveM,
              xlab="lum",ylab="myo",main="3prime methylation level")
points(exon3Specific$lum_aveM,exon3Specific$myo_aveM,col=exon3Specific$expression,pch=19,cex=0.5)
abline(a=0,b=1)
legend("bottomright",c("lum-specific","myo-specific","non-specific"),col=c("black","red","blue"),pch=19)
dev.off()

mean(exon5prime[exon5prime$lumExon==1,]$myo_lum)
t.test(exon5prime[exon5prime$lumExon==1,]$myo_lum)
mean(exon5prime[exon5prime$myoExon==1,]$myo_lum)
t.test(exon5prime[exon5prime$myoExon==1,]$myo_lum)
mean(exon3prime[exon3prime$lumExon==1,]$myo_lum)
t.test(exon3prime[exon3prime$lumExon==1,]$myo_lum)
mean(exon3prime[exon3prime$myoExon==1,]$myo_lum)
t.test(exon3prime[exon3prime$myoExon==1,]$myo_lum)
t.test(exon5prime[exon5prime$lumExon==1,]$lum_aveM,exon5prime[exon5prime$lumExon==1,]$myo_aveM)
t.test(exon5prime[exon5prime$myoExon==1,]$lum_aveM,exon5prime[exon5prime$myoExon==1,]$myo_aveM)
t.test(exon3prime[exon3prime$lumExon==1,]$lum_aveM,exon3prime[exon3prime$lumExon==1,]$myo_aveM)
t.test(exon3prime[exon3prime$myoExon==1,]$lum_aveM,exon3prime[exon3prime$myoExon==1,]$myo_aveM)

# heatmaps 
library(gplots)
#heatmap(as.matrix(cbind(exon3prime$myo_aveM,exon3prime$lum_aveM)))
exon3Express<-as.matrix(data.frame(myo=exon3prime[(exon3prime$lumExpress+exon3prime$myoExpress)>0,]$myo_aveM,
                                   lum=exon3prime[(exon3prime$lumExpress+exon3prime$myoExpress)>0,]$lum_aveM))
#heatmap(exon3Express)
exon5Express<-as.matrix(data.frame(myo=exon5prime[(exon5prime$lumExpress+exon5prime$myoExpress)>0,]$myo_aveM,
                                   lum=exon5prime[(exon5prime$lumExpress+exon5prime$myoExpress)>0,]$lum_aveM))
#heatmap(exon5Express)

pdf("heatmap_CellTypeSpecificExons_3prime.pdf",width=12,height=16)
exon3CellSpecific<-as.matrix(data.frame(myo=exon3Specific$myo_aveM,lum=exon3Specific$lum_aveM))
note3<-as.matrix(data.frame(myo=as.character(exon3Specific$myoExon),lum=as.character(exon3Specific$lumExon)))
par(mar=c(6,6,6,10))
heatmap.2(exon3CellSpecific,Colv=NULL,trace="none",density="none",keysize=0.5,cellnote=note3,notecex=0.5,notecol=1,cexRow=0.5,cexCol=1,
          labRow=exon3prime[exon3prime$lumExon|exon3prime$myoExon,]$ID,
          main="Cell type specific exons 3prime methylation")
dev.off()
exon5Specific<-exon5prime[exon5prime$lumExon|exon5prime$myoExon,]
pdf("heatmap_CellTypeSpecificExons_5prime.pdf",width=12,height=16)
exon5CellSpecific<-as.matrix(data.frame(myo=exon5Specific$myo_aveM,lum=exon5Specific$lum_aveM))
note5<-as.matrix(data.frame(myo=as.character(exon5Specific$myoExon),lum=as.character(exon5Specific$lumExon)))
par(mar=c(6,6,6,10))
heatmap.2(exon5CellSpecific,Colv=NULL,trace="none",density="none",keysize=0.5,cellnote=note5,notecex=0.5,notecol=1,cexRow=0.5,cexCol=1,
          labRow=exon5prime[exon5prime$lumExon|exon5prime$myoExon,]$ID,
          main="Cell type specific exons 5prime methylation")
dev.off()


