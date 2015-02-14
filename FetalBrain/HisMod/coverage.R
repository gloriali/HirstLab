setwd("/Users/epigenomics_lab_02/FetalBrain/chipseq/coverage/")
files<-list.files(pattern="\\.coverage.txt$")
pdf("coverage.pdf")
for (f in files){
  file<-read.table(f)
  name<-unlist(strsplit(f,"[.]"))[1]
  plot(file$V1[file$V2!=0],log2(file$V2[file$V2!=0]),main=name,xlab="coverage",ylab="log2 bases",pch=20,cex=0.5)
}
dev.off()
