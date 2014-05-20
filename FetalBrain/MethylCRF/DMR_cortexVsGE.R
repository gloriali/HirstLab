# differential methylation analysis between neurosphere cortex derived and neurosphere GE derived
# datasets available: MeDIP/MRE (methylCRF) from HuFNSC01 and HuFNSC02, treat the two individuals as biological replicates

getwd()
# read binned methylCRF data (methylation level for 1kb bins and CpG counts per bin)
crf <- read.csv(file="~/FetalBrain/MeDIPMRE/CRF/crf.csv", as.is=T) # sum & #CpGs for each bin
# average methylation level 
crf$brain01 <- crf$brain01/crf$count
crf$brain02 <- crf$brain02/crf$count
crf$cortex01 <- crf$cortex01/crf$count
crf$cortex02 <- crf$cortex02/crf$count
crf$ge01 <- crf$ge01/crf$count
crf$ge02 <- crf$ge02/crf$count
crf$ID <- paste0(crf$chr,":",crf$start,"-",crf$end)
rownames(crf) <- crf$ID
nrow(crf) # Total No. of bins

# identify DMRs between cortex and GE: fold change>2  &  at least one sample CRF>0.5
# calculate log2(fold change)
crf$cortexVsGE01<-log2(crf$cortex01/crf$ge01)
crf$cortexVsGE02<-log2(crf$cortex02/crf$ge02)
# DMRs
DMR_cortexGE01 <- crf[(abs(crf$cortexVsGE01)>1) & ((crf$cortex01>0.5) | (crf$ge01>0.5)),]
nrow(DMR_cortexGE01) # No. of DMRs - cortex01 vs GE01
DMR_cortexGE02 <- crf[(abs(crf$cortexVsGE02)>1) & ((crf$cortex02>0.5) | (crf$ge02>0.5)),]
nrow(DMR_cortexGE02) # No. of DMRs - cortex02 vs GE02
DMR_cortexGE <- crf[intersect(DMR_cortexGE01$ID, DMR_cortexGE02$ID),]
nrow(DMR_cortexGE) # No. of DMRs - cortex vs GE in both individuals
write.csv(DMR_cortexGE, file="DMR_cortexGE.csv", row.names=F, quote=F)
write.csv(DMR_cortexGE01, file="DMR_cortexGE01.csv", row.names=F, quote=F)
write.csv(DMR_cortexGE02, file="DMR_cortexGE02.csv", row.names=F, quote=F)
write.table(cbind(DMR_cortexGE$chr, DMR_cortexGE$start, DMR_cortexGE$end, DMR_cortexGE$ID), 
            file="DMR_cortexGE.bed", quote=F, sep="\t", col.names=F, row.names=F) # input for GREAT
write.table(cbind(DMR_cortexGE01$chr, DMR_cortexGE01$start, DMR_cortexGE01$end, DMR_cortexGE01$ID), 
            file="DMR_cortexGE01.bed", quote=F, sep="\t", col.names=F, row.names=F) # input for GREAT
write.table(cbind(DMR_cortexGE02$chr, DMR_cortexGE02$start, DMR_cortexGE02$end, DMR_cortexGE02$ID), 
            file="DMR_cortexGE02.bed", quote=F, sep="\t", col.names=F, row.names=F) # input for GREAT
# Fisher's exact test for overlapping
fisher.test(matrix(c(nrow(DMR_cortexGE),
                     nrow(DMR_cortexGE01)-nrow(DMR_cortexGE),
                     nrow(DMR_cortexGE02)-nrow(DMR_cortexGE),
                     nrow(crf)-nrow(DMR_cortexGE01)-nrow(DMR_cortexGE02)+nrow(DMR_cortexGE)),ncol=2))

par(mfrow=c(2,2))
xrange<-range(DMR_cortexGE01$cortexVsGE01, DMR_cortexGE02$cortexVsGE02)
yrange<-range(density(DMR_cortexGE01$cortexVsGE01)$y, density(DMR_cortexGE02$cortexVsGE02)$y)
plot(xrange, yrange, type="n", main="DMRs between cortex and GE", xlab="log2(MethylCRF-cortex/MethylCRF-GE)", ylab="density")
lines(density(DMR_cortexGE01$cortexVsGE01), col=1, lty=1, lwd=3)
lines(density(DMR_cortexGE02$cortexVsGE02), col=2, lty=1, lwd=3)
legend("topright", c("cortex01 vs GE01", "cortex02 vs GE02"), col=c(1:2), lty=1, lwd=5, cex=0.5)

xrange<-range(DMR_cortexGE$cortexVsGE01, DMR_cortexGE$cortexVsGE02)
yrange<-range(density(DMR_cortexGE$cortexVsGE01)$y, density(DMR_cortexGE$cortexVsGE02)$y)
plot(xrange, yrange, type="n", main="common DMRs in both individuals between cortex and GE", xlab="log2(MethylCRF-cortex/MethylCRF-GE)", ylab="density")
lines(density(DMR_cortexGE$cortexVsGE01), col=1, lty=1, lwd=3)
lines(density(DMR_cortexGE$cortexVsGE02), col=2, lty=1, lwd=3)
legend("topright", c("cortex01 vs GE01","cortex02 vs GE02"), col=c(1:2), lty=1, lwd=5, cex=0.5)

plot(ecdf(crf$count), xlim=c(0,30), main="Ecdf of #CpGs per bin")
lines(ecdf(DMR_cortexGE$count), xlim=c(0,30), col=2)
abline(v=5)
legend("bottomright", c("all bins", "DMRs - cortex vs GE"), col=c(1:2), lty=1, lwd=5, cex=0.5)


# Annotating DMRs with promoter regions
tss <- read.table("~/hg19/hg19v65_genes_TSS_1500.txt", head=F, as.is=T)
write.table(cbind(tss$V1, tss$V2, tss$V3, tss$V5), file="~/hg19/hg19v65_genes_TSS_1500.bed", row.names=F, col.names=F, quote=F, sep="\t")

# /Applications/bedtools-2.17.0/bin/intersectBed -a ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE.bed -b ~/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_promoter_cortexGE.txt
# /Applications/bedtools-2.17.0/bin/intersectBed -a ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE01.bed -b ~/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_promoter_cortexGE01.txt
# /Applications/bedtools-2.17.0/bin/intersectBed -a ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE02.bed -b ~/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_promoter_cortexGE02.txt

cortexGE <- read.delim("DMR_promoter_cortexGE.txt", head=F, as.is=T); nrow(cortexGE)  # No. of promoter DMRs
cortexGE$Ensembl <- gsub("_[a-zA-Z_]+", "", cortexGE$V8)
cortexGE01 <- read.delim("DMR_promoter_cortexGE01.txt", head=F, as.is=T); nrow(cortexGE01)
cortexGE01$Ensembl <- gsub("_[a-zA-Z_]+", "", cortexGE01$V8)
cortexGE02 <- read.delim("DMR_promoter_cortexGE02.txt", head=F, as.is=T); nrow(cortexGE02)
cortexGE02$Ensembl <- gsub("_[a-zA-Z_]+", "", cortexGE02$V8)

# calculate enrichment
# total length of overlapping with promoter regions: DMR overlaps with TSS+-1.5kb for at least 1bp
sum((tss$V3+500) - (tss$V2-500))
# total lenght of all genome
genome <- 3153508730
# expected percentage of DMR overlapping promoters by chance
e <- sum((tss$V3+500) - (tss$V2-500))/genome
# enrichment: observed/expected by chance
(nrow(cortexGE)/1215)/e    # shared DMRs
(nrow(cortexGE01)/3638)/e  # DMR in HuFNSC01
(nrow(cortexGE02)/4712)/e  # DMR in HuFNSC01

# separate protein_coding and non_coding
pc_cortexGE <- cortexGE[grepl("protein_coding",cortexGE$V8),]; nrow(pc_cortexGE)
nc_cortexGE <- cortexGE[!grepl("protein_coding",cortexGE$V8),]; nrow(nc_cortexGE)
pc_cortexGE01 <- cortexGE01[grepl("protein_coding",cortexGE01$V8),]; nrow(pc_cortexGE01)
nc_cortexGE01 <- cortexGE01[!grepl("protein_coding",cortexGE01$V8),]; nrow(nc_cortexGE01)
pc_cortexGE02 <- cortexGE02[grepl("protein_coding",cortexGE02$V8),]; nrow(pc_cortexGE02)
nc_cortexGE02 <- cortexGE02[!grepl("protein_coding",cortexGE02$V8),]; nrow(nc_cortexGE02)

write.csv(cortexGE, file="DMR_promoter_cortexGE.csv", quote=F, row.names=F)
write.table(cortexGE[,1:4], file="DMR_promoter_cortexGE.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.csv(pc_cortexGE, file="DMR_pcpromoter_cortexGE.csv", quote=F, row.names=F)
write.table(pc_cortexGE[,1:4], file="DMR_pcpromoter_cortexGE.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.csv(nc_cortexGE, file="DMR_ncpromoter_cortexGE.csv", quote=F, row.names=F)
write.table(nc_cortexGE[,1:4], file="DMR_ncpromoter_cortexGE.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.csv(pc_cortexGE01, file="DMR_pcpromoter_cortexGE01.csv", quote=F, row.names=F)
write.table(pc_cortexGE01[,1:4], file="DMR_pcpromoter_cortexGE01.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.csv(nc_cortexGE01, file="DMR_ncpromoter_cortexGE01.csv", quote=F, row.names=F)
write.table(nc_cortexGE01[,1:4], file="DMR_ncpromoter_cortexGE01.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.csv(pc_cortexGE02, file="DMR_pcpromoter_cortexGE02.csv", quote=F, row.names=F)
write.table(pc_cortexGE02[,1:4], file="DMR_pcpromoter_cortexGE02.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.csv(nc_cortexGE02, file="DMR_ncpromoter_cortexGE02.csv", quote=F, row.names=F)
write.table(nc_cortexGE02[,1:4], file="DMR_ncpromoter_cortexGE02.bed", sep="\t", quote=F, row.names=F, col.names=F)

