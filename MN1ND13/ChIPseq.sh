#!/bin/sh

# enhancer homer motif annotation
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn=/projects/epigenomics2/PING/ChIPAnalysis/
/home/lli/bin/homer/bin/annotatePeaks.pl $dirIn/Unique_H3K4me1_DN hg19 -m $dirIn/H3K4me1_Motif/DN/knownResults/known2.motif > $dirIn/H3K4me1_Motif/DN/PU1.annotate
/home/lli/bin/homer/bin/annotatePeaks.pl $dirIn/Unique_H3K4me1_DP hg19 -m $dirIn/H3K4me1_Motif/DP/knownResults/known16.motif > $dirIn/H3K4me1_Motif/DP/RUNX1.annotate
less $dirIn/H3K4me1_Motif/DN/PU1.annotate | awk -F'\t' 'NR>1{if($22!=""){gsub("chr", ""); print $2"\t"$3"\t"$4"\t"$1}}' > $dirIn/H3K4me1_Motif/DN/PU1.annotate.only.bed
less $dirIn/H3K4me1_Motif/DP/RUNX1.annotate | awk -F'\t' 'NR>1{if($22!=""){gsub("chr", ""); print $2"\t"$3"\t"$4"\t"$1}}' > $dirIn/H3K4me1_Motif/DP/RUNX1.annotate.only.bed
## 5mC at enhancer
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dir5mC=/projects/epigenomics2/PING/DMR/
$BEDTOOLS/intersectBed -a $dir5mC/A04_MN1ND13_DP.bed.CpG.txt.gz.combine.5mC.CpG -b $dirIn/H3K4me1_Motif/DN/PU1.annotate.only.bed -wa -wb | awk '{c[$10]=c[$10]+$5; t[$10]=t[$10]+$4; chr[$10]=$7; start[$10]=$8; end[$10]=$9}END{for(i in chr){if(c[i]+t[i]>0){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"c[i]/(c[i]+t[i])}}}' > $dirIn/H3K4me1_Motif/DN/PU1.annotate.only.DP5mC.bed
$BEDTOOLS/intersectBed -a $dir5mC/A05_MN1ND13_DN.bed.CpG.txt.gz.combine.5mC.CpG -b $dirIn/H3K4me1_Motif/DN/PU1.annotate.only.DP5mC.bed -wa -wb | awk '{c[$10]=c[$10]+$5; t[$10]=t[$10]+$4; chr[$10]=$7; start[$10]=$8; end[$10]=$9; dp[$10]=$11}END{for(i in chr){if(c[i]+t[i]>0){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"dp[i]"\t"c[i]/(c[i]+t[i])"\t"dp[i]-c[i]/(c[i]+t[i])}}}' | sort -k1,1 -k 2,2n -T /projects/epigenomics/temp/ > $dirIn/H3K4me1_Motif/DN/PU1.annotate.only.5mC.bed
$BEDTOOLS/intersectBed -a $dir5mC/A04_MN1ND13_DP.bed.CpG.txt.gz.combine.5mC.CpG -b $dirIn/H3K4me1_Motif/DP/RUNX1.annotate.only.bed -wa -wb | awk '{c[$10]=c[$10]+$5; t[$10]=t[$10]+$4; chr[$10]=$7; start[$10]=$8; end[$10]=$9}END{for(i in chr){if(c[i]+t[i]>0){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"c[i]/(c[i]+t[i])}}}' > $dirIn/H3K4me1_Motif/DP/RUNX1.annotate.only.DP5mC.bed
$BEDTOOLS/intersectBed -a $dir5mC/A05_MN1ND13_DN.bed.CpG.txt.gz.combine.5mC.CpG -b $dirIn/H3K4me1_Motif/DP/RUNX1.annotate.only.DP5mC.bed -wa -wb | awk '{c[$10]=c[$10]+$5; t[$10]=t[$10]+$4; chr[$10]=$7; start[$10]=$8; end[$10]=$9; dp[$10]=$11}END{for(i in chr){if(c[i]+t[i]>0){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"dp[i]"\t"c[i]/(c[i]+t[i])"\t"dp[i]-c[i]/(c[i]+t[i])}}}' | sort -k1,1 -k 2,2n -T /projects/epigenomics/temp/ > $dirIn/H3K4me1_Motif/DP/RUNX1.annotate.only.5mC.bed
