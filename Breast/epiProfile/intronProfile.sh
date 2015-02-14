#!/bin/sh

cd /home/lli/REMC/epiProfile/IR/
#introns3p=/projects/epigenomics/acarles/intron/breast/4GloriaAnalysis/hg19v65_introns_for_genes.Id4Gloria.bothStrands.3prime
#introns5p=/projects/epigenomics/acarles/intron/breast/4GloriaAnalysis/hg19v65_introns_for_genes.Id4Gloria.bothStrands.5prime
#less $introns3p | awk '{id=$4; gsub("ENSG[0-9chrXYM:_-]+<", "", $4); print $1"\t"$2"\t"$3"\t"$4"\t"id}' > /home/lli/hg19/hg19v65_introns_for_genes.3prime_200
#less $introns5p | awk '{id=$4; gsub("ENSG[0-9chrXYM:_-]+<", "", $4); print $1"\t"$2"\t"$3"\t"$4"\t"id}' > /home/lli/hg19/hg19v65_introns_for_genes.5prime_200
introns3p=/home/lli/hg19/hg19v65_introns_for_genes.3prime_200
introns5p=/home/lli/hg19/hg19v65_introns_for_genes.5prime_200
#introns=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_introns_for_genes
#less $introns | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"_chr"$1":"$2"-"$3"<"$5}' > /home/lli/hg19/hg19v65_introns_for_genes
introns=/home/lli/hg19/hg19v65_introns_for_genes

#############################################################################################
## WGBS @ intron boundaries
##bed=/projects/edcc_prj2/bs-seq/a22478/bismark/A22478_4_lanes_dupsFlagged.fractional.bedGraph.gz; name=lumRM066_bismark;
#bed=/projects/edcc_prj2/bs-seq/a18473/A18473_5_lanes_dupsFlagged.bam.fractional.bedGraph.gz; name=myoRM045_bismark;
#out=/home/lli/REMC/epiProfile/IR/introns3p_200/; reg=$introns3p
##out=/home/lli/REMC/epiProfile/IR/introns5p_200/; reg=$introns5p
#/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromBEDCalculator.jar -b $bed -r $reg -o $out -s hg19 -n $name -bin 20 -t Y -trueBED

############################################################################################
# MeDIP @ intron boundaries
#wig=/projects/mbilenky/REMC/breast/hg19/MeDIP/wigs/HS1393.MeDIP.RM035.lum.bam.q5.F1028.SET_153.wig.gz; name=lumRM035_MeDIP;
wig=/projects/mbilenky/REMC/breast/hg19/MeDIP/wigs/HS1394.MeDIP.RM035.myo.bam.q5.F1028.SET_120.wig.gz; name=myoRM035_MeDIP;
#out=/home/lli/REMC/epiProfile/IR/introns3p_200/; reg=$introns3p
out=/home/lli/REMC/epiProfile/IR/introns5p_200/; reg=$introns5p
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

#############################################################################################
## H3K4me3 @ intron boundaries: no lum libraries for this mark 
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2760.H3K4me3.myo.RM080.bam.q5.F1028.SET_154.wig.gz; name=myoRM080_H3K4me3;
##out=/home/lli/REMC/epiProfile/IR/introns3p_200/; reg=$introns3p
#out=/home/lli/REMC/epiProfile/IR/introns5p_200/; reg=$introns5p
#/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

############################################################################################
# H3K4me1 @ intron boundaries
wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2791.H3K4me1.lum.RM080.bam.q5.F1028.SET_124.wig.gz; name=lumRM080_H3K4me1;
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2759.H3K4me1.myo.RM080.bam.q5.F1028.SET_159.wig.gz; name=myoRM080_H3K4me1;
out=/home/lli/REMC/epiProfile/IR/introns3p_200/; reg=$introns3p
#out=/home/lli/REMC/epiProfile/IR/introns5p_200/; reg=$introns5p
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

############################################################################################
# H3K9me3 @ intron boundaries
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2793.H3K9me3.lum.RM080.bam.q5.F1028.SET_144.wig.gz; name=lumRM080_H3K9me3;
wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2761.H3K9me3.myo.RM080.bam.q5.F1028.SET_158.wig.gz; name=myoRM080_H3K9me3;
out=/home/lli/REMC/epiProfile/IR/introns3p_200/; reg=$introns3p
#out=/home/lli/REMC/epiProfile/IR/introns5p_200/; reg=$introns5p
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

############################################################################################
# H3K27me3 @ intron boundaries
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2794.H3K27me3.lum.RM080.bam.q5.F1028.SET_144.wig.gz; name=lumRM080_H3K27me3;
wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2762.H3K27me3.myo.RM080.bam.q5.F1028.SET_144.wig.gz; name=myoRM080_H3K27me3;
out=/home/lli/REMC/epiProfile/IR/introns3p_200/; reg=$introns3p
#out=/home/lli/REMC/epiProfile/IR/introns5p_200/; reg=$introns5p
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

#############################################################################################
## H3K36me3 signal @ introns
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2795.H3K36me3.lum.RM080.bam.q5.F1028.SET_174.wig.gz; name=lumRM080_H3K36me3;
##wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2763.H3K36me3.myo.RM080.bam.q5.F1028.SET_139.wig.gz; name=myoRM080_H3K36me3;
#out=/home/lli/REMC/epiProfile/IR/introns/; reg=$introns
#chr=/projects/epigenomics/resources/UCSC_chr/hg19.chrom.sizes
#/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $out -s $chr -n $name
