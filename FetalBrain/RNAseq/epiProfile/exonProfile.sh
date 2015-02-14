#!/bin/sh

cd ~/FetalBrain/RNAseq/epiProfile/
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique
exons=/home/lli/hg19/hg19v65_exons_for_genes
chr=/projects/epigenomics/resources/UCSC_chr/hg19.chrom.sizes

############################################################################################
# WGBS @ exon boundaries
#bed=/projects/edcc_prj2/bs-seq/a22475/A22475_4_lanes_dupsFlagged.fractional.bedGraph.gz; name=cortexHuFNSC02_WGBS;
#bed=/projects/edcc_prj2/bs-seq/joc163/A17784-3_A13819-1_dupsFlagged.fractional.bedGraph.gz; name=GEHuFNSC02_WGBS; 
bed=/projects/edcc_prj2/bs-seq/a22476/A22476_4_lanes_dupsFlagged.fractional.bedGraph.gz; name=GEHuFNSC04_WGBS;
#bed=/projects/edcc_prj2/bs-seq/a22477/A22477_4_lanes_dupsFlagged.fractional.bedGraph.gz; name=cortexHuFNSC04_WGBS; 
out=/home/lli/FetalBrain/RNAseq/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/FetalBrain/RNAseq/epiProfile/exons5p_200/; reg=$exons5p
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromBEDCalculator.jar -b $bed -r $reg -o $out -s hg19 -n $name -bin 20 -t Y -trueBED

############################################################################################
# MeDIP @ exon boundaries
#wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2779.bam.q5.F1028.SET_174.wig.gz; name=cortexHuFNSC02_MeDIP;
#wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2781.bam.q5.F1028.SET_166.wig.gz; name=GEHuFNSC02_MeDIP;
wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2775.bam.q5.F1028.SET_174.wig.gz; name=cortexHuFNSC01_MeDIP;
#wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2777.bam.q5.F1028.SET_157.wig.gz; name=GEHuFNSC01_MeDIP;
#wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2788.bam.q5.F1028.SET_174.wig.gz; name=brainHuFNSC01_MeDIP;
#wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2790.bam.q5.F1028.SET_174.wig.gz; name=brainHuFNSC02_MeDIP;
out=/home/lli/FetalBrain/RNAseq/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/FetalBrain/RNAseq/epiProfile/exons5p_200/; reg=$exons5p
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

###########################################################################################
# H3K36me3 signal @ exons
#wig=/home/lli/FetalBrain/HisMod/wigs/A03285.bam.q5.F1028.SET_238.wig.gz; name=cortexHuFNSC02_H3K36me3;
#wig=/home/lli/FetalBrain/HisMod/wigs/A03481.bam.q5.F1028.SET_215.wig.gz; name=GEHuFNSC02_H3K36me3;
wig=/home/lli/FetalBrain/HisMod/wigs/A03273.bam.q5.F1028.SET_219.wig.gz; name=cortexHuFNSC01_H3K36me3;
#wig=/home/lli/FetalBrain/HisMod/wigs/A03279.bam.q5.F1028.SET_229.wig.gz; name=GEHuFNSC01_H3K36me3;
#wig=/home/lli/FetalBrain/HisMod/wigs/A03489.bam.q5.F1028.SET_175.wig.gz; name=brainHuFNSC01_H3K36me3;
#wig=/home/lli/FetalBrain/HisMod/wigs/A03497.bam.q5.F1028.SET_233.wig.gz; name=brainHuFNSC02_H3K36me3;
out=/home/lli/FetalBrain/RNAseq/epiProfile/exons/; reg=$exons
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $out -s $chr -n $name

#############################################################################################
#############################################################################################
## H3K4me1 @ exon boundaries
##wig=/home/lli/FetalBrain/HisMod/wigs/A03281.bam.q5.F1028.SET_240.wig.gz; name=cortexHuFNSC02_H3K4me1;
#wig=/home/lli/FetalBrain/HisMod/wigs/A03477.bam.q5.F1028.SET_232.wig.gz; name=GEHuFNSC02_H3K4me1;
##out=/home/lli/FetalBrain/RNAseq/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/FetalBrain/RNAseq/epiProfile/exons5p_200/; reg=$exons5p
#/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 
#
#############################################################################################
## H3K4me3 @ exon boundaries
##wig=/home/lli/FetalBrain/HisMod/wigs/A03282.bam.q5.F1028.SET_207.wig.gz; name=cortexHuFNSC02_H3K4me3;
#wig=/home/lli/FetalBrain/HisMod/wigs/A03478.bam.q5.F1028.SET_169.wig.gz; name=GEHuFNSC02_H3K4me3;
##out=/home/lli/FetalBrain/RNAseq/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/FetalBrain/RNAseq/epiProfile/exons5p_200/; reg=$exons5p
#/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 
#
#############################################################################################
## H3K9me3 @ exon boundaries
##wig=/home/lli/FetalBrain/HisMod/wigs/A03283.bam.q5.F1028.SET_224.wig.gz; name=cortexHuFNSC02_H3K9me3;
#wig=/home/lli/FetalBrain/HisMod/wigs/A03479.bam.q5.F1028.SET_206.wig.gz; name=GEHuFNSC02_H3K9me3;
##out=/home/lli/FetalBrain/RNAseq/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/FetalBrain/RNAseq/epiProfile/exons5p_200/; reg=$exons5p
#/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 
#
#############################################################################################
## H3K27me3 @ exon boundaries
#wig=/home/lli/FetalBrain/HisMod/wigs/A03284.bam.q5.F1028.SET_194.wig.gz; name=cortexHuFNSC02_H3K27me3;
##wig=/home/lli/FetalBrain/HisMod/wigs/A03480.bam.q5.F1028.SET_211.wig.gz; name=GEHuFNSC02_H3K27me3;
##out=/home/lli/FetalBrain/RNAseq/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/FetalBrain/RNAseq/epiProfile/exons5p_200/; reg=$exons5p
#/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

