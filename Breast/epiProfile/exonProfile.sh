#!/bin/sh

cd /home/lli/REMC/epiProfile

############################################################################################
# WGBS @ exon boundaries
#bed=/projects/edcc_prj2/bs-seq/a22478/A22478.Cmethyl.cons.bed.gz; name=lumRM066_novoalign;

#bed=/projects/edcc_prj2/bs-seq/a22478/bismark/A22478_4_lanes_dupsFlagged.fractional.bedGraph.gz; name=lumRM066_WGBS;
#bed=/projects/edcc_prj2/bs-seq/a18473/A18473_5_lanes_dupsFlagged.bam.fractional.bedGraph.gz; name=myoRM045_WGBS;
bed=/projects/mbilenky/REMC/H1/hg19/WGBS/h1.SRX006789.SRX006782.SRX006241.SRX006240.SRX006239.dup.sorted.fractional.bedGraph.gz; name=H1_bismark; 
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique

out=/home/lli/REMC/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/REMC/epiProfile/exons5p_200/; reg=$exons5p

/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromBEDCalculator.jar -b $bed -r $reg -o $out -s hg19 -n $name -bin 20 -t Y -trueBED

############################################################################################
# MeDIP @ exon boundaries
wig=/projects/mbilenky/REMC/breast/hg19/MeDIP/wigs/HS1393.MeDIP.RM035.lum.bam.q5.F1028.SET_153.wig.gz; name=lumRM035_MeDIP;
#wig=/projects/mbilenky/REMC/breast/hg19/MeDIP/wigs/HS1394.MeDIP.RM035.myo.bam.q5.F1028.SET_120.wig.gz; name=myoRM035_MeDIP;
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique

#out=/home/lli/REMC/epiProfile/exons3p_200/; reg=$exons3p
out=/home/lli/REMC/epiProfile/exons5p_200/; reg=$exons5p

/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

############################################################################################
# H3K4me3 @ exon boundaries: no lum libraries for this mark 
wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2760.H3K4me3.myo.RM080.bam.q5.F1028.SET_154.wig.gz; name=myoRM080_H3K4me3;
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique

out=/home/lli/REMC/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/REMC/epiProfile/exons5p_200/; reg=$exons5p

/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

############################################################################################
# H3K4me1 @ exon boundaries
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2791.H3K4me1.lum.RM080.bam.q5.F1028.SET_124.wig.gz; name=lumRM080_H3K4me1;
wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2759.H3K4me1.myo.RM080.bam.q5.F1028.SET_159.wig.gz; name=myoRM080_H3K4me1;
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique

out=/home/lli/REMC/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/REMC/epiProfile/exons5p_200/; reg=$exons5p

/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

############################################################################################
# H3K9me3 @ exon boundaries
wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2793.H3K9me3.lum.RM080.bam.q5.F1028.SET_144.wig.gz; name=lumRM080_H3K9me3;
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2761.H3K9me3.myo.RM080.bam.q5.F1028.SET_158.wig.gz; name=myoRM080_H3K9me3;
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique

out=/home/lli/REMC/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/REMC/epiProfile/exons5p_200/; reg=$exons5p

/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

############################################################################################
# H3K27me3 @ exon boundaries
wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2794.H3K27me3.lum.RM080.bam.q5.F1028.SET_144.wig.gz; name=lumRM080_H3K27me3;
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2762.H3K27me3.myo.RM080.bam.q5.F1028.SET_144.wig.gz; name=myoRM080_H3K27me3;
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique

out=/home/lli/REMC/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/REMC/epiProfile/exons5p_200/; reg=$exons5p

/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

############################################################################################
# H3K36me3 @ exon boundaries
wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2795.H3K36me3.lum.RM080.bam.q5.F1028.SET_174.wig.gz; name=lumRM080_H3K36me3;
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2763.H3K36me3.myo.RM080.bam.q5.F1028.SET_139.wig.gz; name=myoRM080_H3K36me3;
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique

out=/home/lli/REMC/epiProfile/exons3p_200/; reg=$exons3p
#out=/home/lli/REMC/epiProfile/exons5p_200/; reg=$exons5p

/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 

############################################################################################
# WGBS along exons
bed=/projects/edcc_prj2/bs-seq/a22478/bismark/A22478_4_lanes_dupsFlagged.fractional.bedGraph.gz; name=lumRM066_bismark;
#bed=/projects/edcc_prj2/bs-seq/a18473/A18473_5_lanes_dupsFlagged.bam.fractional.bedGraph.gz; name=myoRM045_bismark;
exons=/home/lli/hg19/hg19v65_exons_for_genes
out=/home/lli/REMC/epiProfile/exons/; reg=$exons

/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromBEDCalculator.jar -b $bed -r $reg -o $out -s hg19 -n $name -bin 20 -t Y -trueBED

############################################################################################
# H3K36me3 signal @ exons
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2795.H3K36me3.lum.RM080.bam.q5.F1028.SET_174.wig.gz; name=lumRM080_H3K36me3;
#wig=/projects/mbilenky/REMC/breast/hg19/marks/wigs/HS2763.H3K36me3.myo.RM080.bam.q5.F1028.SET_139.wig.gz; name=myoRM080_H3K36me3;
wig=/projects/mbilenky/REMC/H1/UCSF_UBC/wigs_175/HS1032_HS1347.H3K36me3.H1_Rep12.q5.F1028.SET_175.wig.gz; name=H1_H3K36me3;
exons=/home/lli/hg19/hg19v65_exons_for_genes
out=/home/lli/REMC/epiProfile/exons/; reg=$exons
chr=/projects/epigenomics/resources/UCSC_chr/hg19.chrom.sizes

/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $out -s $chr -n $name

############################################################################################
# H3K36me3 signal @ exons
#bed=/projects/mbilenky/REMC/breast/hg19/marks/FindER/H3K36me3/HS2795.H3K36me3.FindER_scan.bin20.step20.signal.bedGraph.gz; name=lumRM080_H3K36me3_FindER;
bed=/projects/mbilenky/REMC/breast/hg19/marks/FindER/H3K36me3/HS2763.H3K36me3.FindER_scan.bin20.step20.signal.bedGraph.gz; name=myoRM080_H3K36me3_FindER;
exons=/home/lli/hg19/hg19v65_exons_for_genes
out=/home/lli/REMC/epiProfile/exons/; reg=$exons
chr=/projects/epigenomics/resources/UCSC_chr/hg19.chrom.sizes

/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromBEDCalculator.jar -b $bed -r $reg -o $out -s $chr -n $name -trueBED

