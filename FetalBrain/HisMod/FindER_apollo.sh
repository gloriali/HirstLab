#!/bin/sh
#$ -S /bin/sh
#$ -m e
#$ -M lli@bcgsc.ca
#$ -N FindER_FetalBrain
#$ -V

chr=/home/mbilenky/UCSC_chr/hg19.chrom.sizes
map=/projects/epigenomics/resources/UCSC_hg19/mappability/mappability_lt_1
 
inputWig=/home/lli/FetalBrain/HisMod/wigs/A03486.bam.q5.F1028.SET_195.wig.gz
out=/home/lli/FetalBrain/HisMod/FindER/
name=A03486_brain01_H3K4me3
java -jar -Xmx30G /home/mbilenky/bin/Solexa_Java/FindER.jar -chr $chr -i $inputWig -o $out -n $name -v -m $map -print

