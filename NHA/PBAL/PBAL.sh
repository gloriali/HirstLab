#!/bin/sh

bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/bam/
for bam in /projects/analysis/analysis30/PX0682/HCTVVCCXY_5/PX0682_*/150nt/hg19a/novo-3.04.06-k17-s2/*.bam; do
    name=$(basename $bam | sed 's/.bam//g')
    echo $name
    $bamstats -g 2864785220 -q 10 -b $bam > $dirOut/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $name $dirOut/$name.bamstats
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut
for fq in /home/aldente/private/Projects/Martin_Hirst/PX0682/AnalyzedData/PX0682-1..1/Solexa/Data/current/BaseCalls/Lane_5/*.fastq.gz; do
    /projects/epigenomics/software/FastQC/fastqc -j $java -o $dirOut -t 10 $fq 
done
for bam in /projects/analysis/analysis30/PX0682/HCTVVCCXY_5/PX0682_*/150nt/hg19a/novo-3.04.06-k17-s2/*.bam; do
    /projects/epigenomics/software/FastQC/fastqc -j $java -o $dirOut -t 10 $bam 
done

