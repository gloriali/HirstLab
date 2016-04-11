#!/bin/sh

# Confirm IDH1 R132H heterozygous mutation (chr2:209113112 C->T)
cd /projects/epigenomics/users/smcconnell/glioma/WGS/bam/
for file in *.bam; do
    echo $file
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools index $file
done
> /projects/epigenomics2/users/lli/glioma/WGS/IDH.txt
for file in *.bam; do
    echo $file >> /projects/epigenomics2/users/lli/glioma/WGS/IDH.txt
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools mpileup -f /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -r 2:209100951-209130798 $file >> /projects/epigenomics2/users/lli/glioma/WGS/IDH.txt
done


