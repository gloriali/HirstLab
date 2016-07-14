#!/bin/sh

## QC 
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirIn=/projects/epigenomics2/Brain/NJabado/ChIPseq/bam/
dirOut=/projects/epigenomics2/Brain/NJabado/ChIPseq/QC/
mkdir -p $dirOut
cd $dirIn
for bam in *.bam; do
    name=$(echo $bam | sed -e 's/.bam//g' | sed -e 's/C9E81ANXX_[3-6]_//g')
    echo $name
    $bamstats -g 2864785220 -q 10 -b $dirIn/$bam  > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $name $dirIn/$name.bamstats
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut
dirIn=/projects/epigenomics2/Brain/NJabado/ChIPseq/bam/
dirOut=/projects/epigenomics2/Brain/NJabado/ChIPseq/QC/PETlen/
mkdir -p $dirOut
cd $dirIn
for bam in *.bam; do
    name=$(echo $bam | sed -e 's/.bam//g' | sed -e 's/C9E81ANXX_[3-6]_//g')
    echo $name
    /home/mbilenky/bin/PETLengthDist.sh $dirIn/$bam 5 $dirOut 10
done
