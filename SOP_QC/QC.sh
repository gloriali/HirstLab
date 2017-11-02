#!/bin/sh

# MeDIP Pos and Neg qPCR primer QC
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/H3K27ac/others/
cd /projects/epigenomics2/users/lli/
echo -e "chr\tstart\tend\tID\tfractional\tn\tsample" > meDIP_qPCR_primers_hg19.5mC.bed
for file in $dir5mC/CEMT*.5mC.CpG; do
    sample=$(basename $file | sed 's/.5mC.CpG//g')
    echo $sample
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$2+1"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a meDIP_qPCR_primers_hg19.bed -b stdin -wa -wb | awk '{if($8+$9 >= 3){n[$4]=n[$4]+1; t[$4]=t[$4]+$8; c[$4]=c[$4]+$9; chr[$4]=$1; start[$4]=$2; end[$4]=$3}} END{for(i in chr){if(t[i]+c[i]>0){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"c[i]/(c[i]+t[i])"\t"n[i]"\t""'$sample'"}}}' >> meDIP_qPCR_primers_hg19.5mC.bed
done
$BEDTOOLS/intersectBed -a meDIP_qPCR_primers_hg19.bed -b /home/lli/hg19/CG.BED -c | awk '{print $0"\t"$5/($3-$2)*100}' > meDIP_qPCR_primers_hg19.CpGdensity.bed

# alignment
java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
genome=/home/lli/hg19/GRCh37-lite.fa
bwa=/home/pubseq/BioSw/bwa/bwa-0.7.5a/bwa
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirIn=/projects/epigenomics3/UBC_miseq/meDIPQC_01NOV2017/
dirOut=/projects/epigenomics3/bams/UBC_meDIPQC_01NOV2017/
mkdir -p $dirOut
cd $dirIn
for f1 in IP*3X*R1*fastq.gz; do
    f2=$(echo $f1 | sed 's/_R1_/_R2_/g'); name=$(echo $f1 | sed 's/_L001.*//g');
    echo $name
    $bwa mem -M -t 12 $genome $dirIn/$f1 $dirIn/$f2 > $dirOut/$name.sam
    $samtools view -Sb $dirOut/$name.sam > $dirOut/$name.bam
    $samtools sort $dirOut/$name.bam $dirOut/$name.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $bamstats -g 2864785220 -q 10 -b $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.bamstats
done 
cd $dirOut
for file in *.bamstats; do
    name=$(echo $file | sed 's/.bamstats//g')
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $name $dirOut/$file
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut

