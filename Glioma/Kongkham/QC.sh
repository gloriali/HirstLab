#!/bin/sh

# QC alignment
## QC original alignment - unsorted bam
java=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirOut
for bam in *_trimmed.bam; do
    name=$(echo $bam | sed -e 's/.bam//g')
    echo $name
    $samtools sort $dirOut/$bam $dirOut/$name.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $samtools flagstat $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $dirOut/$name.sorted.dupsFlagged.bam  > $dirOut/$name.bamstats
done
## re-align a couple of libraries with Fetal Brain pipeline and our current pipeline
### current pipeline
java=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
genome=/home/lli/hg19/GRCh37-lite.fa
bwa=/home/pubseq/BioSw/bwa/bwa-0.7.5a/bwa
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirIn=/projects/epigenomics2/users/lli/glioma/Kongkham/fq/
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirIn
for file in *.fq; do
    name=$(echo $file | sed -e 's/.fastq.trimmed.fq/_current/g')
    echo $name
    $bwa mem -M -t 12 $genome $dirIn/$file  > $dirOut/$name.sam
    $samtools view -Sb $dirOut/$name.sam > $dirOut/$name.bam
    $samtools sort $dirOut/$name.bam $dirOut/$name.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $samtools flagstat $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.bamstats
done 
### Fetal Brain pipeline
java=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
genome=/home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa
bwa=/home/pubseq/BioSw/bwa/bwa-0.5.7/bwa-0.5.7/bwa
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirIn=/projects/epigenomics2/users/lli/glioma/Kongkham/fq/
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirIn
for file in *.fq; do
    name=$(echo $file | sed -e 's/.fastq.trimmed.fq/_FB/g')
    echo $name
    $bwa aln $genome $file > $dirOut/$name.sai
    $bwa samse $genome $dirOut/$name.sai $file > $dirOut/$name.sam
    $samtools view -Sb $dirOut/$name.sam > $dirOut/$name.bam
    $samtools sort $dirOut/$name.bam $dirOut/$name.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $samtools flagstat $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.flagstat   
    $bamstats -g 2864785220 -q 10 -b $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.bamstats
done
## QC reports
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirOut
for file in *.bamstats; do
    lib=$(echo $file | sed -e 's/.bamstats//g')
    echo $lib
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $lib $file
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut

# sort mark dups and QC 
java=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirOut=/projects/epigenomics2/users/lli/glioma/Kongkham/bam/
cd $dirOut
for bam in *_trimmed.bam; do
    name=$(echo $bam | sed -e 's/.bam//g')
    echo $name
    $samtools sort $dirOut/$bam $dirOut/$name.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $samtools flagstat $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $dirOut/$name.sorted.dupsFlagged.bam  > $dirOut/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $name $dirOut/$name.bamstats
    rm $dirOut/$bam $dirOut/$name.sorted.bam
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut

