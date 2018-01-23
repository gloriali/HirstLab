#!/bin/sh

# link bam files
dirIn=/projects/analysis/analysis30/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/
mkdir -p $dirOut/bam/
less $dirOut/sample_info.txt | awk '{"echo $(ls ""'$dirIn'"$4"/*/"$4"_"$3"/125nt/hg19a/bwa-0.5.7/*.bam)" | getline bam; print $1"\t"$4"\t"$3"\t"bam}' > $dirOut/sample_info1.txt
mv $dirOut/sample_info1.txt $dirOut/sample_info.txt
less $dirOut/sample_info.txt | awk '{system("ln -s "$4" ""'$dirOut'""/bam/"$1".bam")}'

# QC
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/bam/
dirWig=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/wig/
dirBW=/gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/hMeDIP/hg19/
mkdir -p $dirWig
mkdir -p $dirBW
cp /gsc/www/bcgsc.ca/downloads/mb/BrainHubs/HistoneHub/genomes.txt /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/hMeDIP/
echo -e "hub VitC_gliomaHub_hMeDIP
shortLabel VitC_glioma Hub (hMeDIP)
longLabel Hub to display VitC glioma data at UCSC (hMeDIP)
genomesFile genomes.txt
email lli@bcgsc.ca" > /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/hMeDIP/hub.txt
> $dirBW/trackDb.txt
function qc {
    java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
    samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
    bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
    chr=/projects/epigenomics/resources/UCSC_chr/hg19.bwa2ucsc.names
    chrsize=/home/lli/hg19/hg19.chrom.sizes
    dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/bam/
    dirWig=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/wig/
    dirBW=/gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/hMeDIP/hg19/
    bam=$1
    name=$(basename $bam | sed 's/.bam//g')
    echo $name 
    $samtools index $bam
    $samtools flagstat $bam > $dirIn/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $bam > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn $name $dirIn/$name.bamstats
    /home/jyzhu/anaconda2/bin/bamCoverage -b $bam -o $dirBW/$name.bw ÐnormalizeUsingRPKM ÐsamFlagExclude 4 ÐminMappingQuality 5 ÐbinSize 10 ÐextendReads ÐignoreDuplicates
    /home/mbilenky/bin/PETLengthDist.sh $bam 5 $dirIn 10
    if [[ "$name" =~ "control" ]]; then
        color="255,0,0"
    else
        color="0,0,255"
    fi
    echo -e "
track $name
shortLabel $name
longLabel hMeDIP $name
type bigWig
visibility full
maxHeightPixels 70:70:32
configurable on
autoScale on
alwaysZero on
priority 0.1
bigDataUrl $name.bw
color $color
" >> $dirBW/trackDb.txt
}
export -f qc
for bam in $dirIn/*.bam; do
    qc $bam
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn
## correlation of MGG replicates
java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
chrom=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
reg=/home/lli/hg19/hg19.chrlen.autoXY.1KB.bed
for wig in $dirWig/MGG*.wig.gz; do
    name=$(basename $wig | cut -d'.' -f1)
    echo "Processing $name"
    $java -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirWig -c $chrom -n $name > $dirWig/$name.coverage.log
done
join $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_control1.coverage $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_control2.coverage -1 4 -2 4 | awk -F' ' '{if($2=="chr1"){print $1"\t"$5"\t"$10}}' > $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_control.coverage
join $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_vitc1.coverage $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_vitc2.coverage -1 4 -2 4 | awk -F' ' '{if($2=="chr1"){print $1"\t"$5"\t"$10}}' > $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_vitc.coverage
join $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_control1.coverage $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_vitc1.coverage -1 4 -2 4 | awk -F' ' '{if($2=="chr1"){print $1"\t"$5"\t"$10}}' > $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_control1_vitc1.coverage
$samtools merge $dirIn/MGG_control.bam $dirIn/MGG_control1.bam $dirIn/MGG_control2.bam 
$samtools merge $dirIn/MGG_vitc.bam $dirIn/MGG_vitc1.bam $dirIn/MGG_vitc2.bam 
## mappability problem: trim adapter and low quality end and re-align
java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
#skewer=/projects/epigenomics/software/skewer/skewer-0.1.127-linux-x86_64
trim_galore=/projects/epigenomics/software/trim_galore/trim_galore
cutadapt=/gsc/software/linux-x86_64-centos5/python-2.7.5/bin/cutadapt
bwa=/home/pubseq/BioSw/bwa/bwa-0.7.5a/bwa
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
genome=/home/lli/hg19/GRCh37-lite.fa
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/fq/
mkdir -p $dirOut
for file in $dirIn/*control.bam $dirIn/*vitc.bam; do
    sample=$(basename $file | sed 's/.bam//g')
    echo $sample
    $samtools sort $file $dirIn/$sample.nsorted -n
    $BEDTOOLS/bamToFastq -i $dirIn/$sample.nsorted.bam -fq $dirOut/$sample.1.fq -fq2 $dirOut/$sample.2.fq
done
rm $dirIn/*.nsorted.bam
for file in $dirIn/*.bam $dirOut/*.fq; do
    /projects/epigenomics/software/FastQC/fastqc -j $java -o $dirOut -t 6 $file 
done
for fq1 in $dirOut/*.1.fq; do
    name=$(basename $fq1 | sed 's/.1.fq//g'); fq2=$dirOut/$name.2.fq
    echo $name $fq1 $fq2
    $trim_galore $fq1 $fq2 -q 30 -o $dirOut --paired --path_to_cutadapt $cutadapt > $dirOut/$name.trim.log
#   $skewer $fq1 $fq2 -o $dirOut/$name -x 'AGATCGGAAGAGCGGTTCAGCAGGAAT' -y 'AGATCGGAAGAGCGTCGTGTAGGGAAA' -l 80 -q 30 -t 8
    $bwa mem -M -t 10 $genome $dirOut/$name-trimmed-pair1.fastq $dirOut/$name-trimmed-pair2.fastq > $dirIn/$name.trim.sam
    $samtools view -Sb $dirIn/$name.trim.sam > $dirIn/$name.trim.bam
    $samtools sort $dirIn/$name.trim.bam $dirIn/$name.trim.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirIn/$name.trim.sorted.bam O=$dirIn/$name.trim.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    rm $dirIn/$name.trim.sam $dirIn/$name.trim.bam $dirIn/$name.trim.sorted.bam
    mv $dirIn/$name.trim.sorted.dupsFlagged.bam $dirIn/$name.trim.bam
done
for bam in $dirIn/*trim.bam; do
    qc $bam
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn
for fq1 in $dirOut/*.1.fq; do
    name=$(basename $fq1 | sed 's/.1.fq//g'); fq2=$dirOut/$name.2.fq
    echo $name $fq1 $fq2
    $bwa mem -M -t 10 $genome $fq1 $fq2 > $dirIn/$name.realign.sam
    $samtools view -Sb $dirIn/$name.realign.sam > $dirIn/$name.realign.bam
    $samtools sort $dirIn/$name.realign.bam $dirIn/$name.realign.sorted
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirIn/$name.realign.sorted.bam O=$dirIn/$name.realign.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    rm $dirIn/$name.realign.sam $dirIn/$name.realign.bam $dirIn/$name.realign.sorted.bam
    mv $dirIn/$name.realign.sorted.dupsFlagged.bam $dirIn/$name.realign.bam
done
for bam in $dirIn/*realign.bam; do
    qc $bam
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn

# MACS2
export PATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/bin:$PATH
export PYTHONPATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/lib/python2.7/site-packages:$PYTHONPATH
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/MACS2/
mkdir -p $dirOut
echo -e "Sample\tN_region\tTotal_length" > $dirOut/ER_summary.txt
for file in $dirIn/*trim.bam; do
    sample=$(basename $file | sed 's/.bam//g')
    echo $sample
    macs2 callpeak -f BAMPE -g hs -t $file -q 0.01 -n $sample --outdir $dirOut
    echo -e $sample"\t"$(less $dirOut/$sample"_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample"_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_summary.txt
done
echo -e "Sample1\tSample2\tunique\tN_region\tTotal_length" > $dirOut/ER_unique_summary.txt
for file in $dirIn/*control.trim.bam; do
    sample1=$(basename $file | sed 's/.bam//g')
    sample2=$(echo $sample1 | sed 's/control/vitc/g')
    echo $sample1 $sample2
    macs2 callpeak -f BAMPE -g hs -t $file -c $dirIn/$sample2.bam -q 0.01 -n $sample1"_"$sample2"."$sample1"_unique" --outdir $dirOut
    echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1"_unique_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1"_unique_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
    macs2 callpeak -f BAMPE -g hs -t $dirIn/$sample2.bam -c $file -q 0.01 -n $sample1"_"$sample2"."$sample2"_unique" --outdir $dirOut
    echo -e $sample1"\t"$sample2"\t"$sample2"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2"_unique_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2"_unique_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
done

## unique ERs
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
enhancer=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/H3K27ac/H3K27ac.union.bed
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/MACS2/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique/
mkdir -p $dirOut
echo -e "Sample1\tSample2\tunique\tN_region\tTotal_length" > $dirOut/ER_unique_summary.txt
for file in $dirIn/*vitc.trim_peaks.narrowPeak; do
    sample1=$(basename $file | sed 's/.trim_peaks.narrowPeak//g')
    sample2=$(echo $sample1 | sed 's/vitc/control/g')
    echo $sample1 $sample2
    $BEDTOOLS/intersectBed -a $file -b $dirIn/$sample2.trim_peaks.narrowPeak -v | awk '{print "chr"$0}' > $dirOut/$sample1"_"$sample2"."$sample1".unique.bed"
    echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1".unique.bed" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1".unique.bed" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
    $BEDTOOLS/intersectBed -a $dirIn/$sample2.trim_peaks.narrowPeak -b $file -v | awk '{print "chr"$0}' > $dirOut/$sample1"_"$sample2"."$sample2".unique.bed"
    echo -e $sample1"\t"$sample2"\t"$sample2"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2".unique.bed" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2".unique.bed" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
done
sample1=NHAR_control; sample2=NHA_control
$BEDTOOLS/intersectBed -a $file -b $dirIn/$sample2.trim_peaks.narrowPeak -v | awk '{print "chr"$0}' > $dirOut/$sample1"_"$sample2"."$sample1".unique.bed"
echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1".unique.bed" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1".unique.bed" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
$BEDTOOLS/intersectBed -a $dirIn/$sample2.trim_peaks.narrowPeak -b $file -v | awk '{print "chr"$0}' > $dirOut/$sample1"_"$sample2"."$sample2".unique.bed"
echo -e $sample1"\t"$sample2"\t"$sample2"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2".unique.bed" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2".unique.bed" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
/home/lli/HirstLab/Pipeline/shell/region.intersect.sh -d $dirOut -r $enhancer -n "enhancer"
### Homer
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique/
dirOut=$dirIn/homer/
mkdir -p $dirOut
for file in $dirIn/*.bed; do
     name=$(basename $file | sed 's/.unique.bed//g')
     echo $name
     mkdir -p $dirOut/$name/
     /home/lli/bin/homer/bin/findMotifsGenome.pl $file hg19 $dirOut/$name/ -size 200 -len 8 
done

