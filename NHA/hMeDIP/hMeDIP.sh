#!/bin/sh

# link bam files
dirIn=/projects/analysis/analysis30/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/
mkdir -p $dirOut/bam/
less $dirOut/sample_info.txt | awk '{"echo $(ls ""'$dirIn'"$4"/*/"$4"_"$3"/125nt/hg19a/bwa-0.5.7/*.bam)" | getline bam; print $1"\t"$4"\t"$3"\t"bam}' > $dirOut/sample_info1.txt
mv $dirOut/sample_info1.txt $dirOut/sample_info.txt
less $dirOut/sample_info.txt | awk '{system("ln -s "$4" ""'$dirOut'""/bam/"$1".bam")}'

# QC
java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
chr=/projects/epigenomics/resources/UCSC_chr/hg19.bwa2ucsc.names
chrsize=/home/lli/hg19/hg19.chrom.sizes
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
for bam in $dirIn/*.bam; do
    name=$(basename $bam | sed 's/.bam//g')
    echo $name 
    $samtools index $bam
    $samtools flagstat $bam > $dirIn/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $bam > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn $name $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirWig -F:1028,-q:5,-n:$name,-chr:$chr,-cp
    /home/lli/HirstLab/Pipeline/UCSC/wigToBigWig $dirWig/$name.q5.F1028.PET.wig.gz $chrsize $dirBW/$name.q5.F1028.PET.bw
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
bigDataUrl $name.q5.F1028.PET.bw
color $color
" >> $dirBW/trackDb.txt
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn
## correlation of MGG replicates
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
for bam in $dirIn/MGG_control.bam $dirIn/MGG_vitc.bam; do
    name=$(basename $bam | sed 's/.bam//g')
    echo $name 
    $samtools index $bam
    $samtools flagstat $bam > $dirIn/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $bam > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn $name $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirWig -F:1028,-q:5,-n:$name,-chr:$chr,-cp
    /home/lli/HirstLab/Pipeline/UCSC/wigToBigWig $dirWig/$name.q5.F1028.PET.wig.gz $chrsize $dirBW/$name.q5.F1028.PET.bw
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
bigDataUrl $name.q5.F1028.PET.bw
color $color
" >> $dirBW/trackDb.txt
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn
## mappability problem
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/fq/
mkdir -p $dirOut
for bam in $dirIn/*.bam; do
    name=$(basename $bam | sed 's/.bam//g')
    echo $name 
    $samtools view -b -f 4 $bam > $dirIn/$name.unmapped.bam
    $samtools view $dirIn/$name.unmapped.bam | awk '{print ">"$1"\n"$10}' > $dirIn/$name.unmapped.fa
done
for bam in $dirIn/*.bam; do
    /projects/epigenomics/software/FastQC/fastqc -j $java -o $dirIn -t 6 $bam 
done
for file in $dirIn/*control.bam $dirIn/*vitc.bam; do
    sample=$(basename $file | sed 's/.bam//g')
    echo $sample
    $samtools sort $file $dirIn/$sample.nsorted -n
    $BEDTOOLS/bamToFastq -i $dirIn/$sample.nsorted.bam -fq $dirOut/$sample.1.fq -fq2 $dirOut/$sample.2.fq
done

# MACS2
export PATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/bin:$PATH
export PYTHONPATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/lib/python2.7/site-packages:$PYTHONPATH
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/MACS2/
mkdir -p $dirOut
echo -e "Sample\tN_region\tTotal_length" > $dirOut/ER_summary.txt
for file in $dirIn/*control.bam $dirIn/*vitc.bam; do
    sample=$(basename $file | sed 's/.bam//g')
    echo $sample
    macs2 callpeak -f BAMPE -g hs -t $file -q 0.01 -n $sample --outdir $dirOut
    echo -e $sample"\t"$(less $dirOut/$sample"_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample"_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_summary.txt
done
echo -e "Sample1\tSample2\tunique\tN_region\tTotal_length" > $dirOut/ER_unique_summary.txt
for file in $dirIn/*control.bam; do
    sample1=$(basename $file | sed 's/.bam//g')
    sample2=$(echo $sample1 | sed 's/control/vitc/g')
    echo $sample1 $sample2
    macs2 callpeak -f BAMPE -g hs -t $file -c $dirIn/$sample2.bam -q 0.01 -n $sample1"_"$sample2"."$sample1"_unique" --outdir $dirOut
    echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1"_unique_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1"_unique_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
    macs2 callpeak -f BAMPE -g hs -t $dirIn/$sample2.bam -c $file -q 0.01 -n $sample1"_"$sample2"."$sample2"_unique" --outdir $dirOut
    echo -e $sample1"\t"$sample2"\t"$sample2"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2"_unique_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2"_unique_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
done


