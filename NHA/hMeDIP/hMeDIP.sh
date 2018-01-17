#!/bin/sh

# link bam files
dirIn=/projects/analysis/analysis30/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/
mkdir -p $dirOut/bam/
less $dirOut/sample_info.txt | awk '{"echo $(ls ""'$dirIn'"$4"/*/"$4"_"$3"/125nt/hg19a/bwa-0.5.7/*.bam)" | getline bam; print $1"\t"$4"\t"$3"\t"bam}' > $dirOut/sample_info1.txt
mv $dirOut/sample_info1.txt $dirOut/sample_info.txt
less $dirOut/sample_info.txt | awk '{system("ln -s "$4" ""'$dirOut'""/bam/"$1".bam")}'

# QC
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
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
reg=/home/lli/hg19/hg19.chrlen.autoXY.1KB.bed
for wig in $dirWig/MGG*.wig.gz; do
    name=$(basename $wig | cut -d'.' -f1)
    echo "Processing $name"
    $JAVA -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirWig -c $chr -n $name > $dirWig/$name.coverage.log
done


