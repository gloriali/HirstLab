#!/bin/sh

# link bam files
dirIn=/projects/analysis/analysis30/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/
less $dirOut/sample_info.txt | awk '{"echo $(ls ""'$dirIn'"$2"/*/"$2"_"$3"/75nt/hg19a/bwa-0.5.7/*.bam)" | getline bam; print $0"\t"bam}' > $dirOut/sample_info1.txt
mv $dirOut/sample_info1.txt $dirOut/sample_info.txt
less $dirOut/sample_info.txt | awk '{system("ln -s "$6" ""'$dirOut'""/bam/"$4"/"$4"_"$5".bam")}'

# QC
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
chr=/projects/epigenomics/resources/UCSC_chr/hg19.bwa2ucsc.names
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/bam/
dirWig=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/wig/
dirQC=/projects/edcc_new/reference_epigenomes/housekeeping/EDCCProd/resources/ChIPQC1/Brads/
echo -e "Library\tNumber_Target_Regions" > $dirIn/QC.target_regions.txt
for bam in $dirIn/*/*.bam; do
    name=$(basename $bam | sed 's/.bam//g')
    mark=$(echo $name | cut -d'_' -f1)
    echo $name $mark
    mkdir -p $dirIn/$mark/
    mkdir -p $dirWig/$mark/
    $samtools index $bam
    $samtools flagstat $bam > $dirIn/$mark/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $bam  > $dirIn/$mark/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn/$mark/ $name $dirIn/$mark/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirWig/$mark/ -F:1028,-q:5,-n:$name,-chr:$chr,-cp
    /home/mbilenky/bin/PETLengthDist.sh $bam 5 $dirIn/$mark/ 10
    if [ $mark == "input" ]; then
        n=0
    else
        if [ "$mark" == "H3K4me3" ] || [ "$mark" == "H3K27ac" ]; then
            region=$dirQC/ensembl_TSS.uniq.sorted.bed.uniq
        fi
        if [ "$mark" == "H3K4me1" ]; then
            region=$dirQC/encode_ChromHMM_enhancer_state_7.sorted.merged.bed
        fi
        if [ "$mark" == "H3K9me3" ]; then
            region=$dirQC/ensembl_Znf.uniq.sorted.bed
        fi
        if [ "$mark" == "H3K27me3" ]; then
            region=$dirQC/HOX_clusters.sorted.bed
        fi
        if [ "$mark" == "H3K36me3" ]; then
            region=$dirQC/ensembl_genenames.uniq.sorted.bed
        fi
        echo $region
        n=$($samtools view -q 5 -F 1028 -L $region $bam | wc -l)
    fi
    echo -e "$name\t$n" >> $dirIn/QC.target_regions.txt
done
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 input; do
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn/$mark/ $dirIn/$mark/
done
cut -f2- $dirIn/H3K27ac/summary.xls > $dirIn/summary.txt
for mark in H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 input; do
    join $dirIn/summary.txt <(cut -f2- $dirIn/$mark/summary.xls) > $dirIn/a
    mv $dirIn/a $dirIn/summary.txt
done

## FindER
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
FindER=/home/mbilenky/bin/Solexa_Java/FindER.1.0.1e.jar
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/
mkdir -p $dirOut/
echo -e "Mark\tSample\tN_region\tTotal_length" > $dirOut/ER_summary.txt
cd $dirIn
for file in H*/*.bam; do
    mark=$(basename $file | sed 's/_.*//g')
    file=$(basename $file)
    sample=$(echo $file | sed 's/.bam//g' | cut -d'_' -f2- | sed 's/\t/_/g')
    echo $mark $sample
    mkdir -p $dirOut/$mark/
    $JAVA -jar -Xmx12G $FindER -signalBam $dirIn/$mark/$file -inputBam $dirIn/input/input_$sample.bam -out $dirOut/$mark/ > $dirOut/$mark/"$mark"_"$sample".log
    echo -e $mark"\t"$sample"\t"$(less $dirOut/$mark/"$mark"_"$sample".vs.input_"$sample".FDR_0.05.FindER.bed.gz | wc -l)"\t"$(less $dirOut/$mark/"$mark"_"$sample".vs.input_"$sample".FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_summary.txt
done

## MACS2
export PATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/bin:$PATH
export PYTHONPATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/lib/python2.7/site-packages:$PYTHONPATH
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/MACS2/
mkdir -p $dirOut/
echo -e "Mark\tSample\tN_region\tTotal_length" > $dirOut/ER_summary.txt
for mark in H3K27ac H3K4me3; do
    cd $dirIn/$mark
    for file in *.bam; do
        sample=$(echo $file | sed 's/.bam//g' | cut -d'_' -f2- | sed 's/\t/_/g')
        echo $mark $sample
        mkdir -p $dirOut/$mark/
        macs2 callpeak -f BAMPE -g hs -t $file -c $dirIn/input/input_$sample.bam -q 0.01 -n $mark"_"$sample --outdir $dirOut/$mark/
        echo -e $mark"\t"$sample"\t"$(less $dirOut/$mark/"$mark"_"$sample"_peaks.narrowPeak | wc -l)"\t"$(less $dirOut/$mark/"$mark"_"$sample"_peaks.narrowPeak | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_summary.txt
    done
done
for mark in H3K27me3 H3K9me3 H3K36me3 H3K4me1; do
    cd $dirIn/$mark
    for file in *.bam; do
        sample=$(echo $file | sed 's/.bam//g' | cut -d'_' -f2- | sed 's/\t/_/g')
        echo $mark $sample
        mkdir -p $dirOut/$mark/
        macs2 callpeak -f BAMPE -g hs -t $file -c $dirIn/input/input_$sample.bam --broad --broad-cutoff 0.01 -n $mark"_"$sample --outdir $dirOut/$mark/
        echo -e $mark"\t"$sample"\t"$(less $dirOut/$mark/"$mark"_"$sample"_peaks.broadPeak | wc -l)"\t"$(less $dirOut/$mark/"$mark"_"$sample"_peaks.broadPeak | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_summary.txt
    done
done
