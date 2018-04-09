#!/bin/sh

## link bam files
dirIn=/projects/analysis/analysis30/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/
less $dirOut/sample_info.txt | awk '{"echo $(ls ""'$dirIn'"$1"/*/"$1"_"$2"/75nt/hg19a/bwa-0.5.7/*.bam)" | getline bam; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"bam}' > $dirOut/sample_info1.txt
mv $dirOut/sample_info1.txt $dirOut/sample_info.txt
less $dirOut/sample_info.txt | awk '{system("ln -s "$6" ""'$dirOut'""/bam/"$3"/"$3"_"$4".bam")}'
## recover bams
function recover {
    file=$1
    dirOut=/projects/epigenomics3/bams/hg19/
    PX=$(echo $file | sed 's/.*Hirst\///' | sed 's/\/Analyzed.*//')
    idx=$(basename $file | cut -d'_' -f3)
    mkdir -p $dirOut/$PX/
    /gsc/software/linux-x86_64-centos6/spec-1.3.2/spec2bam --threads 2 --in $file --out $dirOut/$PX/$PX"_"$idx.bam --ref /projects/sbs_archive2/spec_ref/9606/hg19/1000genomes/GRCh37-lite.fa.spec.ref
}
export -f recover
cd /projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/
less sample_info.txt | awk '{system("echo $(ls /home/aldente/private/Projects/Martin_Hirst/"$1"/AnalyzedData/"$1"*/Solexa/Data/current/BaseCalls/BWA_*/compressed_bams/*_"$2"_*.sorted.bam.spec)")}' > List.txt
cat List.txt | parallel --gnu -j 8 recover
dirIn=/projects/epigenomics3/bams/hg19/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/
rm $dirOut/bam/*/*.bam
less $dirOut/sample_info.txt | awk '{system("ln -s ""'$dirIn'"$1"/"$1"_"$2".bam ""'$dirOut'""bam/"$3"/"$3"_"$4".bam")}'
ls $dirOut/bam/*/*.bam | parallel --gnu -j 5 $sambamba index -t 6

## QC
export PATH=/home/rislam/anaconda2/bin/:$PATH
export PYTHONPATH=/home/rislam/anaconda2/lib/python2.7/site-packages
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/sambamba-bamStats
chr=/projects/epigenomics/resources/UCSC_chr/hg19.bwa2ucsc.names
chrsize=/home/lli/hg19/hg19.chrom.sizes
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/bam/
dirQC=/projects/edcc_new/reference_epigenomes/housekeeping/EDCCProd/resources/ChIPQC1/Brads/
dirWig=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/wig/
dirBW=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/BW/
dirHub=/gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/HistoneHub/hg19/
mkdir -p $dirBW; mkdir -p $dirWig; mkdir -p $dirHub
cp /gsc/www/bcgsc.ca/downloads/mb/BrainHubs/HistoneHub/genomes.txt /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/HistoneHub/
echo -e "hub VitC_gliomaHub_HistoneMods
shortLabel VitC_glioma Hub (Histone)
longLabel Hub to display VitC glioma data at UCSC (HisoneMods)
genomesFile genomes.txt
email lli@bcgsc.ca" > /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/HistoneHub/hub.txt
> $dirHub/trackDb.txt
echo -e "Library\tNumber_Target_Regions\tNumber_total\tPercent_Target_Regions" > $dirIn/QC.target_regions.txt
for bam in $dirIn/*/*.bam; do
    name=$(basename $bam | sed 's/.bam//g')
    mark=$(echo $name | cut -d'_' -f1)
    sample=$(echo $name | sed 's/.*_NHA/NHA/g' | sed 's/.*_MGG/MGG/g' | sed 's/MGG119/MGG/g')
    echo $name $mark $sample
    mkdir -p $dirIn/$mark/; mkdir -p $dirWig/$mark/; mkdir -p $dirBW/$mark/; 
    $sambamba index $bam -t 8
    $sambamba flagstat $bam -t 8 > $dirIn/$mark/$name.flagstat
    $bamstats -g 2864785220 -t 8 $bam > $dirIn/$mark/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn/$mark/ $name $dirIn/$mark/$name.bamstats
    /home/mbilenky/bin/PETLengthDist.sh $bam 5 $dirIn/$mark/ 10
    bamCoverage -b $bam -o $dirBW/$mark/$name.bw --normalizeUsingRPKM --ignoreDuplicates --samFlagExclude 1028 --minMappingQuality 5 --binSize 20 --extendReads
    /home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirWig/$mark -F:1028,-q:5,-n:$name,-chr:$chr,-cp
    /home/lli/HirstLab/Pipeline/UCSC/wigToBigWig $dirWig/$mark/$name.q5.F1028.PET.wig.gz $chrsize $dirHub/$name.q5.F1028.PET.bw
    if [ "$mark" == "input" ]; then
        n=0
        color="0,0,0"
    else
        if [ "$mark" == "H3K4me3" ]; then
            region=$dirQC/ensembl_TSS.uniq.sorted.bed.uniq
            color="255,0,0"
        fi
        if [ "$mark" == "H3K27ac" ]; then
            region=$dirQC/ensembl_TSS.uniq.sorted.bed.uniq
            color="51,102,255"
        fi
        if [ "$mark" == "H3K4me1" ]; then
            region=$dirQC/encode_ChromHMM_enhancer_state_7.sorted.merged.bed
            color="255,102,51"
        fi
        if [ "$mark" == "H3K9me3" ]; then
            region=$dirQC/ensembl_Znf.uniq.sorted.bed
            color="0,0,102"
        fi
        if [ "$mark" == "H3K27me3" ]; then
            region=$dirQC/HOX_clusters.sorted.bed
            color="102,51,0"
        fi
        if [ "$mark" == "H3K36me3" ]; then
            region=$dirQC/ensembl_genenames.uniq.sorted.bed
            color="153,0,153"
        fi
        echo $region $color
        n_all=$($sambamba view $bam -c -F "not (unmapped or duplicate) and mapping_quality >= 5")
        n_region=$($sambamba view $bam -c -F "not (unmapped or duplicate) and mapping_quality >= 5" -L $region)
    fi
    echo -e "$name\t$n_region\t$n_all" | awk '{print $0"\t"$2/$3}' >> $dirIn/QC.target_regions.txt
    echo -e "
track $name
shortLabel $name
longLabel $mark $sample
type bigWig
visibility full
maxHeightPixels 70:70:32
configurable on
autoScale on
alwaysZero on
priority 0.1
bigDataUrl $name.q5.F1028.PET.bw
color $color
" >> $dirHub/trackDb.txt
done
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 input; do
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn/$mark/ $dirIn/$mark/
done
cut -f2- $dirIn/H3K27ac/summary.xls > $dirIn/summary.txt
for mark in H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 input; do
    join $dirIn/summary.txt <(cut -f2- $dirIn/$mark/summary.xls) > $dirIn/a
    mv $dirIn/a $dirIn/summary.txt
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
    for file in $dirIn/$mark/*.bam; do
        sample=$(basename $file | sed 's/.bam//g' | cut -d'_' -f2- | sed 's/\t/_/g')
        echo $mark $sample
        mkdir -p $dirOut/$mark/
        macs2 callpeak -f BAMPE -g hs -t $file -c $dirIn/input/input_$sample.bam -q 0.01 -n $mark"_"$sample --outdir $dirOut/$mark/
        echo -e $mark"\t"$sample"\t"$(less $dirOut/$mark/"$mark"_"$sample"_peaks.narrowPeak | wc -l)"\t"$(less $dirOut/$mark/"$mark"_"$sample"_peaks.narrowPeak | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_summary.txt
    done
done
for mark in H3K27me3 H3K9me3 H3K36me3 H3K4me1; do
    cd $dirIn/$mark
    for file in $dirIn/$mark/*.bam; do
        sample=$(basename $file | sed 's/.bam//g' | cut -d'_' -f2- | sed 's/\t/_/g')
        echo $mark $sample
        mkdir -p $dirOut/$mark/
        macs2 callpeak -f BAMPE -g hs -t $file -c $dirIn/input/input_$sample.bam --broad --broad-cutoff 0.01 -n $mark"_"$sample --outdir $dirOut/$mark/
        echo -e $mark"\t"$sample"\t"$(less $dirOut/$mark/"$mark"_"$sample"_peaks.broadPeak | wc -l)"\t"$(less $dirOut/$mark/"$mark"_"$sample"_peaks.broadPeak | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_summary.txt
    done
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
### H3K9me3 NHAR_vitC
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
Reg=/home/lli/hg19/hg19.chrlen.autoXY.1KB.bed
$BEDTOOLS/windowMaker -b <(less ~/hg19/hg19.chrlen.autoXY.bed | awk '{print $0"\t"$1}') -w 1000 -i srcwinnum > /home/lli/hg19/hg19.chrlen.autoXY.1KB.bed
cd /projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/bam/
> NHAR.H3K9me3.H3K4me3.bed
for bam in H3K9me3/H3K9me3_NHAR_vitc.bam H3K9me3/H3K9me3_NHAR_control.bam H3K4me3/H3K4me3_NHAR_vitc.bam; do
    sample=$(basename $bam | sed 's/.bam//g')
    echo $sample
    depth=$($samtools view -q 5 -F 1028 $bam | wc -l)
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a $Reg -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$5*10^9/($3-$2)/"'$depth'"}' >> NHAR.H3K9me3.H3K4me3.bed
done
### K27ac signal for super enhancer
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
RegCov=/home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
dirWig=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/wig/H3K27ac/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/H3K27ac/
dirOut=$dirIn/signal/; mkdir -p $dirOut
cd $dirIn
for file in *.FindER.bed.gz; do
    sample=$(echo $file | sed 's/.vs.*//g');
    echo $sample
    $JAVA -jar -Xmx15G $RegCov -w $dirWig/$sample.q5.F1028.PET.wig.gz -r $file -o $dirOut -c $chr -n $sample.FindER.signal > $dirOut/$sample.log
done

## unique ER
### unique ER -- strict
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/unique/
dirWig=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/wig/
mkdir -p $dirOut
> $dirOut/ER.unique.log
echo -e "Sample1\tSample2\tMark\tN1\tlen1\tN2\tlen2\tN_unique1\tlen_unique1\tN_unique2\tlen_unique2" > $dirOut/ER.unique.summary
s1=("NHAR_vitc" "NHAR_control" "NHA_vitc" "NHAR_vitc")
s2=("NHAR_control" "NHA_control" "NHA_control" "NHA_control")
for ((i=0; i<4; i++)); do
    sample1=${s1[i]}; sample2=${s2[i]}; name=$sample1"."$sample2
    for mark in H3K27ac H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3; do
        echo $name $mark
        echo -e "\n\n"$name $mark >> $dirOut/ER.unique.log
        mkdir -p $dirOut/$name/$mark/
        /home/lli/HirstLab/Pipeline/shell/ER.unique.sh -r $dirIn/$mark/$mark"_"$sample1.vs.input_$sample1.FDR_0.05.FindER.bed.gz -w $dirWig/$mark/$mark"_"$sample2.q5.F1028.PET.wig.gz -excl $dirIn/$mark/$mark"_"$sample2.vs.input_$sample2.FDR_0.05.FindER.bed.gz -o $dirOut/$name/$mark/ -n $mark.$name.$sample1 >> $dirOut/ER.unique.log
        /home/lli/HirstLab/Pipeline/shell/ER.unique.sh -excl $dirIn/$mark/$mark"_"$sample1.vs.input_$sample1.FDR_0.05.FindER.bed.gz -w $dirWig/$mark/$mark"_"$sample1.q5.F1028.PET.wig.gz -r $dirIn/$mark/$mark"_"$sample2.vs.input_$sample2.FDR_0.05.FindER.bed.gz -o $dirOut/$name/$mark/ -n $mark.$name.$sample2 >> $dirOut/ER.unique.log
        N1=$(less $dirIn/$mark/$mark"_"$sample1.vs.input_$sample1.FDR_0.05.FindER.bed.gz | wc -l)
        len1=$(less $dirIn/$mark/$mark"_"$sample1.vs.input_$sample1.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
        N2=$(less $dirIn/$mark/$mark"_"$sample2.vs.input_$sample2.FDR_0.05.FindER.bed.gz | wc -l)
        len2=$(less $dirIn/$mark/$mark"_"$sample2.vs.input_$sample2.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
        N_unique1=$(less $dirOut/$name/$mark/$mark.$name.$sample1.unique | wc -l)
        len_unique1=$(less $dirOut/$name/$mark/$mark.$name.$sample1.unique | awk '{s=s+$3-$2}END{print s}')
        N_unique2=$(less $dirOut/$name/$mark/$mark.$name.$sample2.unique | wc -l)
        len_unique2=$(less $dirOut/$name/$mark/$mark.$name.$sample2.unique | awk '{s=s+$3-$2}END{print s}')
        echo -e "$sample1\t$sample2\t$mark\t$N1\t$len1\t$N2\t$len2\t$N_unique1\t$len_unique1\t$N_unique2\t$len_unique2" >> $dirOut/ER.unique.summary
    done
done
### unique ER -- non-overlapping
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/unique2/
mkdir -p $dirOut
echo -e "Sample1\tSample2\tMark\tN1\tlen1\tN2\tlen2\tN_unique1\tlen_unique1\tN_unique2\tlen_unique2" > $dirOut/ER.unique.summary
s1=("NHAR_vitc" "NHAR_control" "NHA_vitc" "NHAR_vitc" "MGG_vitc")
s2=("NHAR_control" "NHA_control" "NHA_control" "NHA_control" "MGG119_control")
for ((i=0; i<5; i++)); do
    sample1=${s1[i]}; sample2=${s2[i]}; name=$sample1"."$sample2
    for mark in H3K27ac H3K4me1 H3K4me3 H3K27me3 H3K36me3; do
        echo $name $mark
        mkdir -p $dirOut/$name/$mark/
        $BEDTOOLS/intersectBed -a <(less $dirIn/$mark/$mark"_"$sample1.vs.input_$sample1.FDR_0.05.FindER.bed.gz) -b <(less $dirIn/$mark/$mark"_"$sample2.vs.input_$sample2.FDR_0.05.FindER.bed.gz) -v | awk '{print "chr"$0}' > $dirOut/$name/$mark/$mark.$name.$sample1.unique
        $BEDTOOLS/intersectBed -a <(less $dirIn/$mark/$mark"_"$sample2.vs.input_$sample2.FDR_0.05.FindER.bed.gz) -b <(less $dirIn/$mark/$mark"_"$sample1.vs.input_$sample1.FDR_0.05.FindER.bed.gz) -v | awk '{print "chr"$0}' > $dirOut/$name/$mark/$mark.$name.$sample2.unique
        N1=$(less $dirIn/$mark/$mark"_"$sample1.vs.input_$sample1.FDR_0.05.FindER.bed.gz | wc -l)
        len1=$(less $dirIn/$mark/$mark"_"$sample1.vs.input_$sample1.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
        N2=$(less $dirIn/$mark/$mark"_"$sample2.vs.input_$sample2.FDR_0.05.FindER.bed.gz | wc -l)
        len2=$(less $dirIn/$mark/$mark"_"$sample2.vs.input_$sample2.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
        N_unique1=$(less $dirOut/$name/$mark/$mark.$name.$sample1.unique | wc -l)
        len_unique1=$(less $dirOut/$name/$mark/$mark.$name.$sample1.unique | awk '{s=s+$3-$2}END{print s}')
        N_unique2=$(less $dirOut/$name/$mark/$mark.$name.$sample2.unique | wc -l)
        len_unique2=$(less $dirOut/$name/$mark/$mark.$name.$sample2.unique | awk '{s=s+$3-$2}END{print s}')
        echo -e "$sample1\t$sample2\t$mark\t$N1\t$len1\t$N2\t$len2\t$N_unique1\t$len_unique1\t$N_unique2\t$len_unique2" >> $dirOut/ER.unique.summary
    done
done

### Homer
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
for mark in H3K27ac H3K4me1; do
    for compare in NHAR_vitc.NHAR_control NHAR_control.NHA_control NHAR_vitc.NHA_control NHA_vitc.NHA_control MGG_vitc.MGG119_control; do
        s1=$(echo $compare | sed 's/\..*//g'); s2=$(echo $compare | sed 's/.*\.//g');
        echo $mark $compare $s1 $s2
        dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/unique2/$compare/$mark/
        dirOut=$dirIn/Homer/
        mkdir -p $dirOut
        cd $dirIn
        for file in *.unique; do
            name=$(echo $file | cut -d'.' -f4)
            mkdir -p $dirOut/$name/
            /home/lli/bin/homer/bin/findMotifsGenome.pl $file hg19 $dirOut/$name/ -size 200 -len 8 
        done
        dirOut=$dirIn/Homer2/
        mkdir -p $dirOut; mkdir -p $dirOut/$s1/; mkdir -p $dirOut/$s2/
        /home/lli/bin/homer/bin/findMotifsGenome.pl $mark.$compare.$s1.unique hg19 $dirOut/$s1/ -bg $mark.$compare.$s2.unique -size 200 -len 8
        /home/lli/bin/homer/bin/findMotifsGenome.pl $mark.$compare.$s2.unique hg19 $dirOut/$s2/ -bg $mark.$compare.$s1.unique -size 200 -len 8 
    done
done

### vitc reversed regions
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/unique2/
for mark in H3K27ac H3K4me1 H3K4me3 H3K27me3 H3K36me3; do
    echo $mark
    $BEDTOOLS/intersectBed -a $dirOut/NHAR_control.NHA_control/$mark/$mark.NHAR_control.NHA_control.NHAR_control.unique -b $dirOut/NHAR_vitc.NHAR_control/$mark/$mark.NHAR_vitc.NHAR_control.NHAR_control.unique | awk '$1 !~ /GL/ {print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/$mark.mut_gain.vitc_loss.bed
    $BEDTOOLS/intersectBed -a $dirOut/NHAR_control.NHA_control/$mark/$mark.NHAR_control.NHA_control.NHA_control.unique -b $dirOut/NHAR_vitc.NHAR_control/$mark/$mark.NHAR_vitc.NHAR_control.NHAR_vitc.unique | awk '$1 !~ /GL/ {print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/$mark.mut_loss.vitc_gain.bed
done
/home/lli/HirstLab/Pipeline/shell/region.intersect.sh -d $dirOut -r /projects/epigenomics/resources/UCSC_hg19/rmsk/LTR.bed -n LTR
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
mkdir -p $dirOut/homer/
cd $dirOut
for mark in H3K27ac H3K4me1; do
    for file in $mark.*.bed; do
        name=$(echo $file | sed 's/.bed//g')
        echo $mark $name
        mkdir -p $dirOut/homer/$name/
        /home/lli/bin/homer/bin/findMotifsGenome.pl $file hg19 $dirOut/homer/$name/ -size 200 -len 8 
    done
done
/home/lli/bin/homer/bin/annotatePeaks.pl $dirOut/H3K27ac.mut_gain.vitc_loss.bed hg19 -m $dirOut/homer/H3K27ac.mut_gain.vitc_loss/knownResults/known11.motif -nmotifs -mbed $dirOut/H3K27ac.mut_gain.vitc_loss.GATA3.BS.bed > $dirOut/H3K27ac.mut_gain.vitc_loss.GATA3.annotate
/home/lli/bin/homer/bin/annotatePeaks.pl $dirOut/H3K27ac.mut_gain.vitc_loss.bed hg19 -m $dirOut/homer/H3K27ac.mut_gain.vitc_loss/knownResults/known26.motif -nmotifs -mbed $dirOut/H3K27ac.mut_gain.vitc_loss.Foxo1.BS.bed > $dirOut/H3K27ac.mut_gain.vitc_loss.Foxo1.annotate
/home/lli/bin/homer/bin/annotatePeaks.pl $dirOut/H3K27ac.mut_gain.vitc_loss.bed hg19 -m $dirOut/homer/H3K27ac.mut_gain.vitc_loss/knownResults/known29.motif -nmotifs -mbed $dirOut/H3K27ac.mut_gain.vitc_loss.BMYB.BS.bed > $dirOut/H3K27ac.mut_gain.vitc_loss.BMYB.annotate
/home/lli/bin/homer/bin/annotatePeaks.pl $dirOut/H3K27ac.mut_loss.vitc_gain.bed hg19 -m $dirOut/homer/H3K27ac.mut_loss.vitc_gain/knownResults/known28.motif -nmotifs -mbed $dirOut/H3K27ac.mut_loss.vitc_gain.Bcl6.BS.bed > $dirOut/H3K27ac.mut_loss.vitc_gain.Bcl6.annotate
for file in *.BS.bed; do
    name=$(echo $file | sed 's/.BS.bed//g'); echo $name
    echo -e "chr\tstart\tend\tID\tENSG\tdistance" > $dirOut/$name.closest.gene
    less $file | awk 'NR>1{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' | sort -k1,1 -k2,2n | $BEDTOOLS/closestBed -a stdin -b /home/lli/hg19/hg19v69_genes.TSS.pc.bed -D a | awk '{if($9>=-20000 && $9<=20000){print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$9}}' | sort -k5,5 >> $dirOut/$name.closest.gene
    join $dirOut/$name.closest.gene $dirOut/../../RNAseq/RPKM/vitc.RPKM -1 5 -2 1 | sed 's/ /\t/g' > $dirOut/$name.closest.gene.RPKM
done
cd /projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/unique2/
file=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/H3K27ac/H3K27ac_NHAR_control.vs.input_NHAR_control.FDR_0.05.FindER.bed.gz
$BEDTOOLS/intersectBed -a NHAR_control.NHA_control/H3K27ac/H3K27ac.NHAR_control.NHA_control.NHAR_control.unique -b NHAR_vitc.NHAR_control/H3K27ac/H3K27ac.NHAR_vitc.NHAR_control.NHAR_control.unique -u | $BEDTOOLS/intersectBed -a stdin -b MGG_vitc.MGG119_control/H3K27ac/H3K27ac.MGG_vitc.MGG119_control.MGG119_control.unique -u | $BEDTOOLS/intersectBed -a <(less $file | awk '{print "chr"$0}') -b stdin -u | awk '{print $0"\tMGGintersect"}' > H3K27ac_NHAR_control.category.1
$BEDTOOLS/intersectBed -a NHAR_control.NHA_control/H3K27ac/H3K27ac.NHAR_control.NHA_control.NHAR_control.unique -b NHAR_vitc.NHAR_control/H3K27ac/H3K27ac.NHAR_vitc.NHAR_control.NHAR_control.unique -u | $BEDTOOLS/intersectBed -a stdin -b MGG_vitc.MGG119_control/H3K27ac/H3K27ac.MGG_vitc.MGG119_control.MGG119_control.unique -v | $BEDTOOLS/intersectBed -a <(less $file | awk '{print "chr"$0}') -b stdin -u | awk '{print $0"\tVitCresponse"}' > H3K27ac_NHAR_control.category.2
$BEDTOOLS/intersectBed -a NHAR_control.NHA_control/H3K27ac/H3K27ac.NHAR_control.NHA_control.NHAR_control.unique -b NHAR_vitc.NHAR_control/H3K27ac/H3K27ac.NHAR_vitc.NHAR_control.NHAR_control.unique -u | $BEDTOOLS/intersectBed -a <(less $file | awk '{print "chr"$0}') -b stdin -v | awk '{print $0"\tnon-response"}' > H3K27ac_NHAR_control.category.3
cat H3K27ac_NHAR_control.category.1 H3K27ac_NHAR_control.category.2 H3K27ac_NHAR_control.category.3 > H3K27ac_NHAR_control.category

### use NPC as an outgroup
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/H3K27ac/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/unique2/
file1=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/A19308.NPC_GE04.vs.A19309.NPC_GE04.FDR_0.05.FindER.bed.gz
for sample in MGG119_control NHAR_control MGG_vitc NHAR_vitc; do
    name=$sample.NPC_GE04
    echo $name
    mkdir -p $dirOut/$name/H3K27ac/
    $BEDTOOLS/intersectBed -a <(less $dirIn/H3K27ac_$sample.vs.input_$sample.FDR_0.05.FindER.bed.gz | awk '{print "chr"$0}') -b $file1 -v > $dirOut/$name/H3K27ac/H3K27ac.$name.$sample.unique
    $BEDTOOLS/intersectBed -b <(less $dirIn/H3K27ac_$sample.vs.input_$sample.FDR_0.05.FindER.bed.gz | awk '{print "chr"$0}') -a $file1 -v > $dirOut/$name/H3K27ac/H3K27ac.$name.NPC_GE04.unique
done

### categorize promoters and enhancers
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/K4me3_K27me3/
mkdir -p $dirOut
echo -e "Sample\tBivalent\tH3K4me3\tH3K27me3\tUnmarked" > $dirOut/promoter.K4me3_K27me3.summary
for file in $dirIn/H3K4me3/*.FindER.bed.gz; do
    sample=$(basename $file | cut -d'.' -f1 | sed 's/H3K4me3_//g' | sed 's/MGG119/MGG/g'); file2=$(echo $file | sed 's/H3K4me3/H3K27me3/g')
    echo $sample
    $BEDTOOLS/intersectBed -a $file -b $file2 | $BEDTOOLS/intersectBed -a $promoter -b stdin -u | awk '{print $0"\tbivalent"}' > $dirOut/$sample.promoter.K4me3_K27me3.bed
    bivalent=$($BEDTOOLS/intersectBed -a $file -b $file2 | $BEDTOOLS/intersectBed -a $promoter -b stdin -u | wc -l)
    $BEDTOOLS/intersectBed -a $promoter -b $file -u | $BEDTOOLS/intersectBed -a stdin -b $dirOut/$sample.promoter.K4me3_K27me3.bed -v -f 1 | awk '{print $0"\tH3K4me3"}' >> $dirOut/$sample.promoter.K4me3_K27me3.bed
    H3K4me3=$(expr $(less $dirOut/$sample.promoter.K4me3_K27me3.bed | wc -l) - $bivalent)
    $BEDTOOLS/intersectBed -a $promoter -b $file2 -u | $BEDTOOLS/intersectBed -a stdin -b $dirOut/$sample.promoter.K4me3_K27me3.bed -v -f 1 | awk '{print $0"\tH3K27me3"}' >> $dirOut/$sample.promoter.K4me3_K27me3.bed
    H3K27me3=$(expr $(less $dirOut/$sample.promoter.K4me3_K27me3.bed | wc -l) - $((bivalent+H3K4me3)))
    $BEDTOOLS/intersectBed -a $promoter -b $dirOut/$sample.promoter.K4me3_K27me3.bed -v -f 1 | awk '{print $0"\tUnmarked"}' >> $dirOut/$sample.promoter.K4me3_K27me3.bed
    unmarked=$(expr $(less $dirOut/$sample.promoter.K4me3_K27me3.bed | wc -l) - $((bivalent+H3K4me3+H3K27me3)))
    echo -e "$sample\t$bivalent\t$H3K4me3\t$H3K27me3\t$unmarked" >> $dirOut/promoter.K4me3_K27me3.summary
done
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/K27ac_K4me1/
mkdir -p $dirOut
echo -e "Sample\tH3K27ac\tH3K4me1\tActive\tWeak\tLength_active\tLength_weak" > $dirOut/enhancer.K27ac_K4me1.summary
for file in $dirIn/H3K27ac/*.FindER.bed.gz; do
    sample=$(basename $file | cut -d'.' -f1 | sed 's/H3K27ac_//g' | sed 's/MGG119/MGG/g'); file2=$(echo $file | sed 's/H3K27ac/H3K4me1/g')
    echo $sample
    $BEDTOOLS/intersectBed -a $file -b $file2 -u | $BEDTOOLS/intersectBed -a stdin -b $promoter -v > $dirOut/$sample.active.enhancer.bed
    $BEDTOOLS/intersectBed -a $file2 -b $file -v | $BEDTOOLS/intersectBed -a stdin -b $promoter -v > $dirOut/$sample.weak.enhancer.bed
    echo -e $sample"\t"$(less $file | wc -l)"\t"$(less $file2 | wc -l)"\t"$(less $dirOut/$sample.active.enhancer.bed | wc -l)"\t"$(less $dirOut/$sample.weak.enhancer.bed | wc -l)"\t"$(less $dirOut/$sample.active.enhancer.bed | awk '{s=s+$3-$2}END{print s}')"\t"$(less $dirOut/$sample.weak.enhancer.bed | awk '{s=s+$3-$2}END{print s}') >> $dirOut/enhancer.K27ac_K4me1.summary
done

## FindER 2
java=/gsc/software/linux-x86_64-centos6/jdk1.8.0_162/jre/bin/java
finder2=/home/mbilenky/bin/FindER2/finder2.jar
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER2/
mkdir -p $dirOut
for inp in $dirIn/input/*.bam; do
    sample=$(basename $inp | cut -d'.' -f1 | sed 's/input_//');
    echo $sample
    sig1=$(echo $dirIn/H3K27me3/H3K27me3_$sample.bam)
    sig2=$(echo $dirIn/H3K4me3/H3K4me3_$sample.bam)
    sig3=$(echo $dirIn/H3K4me1/H3K4me1_$sample.bam)
    sig4=$(echo $dirIn/H3K27ac/H3K27ac_$sample.bam)
    sig5=$(echo $dirIn/H3K36me3/H3K36me3_$sample.bam)
    sig=${sig1},${sig2},${sig3},${sig4},${sig5}
    $java -jar -Xmx25G $finder2 inputBam:$inp signalBam:$sig outDir:$dirOut acgtDir:/projects/epigenomics2/users/mbilenky/resources/hg19/ACGT 
done
echo -e "Mark\tSample\tN_region\tTotal_length\tAverage_length" > $dirOut/ER_summary.txt
for file in $dirOut/*.bed; do
    mark=$(basename $file | cut -d'.' -f1 | cut -d'_' -f1)
    sample=$(basename $file | cut -d'.' -f1 | cut -d'_' -f2,3)
    echo $mark $sample
    echo -e $mark"\t"$sample"\t"$(less $file | wc -l)"\t"$(less $file | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"$4/$3}' >> $dirOut/ER_summary.txt
done

