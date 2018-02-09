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
dirHub=/gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/hMeDIP/hg19/
mkdir -p $dirWig
mkdir -p $dirHub
cp /gsc/www/bcgsc.ca/downloads/mb/BrainHubs/HistoneHub/genomes.txt /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/hMeDIP/
echo -e "hub VitC_gliomaHub_hMeDIP
shortLabel VitC_glioma Hub (hMeDIP)
longLabel Hub to display VitC glioma data at UCSC (hMeDIP)
genomesFile genomes.txt
email lli@bcgsc.ca" > /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/hMeDIP/hub.txt
> $dirHub/trackDb.txt
function qc {
    export PATH=/home/rislam/anaconda2/bin/:$PATH
    export PYTHONPATH=/home/rislam/anaconda2/lib/python2.7/site-packages
    sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
    bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/sambamba-bamStats
    chr=/projects/epigenomics/resources/UCSC_chr/hg19.bwa2ucsc.names
    chrsize=/home/lli/hg19/hg19.chrom.sizes
    dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/bam/
    dirBW=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/bw/
    dirWig=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/wig/
    dirHub=/gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/hMeDIP/hg19/
    mkdir -p $dirBW; mkdir -p $dirWig; mkdir -p $dirHub
    bam=$1
    name=$(basename $bam | sed 's/.bam//g')
    echo $name 
    $sambamba index $bam -t 8
    $sambamba flagstat $bam -t 8 > $dirIn/$name.flagstat
    $bamstats -g 2864785220 -t 8 $bam > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn $name $dirIn/$name.bamstats
    /home/mbilenky/bin/PETLengthDist.sh $bam 5 $dirIn 10
    bamCoverage -b $bam -o $dirBW/$name.bw --normalizeUsingRPKM --ignoreDuplicates --samFlagExclude 1028 --minMappingQuality 5 --binSize 20 --extendReads
    /home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirWig -F:1028,-q:5,-n:$name,-chr:$chr,-cp
    /home/lli/HirstLab/Pipeline/UCSC/wigToBigWig $dirWig/$name.q5.F1028.PET.wig.gz $chrsize $dirHub/$name.q5.F1028.PET.bw
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
" >> $dirHub/trackDb.txt
}
export -f qc
for bam in $dirIn/*.bam; do
    qc $bam
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn
## correlation of MGG replicates
java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
chrom=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
reg=/home/lli/hg19/hg19.chrlen.autoXY.1KB.bed
dirWig=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/wig/
for wig in $dirWig/MGG*.wig.gz; do
    name=$(basename $wig | cut -d'.' -f1)
    echo "Processing $name"
    $java -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirWig -c $chrom -n $name > $dirWig/$name.coverage.log
done
join $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_control1.coverage $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_control2.coverage -1 4 -2 4 | awk -F' ' '{if($2=="chr1"){print $1"\t"$5"\t"$10}}' > $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_control.coverage
join $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_vitc1.coverage $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_vitc2.coverage -1 4 -2 4 | awk -F' ' '{if($2=="chr1"){print $1"\t"$5"\t"$10}}' > $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_vitc.coverage
join $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_control1.coverage $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_vitc1.coverage -1 4 -2 4 | awk -F' ' '{if($2=="chr1"){print $1"\t"$5"\t"$10}}' > $dirWig/hg19.chrlen.autoXY.1KB.bed.MGG_control1_vitc1.coverage
$sambamba merge $dirIn/MGG_control.bam $dirIn/MGG_control1.bam $dirIn/MGG_control2.bam -t 8
$sambamba merge $dirIn/MGG_vitc.bam $dirIn/MGG_vitc1.bam $dirIn/MGG_vitc2.bam -t 8
## mappability problem: re-align with bwa mem instead of bwa aln
java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
fastqc=/projects/epigenomics/software/FastQC/fastqc
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
bwa=/home/pubseq/BioSw/bwa/bwa-0.7.5a/bwa
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
genome=/home/lli/hg19/GRCh37-lite.fa
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/fq/
mkdir -p $dirOut
for file in $dirIn/*control.bam $dirIn/*vitc.bam; do
    sample=$(basename $file | sed 's/.bam//g')
    echo $sample
    $fastqc -j $java -o $dirIn -t 8 $file 
    $sambamba sort $file -o $dirIn/$sample.nsorted.bam -n --tmpdir /projects/epigenomics3/temp/ -t 8
    $BEDTOOLS/bamToFastq -i $dirIn/$sample.nsorted.bam -fq $dirOut/$sample.1.fq -fq2 $dirOut/$sample.2.fq
done
rm $dirIn/*.nsorted.bam
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
> $dirHub/trackDb.txt
for bam in $dirIn/*realign.bam; do
    qc $bam
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn
export PATH=/home/rislam/anaconda2/bin/:$PATH
export PYTHONPATH=/home/rislam/anaconda2/lib/python2.7/site-packages
cd /projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/bw/
computeMatrix scale-regions -S MGG_control.realign.bw MGG_vitc.realign.bw NHA_control.realign.bw NHAR_control.realign.bw NHAR_vitc.realign.bw NHA_vitc.realign.bw -R /home/lli/hg19/hg19.chrlen.autoXY.200.bed -out hMeDIP.200.matrix.txt --outFileNameMatrix hMeDIP.200.matrix.txt --regionBodyLength 200 --binSize 200 --sortRegions keep --skipZeros

# FindER
java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
FindER=/home/mbilenky/bin/Solexa_Java/FindER.0.9.3b.jar
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/PET_200/75bp/map/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/FindER/
mkdir -p $dirOut
echo -e "Sample\tN_region\tTotal_length\tTotal_reads\tReads_in_peaks\tAverage_length\tPercent_reads_in_peaks" > $dirOut/ER_summary.txt
for bam in $dirIn/*realign.bam; do
    sample=$(basename $bam | sed 's/.realign.bam//g')
    echo $sample
    $java -jar -Xmx10G $FindER -i $bam -r $reg -o $dirOut -v -m $map -minER 200 -info -cp > $dirOut/FindER_scan.$sample.log
    n_all=$($sambamba view $bam -c -F "not (unmapped or duplicate) and mapping_quality >= 5")
    n_peak=$($sambamba view $bam -c -F "not (unmapped or duplicate) and mapping_quality >= 5" -L <(less $dirOut/FindER_scan.$sample.realign.pctl_0.1.FDR_0.05.bed | sed 's/chr//g'))
    echo -e $sample"\t"$(less $dirOut/FindER_scan.$sample.realign.pctl_0.1.FDR_0.05.bed | wc -l)"\t"$(less $dirOut/FindER_scan.$sample.realign.pctl_0.1.FDR_0.05.bed | awk '{s=s+$3-$2}END{print s}')"\t"$n_all"\t"$n_peak | awk '{print $0"\t"int($3/$2)"\t"$5/$4}' >> $dirOut/ER_summary.txt
done
FindER=/home/mbilenky/bin/Solexa_Java/FindER.1.0.0b.jar
echo -e "Sample1\tSample2\tunique\tN_region\tTotal_length\tAverage_length" > $dirOut/ER_unique_summary.txt
for file in $dirIn/*vitc.realign.bam; do
    sample1=$(basename $file | sed 's/.realign.bam//g')
    sample2=$(echo $sample1 | sed 's/vitc/control/g')
    echo $sample1 $sample2
    $java -jar -Xmx12G $FindER -signalBam $dirIn/$sample1.realign.bam -inputBam $dirIn/$sample2.realign.bam -out $dirOut > $dirOut/"$sample1"_"$sample2".log
    echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirOut/$sample1.realign.vs.$sample2.realign.FDR_0.05.FindER.bed.gz | wc -l)"\t"$(less $dirOut/$sample1.realign.vs.$sample2.realign.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"int($5/$4)}' >> $dirOut/ER_unique_summary.txt
    $java -jar -Xmx12G $FindER -signalBam $dirIn/$sample2.realign.bam -inputBam $dirIn/$sample1.realign.bam -out $dirOut > $dirOut/"$sample2"_"$sample1".log
    echo -e $sample1"\t"$sample2"\t"$sample2"\t"$(less $dirOut/$sample2.realign.vs.$sample1.realign.FDR_0.05.FindER.bed.gz | wc -l)"\t"$(less $dirOut/$sample2.realign.vs.$sample1.realign.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"int($5/$4)}' >> $dirOut/ER_unique_summary.txt
done
sample1=NHAR_control; sample2=NHA_control
echo $sample1 $sample2
$java -jar -Xmx12G $FindER -signalBam $dirIn/$sample1.realign.bam -inputBam $dirIn/$sample2.realign.bam -out $dirOut > $dirOut/"$sample1"_"$sample2".log
echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirOut/$sample1.realign.vs.$sample2.realign.FDR_0.05.FindER.bed.gz | wc -l)"\t"$(less $dirOut/$sample1.realign.vs.$sample2.realign.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"int($5/$4)}' >> $dirOut/ER_unique_summary.txt
$java -jar -Xmx12G $FindER -signalBam $dirIn/$sample2.realign.bam -inputBam $dirIn/$sample1.realign.bam -out $dirOut > $dirOut/"$sample2"_"$sample1".log
echo -e $sample1"\t"$sample2"\t"$sample2"\t"$(less $dirOut/$sample2.realign.vs.$sample1.realign.FDR_0.05.FindER.bed.gz | wc -l)"\t"$(less $dirOut/$sample2.realign.vs.$sample1.realign.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"int($5/$4)}' >> $dirOut/ER_unique_summary.txt

# MACS2
export PATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/bin:$PATH
export PYTHONPATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/lib/python2.7/site-packages:$PYTHONPATH
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/MACS2/
mkdir -p $dirOut
echo -e "Sample\tN_region\tTotal_length\tTotal_reads\tReads_in_peaks\tAverage_length\tPercent_reads_in_peaks" > $dirOut/ER_summary.txt
for file in $dirIn/*realign.bam; do
    sample=$(basename $file | sed 's/.realign.bam//g')
    echo $sample
    macs2 callpeak -f BAMPE -g hs -t $file -q 0.01 -n $sample --outdir $dirOut
    n_all=$($sambamba view $file -c -F "not (unmapped or duplicate) and mapping_quality >= 5")
    n_peak=$($sambamba view $file -c -F "not (unmapped or duplicate) and mapping_quality >= 5" -L $dirOut/$sample"_peaks.narrowPeak")
    echo -e $sample"\t"$(less $dirOut/$sample"_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample"_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}')"\t"$n_all"\t"$n_peak | awk '{print $0"\t"int($3/$2)"\t"$5/$4}' >> $dirOut/ER_summary.txt
done
echo -e "Sample1\tSample2\tunique\tN_region\tTotal_length\tAverage_length" > $dirOut/ER_unique_summary.txt
for file in $dirIn/*vitc.realign.bam; do
    sample1=$(basename $file | sed 's/.realign.bam//g')
    sample2=$(echo $sample1 | sed 's/vitc/control/g')
    echo $sample1 $sample2
    macs2 callpeak -f BAMPE -g hs -t $file -c $dirIn/$sample2.realign.bam -q 0.05 -n $sample1"_"$sample2"."$sample1"_unique" --outdir $dirOut
    echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1"_unique_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1"_unique_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"$5/$4}' >> $dirOut/ER_unique_summary.txt
    macs2 callpeak -f BAMPE -g hs -t $dirIn/$sample2.realign.bam -c $file -q 0.05 -n $sample1"_"$sample2"."$sample2"_unique" --outdir $dirOut
    echo -e $sample1"\t"$sample2"\t"$sample2"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2"_unique_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2"_unique_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"$5/$4}' >> $dirOut/ER_unique_summary.txt
done
sample1=NHAR_control; sample2=NHA_control
macs2 callpeak -f BAMPE -g hs -t $dirIn/$sample1.realign.bam -c $dirIn/$sample2.realign.bam -q 0.05 -n $sample1"_"$sample2"."$sample1"_unique" --outdir $dirOut
echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1"_unique_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1"_unique_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}')"\t"$n_all"\t"$n_peak | awk '{print $0"\t"int($5/$4)"\t"$7/$6}' >> $dirOut/ER_unique_summary.txt
macs2 callpeak -f BAMPE -g hs -t $dirIn/$sample2.realign.bam -c $dirIn/$sample1.realign.bam -q 0.05 -n $sample1"_"$sample2"."$sample2"_unique" --outdir $dirOut
echo -e $sample1"\t"$sample2"\t"$sample2"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2"_unique_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2"_unique_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}')"\t"$n_all"\t"$n_peak | awk '{print $0"\t"int($5/$4)"\t"$7/$6}' >> $dirOut/ER_unique_summary.txt

## unique ERs
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
enhancer=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/H3K27ac/H3K27ac.union.bed
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/MACS2/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique/
mkdir -p $dirOut
echo -e "Sample1\tSample2\tunique\tN_region\tTotal_length" > $dirOut/ER_unique_summary.txt
for file in $dirIn/*vitc_peaks.narrowPeak; do
    sample1=$(basename $file | sed 's/_peaks.narrowPeak//g')
    sample2=$(echo $sample1 | sed 's/vitc/control/g')
    echo $sample1 $sample2
    $BEDTOOLS/intersectBed -a $file -b $dirIn/$sample2'_peaks.narrowPeak' -v | awk '{print "chr"$0}' > $dirOut/$sample1"_"$sample2"."$sample1".unique.bed"
    echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1".unique.bed" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1".unique.bed" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
    $BEDTOOLS/intersectBed -a $dirIn/$sample2'_peaks.narrowPeak' -b $file -v | awk '{print "chr"$0}' > $dirOut/$sample1"_"$sample2"."$sample2".unique.bed"
    echo -e $sample1"\t"$sample2"\t"$sample2"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2".unique.bed" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample2".unique.bed" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
done
sample1=NHAR_control; sample2=NHA_control
$BEDTOOLS/intersectBed -a $dirIn/$sample1'_peaks.narrowPeak' -b $dirIn/$sample2'_peaks.narrowPeak' -v | awk '{print "chr"$0}' > $dirOut/$sample1"_"$sample2"."$sample1".unique.bed"
echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1".unique.bed" | wc -l)"\t"$(less $dirOut/$sample1"_"$sample2"."$sample1".unique.bed" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
$BEDTOOLS/intersectBed -a $dirIn/$sample2'_peaks.narrowPeak' -b $dirIn/$sample1'_peaks.narrowPeak' -v | awk '{print "chr"$0}' > $dirOut/$sample1"_"$sample2"."$sample2".unique.bed"
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
### intersect with enhancer
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique/
dirK27ac=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/H3K27ac/
dirOut=$dirIn/enhancer/
mkdir -p $dirOut
echo -e "Name\tN_total\tlength_total\tN_enhancer\tlength_enhancer\tpercent" > $dirOut/ER_enhancer_summary.txt
for file in $dirIn/*.bed; do
    name=$(basename $file | sed 's/.unique.bed//g')
    sample=$(basename $file | sed 's/.unique.bed//g' | sed 's/.*\.//g' | sed 's/MGG_control/MGG119_control/')
    echo $name $sample
    less $dirK27ac/H3K27ac_"$sample".vs.input_"$sample".FDR_0.05.FindER.bed.gz | awk '{print "chr"$0}' | $BEDTOOLS/intersectBed -a $file -b stdin -u > $dirOut/$name.enhancer.bed
    echo -e $name"\t"$(less $file | wc -l)"\t"$(less $file | awk '{s=s+$3-$2}END{print s}')"\t"$(less $dirOut/$name.enhancer.bed | wc -l)"\t"$(less $dirOut/$name.enhancer.bed | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"$5/$3}' >> $dirOut/ER_enhancer_summary.txt
done
#### Homer
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique/enhancer/
dirOut=$dirIn/homer/
mkdir -p $dirOut
for file in $dirIn/*.bed; do
     name=$(basename $file | sed 's/.enhancer.bed//g')
     echo $name
     mkdir -p $dirOut/$name/
     /home/lli/bin/homer/bin/findMotifsGenome.pl $file hg19 $dirOut/$name/ -size 200 -len 8 
done

## unique ER 2: include magnitude change: pairwise union of ERs -> 2-FC in RPKM & higher one RPKM >20
function uniqueER {
    export PATH=/home/rislam/anaconda2/bin/:$PATH
    export PYTHONPATH=/home/rislam/anaconda2/lib/python2.7/site-packages
    BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
    dirER=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/MACS2/
    dirBW=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/bw/
    dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique2/
    sample1=$1; sample2=$2; name=$sample1'_'$sample2
    echo $sample1 $sample2 $name
    cat $dirER/$sample1'_peaks.narrowPeak' $dirER/$sample2'_peaks.narrowPeak' | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -c 4 -o collapse > $dirER/$name'.union.narrowPeak'
    multiBigwigSummary BED-file -b $dirBW/$sample1.realign.bw $dirBW/$sample2.realign.bw --BED $dirER/$name'.union.narrowPeak' --labels $sample1 $sample2 -out $dirOut/$name'.union.npz' --outRawCounts $dirOut/$name'.union.RPKM'
    less $dirOut/$name'.union.RPKM' | awk 'NR>1 && $1 !~ /GL/ {fc=($4+0.001)/($5+0.001); if(fc>2&&$4>20){print "chr"$1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4"\t"$5"\t"fc}}' | sort -k1,1 -k2,2n > $dirOut/$name"."$sample1".unique.bed"
    echo -e $sample1"\t"$sample2"\t"$sample1"\t"$(less $dirER/$sample1'_peaks.narrowPeak' | wc -l)"\t"$(less $dirER/$sample2'_peaks.narrowPeak' | wc -l)"\t"$(less $dirOut/$name"."$sample1".unique.bed" | wc -l)"\t"$(less $dirOut/$name"."$sample1".unique.bed" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
    less $dirOut/$name'.union.RPKM' | awk 'NR>1 && $1 !~ /GL/ {fc=($5+0.001)/($4+0.001); if(fc>2&&$5>20){print "chr"$1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4"\t"$5"\t"fc}}' | sort -k1,1 -k2,2n > $dirOut/$name"."$sample2".unique.bed"
    echo -e $sample1"\t"$sample2"\t"$sample2"\t"$(less $dirER/$sample1'_peaks.narrowPeak' | wc -l)"\t"$(less $dirER/$sample2'_peaks.narrowPeak' | wc -l)"\t"$(less $dirOut/$name"."$sample2".unique.bed" | wc -l)"\t"$(less $dirOut/$name"."$sample2".unique.bed" | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER_unique_summary.txt
}
export uniqueER
enhancer=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/H3K27ac/H3K27ac.union.bed
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique2/
mkdir -p $dirOut
echo -e "Sample1\tSample2\tunique\tN_ER1\tN_ER2\tN_unique\tlength_unique" > $dirOut/ER_unique_summary.txt
uniqueER MGG_vitc MGG_control
uniqueER NHAR_control NHA_control
uniqueER NHAR_vitc NHAR_control
uniqueER NHA_vitc NHA_control
/home/lli/HirstLab/Pipeline/shell/region.intersect.sh -d $dirOut -r $enhancer -n "enhancer"
### Homer
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique2/
dirOut=$dirIn/homer/
mkdir -p $dirOut
for file in $dirIn/*.bed; do
     name=$(basename $file | sed 's/.unique.bed//g')
     echo $name
     mkdir -p $dirOut/$name/
     /home/lli/bin/homer/bin/findMotifsGenome.pl $file hg19 $dirOut/$name/ -size 200 -len 8 
done
### intersect with enhancer
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique2/
dirK27ac=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/H3K27ac/
dirOut=$dirIn/enhancer/
mkdir -p $dirOut
echo -e "Name\tN_total\tlength_total\tN_enhancer\tlength_enhancer\tpercent" > $dirOut/ER_enhancer_summary.txt
for file in $dirIn/*.bed; do
    name=$(basename $file | sed 's/.unique.bed//g')
    sample=$(basename $file | sed 's/.unique.bed//g' | sed 's/.*\.//g' | sed 's/MGG_control/MGG119_control/')
    echo $name $sample
    less $dirK27ac/H3K27ac_"$sample".vs.input_"$sample".FDR_0.05.FindER.bed.gz | awk '{print "chr"$0}' | $BEDTOOLS/intersectBed -a $file -b stdin -u > $dirOut/$name.enhancer.bed
    echo -e $name"\t"$(less $file | wc -l)"\t"$(less $file | awk '{s=s+$3-$2}END{print s}')"\t"$(less $dirOut/$name.enhancer.bed | wc -l)"\t"$(less $dirOut/$name.enhancer.bed | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"$5/$3}' >> $dirOut/ER_enhancer_summary.txt
done
#### Homer
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/unique2/enhancer/
dirOut=$dirIn/homer/
mkdir -p $dirOut
for file in $dirIn/*.bed; do
     name=$(basename $file | sed 's/.enhancer.bed//g')
     echo $name
     mkdir -p $dirOut/$name/
     /home/lli/bin/homer/bin/findMotifsGenome.pl $file hg19 $dirOut/$name/ -size 200 -len 8 
done
