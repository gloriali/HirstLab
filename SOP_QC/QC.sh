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
human=/home/lli/hg19/GRCh37-lite.fa
Lambda=/projects/epigenomics/users/mbilenky/T7_Phage/NC_001416.1.genome.fasta
T7=/projects/epigenomics/users/mbilenky/T7_Phage/NC_001604.1.genome.fasta
M13=/projects/epigenomics/users/mbilenky/T7_Phage/NC_003287.2.genome.fasta
bwa=/home/pubseq/BioSw/bwa/bwa-0.7.5a/bwa
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/sambamba-bamStats
picard=/home/pubseq/BioSw/picard/picard-tools-1.52/
dirIn=/projects/epigenomics3/UBC_miseq/meDIPQC_01NOV2017/
dirOut=/projects/epigenomics3/bams/hg19/UBC_meDIPQC_01NOV2017/
mkdir -p $dirOut
cd $dirIn
for f1 in IP*R1*fastq.gz meDIP*R1*fastq.gz; do
    f2=$(echo $f1 | sed 's/_R1_/_R2_/g'); name=$(echo $f1 | sed 's/_L001.*//g');
    echo $name
    $bwa mem -M -t 15 $human $dirIn/$f1 $dirIn/$f2 > $dirOut/$name.sam
    $sambamba view -S -f bam $dirOut/$name.sam -o $dirOut/$name.bam -t 10
    $sambamba sort $dirOut/$name.bam -o $dirOut/$name.sorted.bam --tmpdir /projects/epigenomics3/temp/ -t 10 -m 10G
    $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
    $bamstats -g 2864785220 -t 8 $dirOut/$name.sorted.dupsFlagged.bam > $dirOut/$name.bamstats
done
## align spike-ins
echo -e "Library\tGenome\tCoverage" > $dirOut/spike_in.coverage
bwa=/home/pubseq/BioSw/bwa/bwa-0.5.6/bwa
for genome in $Lambda $T7 $M13; do
    for f1 in IP*R1*fastq.gz meDIP*R1*fastq.gz; do
        f2=$(echo $f1 | sed 's/_R1_/_R2_/g'); name=$(echo $f1 | sed 's/_L001.*//g')'__'$(basename $genome | sed 's/.genome.fasta//g');
        echo $name
        $bwa aln -t 12 $genome $f1 > $dirOut/$name.R1.sai
        $bwa aln -t 12 $genome $f2 > $dirOut/$name.R2.sai
        $bwa sampe $genome $dirOut/$name.R1.sai $dirOut/$name.R2.sai $f1 $f2 > $dirOut/$name.sam
        $sambamba view -S -f bam $dirOut/$name.sam -o $dirOut/$name.bam -t 10
        $sambamba sort $dirOut/$name.bam -o $dirOut/$name.sorted.bam --tmpdir /projects/epigenomics3/temp/ -t 10 -m 10G
        $java -jar -Xmx10G $picard/MarkDuplicates.jar I=$dirOut/$name.sorted.bam O=$dirOut/$name.sorted.dupsFlagged.bam M=dups AS=true VALIDATION_STRINGENCY=LENIENT QUIET=true
        mv $dirOut/$name.sorted.dupsFlagged.bam $dirOut/$name.bam
        echo -e $name"\t"$($sambamba view $dirOut/$name.bam -c -F "not (unmapped or duplicate) and mapping_quality >= 5") | sed 's/__/\t/' >> $dirOut/spike_in.coverage
    done
done
dirIn=/projects/epigenomics3/bams/hg19/UBC_meDIPQC_01NOV2017/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/SOP_QC/
mkdir -p $dirOut
cd $dirIn
for file in *.bamstats; do
    name=$(echo $file | sed 's/.bamstats//g')
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $name $dirIn/$file
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut
rm *.sam *sai report* *sorted.bam*

# coverage vs GC content
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
dirOut=/projects/epigenomics3/bams/hg19/UBC_meDIPQC_01NOV2017/
region=/home/lli/hg19/hg19.chrlen.autoXY.1KB.bed
$BEDTOOLS/nucBed -fi /home/lli/hg19/GRCh37-lite.fa -bed $region | awk 'NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$6}' > $region.GC
cd $dirOut
echo -e "chr\tstart\tend\tID\tGC\tcoverage\tsample" > $dirOut/GCcontent.coverage
echo -e "GC\ttotal\tcount\taverage\tsample" > $dirOut/GCcontent.coverage.summary
for bam in *.sorted.dupsFlagged.bam; do
    sample=$(echo $bam | sed 's/.sorted.dupsFlagged.bam//g')
    echo $sample
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a $region.GC -b stdin -counts | awk '{print $0"\t""'$sample'"}' >> $dirOut/GCcontent.coverage
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a $region.GC -b stdin -counts | awk '{s[int($5*10)]=s[int($5*10)]+$6; c[int($5*10)]++}END{for(i in s){print i/10"\t"s[i]"\t"c[i]"\t"s[i]/c[i]"\t""'$sample'"}}' | sort -k1,1n >> $dirOut/GCcontent.coverage.summary
done

# generate BW hub
export PATH=/home/rislam/anaconda2/bin/:$PATH
export PYTHONPATH=/home/rislam/anaconda2/lib/python2.7/site-packages
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
dirIn=/projects/epigenomics3/bams/hg19/UBC_meDIPQC_01NOV2017/
dirBW=/gsc/www/bcgsc.ca/downloads/mb/SOP_QC/hg19/
mkdir -p $dirBW
cp /gsc/www/bcgsc.ca/downloads/mb/BrainHubs/HistoneHub/genomes.txt /gsc/www/bcgsc.ca/downloads/mb/SOP_QC/
echo -e "hub SOP_QC
shortLabel SOP_QC Hub (hMeDIP & MeDIP)
longLabel Hub to display hMeDIP and MeDIP SOP QC data at UCSC
genomesFile genomes.txt
email lli@bcgsc.ca" > /gsc/www/bcgsc.ca/downloads/mb/SOP_QC/hub.txt
> $dirBW/trackDb.txt
for file in $dirIn/IP-*.sorted.dupsFlagged.bam $dirIn/meDIP*.sorted.dupsFlagged.bam; do
    name=$(basename $file | sed 's/.sorted.dupsFlagged.bam//g' | sed 's/meDIP-HL60-B05_S27/IP-5mC-3X.2_S5/g')
    echo $name
    $sambamba index $file -t 8
    bamCoverage -b $file -o $dirBW/$name.bw --normalizeUsingRPKM --ignoreDuplicates --binSize 20 --extendReads
    if [[ "$name" =~ "5hmC" ]]; then
        color="255,0,0"
    else
        color="0,0,255"
    fi
    echo -e "
track $name
shortLabel $name
longLabel SOP_QC $name
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
done

# MACS2
export PATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/bin:$PATH
export PYTHONPATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/lib/python2.7/site-packages:$PYTHONPATH
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/sambamba-bamStats
dirIn=/projects/epigenomics3/bams/hg19/UBC_meDIPQC_01NOV2017/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/SOP_QC/
$sambamba view $dirIn/meDIP-HL60-B05_S27.sorted.dupsFlagged.bam -f bam -s 0.47 -t 8 -o $dirIn/IP-5mC-3X.2sub_S5.sorted.dupsFlagged.bam
$bamstats -g 2864785220 -t 8 $dirIn/IP-5mC-3X.2sub_S5.sorted.dupsFlagged.bam > $dirIn/IP-5mC-3X.2sub_S5.bamstats
/home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut IP-5mC-3X.2sub_S5 $dirIn/IP-5mC-3X.2sub_S5.bamstats
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut
echo -e "Sample\tN_region\tTotal_length\tTotal_reads\tReads_in_peaks\tAverage_length\tPercent_reads_in_peaks" > $dirOut/ER_summary.txt
for file in $dirIn/IP-5mC*.sorted.dupsFlagged.bam $dirIn/meDIP*.sorted.dupsFlagged.bam; do
    sample=$(basename $file | sed 's/.sorted.dupsFlagged.bam//g' | sed 's/meDIP-HL60-B05_S27/IP-5mC-3X.2_S5/g')
    echo $sample
    macs2 callpeak -f BAMPE -g hs -t $file -q 0.01 -n $sample --outdir $dirOut
    n_all=$($sambamba view $file -c -F "not (unmapped or duplicate) and mapping_quality >= 5")
    n_peak=$($sambamba view $file -c -F "not (unmapped or duplicate) and mapping_quality >= 5" -L $dirOut/$sample"_peaks.narrowPeak")
    echo -e $sample"\t"$(less $dirOut/$sample"_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample"_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}')"\t"$n_all"\t"$n_peak | awk '{print $0"\t"int($3/$2)"\t"$5/$4}' >> $dirOut/ER_summary.txt
done

