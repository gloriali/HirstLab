#!/bin/sh

bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/bam/unmerged/
for bam in /projects/analysis/analysis30/PX0682/HCTVVCCXY_5/PX0682_*/150nt/hg19a/novo-3.04.06-k17-s2/*.bam; do
    name=$(basename $bam | sed 's/.bam//g')
    echo $name
    $bamstats -g 2864785220 -q 10 -b $bam > $dirOut/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $name $dirOut/$name.bamstats
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut
for fq in /home/aldente/private/Projects/Martin_Hirst/PX0682/AnalyzedData/PX0682-1..1/Solexa/Data/current/BaseCalls/Lane_5/*.fastq.gz; do
    /projects/epigenomics/software/FastQC/fastqc -j $java -o $dirOut -t 10 $fq 
done
for bam in /projects/analysis/analysis30/PX0682/HCTVVCCXY_5/PX0682_*/150nt/hg19a/novo-3.04.06-k17-s2/*.bam; do
    /projects/epigenomics/software/FastQC/fastqc -j $java -o $dirOut -t 10 $bam 
done

# confirm IDH mutant status: chr2:209113112 C->T (is a CpG)
cd /projects/epigenomics3/epigenomics3_results/users/lli/NHA/
echo -e "Sample\tAssay\tT\tC\trate" > IDH.mutant.rate.txt
for bam in PBAL/bam/*.bam RNAseq/bam/*.bam ChIPseq/bam/H3K36me3/*.bam; do
    sample=$(basename $bam | cut -d'.' -f1)
    assay=$(echo $bam | sed 's/\/bam.*//')
    echo $sample $assay
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools mpileup -f /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -r 2:209113112-209113112 $bam | awk '{s=toupper($5); print "'$sample'""\t""'$assay'""\t"gsub("T", "",s)"\t"gsub(",", "", s)+gsub(".", "", s)}' | awk '{print $0"\t"$3/($3+$4)}' >> IDH.mutant.rate.txt
done

# QC
dirIn=/projects/edcc_prj2/bs-seq/PX0682/
dirBam=/projects/epigenomics3/bams/hg19/PX0682/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/
dirHub=/gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/PBAL/hg19/
mkdir -p $dirHub; mkdir -p $dirBam
cp /gsc/www/bcgsc.ca/downloads/mb/BrainHubs/HistoneHub/genomes.txt /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/PBAL/
echo -e "hub VitC_gliomaHub_PBAL
shortLabel VitC_glioma Hub (PBAL)
longLabel Hub to display VitC glioma data at UCSC (PBAL)
genomesFile genomes.txt
email lli@bcgsc.ca" > /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/PBAL/hub.txt
> $dirHub/trackDb.txt
less $dirOut/SampleInfo.txt | awk '{system("cp /projects/analysis/analysis30/"$2"_"$3"/merge*/150nt/hg19a/"$2"_"$3"_6_lanes_dupsFlagged.bam"" ""'$dirBam'"$2"_"$3"_6_lanes_dupsFlagged.bam")}'
less $dirOut/SampleInfo.txt | awk '{system("ln -s ""'$dirBam'"$2"_"$3"_6_lanes_dupsFlagged.bam"" ""'$dirOut'""/bam/"$1".6_lanes_dupsFlagged.bam")}'
less $dirOut/SampleInfo.txt | awk '{system("ln -s ""'$dirIn'"$2"_"$3"_6_lanes_dupsFlagged.Cmethyl.cons.bed.CpG.txt.gz"" ""'$dirOut'"$1".Cmethyl.cons.bed.CpG.txt.gz")}'
less $dirOut/SampleInfo.txt | awk '{system("cp ""'$dirIn'"$2"_"$3"_6_lanes_dupsFlagged.coverage.bw"" ""'$dirHub'"$1".coverage.bw")}'
less $dirOut/SampleInfo.txt | awk '{system("cp ""'$dirIn'"$2"_"$3"_6_lanes_dupsFlagged.fractional.bw"" ""'$dirHub'"$1".fractional.bw")}'
function qc {
    sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
    bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/sambamba-bamStats
    dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/bam/
    dirHub=/gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/PBAL/hg19/
    bam=$1
    name=$(basename $bam | cut -d'.' -f1)
    echo $name 
    $sambamba index $bam -t 8
    $sambamba flagstat $bam -t 8 > $dirIn/$name.flagstat
    $bamstats -g 2864785220 -t 8 $bam > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn $name $dirIn/$name.bamstats
    /home/mbilenky/bin/PETLengthDist.sh $bam 5 $dirIn 10
    echo -e "
track $name coverage
shortLabel $name cov
longLabel PBAL $name coverage
type bigWig
visibility full
maxHeightPixels 70:70:32
configurable on
autoScale on
alwaysZero on
priority 0.1
bigDataUrl $name.coverage.bw
color 0,0,0

track $name fractional
shortLabel $name 5mC
longLabel PBAL $name fractional
type bigWig
visibility full
maxHeightPixels 70:70:32
configurable on
autoScale on
alwaysZero on
priority 0.1
bigDataUrl $name.fractional.bw
color 0,255,0
" >> $dirHub/trackDb.txt
}
export -f qc
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/bam/
cd $dirIn
ls *.bam > List.txt
cat List.txt | parallel --gnu qc 
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn

# combine strands
function combine {
    file=$1
    dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/
    lib=$(echo $file | sed -e 's/.Cmethyl.cons.bed.CpG.txt.gz//g')
    /home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i $dirIn -o $dirIn -f $file -n $lib 
}
export -f combine
cd /projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/
ls *.Cmethyl.cons.bed.CpG.txt.gz > List.txt
cat List.txt | parallel --gnu combine 

# check coverage profile and 5mC profile
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/
echo -e "sample\tcoverage\tN" > $dirIn/qc_5mC_coverage.txt
echo -e "sample\ttype\tfractional\tN" > $dirIn/qc_5mC_profile.txt 
echo -e "sample\ttype\tmin\tymin\tlower\tmedian\tupper\tymax\tmax" > $dirIn/qc_5mC_quantile.txt #ymin: 10% quantile; ymax: 90% quantile
cd $dirIn
for file in *.combine.5mC.CpG; do
    lib=$(echo $file | sed -e 's/.combine.5mC.CpG//g')
    echo "Processing" $lib
    less $file | awk '{c = $4 + $5; if(c >= 5000){s[5001]++} else {s[c]++}} END{for(i = 3; i <= 5001; i++){print "'$lib'""\t"i"\t"s[i]}}' >> $dirIn/qc_5mC_coverage.txt
    less $file | awk '{s[int($6*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tgenome\t"i/100"\t"s[i]}}' >> $dirIn/qc_5mC_profile.txt 
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | awk '{s[int($1*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tCGI\t"i/100"\t"s[i]}}' >> $dirIn/qc_5mC_profile.txt 
    less $file | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc_5mC_quantile.txt
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc_5mC_quantile.txt
done

# clustering
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/
cd $dirIn
less /home/lli/hg19/CG.strand | awk '{if(NR%2){print $2}}' | sort > x
less /home/lli/hg19/CGI.forProfiles.BED | awk '{print $4}' | sort > a
header="ID"
for file in *combine.5mC.CpG; do
    sample=$(echo $file | sed -e 's/.combine.5mC.CpG//g')
    header=$header" "$sample
    echo -e $sample
    less $file | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1 | join x - > y
    mv y x
    less $file | $BEDTOOLS/intersectBed -a /home/lli/hg19/CGI.forProfiles.BED -b stdin -wa -wb | awk '{t[$4]=t[$4]+$8; c[$4]=c[$4]+$9} END{for(i in t){if(t[i]+c[i]>0){print i"\t"c[i]/(c[i]+t[i])}}}' | sort -k1,1 | join a - > b
    mv b a
done
echo -e $header | cat - x > matrix_genome.5mC
echo -e $header | cat - a > matrix_CGI.5mC
rm x a

# DMRs 
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
enhancer=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/FindER/H3K27ac/H3K27ac.union.bed
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/
dirOut=$dirIn/DMR/
mkdir -p $dirOut/
echo -e "sample\tp-value\tdelta\tm\ttotal\thyper\thypo" > $dirOut/DM.summary.stats
echo -e "sample\tsize\tcut\tmedian_length\tmedian_N_CpG\ttotal\thyper\thypo" > $dirOut/DMR.summary.stats
pth=0.05; delta=0.25; m=0.5; cov=3; size=500; cut=3
cd $dirIn
for file1 in *vitc.combine.5mC.CpG; do
    lib1=$(echo $file1 | sed 's/.combine.5mC.CpG//g')
    lib2=$(echo $lib1 | sed -e 's/vitc/control/g')
    name=$lib1'_'$lib2
    echo $name
    /home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut -f1 $lib1.combine.5mC.CpG -f2 $lib2.combine.5mC.CpG -n $name -p $pth -d $delta -m $m -c $cov
    /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
done
lib1=NHAR_control; lib2=NHA_control
name=$lib1'_'$lib2
echo $name
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut -f1 $lib1.combine.5mC.CpG -f2 $lib2.combine.5mC.CpG -n $name -p $pth -d $delta -m $m -c $cov
/home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
/home/lli/HirstLab/Pipeline/shell/DMR.intersect.sh -d $dirOut -r $enhancer -n "enhancer"

