#!/bin/sh

dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
# index bam files
cd $dirIn
for file in *.bam
do
    echo "Processing" $file
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools index $file
done

# FindER for H3K4me1
## with three sizes of minimal enriched regions 250, 500, 1000bp
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/map/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/H3K4me1/
mkdir -p $dirOut
cd $dirIn
for bam in *.H3K4me1.*.bam
do
    for minER in 250 500 1000
    do
        echo "Processing FindER on" $bam $minER
        /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/ew/Solexa_Java/jarsToDeploy/FindER2.jar -i $bam -r $reg -o $dirOut -v -m $map -minER $minER -info -cs > $dirOut/FindER_scan.$bam.log
    done
done
### 1000bp minER supported by 500 and 250bp
dirER=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/
mkdir -p $dirER 
cd $dirOut
for file in *.bam.log
do
    name=$(echo $file | sed -e 's/.bam.log//g')
    echo "Integrating" $name;
    less $name.minER_250.FDR_0.05.bed | awk '{print $0"\tL"}' > x
    cat x $name.minER_500.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | sort -k1,1 -k2,2n | uniq > y;
    cat y $name.minER_1000.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > $dirER/$name.multi.bed;
done
rm -rf x y

# FindER for H3K4me3
## with three sizes of minimal enriched regions 200, 300bp
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/smap/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/H3K4me3/
mkdir -p $dirOut
cd $dirIn
for bam in *.H3K4me3.*.bam
do
    for minER in 200 300 
    do
        echo "Processing FindER on" $bam $minER
        /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/ew/Solexa_Java/jarsToDeploy/FindER2.jar -i $bam -r $reg -o $dirOut -v -m $map -minER $minER -info -cs > $dirOut/FindER_scan.$bam.log
    done
done

# FindER for H3K9me3
## with three sizes of minimal enriched regions 250, 500, 1000bp
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/smap/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/H3K9me3/
mkdir -p $dirOut
cd $dirIn
for bam in *.H3K9me3.*.bam
do
    for minER in 250 500 1000
    do
        echo "Processing FindER on" $bam $minER
        /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/ew/Solexa_Java/jarsToDeploy/FindER2.jar -i $bam -r $reg -o $dirOut -v -m $map -minER $minER -info -cs > $dirOut/FindER_scan.$bam.log
    done
done
### 1000bp minER supported by 500 and 250bp
dirER=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K9me3/
mkdir -p $dirER 
cd $dirOut
for file in *.bam.log
do
    name=$(echo $file | sed -e 's/.bam.log//g')
    echo "Integrating" $name;
    less $name.minER_250.FDR_0.05.bed | awk '{print $0"\tL"}' > x
    cat x $name.minER_500.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | sort -k1,1 -k2,2n | uniq > y;
    cat y $name.minER_1000.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > $dirER/$name.multi.bed;
done
rm -rf x y

# FindER for H3K27me3
## with three sizes of minimal enriched regions 250, 500, 1000bp
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/smap/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/H3K27me3/
mkdir -p $dirOut
cd $dirIn
for bam in *.H3K27me3.*.bam
do
    for minER in 250 500 1000
    do
        echo "Processing FindER on" $bam $minER
        /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/ew/Solexa_Java/jarsToDeploy/FindER2.jar -i $bam -r $reg -o $dirOut -v -m $map -minER $minER -info -cs > $dirOut/FindER_scan.$bam.log
    done
done
### 1000bp minER supported by 500 and 250bp
dirER=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K27me3/
mkdir -p $dirER 
cd $dirOut
for file in *.bam.log
do
    name=$(echo $file | sed -e 's/.bam.log//g')
    echo "Integrating" $name;
    less $name.minER_250.FDR_0.05.bed | awk '{print $0"\tL"}' > x
    cat x $name.minER_500.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | sort -k1,1 -k2,2n | uniq > y;
    cat y $name.minER_1000.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > $dirER/$name.multi.bed;
done
rm -rf x y

# FindER for H3K36me3
## with three sizes of minimal enriched regions 250, 500, 1000bp
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/smap/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/H3K36me3/
mkdir -p $dirOut
cd $dirIn
for bam in *.H3K36me3.*.bam
do
    for minER in 250 500 1000
    do
        echo "Processing FindER on" $bam $minER
        /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/ew/Solexa_Java/jarsToDeploy/FindER2.jar -i $bam -r $reg -o $dirOut -v -m $map -minER $minER -info -cs > $dirOut/FindER_scan.$bam.log
    done
done
### 1000bp minER supported by 500 and 250bp
dirER=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K36me3/
mkdir -p $dirER 
cd $dirOut
for file in *.bam.log
do
    name=$(echo $file | sed -e 's/.bam.log//g')
    echo "Integrating" $name;
    less $name.minER_250.FDR_0.05.bed | awk '{print $0"\tL"}' > x
    cat x $name.minER_500.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | sort -k1,1 -k2,2n | uniq > y;
    cat y $name.minER_1000.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > $dirER/$name.multi.bed;
done
rm -rf x y

# FindER for Input
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/smap/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/Input/
mkdir -p $dirOut
cd $dirIn
for bam in *.Input.*.bam
do
    echo "Processing FindER on" $bam
    /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/ew/Solexa_Java/jarsToDeploy/FindER2.jar -i $bam -r $reg -o $dirOut -v -m $map -bin 10 -info -cs > $dirOut/FindER_scan.$bam.log
done

