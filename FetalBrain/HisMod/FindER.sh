#!/bin/sh

dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
# index bam files
cd $dirIn
for file in *.bam
do
    echo "Indexing" $file
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools index $file
done

> /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FindER.summary
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
        /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/bin/Solexa_Java/FindER.jar -i $bam -r $reg -o $dirOut -v -m $map -minER $minER -info -cs > $dirOut/FindER_scan.$bam.log
    done
done
### 1000bp minER supported by 500 and 250bp
dirER=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/
mkdir -p $dirER 
cd $dirOut
> $dirER/H3K4me1.combine.summary
for file in *.bam.log
do
    name=$(echo $file | sed -e 's/.bam.log//g')
    echo "Integrating" $name;
    less $name.minER_250.FDR_0.05.bed | awk '{print $0"\tL"}' > x
    cat x $name.minER_500.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | sort -k1,1 -k2,2n | uniq > y;
    cat y $name.minER_1000.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > $dirER/$name.multi.bed;
    N250=`wc -l $name.minER_250.FDR_0.05.bed | cut -d' ' -f 1`
    N500=`wc -l $name.minER_500.FDR_0.05.bed | cut -d' ' -f 1`
    N1000=`wc -l $name.minER_1000.FDR_0.05.bed | cut -d' ' -f 1`
    Npeak=`wc -l $dirER/$name.multi.bed | cut -d' ' -f 1`
    Nbase=`less $dirER/$name.multi.bed | awk '{s=s+$3-$2}END{print s}'`
    echo -e "$name\tH3K4me1\t$Npeak\t$Nbase" >> /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FindER.summary
    echo -e "$name\t$N250\t$N500\t$N1000\t$Npeak" >> $dirER/H3K4me1.combine.summary
done
rm -rf x y

# FindER for H3K4me3
## with three sizes of minimal enriched regions 200, 300bp
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/map/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/H3K4me3/
mkdir -p $dirOut
cd $dirIn
for bam in *.H3K4me3.*.bam
do
    for minER in 200 300 
    do
        echo "Processing FindER on" $bam $minER
        /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/bin/Solexa_Java/FindER.jar -i $bam -r $reg -o $dirOut -v -m $map -minER $minER -info -cs > $dirOut/FindER_scan.$bam.log
    done
done
### summary
dirER=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me3/
mkdir -p $dirER 
cd $dirOut
> $dirER/H3K4me3.combine.summary
for file in *.bam.log
do
    name=$(echo $file | sed -e 's/.bam.log//g')
    cp $name.minER_200.FDR_0.05.bed $dirER
    N200=`wc -l $name.minER_200.FDR_0.05.bed | cut -d' ' -f 1`
    N300=`wc -l $name.minER_300.FDR_0.05.bed | cut -d' ' -f 1`
    Nbase=`less $name.minER_200.FDR_0.05.bed | awk '{s=s+$3-$2}END{print s}'`
    echo -e "$name\tH3K4me3\t$N200\t$Nbase" >> /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FindER.summary
    echo -e "$name\t$N200\t$N300" >> $dirER/H3K4me3.combine.summary
done


# FindER for H3K9me3
## with three sizes of minimal enriched regions 250, 500, 1000bp
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/map/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/H3K9me3/
mkdir -p $dirOut
cd $dirIn
for bam in *.H3K9me3.*.bam
do
    for minER in 250 500 1000
    do
        echo "Processing FindER on" $bam $minER
        /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/bin/Solexa_Java/FindER.jar -i $bam -r $reg -o $dirOut -v -m $map -minER $minER -info -cs > $dirOut/FindER_scan.$bam.log
    done
done
### 1000bp minER supported by 500 and 250bp
dirER=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K9me3/
mkdir -p $dirER 
cd $dirOut
> $dirER/H3K9me3.combine.summary
for file in *.bam.log
do
    name=$(echo $file | sed -e 's/.bam.log//g')
    echo "Integrating" $name;
    less $name.minER_250.FDR_0.05.bed | awk '{print $0"\tL"}' > x
    cat x $name.minER_500.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | sort -k1,1 -k2,2n | uniq > y;
    cat y $name.minER_1000.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > $dirER/$name.multi.bed;
    N250=`wc -l $name.minER_250.FDR_0.05.bed | cut -d' ' -f 1`
    N500=`wc -l $name.minER_500.FDR_0.05.bed | cut -d' ' -f 1`
    N1000=`wc -l $name.minER_1000.FDR_0.05.bed | cut -d' ' -f 1`
    Npeak=`wc -l $dirER/$name.multi.bed | cut -d' ' -f 1`
    Nbase=`less $dirER/$name.multi.bed | awk '{s=s+$3-$2}END{print s}'`
    echo -e "$name\tH3K9me3\t$Npeak\t$Nbase" >> /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FindER.summary
    echo -e "$name\t$N250\t$N500\t$N1000\t$Npeak" >> $dirER/H3K9me3.combine.summary
done
rm -rf x y

# FindER for H3K27me3
## with three sizes of minimal enriched regions 250, 500, 1000bp
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/map/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/H3K27me3/
mkdir -p $dirOut
cd $dirIn
for bam in *.H3K27me3.*.bam
do
    for minER in 250 500 1000
    do
        echo "Processing FindER on" $bam $minER
        /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/bin/Solexa_Java/FindER.jar -i $bam -r $reg -o $dirOut -v -m $map -minER $minER -info -cs > $dirOut/FindER_scan.$bam.log
    done
done
### 1000bp minER supported by 500 and 250bp
dirER=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K27me3/
mkdir -p $dirER 
cd $dirOut
> $dirER/H3K27me3.combine.summary
for file in *.bam.log
do
    name=$(echo $file | sed -e 's/.bam.log//g')
    echo "Integrating" $name;
    less $name.minER_250.FDR_0.05.bed | awk '{print $0"\tL"}' > x
    cat x $name.minER_500.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | sort -k1,1 -k2,2n | uniq > y;
    cat y $name.minER_1000.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > $dirER/$name.multi.bed;
    N250=`wc -l $name.minER_250.FDR_0.05.bed | cut -d' ' -f 1`
    N500=`wc -l $name.minER_500.FDR_0.05.bed | cut -d' ' -f 1`
    N1000=`wc -l $name.minER_1000.FDR_0.05.bed | cut -d' ' -f 1`
    Npeak=`wc -l $dirER/$name.multi.bed | cut -d' ' -f 1`
    Nbase=`less $dirER/$name.multi.bed | awk '{s=s+$3-$2}END{print s}'`
    echo -e "$name\tH3K27me3\t$Npeak\t$Nbase" >> /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FindER.summary
    echo -e "$name\t$N250\t$N500\t$N1000\t$Npeak" >> $dirER/H3K27me3.combine.summary
done
rm -rf x y

# FindER for H3K36me3
## with three sizes of minimal enriched regions 250, 500, 1000bp
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/map/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/H3K36me3/
mkdir -p $dirOut
cd $dirIn
for bam in *.H3K36me3.*.bam
do
    for minER in 250 500 1000
    do
        echo "Processing FindER on" $bam $minER
        /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/bin/Solexa_Java/FindER.jar -i $bam -r $reg -o $dirOut -v -m $map -minER $minER -info -cs > $dirOut/FindER_scan.$bam.log
    done
done
### 1000bp minER supported by 500 and 250bp
dirER=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K36me3/
mkdir -p $dirER 
cd $dirOut
> $dirER/H3K36me3.combine.summary
for file in *.bam.log
do
    name=$(echo $file | sed -e 's/.bam.log//g')
    echo "Integrating" $name;
    less $name.minER_250.FDR_0.05.bed | awk '{print $0"\tL"}' > x
    cat x $name.minER_500.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | sort -k1,1 -k2,2n | uniq > y;
    cat y $name.minER_1000.FDR_0.05.bed | sort -k1,1 -k2,2n -k4,4 | awk '{if($4!="L") {if(chrt==$1 && st<=$3 && et>=$2){print $1"\t"$2"\t"$3"\tL"} else {chr=$1; s=$2; e=$3}} else {if($1==chr && $2<=e && $3>=s){print chr"\t"s"\t"e"\tL"}; chrt=$1; st=$2; et=$3;}}' | cut -f1-3 | sort -k1,1 -k2,2n | uniq > $dirER/$name.multi.bed;
    N250=`wc -l $name.minER_250.FDR_0.05.bed | cut -d' ' -f 1`
    N500=`wc -l $name.minER_500.FDR_0.05.bed | cut -d' ' -f 1`
    N1000=`wc -l $name.minER_1000.FDR_0.05.bed | cut -d' ' -f 1`
    Npeak=`wc -l $dirER/$name.multi.bed | cut -d' ' -f 1`
    Nbase=`less $dirER/$name.multi.bed | awk '{s=s+$3-$2}END{print s}'`
    echo -e "$name\tH3K36me3\t$Npeak\t$Nbase" >> /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FindER.summary
    echo -e "$name\t$N250\t$N500\t$N1000\t$Npeak" >> $dirER/H3K36me3.combine.summary
done
rm -rf x y

# FindER for Input
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/map/
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/FindER/Input/
mkdir -p $dirOut
cd $dirIn
for bam in *.Input.*.bam
do
    echo "Processing FindER on" $bam
    /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx10G /home/mbilenky/bin/Solexa_Java/FindER.jar -i $bam -r $reg -o $dirOut -v -m $map -bin 10 -info -cs > $dirOut/FindER_scan.$bam.log
done
### summary
dirER=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/Input/
mkdir -p $dirER 
cd $dirOut
> $dirER/Input.combine.summary
for file in *.bam.log
do
    name=$(echo $file | sed -e 's/.bam.log//g')
    cp $name.bin_10.minER_200.FDR_0.05.bed $dirER
    Npeak=`wc -l $name.bin_10.minER_200.FDR_0.05.bed | cut -d' ' -f 1`
    Nbase=`less $name.bin_10.minER_200.FDR_0.05.bed | awk '{s=s+$3-$2}END{print s}'`
    echo -e "$name\tInput\t$Npeak\t$Nbase" >> /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FindER.summary
    echo -e "$name\t$Npeak" >> $dirER/Input.combine.summary
done

# closest genes for enhancers (K4me1): exclude overlapping genes or not? 
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/
cd $dirOut
for file in *.multi.bed
do
    name=$(echo $file | sed -e 's/.multi.bed//g')
    name=$(echo $name | sed -e 's/FindER_scan.//g')
    echo "Processing" $name
    # including overlapping genes
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/closestBed -a $file -b /home/lli/hg19/hg19v65_genes.bed -d | awk '$7 ~ /protein_coding/ {gsub("_protein_coding", "", $7); print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$7"\t"$8}' > $dirOut/$name.closest.gene.pc
    # excluding: nearest non-overlapping genes
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/closestBed -a $file -b /home/lli/hg19/hg19v65_genes.bed -d -io | awk '$7 ~ /protein_coding/ {gsub("_protein_coding", "", $7); print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$7"\t"$8}' > $dirOut/$name.closest.gene.pc.io
done 
