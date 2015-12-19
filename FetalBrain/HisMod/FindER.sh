#!/bin/sh

dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/bam/
# index bam files
cd $dirIn
for file in *.bam
do
    echo "Indexing" $file
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools index $file
done

############################################
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
    cp $name.minER_300.FDR_0.05.bed $dirER
    N200=`wc -l $name.minER_200.FDR_0.05.bed | cut -d' ' -f 1`
    N300=`wc -l $name.minER_300.FDR_0.05.bed | cut -d' ' -f 1`
    Nbase=`less $name.minER_300.FDR_0.05.bed | awk '{s=s+$3-$2}END{print s}'`
    echo -e "$name\tH3K4me3\t$N300\t$Nbase" >> /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/FindER.summary
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

############################################
# intersect
## H3K4me3 with promoter
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me3/
cd $dirOut
for file in *.minER_300.FDR_0.05.bed
do
    name=$(echo $file | sed -e 's/.minER_300.FDR_0.05.bed//g')
    name=$(echo $name | sed -e 's/FindER_scan.//g')
    echo "Processing" $name
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $file -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$7}' > $dirOut/$name.promoter.bed
done 
## H3K27me3 with promoter
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K27me3/
cd $dirOut
for file in *.multi.bed
do
    name=$(echo $file | sed -e 's/.multi.bed//g')
    name=$(echo $name | sed -e 's/FindER_scan.//g')
    echo "Processing" $name
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $file -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$7}' > $dirOut/$name.promoter.bed
done 
## H3K36me3 with genebody
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K36me3/
cd $dirOut
for file in *.multi.bed
do
    name=$(echo $file | sed -e 's/.multi.bed//g')
    name=$(echo $name | sed -e 's/FindER_scan.//g')
    echo "Processing" $name
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $file -b /home/lli/hg19/hg19v65_genes.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$7}' > $dirOut/$name.genebody.bed
done 

# closest genes for enhancers (K4me1): exclude overlapping genes or not? 
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/
mkdir -p $dirOut/closest_gene/
cd $dirOut
for file in *.multi.bed
do
    name=$(echo $file | sed -e 's/.multi.bed//g')
    name=$(echo $name | sed -e 's/FindER_scan.//g')
    echo "Processing" $name
    # including overlapping genes
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/closestBed -a $file -b /home/lli/hg19/hg19v65_genes.bed -d | awk '$7 ~ /protein_coding/ {gsub("_protein_coding", "", $7); print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$7"\t"$8}' > $dirOut/closest_gene/$name.closest.gene.pc
    # excluding: nearest non-overlapping genes
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/closestBed -a $file -b /home/lli/hg19/hg19v65_genes.bed -d -io | awk '$7 ~ /protein_coding/ {gsub("_protein_coding", "", $7); print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$7"\t"$8}' > $dirOut/closest_gene/$name.closest.gene.pc.io
done 

############################################
# core enhancers: overlapping all NPC enhancers
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/
cd $dirOut
mkdir -p $dirOut/core/
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed -b $dirOut/FindER_scan.A03477.H3K4me1.GE02.multi.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/H3K4me1.Cortex02.GE02.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/H3K4me1.Cortex02.GE02.bed -b $dirOut/FindER_scan.A19303.H3K4me1.GE04.multi.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/H3K4me1.Cortex02.GE02.GE04.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/H3K4me1.Cortex02.GE02.GE04.bed -b $dirOut/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/H3K4me1.Cortex02.GE02.GE04.Cortex01.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/H3K4me1.Cortex02.GE02.GE04.Cortex01.bed -b $dirOut/FindER_scan.A03275.H3K4me1.GE01.multi.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/core/core_enhancers.bed
rm H3K4me1.Cortex02.GE02*.bed
## Intersect with GWAS sites
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/core/core_enhancers.bed -b /home/lli/hg19/gwasCatalog_July2014.bed -wa -wb > $dirOut/core/core_enhancers.GWAS.txt
less $dirOut/core/core_enhancers.GWAS.txt | awk -F"\t" '{print $5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$4}' > $dirOut/core/core_enhancers.GWAS.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/core/core_enhancers.GWAS.bed -b /projects/mbilenky/REMC/breast/ENCODE/TFs/wgEncodeRegTfbsClusteredV3.bed.gz -wa -wb > $dirOut/core/core_enhancers.GWAS.TFBS.txt
## Homer for TFBS motifs
PATH=$PATH:/home/acarles/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
/home/acarles/homer/bin/findMotifsGenome.pl $dirOut/core/core_enhancers.bed hg19 $dirOut/core/homer/ -size 200 -len 8
mkdir -p $dirOut/core/homer/annotate/ 
less $dirOut/core/homer/homer_core_enhancer_top.txt | awk 'NR>=2 {system("/home/acarles/homer/bin/annotatePeaks.pl ""'$dirOut'""/core/core_enhancers.bed hg19 -m ""'$dirOut'""/core/homer/"$10" > ""'$dirOut'""/core/homer/annotate/"$11".annotate")}'
## Intersect with WGBS UMRs
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/core/UMR/
mkdir -p $dirOut
enhancer=`wc -l $dirOut/../core_enhancers.bed | cut -d' ' -f 1`
echo -e "Sample\tUMR\tNo.enhancers\tNo.UMR\tNo.enhancerUMR" > $dirOut/WGBS_UMR_enhancers.summary
dirIn=/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/../core_enhancers.bed -b $dirIn/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed -wa -wb > $dirOut/neurosphere02hyper_enhancer.bed
neurosphere02hyper=`wc -l $dirIn/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed | cut -d' ' -f 1`
neurosphere02hyper_enhancer=`wc -l $dirOut/neurosphere02hyper_enhancer.bed | cut -d' ' -f 1`
echo -e "neurosphere02\thyper\t$enhancer\t$neurosphere02hyper\t$neurosphere02hyper_enhancer" >> $dirOut/WGBS_UMR_enhancers.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/../core_enhancers.bed -b $dirIn/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed -wa -wb > $dirOut/neurosphere02hypo_enhancer.bed
neurosphere02hypo=`wc -l $dirIn/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed | cut -d' ' -f 1`
neurosphere02hypo_enhancer=`wc -l $dirOut/neurosphere02hypo_enhancer.bed | cut -d' ' -f 1`
echo -e "neurosphere02\thypo\t$enhancer\t$neurosphere02hypo\t$neurosphere02hypo_enhancer" >> $dirOut/WGBS_UMR_enhancers.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/../core_enhancers.bed -b $dirIn/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed -wa -wb > $dirOut/neurosphere04hyper_enhancer.bed
neurosphere04hyper=`wc -l $dirIn/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed | cut -d' ' -f 1`
neurosphere04hyper_enhancer=`wc -l $dirOut/neurosphere04hyper_enhancer.bed | cut -d' ' -f 1`
echo -e "neurosphere04\thyper\t$enhancer\t$neurosphere04hyper\t$neurosphere04hyper_enhancer" >> $dirOut/WGBS_UMR_enhancers.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/../core_enhancers.bed -b $dirIn/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed -wa -wb > $dirOut/neurosphere04hypo_enhancer.bed
neurosphere04hypo=`wc -l $dirIn/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed | cut -d' ' -f 1`
neurosphere04hypo_enhancer=`wc -l $dirOut/neurosphere04hypo_enhancer.bed | cut -d' ' -f 1`
echo -e "neurosphere04\thypo\t$enhancer\t$neurosphere04hypo\t$neurosphere04hypo_enhancer" >> $dirOut/WGBS_UMR_enhancers.summary
dirIn=/projects/epigenomics/users/lli/FetalBrain/GW/DMR/
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/../core_enhancers.bed -b $dirIn/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed -wa -wb > $dirOut/GW_Cortex_hyper_enhancer.bed
GW_Cortex_hyper=`wc -l $dirIn/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed | cut -d' ' -f 1`
GW_Cortex_hyper_enhancer=`wc -l $dirOut/GW_Cortex_hyper_enhancer.bed | cut -d' ' -f 1`
echo -e "GW_Cortex\thyper\t$enhancer\t$GW_Cortex_hyper\t$GW_Cortex_hyper_enhancer" >> $dirOut/WGBS_UMR_enhancers.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/../core_enhancers.bed -b $dirIn/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed -wa -wb > $dirOut/GW_Cortex_hypo_enhancer.bed
GW_Cortex_hypo=`wc -l $dirIn/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed | cut -d' ' -f 1`
GW_Cortex_hypo_enhancer=`wc -l $dirOut/GW_Cortex_hypo_enhancer.bed | cut -d' ' -f 1`
echo -e "GW_Cortex\thypo\t$enhancer\t$GW_Cortex_hypo\t$GW_Cortex_hypo_enhancer" >> $dirOut/WGBS_UMR_enhancers.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/../core_enhancers.bed -b $dirIn/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed -wa -wb > $dirOut/GW_GE_hyper_enhancer.bed
GW_GE_hyper=`wc -l $dirIn/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed | cut -d' ' -f 1`
GW_GE_hyper_enhancer=`wc -l $dirOut/GW_GE_hyper_enhancer.bed | cut -d' ' -f 1`
echo -e "GW_GE\thyper\t$enhancer\t$GW_GE_hyper\t$GW_GE_hyper_enhancer" >> $dirOut/WGBS_UMR_enhancers.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/../core_enhancers.bed -b $dirIn/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed -wa -wb > $dirOut/GW_GE_hypo_enhancer.bed
GW_GE_hypo=`wc -l $dirIn/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed | cut -d' ' -f 1`
GW_GE_hypo_enhancer=`wc -l $dirOut/GW_GE_hypo_enhancer.bed | cut -d' ' -f 1`
echo -e "GW_GE\thypo\t$enhancer\t$GW_GE_hypo\t$GW_GE_hypo_enhancer" >> $dirOut/WGBS_UMR_enhancers.summary
### input files for GREAT
cd $dirOut
for file in *.bed
do
    less $file | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $file.GREAT.bed
done
### enhancer UMRs - enrichment
dirOut=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/core/UMR/CpG/
mkdir -p $dirOut
cd $dirOut
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/core/core_enhancers.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/core_enhancers.CpG.bed
enhancer=`wc -l $dirOut/core_enhancers.CpG.bed| cut -d' ' -f 1`
all_total=28217448 # wc -l /home/lli/hg19/CG.BED
echo -e "Sample\tUMR\tNo.total\tNo.enhancers\tNo.UMR\tNo.enhancerUMR\tenrich" > $dirOut/WGBS_UMR_enhancers_enrich.summary
dirIn=/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/CpG/
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/core_enhancers.CpG.bed -b $dirIn/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.CpG.bed -wa -wb > $dirOut/neurosphere02hyper_enhancer.CpG.bed
neurosphere02hyper=`wc -l $dirIn/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.CpG.bed | cut -d' ' -f 1`
neurosphere02hyper_enhancer=`wc -l $dirOut/neurosphere02hyper_enhancer.CpG.bed | cut -d' ' -f 1`
echo -e "neurosphere02\thyper\t$all_total\t$enhancer\t$neurosphere02hyper\t$neurosphere02hyper_enhancer" | awk '{print $0"\t"($6/$5)/($4/$3)}' >> $dirOut/WGBS_UMR_enhancers_enrich.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/core_enhancers.CpG.bed -b $dirIn/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.CpG.bed -wa -wb > $dirOut/neurosphere02hypo_enhancer.CpG.bed
neurosphere02hypo=`wc -l $dirIn/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.CpG.bed | cut -d' ' -f 1`
neurosphere02hypo_enhancer=`wc -l $dirOut/neurosphere02hypo_enhancer.CpG.bed | cut -d' ' -f 1`
echo -e "neurosphere02\thypo\t$all_total\t$enhancer\t$neurosphere02hypo\t$neurosphere02hypo_enhancer" | awk '{print $0"\t"($6/$5)/($4/$3)}' >> $dirOut/WGBS_UMR_enhancers_enrich.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/core_enhancers.CpG.bed -b $dirIn/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.CpG.bed -wa -wb > $dirOut/neurosphere04hyper_enhancer.CpG.bed
neurosphere04hyper=`wc -l $dirIn/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.CpG.bed | cut -d' ' -f 1`
neurosphere04hyper_enhancer=`wc -l $dirOut/neurosphere04hyper_enhancer.CpG.bed | cut -d' ' -f 1`
echo -e "neurosphere04\thyper\t$all_total\t$enhancer\t$neurosphere04hyper\t$neurosphere04hyper_enhancer" | awk '{print $0"\t"($6/$5)/($4/$3)}' >> $dirOut/WGBS_UMR_enhancers_enrich.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/core_enhancers.CpG.bed -b $dirIn/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.CpG.bed -wa -wb > $dirOut/neurosphere04hypo_enhancer.CpG.bed
neurosphere04hypo=`wc -l $dirIn/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.CpG.bed | cut -d' ' -f 1`
neurosphere04hypo_enhancer=`wc -l $dirOut/neurosphere04hypo_enhancer.CpG.bed | cut -d' ' -f 1`
echo -e "neurosphere04\thypo\t$all_total\t$enhancer\t$neurosphere04hypo\t$neurosphere04hypo_enhancer" | awk '{print $0"\t"($6/$5)/($4/$3)}' >> $dirOut/WGBS_UMR_enhancers_enrich.summary
dirIn=/projects/epigenomics/users/lli/FetalBrain/GW/DMR/CpG/
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/core_enhancers.CpG.bed -b $dirIn/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.CpG.bed -wa -wb > $dirOut/GW_Cortex_hyper_enhancer.CpG.bed
GW_Cortex_hyper=`wc -l $dirIn/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.CpG.bed | cut -d' ' -f 1`
GW_Cortex_hyper_enhancer=`wc -l $dirOut/GW_Cortex_hyper_enhancer.CpG.bed | cut -d' ' -f 1`
echo -e "GW_Cortex\thyper\t$all_total\t$enhancer\t$GW_Cortex_hyper\t$GW_Cortex_hyper_enhancer" | awk '{print $0"\t"($6/$5)/($4/$3)}' >> $dirOut/WGBS_UMR_enhancers_enrich.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/core_enhancers.CpG.bed -b $dirIn/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.CpG.bed -wa -wb > $dirOut/GW_Cortex_hypo_enhancer.CpG.bed
GW_Cortex_hypo=`wc -l $dirIn/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.CpG.bed | cut -d' ' -f 1`
GW_Cortex_hypo_enhancer=`wc -l $dirOut/GW_Cortex_hypo_enhancer.CpG.bed | cut -d' ' -f 1`
echo -e "GW_Cortex\thypo\t$all_total\t$enhancer\t$GW_Cortex_hypo\t$GW_Cortex_hypo_enhancer" | awk '{print $0"\t"($6/$5)/($4/$3)}' >> $dirOut/WGBS_UMR_enhancers_enrich.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/core_enhancers.CpG.bed -b $dirIn/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.CpG.bed -wa -wb > $dirOut/GW_GE_hyper_enhancer.CpG.bed
GW_GE_hyper=`wc -l $dirIn/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.CpG.bed | cut -d' ' -f 1`
GW_GE_hyper_enhancer=`wc -l $dirOut/GW_GE_hyper_enhancer.CpG.bed | cut -d' ' -f 1`
echo -e "GW_GE\thyper\t$all_total\t$enhancer\t$GW_GE_hyper\t$GW_GE_hyper_enhancer" | awk '{print $0"\t"($6/$5)/($4/$3)}' >> $dirOut/WGBS_UMR_enhancers_enrich.summary
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/core_enhancers.CpG.bed -b $dirIn/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.CpG.bed -wa -wb > $dirOut/GW_GE_hypo_enhancer.CpG.bed
GW_GE_hypo=`wc -l $dirIn/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.CpG.bed | cut -d' ' -f 1`
GW_GE_hypo_enhancer=`wc -l $dirOut/GW_GE_hypo_enhancer.CpG.bed | cut -d' ' -f 1`
echo -e "GW_GE\thypo\t$all_total\t$enhancer\t$GW_GE_hypo\t$GW_GE_hypo_enhancer" | awk '{print $0"\t"($6/$5)/($4/$3)}' >> $dirOut/WGBS_UMR_enhancers_enrich.summary
### Homer for TFBS motifs
PATH=$PATH:/home/acarles/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/core/UMR/'
mkdir -p $dirIn/homer/
cd $dirIn/
for file in *.GREAT.bed
do
    name=$(echo $file | sed -e 's/.GREAT//g' | sed -e 's/.bed//g')
    echo "Processing "$name
    dirOut=$dirIn/homer/$name/
    mkdir -p $dirOut
    /home/acarles/homer/bin/findMotifsGenome.pl $file hg19 $dirOut -size 200 -len 8 
done

############################################
# distance between closest peaks in pairwise comparisons
dirIn=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/
dirOut=$dirIn/dis/
BEDTOOLS=/projects/epigenomics/software/bedtools-2.23.0/bin/
mkdir -p $dirOut
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03485.H3K4me1.Brain01.multi.bed -b $dirIn/FindER_scan.A03493.H3K4me1.Brain02.multi.bed -d | awk '{if($7!=0)print $0}' > $dirOut/H3K4me1_Brain01_Brain02.closest
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03485.H3K4me1.Brain01.multi.bed -b $dirIn/FindER_scan.A03493.H3K4me1.Brain02.multi.bed -d | awk '{if($7==0)print $1"\t"$2"\t"$3}' | $BEDTOOLS/intersectBed -a stdin -b $dirIn/FindER_scan.A03493.H3K4me1.Brain02.multi.bed -wo | awk '{dis=-$7; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"dis}' > $dirOut/H3K4me1_Brain01_Brain02.overlap
cat $dirOut/H3K4me1_Brain01_Brain02.closest $dirOut/H3K4me1_Brain01_Brain02.overlap > $dirOut/H3K4me1_Brain01_Brain02.dis
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed -b $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed -d | awk '{if($7!=0)print $0}' > $dirOut/H3K4me1_Cortex01_Cortex02.closest
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed -b $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed -d | awk '{if($7==0)print $1"\t"$2"\t"$3}' | $BEDTOOLS/intersectBed -a stdin -b $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed -wo | awk '{dis=-$7; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"dis}' > $dirOut/H3K4me1_Cortex01_Cortex02.overlap
cat $dirOut/H3K4me1_Cortex01_Cortex02.closest $dirOut/H3K4me1_Cortex01_Cortex02.overlap > $dirOut/H3K4me1_Cortex01_Cortex02.dis
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -b $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -d | awk '{if($7!=0)print $0}' > $dirOut/H3K4me1_GE01_GE02.closest
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -b $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -d | awk '{if($7==0)print $1"\t"$2"\t"$3}' | $BEDTOOLS/intersectBed -a stdin -b $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -wo | awk '{dis=-$7; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"dis}' > $dirOut/H3K4me1_GE01_GE02.overlap
cat $dirOut/H3K4me1_GE01_GE02.closest $dirOut/H3K4me1_GE01_GE02.overlap > $dirOut/H3K4me1_GE01_GE02.dis
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed -b $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -d | awk '{if($7!=0)print $0}' > $dirOut/H3K4me1_Cortex01_GE01.closest
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed -b $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -d | awk '{if($7==0)print $1"\t"$2"\t"$3}' | $BEDTOOLS/intersectBed -a stdin -b $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -wo | awk '{dis=-$7; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"dis}' > $dirOut/H3K4me1_Cortex01_GE01.overlap
cat $dirOut/H3K4me1_Cortex01_GE01.closest $dirOut/H3K4me1_Cortex01_GE01.overlap > $dirOut/H3K4me1_Cortex01_GE01.dis
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed -b $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -d | awk '{if($7!=0)print $0}' > $dirOut/H3K4me1_Cortex02_GE02.closest
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed -b $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -d | awk '{if($7==0)print $1"\t"$2"\t"$3}' | $BEDTOOLS/intersectBed -a stdin -b $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -wo | awk '{dis=-$7; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"dis}' > $dirOut/H3K4me1_Cortex02_GE02.overlap
cat $dirOut/H3K4me1_Cortex02_GE02.closest $dirOut/H3K4me1_Cortex02_GE02.overlap > $dirOut/H3K4me1_Cortex02_GE02.dis
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -b $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed -d | awk '{if($7!=0)print $0}' > $dirOut/H3K4me1_GE01_GE04.closest
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -b $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed -d | awk '{if($7==0)print $1"\t"$2"\t"$3}' | $BEDTOOLS/intersectBed -a stdin -b $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed -wo | awk '{dis=-$7; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"dis}' > $dirOut/H3K4me1_GE01_GE04.overlap
cat $dirOut/H3K4me1_GE01_GE04.closest $dirOut/H3K4me1_GE01_GE04.overlap > $dirOut/H3K4me1_GE01_GE04.dis
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -b $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed -d | awk '{if($7!=0)print $0}' > $dirOut/H3K4me1_GE02_GE04.closest
$BEDTOOLS/closestBed -a $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -b $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed -d | awk '{if($7==0)print $1"\t"$2"\t"$3}' | $BEDTOOLS/intersectBed -a stdin -b $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed -wo | awk '{dis=-$7; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"dis}' > $dirOut/H3K4me1_GE02_GE04.overlap
cat $dirOut/H3K4me1_GE02_GE04.closest $dirOut/H3K4me1_GE02_GE04.overlap > $dirOut/H3K4me1_GE02_GE04.dis

############################################
# Unique enhancers
dirIn='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/'
dirOut=$dirIn/unique/
BEDTOOLS=/projects/epigenomics/software/bedtools-2.23.0/bin/
mkdir -p $dirOut
cd $dirIn
echo -e "Mark\tComparison\tSamples\tSample1\tSample2\tSample1_unique\tSample2_unique" > $dirOut/unique_enhancer.summary
## Between MZ twins
### Brain
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03485.H3K4me1.Brain01.multi.bed -b $dirIn/FindER_scan.A03493.H3K4me1.Brain02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.MZ.Brain01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03493.H3K4me1.Brain02.multi.bed -b $dirIn/FindER_scan.A03485.H3K4me1.Brain01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.MZ.Brain02.bed
s1=`less $dirIn/FindER_scan.A03485.H3K4me1.Brain01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03493.H3K4me1.Brain02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_enhancer.MZ.Brain01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_enhancer.MZ.Brain02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K4me1\tMZ\tBrain01_Brain02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_enhancer.summary
### Cortex
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed -b $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.MZ.Cortex01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed -b $dirIn/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.MZ.Cortex02.bed
s1=`less $dirIn/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_enhancer.MZ.Cortex01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_enhancer.MZ.Cortex02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K4me1\tMZ\tCortex01_Cortex02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_enhancer.summary
### GE
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -b $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.MZ.GE01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -b $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.MZ.GE02.bed
s1=`less $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_enhancer.MZ.GE01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_enhancer.MZ.GE02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K4me1\tMZ\tGE01_GE02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_enhancer.summary
### Intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_enhancer.MZ.Cortex01.bed -b $dirOut/unique_enhancer.MZ.GE01.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.MZ.HuFNSC01.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_enhancer.MZ.Cortex02.bed -b $dirOut/unique_enhancer.MZ.GE02.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.MZ.HuFNSC02.bed
## Between neurospheres
### HuFNSC01
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed -b $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.Neurospheres.Cortex01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -b $dirIn/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.Neurospheres.GE01.bed
s1=`less $dirIn/FindER_scan.A03269.H3K4me1.Cortex01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_enhancer.Neurospheres.Cortex01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_enhancer.Neurospheres.GE01.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K4me1\tNPCs\tCortex01_GE01\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_enhancer.summary
### HuFNSC02
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed -b $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.Neurospheres.Cortex02.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -b $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.Neurospheres.GE02.bed
s1=`less $dirIn/FindER_scan.A03281.H3K4me1.Cortex02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_enhancer.Neurospheres.Cortex02.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_enhancer.Neurospheres.GE02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K4me1\tNPCs\tCortex02_GE02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_enhancer.summary
### intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_enhancer.Neurospheres.Cortex01.bed -b $dirOut/unique_enhancer.Neurospheres.Cortex02.bed > $dirOut/unique_enhancer.Neurospheres.Cortex.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_enhancer.Neurospheres.GE01.bed -b $dirOut/unique_enhancer.Neurospheres.GE02.bed > $dirOut/unique_enhancer.Neurospheres.GE.bed
## Between GW
### GE01 vs GE04
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -b $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.GW.GW17_01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed -b $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.GW.GW13_01.bed
s1=`less $dirIn/FindER_scan.A03275.H3K4me1.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_enhancer.GW.GW17_01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_enhancer.GW.GW13_01.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K4me1\tGW\tGE01_GE04\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_enhancer.summary
### GE02 vs GE04
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -b $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.GW.GW17_02.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed -b $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_enhancer.GW.GW13_02.bed
s1=`less $dirIn/FindER_scan.A03477.H3K4me1.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A19303.H3K4me1.GE04.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_enhancer.GW.GW17_02.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_enhancer.GW.GW13_02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K4me1\tGW\tGE02_GE04\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_enhancer.summary
### intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_enhancer.GW.GW13_01.bed -b $dirOut/unique_enhancer.GW.GW13_02.bed > $dirOut/unique_enhancer.GW.GW13.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_enhancer.GW.GW17_01.bed -b $dirOut/unique_enhancer.GW.GW17_02.bed > $dirOut/unique_enhancer.GW.GW17.bed
## Intersect with GWAS sites
mkdir -p $dirOut/GWAS/
cd $dirOut/
for file in unique_enhancer.*.bed
do
    name=$(echo $file | sed -e 's/.bed//g')
    echo "Processing "$name
    /projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/$file -b /home/lli/hg19/gwasCatalog_July2014.bed -wa -wb > $dirOut/GWAS/$name.GWAS.bed
done
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/GWAS/GWAS_unique_enhancer_all.bed -b /projects/mbilenky/REMC/breast/ENCODE/TFs/wgEncodeRegTfbsClusteredV3.bed.gz -wa -wb > /projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/GWAS/GWAS_unique_enhancer_all_TFBS.txt
## Homer for TFBS motifs
PATH=$PATH:/home/acarles/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique'
mkdir -p $dirIn/homer/
cd $dirIn/
for file in unique_enhancer.*.bed
do
    name=$(echo $file | sed -e 's/.bed//g' | sed -e 's/unique_enhancer.//g')
    echo "Processing "$name
    dirOut=$dirIn/homer/$name/
    mkdir -p $dirOut
    /home/acarles/homer/bin/findMotifsGenome.pl $file hg19 $dirOut -size 200 -len 8 
done
### GW
PATH=$PATH:/home/acarles/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirOut='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/'
dirMotif='/home/acarles/homer/data/knownTFs/motifs/'
#/home/acarles/homer/bin/annotatePeaks.pl $dirOut/unique_enhancer.GW.GW17.bed hg19 -m $dirMotif/ar-half.motif $dirMotif/ptf1a.motif $dirMotif/pr.motif $dirMotif/foxo1.motif $dirMotif/olig2.motif -mscore > $dirOut/homer/GW/GW17only_motif.annotate
#/home/acarles/homer/bin/findMotifsGenome.pl $dirOut/unique_enhancer.GW.GW17.bed hg19 $dirOut/homer/GW/ -find $dirOut/homer/GW.GW17/knownResults/known21.motif > $dirOut/homer/GW/Olig2.test.txt
#less $dirOut/homer/GW/homer_unique_enhancer_GW_common.txt | awk 'NR==2 {system("/home/acarles/homer/bin/annotatePeaks.pl ""'$dirOut'""/unique_enhancer.GW.GW13.bed hg19 -m ""'$dirOut'""/homer/GW.GW13/"$7" -nmotifs > ""'$dirOut'""/homer/GW/Common_GW13_"$3".annotate")}'
#less $dirOut/homer/GW/homer_unique_enhancer_GW_common.txt | awk 'NR==2 {system("/home/acarles/homer/bin/annotatePeaks.pl ""'$dirOut'""/unique_enhancer.GW.GW17.bed hg19 -m ""'$dirOut'""/homer/GW.GW17/"$12" -nmotifs > ""'$dirOut'""/homer/GW/Common_GW17_"$3".annotate")}'
less $dirOut/homer/GW/homer_unique_enhancer_GW_GW17only.txt | awk 'NR<=6&&NR>=2 {system("/home/acarles/homer/bin/findMotifsGenome.pl ""'$dirOut'""/unique_enhancer.GW.GW17.bed hg19 ""'$dirOut'""/homer/GW/ -find ""'$dirOut'""/homer/GW.GW17/"$4" > ""'$dirOut'""/homer/GW/GW17only_"$2".annotate")}'
cd $dirOut/homer/GW/
for file in GW17only_*.annotate; do
    tf=$(echo -e $file | sed -e 's/GW17only_//g' | sed -e 's/.annotate//g')
    echo -e "Processing "$tf
    less $file | sed -e s/:/"\t"/g | sed -e s/-/"\t"/g | awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"}' | sort -k1,1 -k2,2n -T /projects/epigenomics/temp/ > $file.motif.bed
    /projects/epigenomics/software/bedtools-2.23.0/bin/closestBed -a $file.motif.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -d | awk '{print $0"\t""'$tf'"}' > $file.motif.gene
done
cat GW17only_*.annotate.motif.gene > GW17only_TF_motif.gene

### Neurospheres
PATH=$PATH:/home/acarles/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirOut='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/'
less $dirOut/homer/Neurospheres/homer_unique_enhancer_Neurospheres_common.txt | awk 'NR==3 {system("/home/acarles/homer/bin/annotatePeaks.pl ""'$dirOut'""/unique_enhancer.Neurospheres.Cortex.bed hg19 -m ""'$dirOut'""/homer/Neurospheres.Cortex/"$7" > ""'$dirOut'""/homer/Neurospheres/Common_Cortex_"$3".annotate")}'
less $dirOut/homer/Neurospheres/homer_unique_enhancer_Neurospheres_common.txt | awk 'NR==3 {system("/home/acarles/homer/bin/annotatePeaks.pl ""'$dirOut'""/unique_enhancer.Neurospheres.GE.bed hg19 -m ""'$dirOut'""/homer/Neurospheres.GE/"$12" > ""'$dirOut'""/homer/Neurospheres/Common_GE_"$3".annotate")}'
less $dirOut/homer/Neurospheres/homer_unique_enhancer_Neurospheres_GEonly.txt | awk 'NR==5 {system("/home/acarles/homer/bin/annotatePeaks.pl ""'$dirOut'""/unique_enhancer.Neurospheres.GE.bed hg19 -m ""'$dirOut'""/homer/Neurospheres.GE/"$4" > ""'$dirOut'""/homer/Neurospheres/GEonly_"$2".annotate")}'

## UMR enhancers
dirGW='/projects/epigenomics/users/lli/FetalBrain/GW/DMR/'
dirNPC='/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/'
dirIn='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/'
dirOut=$dirIn'/UMR_enhancer/'
mkdir -p $dirOut
cd $dirOut
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirIn/unique_enhancer.GW.GW13_02.bed -b $dirGW/DMR.HuFNSC02_HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed -wa > $dirOut/UMR_enhancer.GW.GW13_02.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirIn/unique_enhancer.GW.GW17_02.bed -b $dirGW/DMR.HuFNSC02_HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed -wa > $dirOut/UMR_enhancer.GW.GW17_02.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirIn/unique_enhancer.Neurospheres.Cortex02.bed -b $dirNPC/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed -wa > $dirOut/UMR_enhancer.Neurospheres.Cortex02.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirIn/unique_enhancer.Neurospheres.GE02.bed -b $dirNPC/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed -wa > $dirOut/UMR_enhancer.Neurospheres.GE02.bed
for file in UMR_enhancer*.bed; do
    name=$(echo $file | sed -e 's/.bed//g')
    echo "Processing"$name
    /projects/epigenomics/software/bedtools-2.23.0/bin/closestBed -a $file -b /home/lli/hg19/hg19v65_genes.pc.bed -d > $dirOut/$name.closest.gene
done

## Enhancer fractional methylation call
dirr='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique'
dirIn='/projects/epigenomics/users/lli/FetalBrain/WGBS/'
dirOut='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/fractional'
mkdir -p $dirOut
/home/lli/HirstLab/Pipeline/shell/region.mean.sh -i $dirIn -o $dirOut -f A22476.WGBS.NeurospheresGE04.sam.bedGraph -r $dirr/unique_enhancer.GW.GW13.bed -n unique_enhancer.GW.GW13.5mC_GE04 -w 
/home/lli/HirstLab/Pipeline/shell/region.mean.sh -i $dirIn -o $dirOut -f A17784-A13819.WGBS.NeurospheresGE02.sam.bedGraph -r $dirr/unique_enhancer.GW.GW13.bed -n unique_enhancer.GW.GW13.5mC_GE02 -w 
paste $dirOut/unique_enhancer.GW.GW13.5mC_GE02.weighted.mean.bed $dirOut/unique_enhancer.GW.GW13.5mC_GE04.weighted.mean.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10"\t"$5-$10}' > $dirOut/unique_enhancer.GW.GW13.5mC_GE02-04.bed
rm $dirOut/unique_enhancer.GW.GW13.5mC_GE02.weighted.mean.bed $dirOut/unique_enhancer.GW.GW13.5mC_GE04.weighted.mean.bed
/home/lli/HirstLab/Pipeline/shell/region.mean.sh -i $dirIn -o $dirOut -f A22476.WGBS.NeurospheresGE04.sam.bedGraph -r $dirr/unique_enhancer.GW.GW17.bed -n unique_enhancer.GW.GW17.5mC_GE04 -w 
/home/lli/HirstLab/Pipeline/shell/region.mean.sh -i $dirIn -o $dirOut -f A17784-A13819.WGBS.NeurospheresGE02.sam.bedGraph -r $dirr/unique_enhancer.GW.GW17.bed -n unique_enhancer.GW.GW17.5mC_GE02 -w 
paste $dirOut/unique_enhancer.GW.GW17.5mC_GE02.weighted.mean.bed $dirOut/unique_enhancer.GW.GW17.5mC_GE04.weighted.mean.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10"\t"$5-$10}' > $dirOut/unique_enhancer.GW.GW17.5mC_GE02-04.bed
rm $dirOut/unique_enhancer.GW.GW17.5mC_GE02.weighted.mean.bed $dirOut/unique_enhancer.GW.GW17.5mC_GE04.weighted.mean.bed
cd $dirr/homer/GW/
for file in GW17only_*.annotate.motif; do
    echo "Processing "$file
    awk 'NR==FNR {mc[$4]=$7; mc1[$4]=$5; mc2[$4]=$6; next} {print $0"\t"mc1[$1]"\t"mc2[$1]"\t"mc[$1]}' $dirOut/unique_enhancer.GW.GW17.5mC_GE02-04.bed $file > $dirOut/$file.mC
    less $dirOut/$file.mC | awk 'BEGIN {FS="\t"}; NR==1 {print $0} $0~ /protein-coding/ {if($25<-0.2){print $0}}' > $dirr/$file.UMR
done

## Enhancer K4me1 signal
out=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/signal/
mkdir -p $out
cd $out
reg=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me1/unique/unique_enhancer.GW.GW17.bed
name=GE01.K4me1
wig=/home/lli/FetalBrain/HisMod/wigs/A03275.bam.q5.F1028.SET_204.wig.gz
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name > $out$name.coverage.log
name=GE02.K4me1
wig=/home/lli/FetalBrain/HisMod/wigs/A03477.bam.q5.F1028.SET_232.wig.gz
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name > $out$name.coverage.log
name=GE04.K4me1
wig=/home/lli/FetalBrain/HisMod/wigs/A19303.q5.F1028.SET_194.wig.gz
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name > $out$name.coverage.log
# normaliazed to No. of mapped reads (GE04 as reference)
paste unique_enhancer.GW.GW17.bed.GE01.K4me1.coverage unique_enhancer.GW.GW17.bed.GE02.K4me1.coverage unique_enhancer.GW.GW17.bed.GE04.K4me1.coverage | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5*30778793/96165980"\t"$11*30778793/38955337"\t"$17}' > unique_enhancer.GW.GW17.bed.K4me1.coverage

############################################
# Unique H3K4me3
dirIn='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K4me3/'
dirOut=$dirIn/unique/
BEDTOOLS=/projects/epigenomics/software/bedtools-2.23.0/bin/
mkdir -p $dirOut
cd $dirIn
echo -e "Mark\tComparison\tSamples\tSample1\tSample2\tSample1_unique\tSample2_unique" > $dirOut/unique_H3K4me3.summary
## Between MZ twins
### Brain
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03486.H3K4me3.Brain01.minER_300.FDR_0.05.bed -b $dirIn/FindER_scan.A03494.H3K4me3.Brain02.minER_300.FDR_0.05.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K4me3.MZ.Brain01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03494.H3K4me3.Brain02.minER_300.FDR_0.05.bed -b $dirIn/FindER_scan.A03486.H3K4me3.Brain01.minER_300.FDR_0.05.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K4me3.MZ.Brain02.bed
s1=`less $dirIn/FindER_scan.A03486.H3K4me3.Brain01.minER_300.FDR_0.05.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03494.H3K4me3.Brain02.minER_300.FDR_0.05.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K4me3.MZ.Brain01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K4me3.MZ.Brain02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K4me3\tMZ\tBrain01_Brain02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K4me3.summary
## Between neurospheres
### HuFNSC02
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03282.H3K4me3.Cortex02.minER_300.FDR_0.05.bed -b $dirIn/FindER_scan.A03478.H3K4me3.GE02.minER_300.FDR_0.05.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K4me3.Neurospheres.Cortex02.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03478.H3K4me3.GE02.minER_300.FDR_0.05.bed -b $dirIn/FindER_scan.A03282.H3K4me3.Cortex02.minER_300.FDR_0.05.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K4me3.Neurospheres.GE02.bed
s1=`less $dirIn/FindER_scan.A03282.H3K4me3.Cortex02.minER_300.FDR_0.05.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03478.H3K4me3.GE02.minER_300.FDR_0.05.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K4me3.Neurospheres.Cortex02.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K4me3.Neurospheres.GE02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K4me3\tNPCs\tCortex02_GE02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K4me3.summary
## Between GW
### GE02 vs GE04
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03478.H3K4me3.GE02.minER_300.FDR_0.05.bed -b $dirIn/FindER_scan.A19304.H3K4me3.GE04.minER_300.FDR_0.05.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K4me3.GW.GW17_02.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A19304.H3K4me3.GE04.minER_300.FDR_0.05.bed -b $dirIn/FindER_scan.A03478.H3K4me3.GE02.minER_300.FDR_0.05.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K4me3.GW.GW13_02.bed
s1=`less $dirIn/FindER_scan.A03478.H3K4me3.GE02.minER_300.FDR_0.05.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A19304.H3K4me3.GE04.minER_300.FDR_0.05.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K4me3.GW.GW17_02.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K4me3.GW.GW13_02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K4me3\tGW\tGE02_GE04\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K4me3.summary

###################################
# Unique H3K9me3
dirIn='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K9me3/'
dirOut=$dirIn/unique/
BEDTOOLS=/projects/epigenomics/software/bedtools-2.23.0/bin/
mkdir -p $dirOut
cd $dirIn
echo -e "Mark\tComparison\tSamples\tSample1\tSample2\tSample1_unique\tSample2_unique" > $dirOut/unique_H3K9me3.summary
## Between MZ twins
### Brain
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03487.H3K9me3.Brain01.multi.bed -b $dirIn/FindER_scan.A03495.H3K9me3.Brain02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.MZ.Brain01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03495.H3K9me3.Brain02.multi.bed -b $dirIn/FindER_scan.A03487.H3K9me3.Brain01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.MZ.Brain02.bed
s1=`less $dirIn/FindER_scan.A03487.H3K9me3.Brain01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03495.H3K9me3.Brain02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K9me3.MZ.Brain01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K9me3.MZ.Brain02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K9me3\tMZ\tBrain01_Brain02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K9me3.summary
### Cortex
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03271.H3K9me3.Cortex01.multi.bed -b $dirIn/FindER_scan.A03283.H3K9me3.Cortex02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.MZ.Cortex01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03283.H3K9me3.Cortex02.multi.bed -b $dirIn/FindER_scan.A03271.H3K9me3.Cortex01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.MZ.Cortex02.bed
s1=`less $dirIn/FindER_scan.A03271.H3K9me3.Cortex01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03283.H3K9me3.Cortex02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K9me3.MZ.Cortex01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K9me3.MZ.Cortex02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K9me3\tMZ\tCortex01_Cortex02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K9me3.summary
### GE
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03277.H3K9me3.GE01.multi.bed -b $dirIn/FindER_scan.A03479.H3K9me3.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.MZ.GE01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03479.H3K9me3.GE02.multi.bed -b $dirIn/FindER_scan.A03277.H3K9me3.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.MZ.GE02.bed
s1=`less $dirIn/FindER_scan.A03277.H3K9me3.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03479.H3K9me3.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K9me3.MZ.GE01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K9me3.MZ.GE02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K9me3\tMZ\tGE01_GE02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K9me3.summary
### Intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K9me3.MZ.Cortex01.bed -b $dirOut/unique_H3K9me3.MZ.GE01.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.MZ.HuFNSC01.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K9me3.MZ.Cortex02.bed -b $dirOut/unique_H3K9me3.MZ.GE02.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.MZ.HuFNSC02.bed
## Between neurospheres
### HuFNSC01
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03271.H3K9me3.Cortex01.multi.bed -b $dirIn/FindER_scan.A03277.H3K9me3.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.Neurospheres.Cortex01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03277.H3K9me3.GE01.multi.bed -b $dirIn/FindER_scan.A03271.H3K9me3.Cortex01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.Neurospheres.GE01.bed
s1=`less $dirIn/FindER_scan.A03271.H3K9me3.Cortex01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03277.H3K9me3.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K9me3.Neurospheres.Cortex01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K9me3.Neurospheres.GE01.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K9me3\tNPCs\tCortex01_GE01\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K9me3.summary
### HuFNSC02
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03283.H3K9me3.Cortex02.multi.bed -b $dirIn/FindER_scan.A03479.H3K9me3.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.Neurospheres.Cortex02.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03479.H3K9me3.GE02.multi.bed -b $dirIn/FindER_scan.A03283.H3K9me3.Cortex02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.Neurospheres.GE02.bed
s1=`less $dirIn/FindER_scan.A03283.H3K9me3.Cortex02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03479.H3K9me3.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K9me3.Neurospheres.Cortex02.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K9me3.Neurospheres.GE02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K9me3\tNPCs\tCortex02_GE02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K9me3.summary
### intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K9me3.Neurospheres.Cortex01.bed -b $dirOut/unique_H3K9me3.Neurospheres.Cortex02.bed > $dirOut/unique_H3K9me3.Neurospheres.Cortex.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K9me3.Neurospheres.GE01.bed -b $dirOut/unique_H3K9me3.Neurospheres.GE02.bed > $dirOut/unique_H3K9me3.Neurospheres.GE.bed
## Between GW
### GE01 vs GE04
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03277.H3K9me3.GE01.multi.bed -b $dirIn/FindER_scan.A19305.H3K9me3.GE04.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.GW.GW17_01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A19305.H3K9me3.GE04.multi.bed -b $dirIn/FindER_scan.A03277.H3K9me3.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.GW.GW13_01.bed
s1=`less $dirIn/FindER_scan.A03277.H3K9me3.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A19305.H3K9me3.GE04.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K9me3.GW.GW17_01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K9me3.GW.GW13_01.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K9me3\tGW\tGE01_GE04\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K9me3.summary
### GE02 vs GE04
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03479.H3K9me3.GE02.multi.bed -b $dirIn/FindER_scan.A19305.H3K9me3.GE04.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.GW.GW17_02.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A19305.H3K9me3.GE04.multi.bed -b $dirIn/FindER_scan.A03479.H3K9me3.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K9me3.GW.GW13_02.bed
s1=`less $dirIn/FindER_scan.A03479.H3K9me3.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A19305.H3K9me3.GE04.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K9me3.GW.GW17_02.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K9me3.GW.GW13_02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K9me3\tGW\tGE02_GE04\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K9me3.summary
### intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K9me3.GW.GW13_01.bed -b $dirOut/unique_H3K9me3.GW.GW13_02.bed > $dirOut/unique_H3K9me3.GW.GW13.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K9me3.GW.GW17_01.bed -b $dirOut/unique_H3K9me3.GW.GW17_02.bed > $dirOut/unique_H3K9me3.GW.GW17.bed

###################################
# Unique H3K27me3
dirIn='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K27me3/'
dirOut=$dirIn/unique/
BEDTOOLS=/projects/epigenomics/software/bedtools-2.23.0/bin/
mkdir -p $dirOut
cd $dirIn
echo -e "Mark\tComparison\tSamples\tSample1\tSample2\tSample1_unique\tSample2_unique" > $dirOut/unique_H3K27me3.summary
## Between MZ twins
### Brain
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03488.H3K27me3.Brain01.multi.bed -b $dirIn/FindER_scan.A03496.H3K27me3.Brain02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.MZ.Brain01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03496.H3K27me3.Brain02.multi.bed -b $dirIn/FindER_scan.A03488.H3K27me3.Brain01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.MZ.Brain02.bed
s1=`less $dirIn/FindER_scan.A03488.H3K27me3.Brain01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03496.H3K27me3.Brain02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K27me3.MZ.Brain01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K27me3.MZ.Brain02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K27me3\tMZ\tBrain01_Brain02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K27me3.summary
### Cortex
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03272.H3K27me3.Cortex01.multi.bed -b $dirIn/FindER_scan.A03284.H3K27me3.Cortex02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.MZ.Cortex01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03284.H3K27me3.Cortex02.multi.bed -b $dirIn/FindER_scan.A03272.H3K27me3.Cortex01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.MZ.Cortex02.bed
s1=`less $dirIn/FindER_scan.A03272.H3K27me3.Cortex01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03284.H3K27me3.Cortex02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K27me3.MZ.Cortex01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K27me3.MZ.Cortex02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K27me3\tMZ\tCortex01_Cortex02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K27me3.summary
### GE
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03278.H3K27me3.GE01.multi.bed -b $dirIn/FindER_scan.A03480.H3K27me3.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.MZ.GE01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03480.H3K27me3.GE02.multi.bed -b $dirIn/FindER_scan.A03278.H3K27me3.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.MZ.GE02.bed
s1=`less $dirIn/FindER_scan.A03278.H3K27me3.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03480.H3K27me3.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K27me3.MZ.GE01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K27me3.MZ.GE02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K27me3\tMZ\tGE01_GE02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K27me3.summary
### Intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K27me3.MZ.Cortex01.bed -b $dirOut/unique_H3K27me3.MZ.GE01.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.MZ.HuFNSC01.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K27me3.MZ.Cortex02.bed -b $dirOut/unique_H3K27me3.MZ.GE02.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.MZ.HuFNSC02.bed
## Between neurospheres
### HuFNSC01
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03272.H3K27me3.Cortex01.multi.bed -b $dirIn/FindER_scan.A03278.H3K27me3.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.Neurospheres.Cortex01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03278.H3K27me3.GE01.multi.bed -b $dirIn/FindER_scan.A03272.H3K27me3.Cortex01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.Neurospheres.GE01.bed
s1=`less $dirIn/FindER_scan.A03272.H3K27me3.Cortex01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03278.H3K27me3.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K27me3.Neurospheres.Cortex01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K27me3.Neurospheres.GE01.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K27me3\tNPCs\tCortex01_GE01\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K27me3.summary
### HuFNSC02
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03284.H3K27me3.Cortex02.multi.bed -b $dirIn/FindER_scan.A03480.H3K27me3.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.Neurospheres.Cortex02.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03480.H3K27me3.GE02.multi.bed -b $dirIn/FindER_scan.A03284.H3K27me3.Cortex02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.Neurospheres.GE02.bed
s1=`less $dirIn/FindER_scan.A03284.H3K27me3.Cortex02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03480.H3K27me3.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K27me3.Neurospheres.Cortex02.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K27me3.Neurospheres.GE02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K27me3\tNPCs\tCortex02_GE02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K27me3.summary
### intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K27me3.Neurospheres.Cortex01.bed -b $dirOut/unique_H3K27me3.Neurospheres.Cortex02.bed > $dirOut/unique_H3K27me3.Neurospheres.Cortex.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K27me3.Neurospheres.GE01.bed -b $dirOut/unique_H3K27me3.Neurospheres.GE02.bed > $dirOut/unique_H3K27me3.Neurospheres.GE.bed
## Between GW
### GE01 vs GE04
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03278.H3K27me3.GE01.multi.bed -b $dirIn/FindER_scan.A19306.H3K27me3.GE04.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.GW.GW17_01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A19306.H3K27me3.GE04.multi.bed -b $dirIn/FindER_scan.A03278.H3K27me3.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.GW.GW13_01.bed
s1=`less $dirIn/FindER_scan.A03278.H3K27me3.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A19306.H3K27me3.GE04.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K27me3.GW.GW17_01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K27me3.GW.GW13_01.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K27me3\tGW\tGE01_GE04\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K27me3.summary
### GE02 vs GE04
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03480.H3K27me3.GE02.multi.bed -b $dirIn/FindER_scan.A19306.H3K27me3.GE04.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.GW.GW17_02.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A19306.H3K27me3.GE04.multi.bed -b $dirIn/FindER_scan.A03480.H3K27me3.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K27me3.GW.GW13_02.bed
s1=`less $dirIn/FindER_scan.A03480.H3K27me3.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A19306.H3K27me3.GE04.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K27me3.GW.GW17_02.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K27me3.GW.GW13_02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K27me3\tGW\tGE02_GE04\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K27me3.summary
### intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K27me3.GW.GW13_01.bed -b $dirOut/unique_H3K27me3.GW.GW13_02.bed > $dirOut/unique_H3K27me3.GW.GW13.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K27me3.GW.GW17_01.bed -b $dirOut/unique_H3K27me3.GW.GW17_02.bed > $dirOut/unique_H3K27me3.GW.GW17.bed

###################################
# Unique H3K36me3
dirIn='/projects/epigenomics/users/lli/FetalBrain/ChIPseq/ER/H3K36me3/'
dirOut=$dirIn/unique/
BEDTOOLS=/projects/epigenomics/software/bedtools-2.23.0/bin/
mkdir -p $dirOut
cd $dirIn
echo -e "Mark\tComparison\tSamples\tSample1\tSample2\tSample1_unique\tSample2_unique" > $dirOut/unique_H3K36me3.summary
## Between MZ twins
### Brain
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03489.H3K36me3.Brain01.multi.bed -b $dirIn/FindER_scan.A03497.H3K36me3.Brain02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.MZ.Brain01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03497.H3K36me3.Brain02.multi.bed -b $dirIn/FindER_scan.A03489.H3K36me3.Brain01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.MZ.Brain02.bed
s1=`less $dirIn/FindER_scan.A03489.H3K36me3.Brain01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03497.H3K36me3.Brain02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K36me3.MZ.Brain01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K36me3.MZ.Brain02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K36me3\tMZ\tBrain01_Brain02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K36me3.summary
### Cortex
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03273.H3K36me3.Cortex01.multi.bed -b $dirIn/FindER_scan.A03285.H3K36me3.Cortex02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.MZ.Cortex01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03285.H3K36me3.Cortex02.multi.bed -b $dirIn/FindER_scan.A03273.H3K36me3.Cortex01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.MZ.Cortex02.bed
s1=`less $dirIn/FindER_scan.A03273.H3K36me3.Cortex01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03285.H3K36me3.Cortex02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K36me3.MZ.Cortex01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K36me3.MZ.Cortex02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K36me3\tMZ\tCortex01_Cortex02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K36me3.summary
### GE
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03279.H3K36me3.GE01.multi.bed -b $dirIn/FindER_scan.A03481.H3K36me3.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.MZ.GE01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03481.H3K36me3.GE02.multi.bed -b $dirIn/FindER_scan.A03279.H3K36me3.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.MZ.GE02.bed
s1=`less $dirIn/FindER_scan.A03279.H3K36me3.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03481.H3K36me3.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K36me3.MZ.GE01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K36me3.MZ.GE02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K36me3\tMZ\tGE01_GE02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K36me3.summary
### Intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K36me3.MZ.Cortex01.bed -b $dirOut/unique_H3K36me3.MZ.GE01.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.MZ.HuFNSC01.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K36me3.MZ.Cortex02.bed -b $dirOut/unique_H3K36me3.MZ.GE02.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.MZ.HuFNSC02.bed
## Between neurospheres
### HuFNSC01
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03273.H3K36me3.Cortex01.multi.bed -b $dirIn/FindER_scan.A03279.H3K36me3.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.Neurospheres.Cortex01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03279.H3K36me3.GE01.multi.bed -b $dirIn/FindER_scan.A03273.H3K36me3.Cortex01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.Neurospheres.GE01.bed
s1=`less $dirIn/FindER_scan.A03273.H3K36me3.Cortex01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03279.H3K36me3.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K36me3.Neurospheres.Cortex01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K36me3.Neurospheres.GE01.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K36me3\tNPCs\tCortex01_GE01\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K36me3.summary
### HuFNSC02
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03285.H3K36me3.Cortex02.multi.bed -b $dirIn/FindER_scan.A03481.H3K36me3.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.Neurospheres.Cortex02.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03481.H3K36me3.GE02.multi.bed -b $dirIn/FindER_scan.A03285.H3K36me3.Cortex02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.Neurospheres.GE02.bed
s1=`less $dirIn/FindER_scan.A03285.H3K36me3.Cortex02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A03481.H3K36me3.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K36me3.Neurospheres.Cortex02.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K36me3.Neurospheres.GE02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K36me3\tNPCs\tCortex02_GE02\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K36me3.summary
### intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K36me3.Neurospheres.Cortex01.bed -b $dirOut/unique_H3K36me3.Neurospheres.Cortex02.bed > $dirOut/unique_H3K36me3.Neurospheres.Cortex.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K36me3.Neurospheres.GE01.bed -b $dirOut/unique_H3K36me3.Neurospheres.GE02.bed > $dirOut/unique_H3K36me3.Neurospheres.GE.bed
## Between GW
### GE01 vs GE04
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03279.H3K36me3.GE01.multi.bed -b $dirIn/FindER_scan.A19307.H3K36me3.GE04.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.GW.GW17_01.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A19307.H3K36me3.GE04.multi.bed -b $dirIn/FindER_scan.A03279.H3K36me3.GE01.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.GW.GW13_01.bed
s1=`less $dirIn/FindER_scan.A03279.H3K36me3.GE01.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A19307.H3K36me3.GE04.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K36me3.GW.GW17_01.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K36me3.GW.GW13_01.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K36me3\tGW\tGE01_GE04\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K36me3.summary
### GE02 vs GE04
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A03481.H3K36me3.GE02.multi.bed -b $dirIn/FindER_scan.A19307.H3K36me3.GE04.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.GW.GW17_02.bed
$BEDTOOLS/intersectBed -a $dirIn/FindER_scan.A19307.H3K36me3.GE04.multi.bed -b $dirIn/FindER_scan.A03481.H3K36me3.GE02.multi.bed -v | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/unique_H3K36me3.GW.GW13_02.bed
s1=`less $dirIn/FindER_scan.A03481.H3K36me3.GE02.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s2=`less $dirIn/FindER_scan.A19307.H3K36me3.GE04.multi.bed | awk '{s=s+$3-$2} END{print s}'`
s1_unique=`less $dirOut/unique_H3K36me3.GW.GW17_02.bed | awk '{s=s+$3-$2} END{print s}'`
s2_unique=`less $dirOut/unique_H3K36me3.GW.GW13_02.bed | awk '{s=s+$3-$2} END{print s}'`
echo -e "H3K36me3\tGW\tGE02_GE04\t"$s1"\t"$s2"\t"$s1_unique"\t"$s2_unique >> $dirOut/unique_H3K36me3.summary
### intersect
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K36me3.GW.GW13_01.bed -b $dirOut/unique_H3K36me3.GW.GW13_02.bed > $dirOut/unique_H3K36me3.GW.GW13.bed
/projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a $dirOut/unique_H3K36me3.GW.GW17_01.bed -b $dirOut/unique_H3K36me3.GW.GW17_02.bed > $dirOut/unique_H3K36me3.GW.GW17.bed







