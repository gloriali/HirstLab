#!/bin/sh

# QC and RPKM
samtools=/gsc/software/linux-x86_64/samtools-0.1.13/samtools
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/RPKM/
ens=hg19v69
cd $dirIn
for bam in *.bam; do
    name=$(echo $bam | cut -d'.' -f1)
    echo $name 
    rm -rf $dirOut/$name/
    /home/lli/bin/Solexa_Shell/src/RNAseqMaster.sh $(readlink -f $bam) $name $dirOut $ens S 0 "1,1,1,1,1" /projects/epigenomics/resources/ $JAVA $samtools 
done
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/bam/
cd $dirIn
for bam in *.bam; do
    name=$(echo $bam | cut -d'.' -f1)
    echo $name
    $samtools index $bam
    $bamstats -g 2864785220 -q 10 -b $bam  > $name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn $name $name.bamstats
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn

# RPKM matrix
cd /projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/RPKM/
> /projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/RPKM/RPKM.long
echo "ENSG" > vitc.RPKM
less MGG119_control/coverage/MGG119_control.G.A.rpkm.pc | awk '{print $1}' >> vitc.RPKM
for file in */coverage/*.G.A.rpkm.pc; do
    lib=$(echo $file | sed 's/\/.*//g')
    echo -e "ENSG\t$lib" > x
    less $file | awk '{print $1"\t"$3}' >> x
    join vitc.RPKM x | sed 's/ /\t/g' >y
    mv y vitc.RPKM
    less $file | awk '{print "'$lib'""\t"$1"\t"$3}' >> /projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/RPKM/RPKM.long
done
rm x

# DE between glioma and NPCs
## generate matlab code for DEfine
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/RPKM/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/DEfine/
mkdir -p $dirOut
echo -e "%/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab
addpath /home/mbilenky/matlab -end" > $dirOut/DEfine.vitc.m
s1=("NHAR_vitc" "NHAR_control" "NHA_vitc" "NHAR_vitc" "MGG119_vitc")
s2=("NHAR_control" "NHA_control" "NHA_control" "NHA_control" "MGG119_control")
for ((i=0; i<5; i++)); do
    sample1=${s1[i]}; sample2=${s2[i]};
    name1=$sample1; name2=$sample2
    echo $sample1 $sample2
    echo -e "
%%%%%%%%%%%%%%%%%%%%%%%%%% $name1 vs $name2 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='$dirIn';
dirOut='$dirOut';
sample1='$sample1'; name1='$name1';
sample2='$sample2'; name2='$name2';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, sample1,'/coverage/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);
" >> $dirOut/DEfine.vitc.m
done

# Call SNPs
genome=/home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa
SAMTOOLS=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
BCFTOOLS=/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
snpEff=/projects/wtsspipeline/programs/external_programs/snpEff3.3/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/SNP/
mkdir -p $dirOut
for file in $dirIn/*.bam; do
    name=$(basename $file | cut -d'.' -f1)
    echo $name
    $SAMTOOLS mpileup -C50 -uf $genome $file | $BCFTOOLS/bcftools view -bvcg - > $dirOut/$name.bcf
    $BCFTOOLS/bcftools view $dirOut/$name.bcf | $BCFTOOLS/vcfutils.pl varFilter -D100 > $dirOut/$name.vcf
    $JAVA -jar $snpEff/snpEff.jar GRCh37.69 $dirOut/$name.vcf -c $snpEff/snpEff.config > $dirOut/$name.snpEff.vcf
    $JAVA -jar $snpEff/SnpSift.jar annotate -id $snpEff/cosmic_v64.vcf $dirOut/$name.snpEff.vcf > $dirOut/$name.snpEff.COSMIC.vcf
    less $dirOut/$name.snpEff.COSMIC.vcf | awk '$1 !~ /#/ {chr=$1; pos=$2; $1=$2=""; print chr"\t"pos"\t"pos+1"\t"$0}' > $dirOut/$name.snpEff.COSMIC.bed
done

