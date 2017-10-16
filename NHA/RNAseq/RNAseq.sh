#!/bin/sh

## QC and RPKM
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

# DE between glioma and NPCs
## generate matlab code for DEfine
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/RPKM/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/DEfine/
mkdir -p $dirOut
echo -e "%/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab
addpath /home/mbilenky/matlab -end" > $dirOut/DEfine.vitc.m
s1=("NHAR_vitc" "NHAR_control" "NHA_vitc" "NHAR_vitc" "MGG119_control")
s2=("NHAR_control" "NHA_control" "NHA_control" "NHA_control" "MGG119_vitc")
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

