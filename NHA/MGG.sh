#!/bin/sh

export PATH=/home/lli/anaconda2/bin/:$PATH
export PYTHONPATH=/home/lli/anaconda2/lib/python2.7/site-packages
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirNPC=/projects/epigenomics2/users/lli/glioma/
dirRNA=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/
dirHM=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/
dir5mC=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/
dir5hmC=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/
mkdir -p $dirOut

## ChIPseq
java=/gsc/software/linux-x86_64-centos6/jdk1.8.0_162/jre/bin/java
finder2=/home/mbilenky/bin/FindER2/finder2.jar
dirIn=$dirOut/ChIPseq/bam/
mkdir -p $dirOut/ChIPseq/FindER2/
for file in $dirNPC/ChIPseq/bam/*/*NPC*.bam $dirHM/bam/*/*MGG*.bam; do
    mark=$(echo $file | sed 's/.*bam\///' | sed 's/\/.*//')
    sample=$(basename $file | sed 's/.bam//' | sed 's/\./_/' | cut -d'_' -f2,3 | sed 's/119//')
    echo $mark $sample
    mkdir -p $dirIn/$mark/
    ln -s $file $dirIn/$mark/$mark"."$sample.bam
    ln -s $file.bai $dirIn/$mark/$mark"."$sample.bam.bai
done
for inp in $dirIn/input/*.bam; do
    sample=$(basename $inp | cut -d'.' -f2);
    echo $sample
    sig1=$(echo $dirIn/H3K27me3/H3K27me3.$sample.bam)
    sig2=$(echo $dirIn/H3K4me3/H3K4me3.$sample.bam)
    sig3=$(echo $dirIn/H3K4me1/H3K4me1.$sample.bam)
    sig4=$(echo $dirIn/H3K27ac/H3K27ac.$sample.bam)
    sig5=$(echo $dirIn/H3K36me3/H3K36me3.$sample.bam)
    sig=${sig1},${sig2},${sig3},${sig4},${sig5}
    $java -jar -Xmx25G $finder2 inputBam:$inp signalBam:$sig outDir:$dirOut/ChIPseq/FindER2/ acgtDir:/projects/epigenomics2/users/mbilenky/resources/hg19/ACGT 
done
echo -e "Mark\tSample\tN_region\tTotal_length\tAverage_length" > $dirOut/ChIPseq/FindER2/ER_summary.txt
for file in $dirOut/ChIPseq/FindER2/*.bed; do
    mark=$(basename $file | cut -d'.' -f1)
    sample=$(basename $file | cut -d'.' -f2)
    echo $mark $sample
    less $file | cut -f1-3 > $file.1; mv $file.1 $file
    echo -e $mark"\t"$sample"\t"$(less $file | wc -l)"\t"$(less $file | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"$4/$3}' >> $dirOut/ChIPseq/FindER2/ER_summary.txt
done
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    intervene upset -i $dirOut/ChIPseq/FindER2/$mark.*.bed --project $mark -o $dirOut/ChIPseq/FindER2/
done

## 5mC
mkdir -p $dirOut/WGBS/
for file in $dir5mC/MGG*.combine.5mC.CpG $dirNPC/WGBS/NPC_GE04.combine.5mC.CpG; do
    ln -s $file $dirOut/WGBS/
done
echo -e "sample\tcoverage\tN" > $dirOut/WGBS/qc_5mC_coverage.txt
echo -e "sample\ttype\tfractional\tN" > $dirOut/WGBS/qc_5mC_profile.txt 
echo -e "sample\ttype\tmin\tymin\tlower\tmedian\tupper\tymax\tmax" > $dirOut/WGBS/qc_5mC_quantile.txt #ymin: 10% quantile; ymax: 90% quantile
cd $dirOut/WGBS/
for file in *.combine.5mC.CpG; do
    lib=$(echo $file | sed -e 's/.combine.5mC.CpG//g')
    echo "Processing" $lib
    less $file | awk '{c = $4 + $5; if(c >= 5000){s[5001]++} else {s[c]++}} END{for(i = 3; i <= 5001; i++){print "'$lib'""\t"i"\t"s[i]}}' >> $dirOut/WGBS/qc_5mC_coverage.txt
    less $file | awk '{s[int($6*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tgenome\t"i/100"\t"s[i]}}' >> $dirOut/WGBS/qc_5mC_profile.txt 
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | awk '{s[int($1*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tCGI\t"i/100"\t"s[i]}}' >> $dirOut/WGBS/qc_5mC_profile.txt 
    less $file | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirOut/WGBS/qc_5mC_quantile.txt
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirOut/WGBS/qc_5mC_quantile.txt
done
mkdir -p $dirOut/WGBS/DMR/
echo -e "sample\tp-value\tdelta\tm\ttotal\thyper\thypo" > $dirOut/WGBS/DMR/DM.summary.stats
echo -e "sample\tsize\tcut\tmedian_length\tmedian_N_CpG\ttotal\thyper\thypo" > $dirOut/WGBS/DMR/DMR.summary.stats
pth=0.01; delta=0.5; m=0.7; cov=3; size=500; cut=3; lib1=MGG_control
for lib2 in MGG_vitc NPC_GE04; do
    name=$lib1'_'$lib2; echo $name
    /home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirOut/WGBS/ -o $dirOut/WGBS/DMR/ -f1 $lib1.combine.5mC.CpG -f2 $lib2.combine.5mC.CpG -n $name -p $pth -d $delta -m $m -c $cov
    /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut/WGBS/DMR/ -o $dirOut/WGBS/DMR/ -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
done

## hMeDIP

## RNAseq
mkdir -p $dirOut/RNAseq/DEfine/
for file in $dirRNA/RPKM/MGG*/coverage/*.G.A.rpkm.pc $dirNPC/RNAseq/NPC_RPKM/*GE04/coverage/*.G.A.rpkm.pc; do
    name=$(basename $file | sed 's/119//' | sed 's/A15299.GE04/NPC_GE04/')
    ln -s $file $dirOut/RNAseq/$name
done
echo -e "%/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab
addpath /home/mbilenky/matlab -end" > $dirOut/RNAseq/DEfine/DEfine.MGG.m
sample1=MGG_control
for sample2 in MGG_vitc NPC_GE04; do
    name1=$sample1; name2=$sample2
    echo $sample1 $sample2
    echo -e "
%%%%%%%%%%%%%%%%%%%%%%%%%% $name1 vs $name2 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.05; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='$dirOut/RNAseq/';
dirOut='$dirOut/RNAseq/DEfine/';
sample1='$sample1'; name1='$name1';
sample2='$sample2'; name2='$name2';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);
" >> $dirOut/RNAseq/DEfine/DEfine.MGG.m
done

