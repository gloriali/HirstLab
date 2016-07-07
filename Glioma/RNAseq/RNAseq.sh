#!/bin/sh

# link RPKM files
dirIn=/projects/edcc_new/reference_epigenomes/
dirOut=/projects/epigenomics2/users/lli/glioma/RNAseq/RPKM/
mkdir -p $dirOut
for lib in 19 21 22 23 47; do
    echo $lib;
    ln -s $dirIn/CEMT_$lib/bams/RNA-Seq/qca/*/coverage/*.G.A.rpkm.pc $dirOut/CEMT_$lib.G.A.rpkm.pc
done

# IDH1/2 expression
dirOut=/projects/epigenomics2/users/lli/glioma/RNAseq/RPKM/
> $dirOut/IDH.RPKM
for lib in 19 21 22 23 47; do
    echo $lib;
    less $dirOut/CEMT_$lib.G.A.rpkm.pc | awk '$1 ~ /ENSG00000138413/ {print $1"\tIDH1\tCEMT_""'$lib'""\t"$3}' >> $dirOut/IDH.RPKM
    less $dirOut/CEMT_$lib.G.A.rpkm.pc | awk '$1 ~ /ENSG00000182054/ {print $1"\tIDH2\tCEMT_""'$lib'""\t"$3}' >> $dirOut/IDH.RPKM
done

# Update Fetal Brain RPKM to hg19v69
## re-sort bam on read names
samtools=/gsc/software/linux-x86_64/samtools-0.1.13/samtools
dirIn=/projects/epigenomics/users/lli/FetalBrain/RNAseq/bam/
dirOut=/projects/epigenomics2/users/lli/glioma/RNAseq/NPC_fq/
mkdir -p $dirOut
cd $dirIn
for file in *.bam; do
    name=$(echo $file | sed -e 's/.bam//g')
    echo $name
    $samtools sort $dirIn/$file $dirOut/$name.sorted -n 
done
## bam to fastq
cd $dirOut
ls *.bam > $dirOut/BamList.txt
function bam2fq {
    BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
    dirOut=/projects/epigenomics2/users/lli/glioma/RNAseq/NPC_fq/
    file=$1
    name=$(echo $file | sed -e 's/.sorted.bam//g')
    echo $name
    $BEDTOOLS/bamToFastq -i $dirOut/$file -fq $dirOut/$name.1.fq -fq2 $dirOut/$name.2.fq
}
export -f bam2fq
cat BamList.txt | parallel --gnu bam2fq 
cat <(ls *.fq) | parallel --gnu gzip
rm *.sorted.bam
## jaguar alignment
ref=/home/pubseq/genomes/Homo_sapiens/hg19a/jaguar/1.7.5/ens69/bwa_ind/transcriptome/75/ref.fa
ens=hg19_ens69
dirIn=/projects/epigenomics2/users/lli/glioma/RNAseq/NPC_fq/
dirOut=/projects/epigenomics2/users/lli/glioma/RNAseq/NPC_bam/
mkdir -p $dirOut
for f1 in *.1.fq.gz; do
    name=$(echo $f1 | sed -e 's/.1.fq.gz//g')
    f2=$name.2.fq.gz
    echo $f1 $f2 $name
    /home/lli/HirstLab/Pipeline/shell/jaguar.sh -i $dirIn -o $dirOut -f1 $f1 -f2 $f2 -n $name -r $ref -v $ens
done
## QC and RPKM
samtools=/gsc/software/linux-x86_64/samtools-0.1.13/samtools
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
dirIn=/projects/epigenomics2/users/lli/glioma/RNAseq/NPC_bam/
dirOut=/projects/epigenomics2/users/lli/glioma/RNAseq/NPC_RPKM/
ens=hg19v69
name=A03473.Cortex01
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A03474.GE01
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A03475.Cortex02
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A03476.GE02
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A03484.Brain01
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A04599.Cortex03
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A07825.Brain02
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens S 0 "1,1,1,1,1" $JAVA $samtools 
name=A15295.GE03
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens S 0 "1,1,1,1,1" $JAVA $samtools 
name=A15298.Cortex04
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens S 0 "1,1,1,1,1" $JAVA $samtools 
name=A15299.GE04
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens S 0 "1,1,1,1,1" $JAVA $samtools 

# DE between glioma and NPCs
## generate matlab code for DEfine
dirIn='/projects/epigenomics2/users/lli/glioma/RNAseq/';
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine/';
mkdir -p $dirOut
echo -e "%/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab
addpath /home/mbilenky/matlab -end" > $dirOut/DEfine.glioma.m
for name1 in CEMT_19 CEMT_21 CEMT_22 CEMT_23 CEMT_47; do
    mkdir -p $dirOut/$name1/
    for name2 in Cortex02 GE02 Cortex04 GE04; do
        sample1=$name1
        sample2=$(ls $dirIn/NPC_RPKM/*$name2/*.bam | sed -e 's/.*\///g' | sed -e 's/.bam//g')
        echo $sample1 $sample2
        if [ $name2 == 'Cortex02' ]; then
            echo -e "
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% $name1 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
" >> $dirOut/DEfine.glioma.m
        fi
        echo -e "
%%%%%%%%%%%%%%%%%%%%%%%%%% $name1 vs $name2 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.01; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='$dirIn';
dirOut='$dirOut/$name1/';
sample1='$sample1'; name1='$name1';
sample2='$sample2'; name2='$name2';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, 'RPKM/', sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, 'NPC_RPKM/', sample2,'/coverage/', sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);
" >> $dirOut/DEfine.glioma.m
    done
done
