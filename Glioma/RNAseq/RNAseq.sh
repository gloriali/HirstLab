#!/bin/sh

# link files
dirIn=/projects/edcc_new/reference_epigenomes/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/
mkdir -p $dirOut/bam/; mkdir -p $dirOut/RPKM/
for lib in 19 21 22 23 47 73 74 75 76 78 79 81; do
    IDH=$(less $dirOut/../samples.txt | awk '{if($1=="CEMT_""'$lib'")print $2}')
    echo $IDH $lib 
    ln -s $dirIn/CEMT_$lib/bams/RNA-Seq/*.bam $dirOut/bam/$IDH.CEMT_$lib.bam
    ln -s $dirIn/CEMT_$lib/bams/RNA-Seq/*.bam.bai $dirOut/bam/$IDH.CEMT_$lib.bam.bai
    ln -s $dirIn/CEMT_$lib/bams/RNA-Seq/qca*/*/coverage/*.G.A.rpkm.pc $dirOut/RPKM/$IDH.CEMT_$lib.G.A.rpkm.pc
done
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/sambamba-bamStats
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/bam/
for bam in $dirIn/*.bam; do
    name=$(basename $bam | sed 's/\.bam//'); echo $name
    $bamstats -g 2864785220 -t 8 $bam > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn $name $dirIn/$name.bamstats
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn

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
rm -rf $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A03474.GE01
rm -rf $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A03475.Cortex02
rm -rf $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A03476.GE02
rm -rf $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A03484.Brain01
rm -rf $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A04599.Cortex03
rm -rf $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A07825.Brain02
rm -rf $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A15295.GE03
rm -rf $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A15298.Cortex04
rm -rf $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 
name=A15299.GE04
rm -rf $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens R 0 "1,1,1,1,1" $JAVA $samtools 

# RPKM matrix
> /projects/epigenomics2/users/lli/glioma/RNAseq/RPKM.long
cd /projects/epigenomics2/users/lli/glioma/RNAseq/RPKM/
echo "ENSG" > glioma.RPKM
less CEMT_19.G.A.rpkm.pc | awk '{print $1}' >> glioma.RPKM
for file in *.G.A.rpkm.pc; do
    lib=$(echo $file | sed 's/.G.A.rpkm.pc//g')
    echo -e "ENSG\t$lib" > x
    less $file | awk '{print $1"\t"$3}' >> x
    join glioma.RPKM x | sed 's/ /\t/g' >y
    mv y glioma.RPKM
    less $file | awk '{print "'$lib'""."$1"\t"$3}' >> /projects/epigenomics2/users/lli/glioma/RNAseq/RPKM.long
done
rm x
cd /projects/epigenomics2/users/lli/glioma/RNAseq/NPC_RPKM/
echo "ENSG" > NPC.RPKM
less A03473.Cortex01/coverage/A03473.Cortex01.G.A.rpkm.pc | awk '{print $1}' >> NPC.RPKM
for file in */coverage/*.G.A.rpkm.pc; do
    lib=$(echo $file | sed 's/\/.*//g' | sed 's/A[0-9]*.//g')
    echo -e "ENSG\t$lib" > x
    less $file | awk '{print $1"\t"$3}' >> x
    join NPC.RPKM x | sed 's/ /\t/g' >y
    mv y NPC.RPKM
    less $file | awk '{print "NPC_""'$lib'""."$1"\t"$3}' >> /projects/epigenomics2/users/lli/glioma/RNAseq/RPKM.long
done
rm x
join /projects/epigenomics2/users/lli/glioma/RNAseq/RPKM/glioma.RPKM /projects/epigenomics2/users/lli/glioma/RNAseq/NPC_RPKM/NPC.RPKM | sed 's/ /\t/g' > /projects/epigenomics2/users/lli/glioma/RNAseq/RPKM.matrix

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
## intersect DE from NPCs
dirOut='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine/';
annotate='/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes';
echo -e "Sample\tDE\tCortex02\tCortex04\tGE02\tGE04\tintersect" > $dirOut/DE.pc.summary
for lib in CEMT_19 CEMT_21 CEMT_22 CEMT_23 CEMT_47; do
    echo $lib
    cat $dirOut/$lib/UP.$lib* | cut -f 1 | sort | uniq -c | awk '$1 ~ 4 {print $2}' | join - $annotate | sed 's/ /\t/g' > $dirOut/UP.$lib'_NPC.FDR_0.01.rmin_0.005.Nmin_25'
    cat $dirOut/$lib/DN.$lib* | cut -f 1 | sort | uniq -c | awk '$1 ~ 4 {print $2}' | join - $annotate | sed 's/ /\t/g' > $dirOut/DN.$lib'_NPC.FDR_0.01.rmin_0.005.Nmin_25'
    up=$(wc -l $dirOut/$lib/UP.$lib* | awk 'NR==1{s1=$1} NR==2{s2=$1} NR==3{s3=$1} NR==4{s4=$1} END{print s1"\t"s2"\t"s3"\t"s4}')
    up_intersect=$(less $dirOut/UP.$lib'_NPC.FDR_0.01.rmin_0.005.Nmin_25' | wc -l)
    echo -e "$lib\tUP\t$up\t$up_intersect" >> $dirOut/DE.pc.summary
    dn=$(wc -l $dirOut/$lib/DN.$lib* | awk 'NR==1{s1=$1} NR==2{s2=$1} NR==3{s3=$1} NR==4{s4=$1} END{print s1"\t"s2"\t"s3"\t"s4}')
    dn_intersect=$(less $dirOut/DN.$lib'_NPC.FDR_0.01.rmin_0.005.Nmin_25' | wc -l)
    echo -e "$lib\tDN\t$dn\t$dn_intersect" >> $dirOut/DE.pc.summary
done
cd $dirOut
cat UP.CEMT_19_NPC.FDR_0.01.rmin_0.005.Nmin_25 UP.CEMT_22_NPC.FDR_0.01.rmin_0.005.Nmin_25 UP.CEMT_47_NPC.FDR_0.01.rmin_0.005.Nmin_25 | awk '{print $1}' | sort | uniq -d > UP.IDHmut.NPC
cat DN.CEMT_19_NPC.FDR_0.01.rmin_0.005.Nmin_25 DN.CEMT_22_NPC.FDR_0.01.rmin_0.005.Nmin_25 DN.CEMT_47_NPC.FDR_0.01.rmin_0.005.Nmin_25 | awk '{print $1}' | sort | uniq -d > DN.IDHmut.NPC

## gene fusion detection
### bam to fastq
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.19/samtools
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics2/users/lli/glioma/RNAseq/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/fq/
mkdir -p $dirOut
cd $dirIn
for bam in $dirIn/*.bam; do
    name=$(basename $bam | sed 's/.bam//g')
    echo $name
    $samtools sort -n -@ 8 $bam $dirIn/$name.nsort
    $BEDTOOLS/bamToFastq -i $dirIn/$name.nsort.bam -fq $dirOut/$name.1.fq -fq2 $dirOut/$name.2.fq
done
### deFuse
export PATH=$PATH":/home/lli/anaconda2/bin/"
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/fq/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/fusion/
dirRef=/projects/epigenomics3/epigenomics3_results/users/lli/defuse_ref/hg19v69/
mkdir -p $dirOut
mkdir -p $dirRef
defuse_create_ref.pl -d $dirRef
cd $dirIn
for fq1 in *.1.fq; do
    name=$(basename $fq1 | sed 's/.1.fq//g'); fq2=$name.2.fq
    echo $name 
    mkdir -p $dirOut/$name/
    defuse_run.pl -d $dirRef -1 $fq1 -2 $fq2 -o $dirOut/$name/ -p 15 -n $name
done

