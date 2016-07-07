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
dirIn=/projects/epigenomics2/users/lli/glioma/RNAseq/NPC_bam/
echo -e "A03473.Cortex01 R
A03474.GE01 R
A03475.Cortex02 R
A03476.GE02 R
A03484.Brain01 R
A04599.Cortex03 R
A07825.Brain02 S
A15295.GE03 S
A15298.Cortex04 S
A15299.GE04 S" > $dirIn/strand.specificity
function RNAseqQC {
    samtools=/gsc/software/linux-x86_64/samtools-0.1.13/samtools
    JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
    dirIn=/projects/epigenomics2/users/lli/glioma/RNAseq/NPC_bam/
    dirOut=/projects/epigenomics2/users/lli/glioma/RNAseq/NPC_RPKM/
    ens=hg19v69
    name=$1
    s=$2
    mkdir -p $dirOut/$name/
    /home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/$name/$name'_withJunctionsOnGenome_dupsFlagged.bam' $name $dirOut $ens $s 0 "1,1,1,1,1" $JAVA $samtools > $dirOut/$name/qca.log
}
export -f RNAseqQC
cat $dirIn/strand.specificity | parallel --gnu RNAseqQC 

