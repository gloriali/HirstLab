#!/bin/sh

# align to hg19v65
## jaguar alignment
ref=/home/pubseq/genomes/Homo_sapiens/hg19a/jaguar/1.7.5/ens65only/bwa_ind/transcriptome/100/ref.fa
ens=hg19_ens65
dirIn=/projects/sftp/ckleinman/incoming/fetal_brain_RNAseq/
dirOut=/projects/epigenomics3/users/lli/Claudia/RNAseq/bam/
mkdir -p $dirOut
cd $dirIn
for f1 in *_R1.fastq.gz; do
    f2=$(echo $f1 | sed 's/_R1./_R2./')
    id=$(echo $f1 | cut -d'.' -f5)
    age=$(less /projects/epigenomics3/users/lli/Claudia/SampleInfo.txt | awk '{if($1=="'$id'"){print $2}}')
    name=$id.$age
    echo $name
    /home/lli/HirstLab/Pipeline/shell/jaguar.sh -i $dirIn -o $dirOut -f1 $f1 -f2 $f2 -n $name -r $ref -v $ens
done
function JRalign {
    ens=hg19_ens65
    dirOut=/projects/epigenomics3/users/lli/Claudia/RNAseq/bam/
    name=$1
    rm -rf $dirOut/$name/
    /home/mbilenky/bin/Solexa_Shell/RunJR.sh $dirOut/$name".sortedByName.bam" $dirOut/$name $ens &> $dirOut/run/$name.j.log
}
export -f JRalign
dirOut=/projects/epigenomics3/users/lli/Claudia/RNAseq/bam/
cd $dirOut
ls *.sortedByName.bam | awk '{gsub(".sortedByName.bam", ""); print $0}' > $dirOut/List.txt
cat List.txt | parallel --gnu JRalign 
## QC and RPKM
samtools=/gsc/software/linux-x86_64/samtools-0.1.13/samtools
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
dirIn=/projects/epigenomics3/users/lli/Claudia/RNAseq/bam/
dirOut=/projects/epigenomics3/users/lli/Claudia/RNAseq/RPKM/
ens=hg19v65
mkdir -p $dirOut
cd $dirIn
for bam in *_withJunctionsOnGenome_dupsFlagged.bam; do
    name=$(echo $bam | sed 's/_withJunctionsOnGenome_dupsFlagged.bam//g')
    echo $name
    rm -rf $dirOut/$name/
    /home/lli/bin/Solexa_Shell/src/RNAseqMaster.sh $bam $name $dirOut $ens S 0 "1,1,1,1,1" /project/epigenomics/resources/
done

