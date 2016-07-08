#!/bin/sh

## QC and RPKM
samtools=/gsc/software/linux-x86_64/samtools-0.1.13/samtools
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
dirIn=/projects/epigenomics2/Brain/NJabado/RNAseq/bam/
dirOut=/projects/epigenomics2/Brain/NJabado/RNAseq/hg19v69/
ens=hg19v69
name=B07_PC_AK1189_RNAseq
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/PX0380_C9587ANXX_8_CCACGC_withJunctionsOnGenome_dupsFlagged.bam $name $dirOut $ens S 0 "1,1,1,1,1" $JAVA $samtools 
name=B08_PC_OPK69_RNAseq
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/PX0380_C9587ANXX_8_CTATAC_withJunctionsOnGenome_dupsFlagged.bam $name $dirOut $ens S 0 "1,1,1,1,1" $JAVA $samtools 
name=B09_PC_JN27BS+4_RNAseq
mkdir -p $dirOut/$name/
/home/acarles/Solexa_Shell/src/RNAseqMaster.sh $dirIn/PX0380_C9587ANXX_8_GCAAGG_withJunctionsOnGenome_dupsFlagged.bam $name $dirOut $ens S 0 "1,1,1,1,1" $JAVA $samtools 

