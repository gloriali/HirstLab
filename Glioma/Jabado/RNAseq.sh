#!/bin/sh

## QC and RPKM
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirIn=/projects/epigenomics2/Brain/NJabado/RNAseq/bam/
dirOut=/projects/epigenomics2/Brain/NJabado/QC/
mkdir -p $dirOut
cd $dirIn
for bam in *.bam; do
    name=$(echo $bam | sed -e 's/_withJunctionsOnGenome_dupsFlagged.bam//g' | sed -e 's/C9587ANXX_8_//g')
    echo $name
    $bamstats -g 2864785220 -q 10 -b $dirIn/$bam  > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $name $dirIn/$name.bamstats
done
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
join $dirOut/B07_PC_AK1189_RNAseq/coverage/B07_PC_AK1189_RNAseq.G.A.rpkm.pc $dirOut/B08_PC_OPK69_RNAseq/coverage/B08_PC_OPK69_RNAseq.G.A.rpkm.pc | join - $dirOut/B09_PC_JN27BS+4_RNAseq/coverage/B09_PC_JN27BS+4_RNAseq.G.A.rpkm.pc | awk 'BEGIN {print "ID\tB07_PC_AK1189\tB08_PC_OPK69\tB09_PC_JN27BS+4"} {print $1"\t"$3"\t"$8"\t"$13}' > $dirOut/rpkm.pc
join $dirOut/B07_PC_AK1189_RNAseq/coverage/B07_PC_AK1189_RNAseq.G.A.rpkm.nc $dirOut/B08_PC_OPK69_RNAseq/coverage/B08_PC_OPK69_RNAseq.G.A.rpkm.nc | join - $dirOut/B09_PC_JN27BS+4_RNAseq/coverage/B09_PC_JN27BS+4_RNAseq.G.A.rpkm.nc | awk 'BEGIN {print "ID\tB07_PC_AK1189\tB08_PC_OPK69\tB09_PC_JN27BS+4"} {print $1"\t"$3"\t"$8"\t"$13}' > $dirOut/rpkm.nc

