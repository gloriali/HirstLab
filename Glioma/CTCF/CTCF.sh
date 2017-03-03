#!/bin/sh

dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/bam/
mkdir -p $dirOut
cd /projects/glioma_dart/glimoa_ctcf/
for bam in *.bam; do
	ln -s /projects/glioma_dart/glimoa_ctcf/$bam $dirOut/$bam
done
cat mgh59_idh1m_ctcf_10/subfile1.bin mgh59_idh1m_ctcf_10/subfile2.bin > $dirOut/mgh59_idh1m_ctcf_10.bam
cat mgh76_ctcf_19/subfile1.bin mgh76_ctcf_19/subfile2.bin > $dirOut/mgh76_ctcf_19.bam
cat mgh80_ctcf_21/subfile1.bin mgh80_ctcf_21/subfile2.bin > $dirOut/mgh80_ctcf_21.bam
cat mgh81_idh1m_ctcf_22/subfile1\ \(1\).bin mgh81_idh1m_ctcf_22/subfile2.bin > $dirOut/mgh81_idh1m_ctcf_22.bam

## FindER
JAVA=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/PET_200/map_old/
dirIn=/projects/epigenomics2/users/lli/glioma/CTCF/bam/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/
mkdir -p $dirOut
cd $dirIn
for bam in *.bam; do
	lib=$(echo $bam | sed -e 's/.bam//g')
	echo $lib
	$samtools index $bam
	$bamstats -g 2864785220 -q 10 -b $bam > $lib.bamstats
	/home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn $lib $lib.bamstats
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn
for bam in *ctcf*.bam; do
	echo $bam
	$JAVA -jar -Xmx10G /home/mbilenky/bin/Solexa_Java/FindER.0.9.3b.jar -i $bam -r $reg -o $dirOut -v -m $map -info -TF > $dirOut/FindER.$bam.log
done
for bam in *h3k27ac*.bam; do
	echo $bam
	$JAVA -jar -Xmx10G /home/mbilenky/bin/Solexa_Java/FindER.0.9.3b.jar -i $bam -r $reg -o $dirOut -v -m $map -info -minER 300 > $dirOut/FindER.$bam.log
done
