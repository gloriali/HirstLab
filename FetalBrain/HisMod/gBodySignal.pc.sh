#!/bin/sh

cd /home/lli/FetalBrain/HisMod/wigs/
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
reg=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes_body.pc
out=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/signal/

for wig in A03273.bam.q5.F1028.SET_219.wig.gz A03279.bam.q5.F1028.SET_229.wig.gz A03285.bam.q5.F1028.SET_238.wig.gz A03481.bam.q5.F1028.SET_215.wig.gz A03489.bam.q5.F1028.SET_175.wig.gz A03497.bam.q5.F1028.SET_233.wig.gz A19307.q5.F1028.SET_164.wig.gz
do
    name=$(echo $wig | cut -d'.' -f1)
    echo "Processing $name"
    /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $out -s $chr -n $name > $out"pc."$name.log
done
