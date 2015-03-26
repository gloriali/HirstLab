#!/bin/sh

chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
reg=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes_TSS_2000
out=/projects/epigenomics/users/lli/FetalBrain/ChIPseq/signal/
mkdir -p $out
cd /home/lli/FetalBrain/HisMod/wigs/

for wig in *.wig.gz
do
    name=$(echo $wig | cut -d'.' -f1)
    echo "Processing $name"
    /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $out -s $chr -n $name > $out$name.log
done
