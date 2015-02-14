#!/bin/sh

cd /home/lli/FetalBrain/HisMod/signal/
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
reg=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_genes_TSS_2000
out=/home/lli/FetalBrain/HisMod/signal/

for wig in /home/lli/FetalBrain/HisMod/wigs/*.wig.gz
do
    name=$(basename $wig | sed 's/.gz//' | sed 's/.wig//')
    echo "Processing $name"
    /gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $out -s $chr -n $name > $out$name.log
done
