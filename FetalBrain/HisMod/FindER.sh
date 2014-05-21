#!/bin/sh

cd /home/lli/FetalBrain/HisMod/wigs/
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
map=/projects/epigenomics/resources/UCSC_hg19/mappability/mappability_lt_1
out=/home/lli/FetalBrain/HisMod/FindER/

for inputWig in *.wig.gz
do
    echo "Processing $inputWig"
    name=`echo $inputWig | cut -d'.' -f 1`
    java -jar -Xmx30G /home/mbilenky/bin/Solexa_Java/FindER.jar -chr $chr -i $inputWig -o $out -n $name -v -m $map -print > /home/lli/FetalBrain/HisMod/FindER/$name.log
done
