#!/bin/sh

dirIn="/projects/epigenomics/ep50/internal/jqc.1.7.6/"
chrsize="/home/lli/hg19/hg19.chrom.sizes"
dirOut="/gsc/www/bcgsc.ca/downloads/mb/BrainHubs/RNAseqHub/hg19/"
mkdir -p $dirOut
for lib in A03484 A03473 A03474 A07825 A03475 A03476 A04599 A15295 A15298 A15299
do
    cd $dirIn/$lib/wig/
    for file in *.wig.gz
    do
        name=$(echo $file | sed -e 's/.wig.gz//g')
        echo "Processing" $name
        /home/lli/HirstLab/Pipeline/UCSC/wigToBigWig $file $chrsize $dirOut/$name.bw
    done
done

