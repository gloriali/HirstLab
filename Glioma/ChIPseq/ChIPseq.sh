#!/bin/sh

# FindER 1.0.0b results
dirIn=/projects/edcc_new/reference_epigenomes/
dirOut=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/
for lib in 19 21 22 23 47; do
    for mark in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3 H3K27ac; do
        echo $lib $mark
        mkdir -p $dirOut/$mark/
        ln -s $dirIn/CEMT_$lib/bams/ChIP-Seq/$mark/FindER.1.0.0b/*.FDR_0.05.FindER.bed.gz $dirOut/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz   
    done
done

echo -e "Mark\tSample\tN_region\tTotal_length" > $dirOut/ER.summary.stats
for mark in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3 H3K27ac; do
    cd $dirOut/$mark/
    for file in *.bed.gz; do
        lib=$(echo $file | sed -e 's/.FDR_0.05.FindER.bed.gz//g' | sed -e 's/.*Input.//g' | sed -e 's/GE04/NPC_GE04/g')
        echo -e $mark"\t"$lib"\t"$(less $file | wc -l)"\t"$(less $file | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER.summary.stats
    done
done

