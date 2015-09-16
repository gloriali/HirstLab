#!/bin/sh

cd /projects/epigenomics/users/lli/DNM3
for file in *.bam; do
    echo $file
    #/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar /gsc/software/linux-x86_64-centos5/picard-tools-1.92/CollectAlignmentSummaryMetrics.jar INPUT=$file OUTPUT=$file.stats
    t=$(less $file.stats | awk 'NR==8{print $2}')
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools view $file 'chr1:171841498-172418466' -b | /projects/epigenomics/software/bedtools-2.23.0/bin/bamToBed -i stdin | sort -k1,1 -k2,2n -T /projects/epigenomics/temp/ > DNM3.$file.bed
    /projects/epigenomics/software/bedtools-2.23.0/bin/intersectBed -a ./ensembl/DNM3_ENST.bed -b DNM3.$file.bed -wa -c | awk '{print $0"\t""'$t'""\t"$3-$2"\t"($6*10^9/($3-$2)/"'$t'")}' > DNM3.$file.exon.RPKM.bed
    less DNM3.$file.exon.RPKM.bed | awk '{read[$5]=read[$5]+$6; len[$5]=len[$5]+$8; total=$7} END{for(id in read){print id"\t"read[id]"\t"len[id]"\t"total"\t"(read[id]*10^9/len[id]/total)}}' > DNM3.$file.transcript.RPKM
done
