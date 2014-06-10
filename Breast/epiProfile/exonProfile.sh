#!/bin/sh

cd /home/lli/REMC/epiProfile

#bed=/projects/edcc_prj2/bs-seq/a22478/A22478.Cmethyl.cons.bed.gz; name=lumRM066_novoalign;

#bed=/projects/edcc_prj2/bs-seq/a22478/bismark/A22478_4_lanes_dupsFlagged.fractional.bedGraph.gz; name=lumRM066_bismark;
bed=/projects/edcc_prj2/bs-seq/a18473/A18473_5_lanes_dupsFlagged.bam.fractional.bedGraph.gz; name=myoRM045_bismark;
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique

#out=/home/lli/REMC/epiProfile/exons3p_200/; reg=$exons3p
out=/home/lli/REMC/epiProfile/exons5p_200/; reg=$exons5p

java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromBEDCalculator.jar -b $bed -r $reg -o $out -s hg19 -n $name -bin 20 -t Y -trueBED
