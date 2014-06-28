cd /home/lli/REMC/epiProfile/IR/
exons3p=/projects/epigenomics/acarles/intron/breast/4GloriaAnalysis/hg19v65_introns_for_genes.Id4Gloria.bothStrands.3prime
exons5p=/projects/epigenomics/acarles/intron/breast/4GloriaAnalysis/hg19v65_introns_for_genes.Id4Gloria.bothStrands.5prime
less $exons3p | awk '{id=$4; gsub("ENSG[0-9chrXYM:_-]+<", "", $4); print $1"\t"$2"\t"$3"\t"$4"\t"id}' > /home/lli/hg19/hg19v65_introns_for_genes.3prime_200
less $exons5p | awk '{id=$4; gsub("ENSG[0-9chrXYM:_-]+<", "", $4); print $1"\t"$2"\t"$3"\t"$4"\t"id}' > /home/lli/hg19/hg19v65_introns_for_genes.5prime_200
############################################################################################
# WGBS @ exon boundaries
exons3p=/home/lli/hg19/hg19v65_introns_for_genes.3prime_200
exons5p=/home/lli/hg19/hg19v65_introns_for_genes.5prime_200
#bed=/projects/edcc_prj2/bs-seq/a22478/bismark/A22478_4_lanes_dupsFlagged.fractional.bedGraph.gz; name=lumRM066_bismark;
bed=/projects/edcc_prj2/bs-seq/a18473/A18473_5_lanes_dupsFlagged.bam.fractional.bedGraph.gz; name=myoRM045_bismark;
out=/home/lli/REMC/epiProfile/IR/exons3p_200/; reg=$exons3p
#out=/home/lli/REMC/epiProfile/IR/exons5p_200/; reg=$exons5p
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromBEDCalculator.jar -b $bed -r $reg -o $out -s hg19 -n $name -bin 20 -t Y -trueBED

