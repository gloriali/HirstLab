#!/bin/sh

# create Bed file with exon boundaries binned to 20bp bins
less /projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' > /home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique
exons3p=/home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools makewindows -b $exons3p -w 20 -i src > /home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique_bin_20
rm $exons3p
exons3p=/home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique_bin_20
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools nuc -fi /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -bed $exons3p > /home/lli/REMC/epiProfile/exons3p_200/GC.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt
less /home/lli/REMC/epiProfile/exons3p_200/GC.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt | awk 'BEGIN {id=""; record=""} !/#/ {if ($4 == id){record=(record" "$6)} else {print record; id=$4; record=(id" "$6)}}' > /home/lli/REMC/epiProfile/exons3p_200/GC.hg19v65_exons_for_genes.3prime_200.unique
rm /home/lli/REMC/epiProfile/exons3p_200/GC.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt

less /projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' > /home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique
exons5p=/home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools makewindows -b $exons5p -w 20 -i src > /home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique_bin_20
rm $exons5p
exons5p=/home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique_bin_20
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools nuc -fi /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -bed $exons5p > /home/lli/REMC/epiProfile/exons5p_200/GC.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt
less /home/lli/REMC/epiProfile/exons5p_200/GC.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt | awk 'BEGIN {id=""; record=""} !/#/ {if ($4 == id){record=(record" "$6)} else {print record; id=$4; record=(id" "$6)}}' > /home/lli/REMC/epiProfile/exons5p_200/GC.hg19v65_exons_for_genes.5prime_200.unique
rm /home/lli/REMC/epiProfile/exons5p_200/GC.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt
