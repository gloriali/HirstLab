#!/bin/sh

# create Bed file with exon boundaries binned to 20bp bins 
less /projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' > /home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique
exons3p=/home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools makewindows -b $exons3p -w 20 -i src > /home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique_bin_20
rm $exons3p
exons3p=/home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique_bin_20
less /projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' > /home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique
exons5p=/home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools makewindows -b $exons5p -w 20 -i src > /home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique_bin_20
rm $exons5p
exons5p=/home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique_bin_20


# GC content profile @ exon boundaries
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools nuc -fi /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -bed $exons3p > /home/lli/REMC/epiProfile/exons3p_200/GC.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt
less /home/lli/REMC/epiProfile/exons3p_200/GC.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt | awk 'BEGIN {id=""; record=""} !/#/ {if ($4 == id){record=(record" "$6)} else {print record; id=$4; record=(id" "$6)}}' > /home/lli/REMC/epiProfile/exons3p_200/GC.hg19v65_exons_for_genes.3prime_200.unique
rm /home/lli/REMC/epiProfile/exons3p_200/GC.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools nuc -fi /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -bed $exons5p > /home/lli/REMC/epiProfile/exons5p_200/GC.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt
less /home/lli/REMC/epiProfile/exons5p_200/GC.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt | awk 'BEGIN {id=""; record=""} !/#/ {if ($4 == id){record=(record" "$6)} else {print record; id=$4; record=(id" "$6)}}' > /home/lli/REMC/epiProfile/exons5p_200/GC.hg19v65_exons_for_genes.5prime_200.unique
rm /home/lli/REMC/epiProfile/exons5p_200/GC.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt


# No. of CpGs instead of GC% @ exon boundaries
exons3p=/home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique_bin_20
exons5p=/home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique_bin_20
less $exons3p | awk '{print "chr"$0}' > /home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique_bin_20.1
mv -f /home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique_bin_20.1 $exons3p
less $exons5p | awk '{print "chr"$0}' > /home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique_bin_20.1
mv -f /home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique_bin_20.1 $exons5p
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $exons3p -b /home/lli/hg19/CG.BED -c > /home/lli/REMC/epiProfile/exons3p_200/CpG.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt   
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $exons5p -b /home/lli/hg19/CG.BED -c > /home/lli/REMC/epiProfile/exons5p_200/CpG.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt   
less /home/lli/REMC/epiProfile/exons3p_200/CpG.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt | awk 'BEGIN {id=""; record=""} {if ($4 == id){record=(record" "$5)} else {print record; id=$4; record=(id" "$5)}}' > /home/lli/hg19/CpG.hg19v65_exons_for_genes.3prime_200.unique
rm /home/lli/REMC/epiProfile/exons3p_200/CpG.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt
less /home/lli/REMC/epiProfile/exons5p_200/CpG.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt | awk 'BEGIN {id=""; record=""} {if ($4 == id){record=(record" "$5)} else {print record; id=$4; record=(id" "$5)}}' > /home/lli/hg19/CpG.hg19v65_exons_for_genes.5prime_200.unique
rm /home/lli/REMC/epiProfile/exons5p_200/CpG.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt



# No. of CpGs profile @ intron boundaries
introns3p=/home/lli/hg19/hg19v65_introns_for_genes.3prime_200
introns5p=/home/lli/hg19/hg19v65_introns_for_genes.5prime_200
less $introns3p | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' > /home/lli/hg19/hg19v65_introns_for_genes.3prime_200_1
less $introns5p | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' > /home/lli/hg19/hg19v65_introns_for_genes.5prime_200_1
introns3p=/home/lli/hg19/hg19v65_introns_for_genes.3prime_200_1
introns5p=/home/lli/hg19/hg19v65_introns_for_genes.5prime_200_1
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools makewindows -b $introns3p -w 20 -i src > /home/lli/hg19/hg19v65_introns_for_genes.3prime_200_bin_20
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools makewindows -b $introns5p -w 20 -i src > /home/lli/hg19/hg19v65_introns_for_genes.5prime_200_bin_20
rm $introns3p; rm $introns5p;
introns3p=/home/lli/hg19/hg19v65_introns_for_genes.3prime_200_bin_20
introns5p=/home/lli/hg19/hg19v65_introns_for_genes.5prime_200_bin_20
less $introns3p | awk '{print "chr"$0}' > /home/lli/hg19/hg19v65_introns_for_genes.3prime_200_bin_20.1
mv -f /home/lli/hg19/hg19v65_introns_for_genes.3prime_200_bin_20.1 $introns3p
less $introns5p | awk '{print "chr"$0}' > /home/lli/hg19/hg19v65_introns_for_genes.5prime_200_bin_20.1
mv -f /home/lli/hg19/hg19v65_introns_for_genes.5prime_200_bin_20.1 $introns5p
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $introns3p -b /home/lli/hg19/CG.BED -c > /home/lli/REMC/epiProfile/IR/introns3p_200/CpG.hg19v65_introns_for_genes.3prime_200_bin_20.txt   
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $introns5p -b /home/lli/hg19/CG.BED -c > /home/lli/REMC/epiProfile/IR/introns5p_200/CpG.hg19v65_introns_for_genes.5prime_200_bin_20.txt   
less /home/lli/REMC/epiProfile/IR/introns3p_200/CpG.hg19v65_introns_for_genes.3prime_200_bin_20.txt | awk 'BEGIN {id=""; record=""} {if ($4 == id){record=(record" "$5)} else {print record; id=$4; record=(id" "$5)}}' > /home/lli/hg19/CpG.hg19v65_introns_for_genes.3prime_200
rm /home/lli/REMC/epiProfile/IR/introns3p_200/CpG.hg19v65_introns_for_genes.3prime_200_bin_20.txt
less /home/lli/REMC/epiProfile/IR/introns5p_200/CpG.hg19v65_introns_for_genes.5prime_200_bin_20.txt | awk 'BEGIN {id=""; record=""} {if ($4 == id){record=(record" "$5)} else {print record; id=$4; record=(id" "$5)}}' > /home/lli/hg19/CpG.hg19v65_introns_for_genes.5prime_200
rm /home/lli/REMC/epiProfile/IR/introns5p_200/CpG.hg19v65_introns_for_genes.5prime_200_bin_20.txt

