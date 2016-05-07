#!/bin/sh

# Confirm IDH1 R132H heterozygous mutation (chr2:209113112 C->T)
cd /projects/epigenomics/users/smcconnell/glioma/WGS/bam/
for file in *.bam; do
    echo $file
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools index $file
done
> /projects/epigenomics2/users/lli/glioma/WGS/IDH.txt
for file in *.bam; do
    echo $file >> /projects/epigenomics2/users/lli/glioma/WGS/IDH.txt
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools mpileup -f /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -r 2:209100951-209130798 $file >> /projects/epigenomics2/users/lli/glioma/WGS/IDH.txt
done
cd /projects/edcc_new/reference_epigenomes/
for V in ./CEMT_19/bams/WGS/*varFilter.eff.snvs.vcf ./CEMT_2[1-3]/bams/WGS/*varFilter.eff.snvs.vcf ./CEMT_47/bams/WGS/*varFilter.eff.snvs.vcf; do 
    echo $V;
    cat $V | grep -E 'IDH1' 
done

# Confirm common mutations in gliomas
dirVCF='/projects/edcc_new/reference_epigenomes/'
dirOut='/projects/epigenomics2/users/lli/glioma/WGS/VCF/'
BEDTOOLS='/projects/epigenomics/software/bedtools-2.23.0/bin/'
mkdir -p $dirOut
for lib in 19 21 22 23 47; do
    ln -s $dirVCF/CEMT_$lib/bams/WGS/*dbSNP_v137.COSMIC_v64.annotations.vcf $dirOut/CEMT_$lib.vcf
    less $dirOut/CEMT_$lib.vcf | awk '!/^#/ {print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"_"$9"_"$10}' > $dirOut/CEMT_$lib.vcf.bed
    $BEDTOOLS/intersectBed -a $dirOut/glioma_mutations.bed -b $dirOut/CEMT_$lib.vcf.bed -wa -wb > $dirOut/CEMT_$lib.vcf.glioma_mutations.txt
done

