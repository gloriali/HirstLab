#!/bin/sh

# link files
dirIn=/projects/edcc_new/reference_epigenomes/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/glioma/WGS/
for lib in 19 21 22 23 47 73 74 75 76 78 79 81; do
    IDH=$(less $dirOut/../samples.txt | awk '{if($1=="CEMT_""'$lib'")print $2}')
    echo $IDH $lib 
    mkdir -p $dirOut/bam/
    mkdir -p $dirOut/VCF/
    mkdir -p $dirOut/CNV/
    ln -s $dirIn/CEMT_$lib/bams/WGS/*.bam $dirOut/bam/$IDH.CEMT_$lib.bam
    ln -s $dirIn/CEMT_$lib/bams/WGS/*.eff.dbSNP_v137.COSMIC_v64.annotations.vcf $dirOut/VCF/$IDH.CEMT_$lib.eff.dbSNP_v137.COSMIC_v64.annotations.vcf
    ln -s $dirIn/CEMT_$lib/bams/WGS/CNVs/*.bam_CNVs.corr.*list*.sorted $dirOut/CNV/$IDH.CEMT_$lib.CNVs.corr.blacklistFiltered.sorted
done
ln -s /projects/edcc_prj2/upstream_data/P00015_8_lanes_dupsFlagged.bam $dirOut/bam/IDHmut.CEMT_47.bam

# Confirm IDH1 R132H heterozygous mutation (chr2:209113112 C->T)
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/glioma/WGS/bam/
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/sambamba-bamStats
for bam in $dirIn/*.bam; do
    echo $(basename $bam)
    $sambamba index $bam -t 8
    $bamstats -g 2864785220 -t 8 $bam > $(echo $bam | sed 's/\.bam//').bamstats
done
> $dirIn/IDH.txt
for file in $dirIn/*.bam; do
    sample=$(basename $file | sed 's/\.bam//'); echo $sample
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools mpileup -f /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -r 2:209113112-209113112 $file | awk '{s=toupper($5); print "'$sample'""\tWGS\t"gsub("T", "",s)"\t"gsub(",", "", s)+gsub(".", "", s)}' | awk '{print $0"\t"$3/($3+$4)}' >> $dirIn/IDH.txt
done
dirRNA=/projects/epigenomics3/epigenomics3_results/users/lli/glioma/RNAseq/bam/
for file in $dirRNA/*.bam; do
    sample=$(basename $file | sed 's/\.bam//'); echo $sample
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools mpileup -f /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -r 2:209113112-209113112 $file | awk '{s=toupper($5); print "'$sample'""\tRNAseq\t"gsub("T", "",s)"\t"gsub(",", "", s)+gsub(".", "", s)}' | awk '{print $0"\t"$3/($3+$4)}' >> $dirIn/IDH.txt
done
for V in $dirIn/../VCF/*.vcf; do 
    echo $V;
    cat $V | grep -E 'IDH1' 
done

# Confirm common mutations in gliomas
dirVCF='/projects/edcc_new/reference_epigenomes/'
dirOut='/projects/epigenomics2/users/lli/glioma/WGS/VCF/'
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
mkdir -p $dirOut
for lib in 19 21 22 23; do
    ln -s $dirVCF/CEMT_$lib/bams/WGS/*dbSNP_v137.COSMIC_v64.annotations.vcf $dirOut/CEMT_$lib.vcf
    less $dirOut/CEMT_$lib.vcf | awk '!/^#/ {print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"_"$9"_"$10}' > $dirOut/CEMT_$lib.vcf.bed
    $BEDTOOLS/intersectBed -a $dirOut/glioma_mutations.bed -b $dirOut/CEMT_$lib.vcf.bed -wa -wb > $dirOut/CEMT_$lib.vcf.glioma_mutations.txt
done
lib=47
ln -s $dirVCF/CEMT_$lib/bams/WGS/P00015*dbSNP_v137.COSMIC_v64.annotations.vcf $dirOut/CEMT_$lib.vcf
less $dirOut/CEMT_$lib.vcf | awk '!/^#/ {print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"_"$9"_"$10}' > $dirOut/CEMT_$lib.vcf.bed
$BEDTOOLS/intersectBed -a $dirOut/glioma_mutations.bed -b $dirOut/CEMT_$lib.vcf.bed -wa -wb > $dirOut/CEMT_$lib.vcf.glioma_mutations.txt
# intersecting with COSMIC annotation
cd $dirOut
for file in *.glioma_mutations.txt; do
    lib=$(echo $file | sed -e 's/.txt//g')
    less $file | awk '$8~/COSM/ {print $0}' > $dirOut/$lib.COSMIC.txt
done

# Somatic mutations in CEMT_47 (only one with matched blood sample)
dirVCF='/projects/edcc_new/reference_epigenomes/'
dirOut='/projects/epigenomics2/users/lli/glioma/WGS/VCF/'
ln -s $dirVCF/CEMT_47/bams/WGS/P00017.4_lanes_dupsFlagged.varFilter.eff.dbSNP_v137.COSMIC_v64.annotations.vcf $dirOut/CB.CEMT_47.vcf
less $dirOut/CB.CEMT_47.vcf | awk '!/^#/ {print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"_"$9"_"$10}' > $dirOut/CB.CEMT_$lib.vcf.bed
cd /projects/epigenomics2/users/lli/glioma/WGS/VCF/
awk 'NR==FNR {id=$1":"$2; h[id]=1} {i=$1":"$2; if(!(i in h)){print $0}}' CB.CEMT_47.vcf.bed CEMT_47.vcf.bed > CEMT_47.somatic.vcf.bed
less CEMT_47.somatic.vcf.bed | awk '$0~/COSM/ {print $0}' > CEMT_47.somatic.COSMIC.vcf.bed
less CEMT_47.somatic.COSMIC.vcf.bed | awk '$0~/NON_SYNONYMOUS_CODING/ {print $0}' > CEMT_47.somatic.COSMIC.non_synonymous.vcf.bed
## TET2?
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
$BEDTOOLS/intersectBed -a glioma_mutations.bed -b CEMT_47.somatic.vcf.bed -wa -wb | grep 'TET2' # frameshift INDEL @ 4:106190839-106190840
/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools mpileup -f /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -r 4:106190839-106190840 /projects/epigenomics2/users/lli/glioma/WGS/bam/P00015.CEMT_47.bam
/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools mpileup -f /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -r 4:106190839-106190840 /projects/epigenomics2/users/lli/glioma/WGS/bam/P00017.CEMT_47.bam

# CNVs
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirIn='/projects/epigenomics2/users/lli/glioma/WGS/CNV/'
mkdir -p $dirIn
cd $dirIn
for file in *.CNV; do
    echo $file
    less $file | awk '!/no_change/ {print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4"\t"$5}' > $file.bed
    $BEDTOOLS/intersectBed -a $file.bed -b /home/lli/hg19/cytoBand.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10}' > $file.cytoband
    $BEDTOOLS/intersectBed -a $file.bed -b /home/lli/hg19/hg19v69_genes.bed -wa -wb | awk '{gsub("_", "\t"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10"\t"$11}' > $file.gene.tmp
    awk 'NR==FNR {name[$1]=$2; next} {print $0"\t"name[$7]}' /projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.EnsID_sorted.HUGO $file.gene.tmp > $file.gene
    rm $file.gene.tmp
done

