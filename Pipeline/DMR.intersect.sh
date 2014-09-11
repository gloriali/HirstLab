#!/bin/sh

# intersect with genomic regions
while [ $# -gt 0 ]
do
    case "$1" in
        -d) dirIn="$2"; shift;
    esac
    shift
done
dirOut=$dirIn/CpG
mkdir -p $dirOut
> $dirOut/genomic.breakdown.summary
cd $dirIn
echo -e "Intersecting DM CpGs with genomic regions
Output files: $dirOut/DMR.$name.CpG.bed: chr\tstart\tend\tDMR_ID>DM
$dirOut/DMR.$name.CpG_gene.bed: chr\tstart\tend\tDMR_ID>DM\tEnsembl
$dirOut/DMR.$name.CpG_exon.bed: chr\tstart\tend\tDMR_ID>DM\tExon_ID
$dirOut/DMR.$name.CpG_promoter.bed: chr\tstart\tend\tDMR_ID>DM\tEnsembl
$dirOut/DMR.$name.CpG_gene_pc.bed: chr\tstart\tend\tDMR_ID>DM\tEnsembl
$dirOut/DMR.$name.CpG_gene_promoter_pc.bed: chr\tstart\tend\tDMR_ID>DM\tEnsembl
$dirOut/genomic.breakdown.summary: Name\tTotal No. of DM CpGs\tIntergenic\tIntron\tExon\tGene\tPromoter"\n
for dmr in DMR.*.bed
do
    name=`echo $dmr | cut -d'.' -f 2`
    echo "Processing $name" 
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b $dmr -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$7">"$8}' > $dirOut/DMR.$name.CpG.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.CpG_gene.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v65_exons.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5":"$6"-"$7}' > $dirOut/DMR.$name.CpG_exon.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.CpG_promoter.bed
    less $dirOut/DMR.$name.CpG_gene.bed | awk '/protein_coding/ {gsub("_protein_coding", "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $dirOut/DMR.$name.CpG_gene_pc.bed
    less $dirOut/DMR.$name.CpG_promoter.bed | awk '/protein_coding/ {gsub("_protein_coding", "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $dirOut/DMR.$name.CpG_promoter_pc.bed
    total=`wc -l $dirOut/DMR.$name.CpG.bed | cut -d' ' -f 1`
    gene=`wc -l $dirOut/DMR.$name.CpG_gene.bed | cut -d' ' -f 1`
    exon=`wc -l $dirOut/DMR.$name.CpG_exon.bed | cut -d' ' -f 1`
    echo -e "$name\t$total\t$gene\t$exon\t$promoter" | awk '{print $1"\t"$2"\t"($2-$3)/$2"\t"($3-$4)/$2"\t"$4/$2"\t"$3/$2"\t"$5/$2}' >> $dirOut/genomic.breakdown.summary
done

