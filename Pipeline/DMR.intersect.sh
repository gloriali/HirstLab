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
Output files: $dirOut/DMR.$name.CpG.bed: chr\tstart\tend\tDMR_ID
$dirOut/DMR.$name.CpG_gene.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.$name.CpG_exon.bed: chr\tstart\tend\tDMR_IDtExon_ID
$dirOut/DMR.$name.CpG_promoter.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.$name.CpG_gene_pc.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.$name.CpG_promoter_pc.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.$name.CpG_CGI.bed: chr\tstart\tend\tDMR_ID\tCGI_ID
$dirOut/genomic.breakdown.summary: Name\tTotal No. of DM CpGs\tIntergenic\tIntron\tExon\tGene\tPromoter\tCGI"\n
for dmr in DMR.*.bed
do
    name=$(echo $dmr | sed -e s/'DMR.'//g)
    name=$(echo $name | sed -e s/'.bed'//g)
    echo "Processing $name" 
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b $dmr -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$7"_"$1":"$2"-"$3}' > $dirOut/DMR.$name.CpG.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.CpG_gene.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v65_exons.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5":"$6"-"$7}' > $dirOut/DMR.$name.CpG_exon.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.CpG_promoter.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /projects/epigenomics/resources/UCSC_hg19/CGI/CGI.forProfiles.BED -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > $dirOut/DMR.$name.CpG_CGI.bed
    less $dirOut/DMR.$name.CpG_gene.bed | awk '/protein_coding/ {gsub("_protein_coding", "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $dirOut/DMR.$name.CpG_gene_pc.bed
    less $dirOut/DMR.$name.CpG_promoter.bed | awk '/protein_coding/ {gsub("_protein_coding", "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $dirOut/DMR.$name.CpG_promoter_pc.bed
    total=`wc -l $dirOut/DMR.$name.CpG.bed | cut -d' ' -f 1`
    gene=`less $dirOut/DMR.$name.CpG_gene.bed | awk '{if(!($4 in h)){s=s+1; h[$4]=1}} END{print s}'`
    exon=`less $dirOut/DMR.$name.CpG_exon.bed | awk '{if(!($4 in h)){s=s+1; h[$4]=1}} END{print s}'`
    promoter=`less $dirOut/DMR.$name.CpG_promoter.bed | awk '{if(!($4 in h)){s=s+1; h[$4]=1}} END{print s}'`
    CGI=`wc -l $dirOut/DMR.$name.CpG_CGI.bed | cut -d' ' -f 1`
    echo -e "$name\t$total\t$gene\t$exon\t$promoter\t$CGI" | awk '{print $1"\t"$2"\t"($2-$3)/$2"\t"($3-$4)/$2"\t"$4/$2"\t"$3/$2"\t"$5/$2"\t"$6/$2}' >> $dirOut/genomic.breakdown.summary
done

