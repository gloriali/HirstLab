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
echo -e "Intersecting DM CpGs with genomic regions\nOutput directory: $dirOut\nOutput DM CpG percentage to genomic.breakdown.summary\ngenomic.breakdown.summary format: Name\tTotal No. of DM CpGs\tIntergenic\tIntron\tExon\tGene\tPromoter"
for dmr in DMR.*.bed
do
    name=`echo $dmr | cut -d'.' -f 2`
    echo "Processing $name" 
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b $dmr -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$7}' > $dirOut/DMR.$name.CpG.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.CpG_gene.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v65_exons.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5":"$6"-"$7}' > $dirOut/DMR.$name.CpG_exon.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.CpG_promoter.bed
    total=`wc -l $dirOut/DMR.$name.CpG.bed | cut -d' ' -f 1`
    gene=`wc -l $dirOut/DMR.$name.CpG_gene.bed | cut -d' ' -f 1`
    exon=`wc -l $dirOut/DMR.$name.CpG_exon.bed | cut -d' ' -f 1`
    promoter=`wc -l $dirOut/DMR.$name.CpG_promoter.bed | cut -d' ' -f 1`
    echo -e "$name\t$total\t$gene\t$exon\t$promoter" | awk '{print $1"\t"$2"\t"($2-$3)/$2"\t"($3-$4)/$2"\t"$4/$2"\t"$3/$2"\t"$5/$2}' >> $dirOut/genomic.breakdown.summary
done

