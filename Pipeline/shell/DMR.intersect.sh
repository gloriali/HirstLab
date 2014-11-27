#!/bin/sh

# intersect with genomic regions
if [ "$1" == "-h" ] ; then
    echo -e "Intersecting DM CpGs with genomic regions and report fold enrichment for each genomic feature 
Usage: `basename $0` -d <dirIn>
    <dirIn>: input directory with DMR bed files"
    exit 0
fi

while [ $# -gt 0 ]
do
    case "$1" in
        -d) dirIn="$2"; shift;
    esac
    shift
done
dirOut=$dirIn/CpG
mkdir -p $dirOut
cd $dirIn
echo -e "Output files:
$dirOut/DMR.<name>.CpG.bed: chr\tstart\tend\tDMR_ID
$dirOut/DMR.<name>.CpG_gene.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.<name>.CpG_exon.bed: chr\tstart\tend\tDMR_IDtExon_ID
$dirOut/DMR.<name>.CpG_promoter.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.<name>.CpG_gene_pc.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.<name>.CpG_promoter_pc.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.<name>.CpG_CGI.bed: chr\tstart\tend\tDMR_ID\tCGI_ID
$dirOut/genomic.breakdown.summary: Name\tTotal No. of CpGs\tIntergenic\tIntron\tExon\tGene\tPromoter\tCGI\n"
all_total=28217448 # wc -l /home/lli/hg19/CG.BED
all_gene=18270435 # /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b /home/lli/hg19/hg19v65_genes.bed -wa -wb | wc -l
all_exon=7064050 # /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b /home/lli/hg19/hg19v65_exons.bed -wa -wb | wc -l
all_promoter=3727169 # /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb | wc -l
all_CGI=2089538 # /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b /projects/epigenomics/resources/UCSC_hg19/CGI/CGI.forProfiles.BED -wa -wb | wc -l
> $dirOut/genomic.breakdown.summary
echo -e "$all_total\t$all_gene\t$all_exon\t$all_promoter\t$all_CGI" | awk '{print "ExpectedPercent\t"$1"\t"($1-$2)/$1"\t"($2-$3)/$1"\t"$3/$1"\t"$2/$1"\t"$4/$1"\t"$5/$1}' >> $dirOut/genomic.breakdown.summary
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
    echo -e "$name\t$total\t$gene\t$exon\t$promoter\t$CGI\t$all_total\t$all_gene\t$all_exon\t$all_promoter\t$all_CGI" | awk '{print $1"\t"$2"\t"(($2-$3)/$2)/(($7-$8)/$7)"\t"(($3-$4)/$2)/(($8-$9)/$7)"\t"($4/$2)/($9/$7)"\t"($3/$2)/($8/$7)"\t"($5/$2)/($10/$7)"\t"($6/$2)/($11/$7)}' >> $dirOut/genomic.breakdown.summary
done

