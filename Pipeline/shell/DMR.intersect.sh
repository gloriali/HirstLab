#!/bin/sh

# intersect with genomic regions
if [ "$1" == "-h" ] ; then
    echo -e "Intersecting DM CpGs with genomic regions and report fold enrichment for each genomic feature 
Usage: `basename $0` -d <dirIn> -r <region> -n <region_name>
    <dirIn>: input directory with DMR bed files;
    <region>: customized region;
    <region_name>: name of the customized region."
    exit 0
fi

region=''
region_name=''

while [ $# -gt 0 ]
do
    case "$1" in
        -d) dirIn="$2"; shift;;
        -r) region="$2"; shift;;
        -n) region_name="$2"; shift;;
    esac
    shift
done
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirOut=$dirIn/intersect/
mkdir -p $dirOut
cd $dirIn
echo -e "Output files:
$dirOut/DMR.<name>.CpG.bed: chr\tstart\tend\tDMR_ID
$dirOut/DMR.<name>.gene.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.<name>.exon.bed: chr\tstart\tend\tDMR_IDtExon_ID
$dirOut/DMR.<name>.promoter.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.<name>.gene_pc.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.<name>.CGI.bed: chr\tstart\tend\tDMR_ID\tCGI_ID
$dirOut/DMR.<name>.CGI_shore.bed: chr\tstart\tend\tDMR_ID\tCGI_shore_ID
$dirOut/DMR.<name>.$region_name.bed: chr\tstart\tend\tDMR_ID\tregion_ID
$dirOut/genomic.breakdown.summary (fold enrichment): Name\tTotal No. of CpGs\tIntergenic\tIntron\tExon\tGene\tPromoter\tCGI\tCGI_shore\t$region_name\n"
all_total=28217448 # wc -l /home/lli/hg19/CG.BED
all_gene=16755181 # $BEDTOOLS/intersectBed -a /home/lli/hg19/CG.BED -b /home/lli/hg19/hg19v69_genes.bed -u | wc -l
all_exon=2707705 # $BEDTOOLS/intersectBed -a /home/lli/hg19/CG.BED -b /home/lli/hg19/hg19v69_exons_for_genes.bed -u | wc -l
all_promoter=3888980 # $BEDTOOLS/intersectBed -a /home/lli/hg19/CG.BED -b /home/lli/hg19/hg19v69_genes_TSS_2000.bed -u | wc -l
all_CGI=2089538 # $BEDTOOLS/intersectBed -a /home/lli/hg19/CG.BED -b /home/lli/hg19/CGI.forProfiles.BED -u | wc -l
all_CGI_shore=2288095 # $BEDTOOLS/intersectBed -a /home/lli/hg19/CG.BED -b /home/lli/hg19/CGI.2000shores.BED -u |wc -l
if [ -f "$region" ]; then
    all_region=$(less $region | sed 's/chr//g' | $BEDTOOLS/intersectBed -a /home/lli/hg19/CG.BED -b stdin -u | wc -l)
fi
echo -e "$all_total\t$all_gene\t$all_exon\t$all_promoter\t$all_CGI\t$all_CGI_shore\t$all_region" | awk '{print "Total No. of CpGs: "$1"\nExpectedPercentage:\n-Intergenic: "($1-$2)/$1"\n-Intron: "($2-$3)/$1"\n-Exon: "$3/$1"\n-Gene: "$2/$1"\n-Promoter: "$4/$1"\n-CGI: "$5/$1"\n-CGI shore: "$6/$1"\n-""'$region_name'"": "$7/$1}' 
echo -e "Name\tNCpG\tIntergenic\tIntron\tExon\tGene\tPromoter\tCGI\tCGI_shore\t$region_name" > $dirOut/genomic.breakdown.summary
for dmr in DMR.*.bed; do
    name=$(echo $dmr | sed -e s/'DMR.'//g | sed -e s/'.bed'//g)
    echo "Processing $name" 
    less $dmr | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dmr.tmp
    $BEDTOOLS/intersectBed -a $dmr.tmp -b /home/lli/hg19/hg19v69_genes.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.gene.bed
    $BEDTOOLS/intersectBed -a $dmr.tmp -b /home/lli/hg19/hg19v69_exons_for_genes.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.exon.bed
    $BEDTOOLS/intersectBed -a $dmr.tmp -b /home/lli/hg19/hg19v69_genes_TSS_2000.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.promoter.bed
    $BEDTOOLS/intersectBed -a $dmr.tmp -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.CGI.bed
    $BEDTOOLS/intersectBed -a $dmr.tmp -b /home/lli/hg19/CGI.2000shores.BED -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/DMR.$name.CGI_shore.bed
    less $dirOut/DMR.$name.gene.bed | awk '/protein_coding/ {gsub("_protein_coding", "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $dirOut/DMR.$name.gene_pc.bed
    $BEDTOOLS/intersectBed -a /home/lli/hg19/CG.BED -b $dmr.tmp -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$8"_"$1":"$2"-"$3}' > $dirOut/DMR.$name.CpG.bed
    total=$(less $dirOut/DMR.$name.CpG.bed | wc -l)
    gene=$($BEDTOOLS/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v69_genes.bed -u | wc -l)
    exon=$($BEDTOOLS/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v69_exons_for_genes.bed -u | wc -l)
    promoter=$($BEDTOOLS/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/hg19v69_genes_TSS_2000.bed -u | wc -l)
    CGI=$($BEDTOOLS/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/CGI.forProfiles.BED -u | wc -l)
    CGI_shore=$($BEDTOOLS/intersectBed -a $dirOut/DMR.$name.CpG.bed -b /home/lli/hg19/CGI.2000shores.BED -u | wc -l)
    if [ -f "$region" ]; then
        less $region | sed 's/chr//g'| $BEDTOOLS/intersectBed -a $dmr.tmp -b stdin -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5":"$6"-"$7}' > $dirOut/DMR.$name.$region_name.bed
        DMR_region=$(less $region | sed 's/chr//g' | $BEDTOOLS/intersectBed -a $dirOut/DMR.$name.CpG.bed -b stdin -u | wc -l)
        echo -e "$name\t$total\t$gene\t$exon\t$promoter\t$CGI\t$CGI_shore\t$DMR_region\t$all_total\t$all_gene\t$all_exon\t$all_promoter\t$all_CGI\t$all_CGI_shore\t$all_region" | awk '{print $1"\t"$2"\t"(($2-$3)/$2)/(($9-$10)/$9)"\t"(($3-$4)/$2)/(($10-$11)/$9)"\t"($4/$2)/($11/$9)"\t"($3/$2)/($10/$9)"\t"($5/$2)/($12/$9)"\t"($6/$2)/($13/$9)"\t"($7/$2)/($14/$9)"\t"($8/$2)/($15/$9)}' >> $dirOut/genomic.breakdown.summary
    else
        echo -e "$name\t$total\t$gene\t$exon\t$promoter\t$CGI\t$CGI_shore\t$all_total\t$all_gene\t$all_exon\t$all_promoter\t$all_CGI\t$all_CGI_shore" | awk '{print $1"\t"$2"\t"(($2-$3)/$2)/(($8-$9)/$8)"\t"(($3-$4)/$2)/(($9-$10)/$8)"\t"($4/$2)/($10/$8)"\t"($3/$2)/($9/$8)"\t"($5/$2)/($11/$8)"\t"($6/$2)/($12/$8)"\t"($7/$2)/($13/$8)}' >> $dirOut/genomic.breakdown.summary        
    fi
    rm $dmr.tmp
done

