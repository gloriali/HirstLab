#!/bin/sh

# intersect with genomic regions
## instruction message
if [ "$1" == "-h" ] ; then
    echo -e "Intersecting regions of interest with genomic regions and report fold enrichment for each genomic feature 
Usage: `basename $0` -d <dirIn> -r <region> -n <region_name>
    <dirIn>: input directory with DMR bed files;
    <region>: customized region, e.g. enhancers;
    <region_name>: name of the customized region."
    exit 0
fi

## command line arguments
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

## path for reference files
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
gene=/home/lli/hg19/hg19v69_genes.bed
exon=/home/lli/hg19/hg19v69_exons_for_genes.bed
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
CGI_promoter=/home/lli/hg19/CGI.promoter.BED # $BEDTOOLS/intersectBed -a /home/lli/hg19/CGI.forProfiles.BED -b /home/lli/hg19/hg19v69_genes_TSS_2000.bed -u > /home/lli/hg19/CGI.promoter.BED
CGI_genebody=/home/lli/hg19/CGI.genebody.BED # $BEDTOOLS/intersectBed -a /home/lli/hg19/CGI.forProfiles.BED -b /home/lli/hg19/hg19v69_genes.bed -u > /home/lli/hg19/CGI.genebody.BED
CGI_intergenic=/home/lli/hg19/CGI.intergenic.BED # $BEDTOOLS/intersectBed -a /home/lli/hg19/CGI.forProfiles.BED -b /home/lli/hg19/CGI.promoter.BED /home/lli/hg19/CGI.genebody.BED -v > /home/lli/hg19/CGI.intergenic.BED
CGI_shore=/home/lli/hg19/CGI.2000shores.BED
dirOut=$dirIn/intersect/
mkdir -p $dirOut
cd $dirIn
echo -e "Processing all DMR.<name>.bed files in the input directory. Reference: hg19, Ensembl hg19v69. 
Output files:
$dirOut/genomic.breakdown.summary (fold enrichment for plotting): Name\tTotal length\tIntergenic\tIntron\tExon\tGene\tPromoter\tCGI_promoter\tCGI_genebody\tCGI_intergenic\tCGI_shore\t$region_name
$dirOut/DMR.<name>.gene.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.<name>.exon.bed: chr\tstart\tend\tDMR_ID\tExon_ID
$dirOut/DMR.<name>.promoter.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.<name>.gene_pc.bed: chr\tstart\tend\tDMR_ID\tEnsembl
$dirOut/DMR.<name>.CGI_promoter.bed: chr\tstart\tend\tDMR_ID\tCGI_ID
$dirOut/DMR.<name>.CGI_genebody.bed: chr\tstart\tend\tDMR_ID\tCGI_ID
$dirOut/DMR.<name>.CGI_intergenic.bed: chr\tstart\tend\tDMR_ID\tCGI_ID
$dirOut/DMR.<name>.CGI_shore.bed: chr\tstart\tend\tDMR_ID\tCGI_shore_ID
$dirOut/DMR.<name>.$region_name.bed: chr\tstart\tend\tDMR_ID\tregion_ID\n"

## genomic features region size
all_total=3095677412 # less /home/lli/hg19/hg19.chrom.len.autoXY | awk '{s=s+$2}END{print s}'
all_gene=1712708126 # less $gene | awk '{s=s+$3-$2}END{print s}'
all_exon=120342490 # less $exon | awk '{s=s+$3-$2}END{print s}'
all_promoter=223792038 # less $promoter | awk '{s=s+$3-$2}END{print s}'
all_CGI_promoter=14298323 # less $CGI_promoter | awk '{s=s+$3-$2}END{print s}'
all_CGI_genebody=18842834 # less $CGI_genebody | awk '{s=s+$3-$2}END{print s}'
all_CGI_intergenic=2624167 # less $CGI_intergenic | awk '{s=s+$3-$2}END{print s}'
all_CGI_shore=114702620 # less $CGI_shore | awk '{s=s+$3-$2}END{print s}'
if [ -f "$region" ]; then
    all_region=$(less $region | awk '{s=s+$3-$2}END{print s}')
fi

## intersect DMRs with regions and calculate fold enrichment
echo -e "$all_total\t$all_gene\t$all_exon\t$all_promoter\t$all_CGI_promoter\t$all_CGI_genebody\t$all_CGI_intergenic\t$all_CGI_shore\t$all_region" | awk '{print "Total No. of CpGs: "$1"\nExpectedPercentage:\n-Intergenic: "($1-$2)/$1"\n-Intron: "($2-$3)/$1"\n-Exon: "$3/$1"\n-Gene: "$2/$1"\n-Promoter: "$4/$1"\n-CGI_promoter: "$5/$1"\n-CGI_genebody: "$6/$1"\n-CGI_intergenic: "$7/$1"\n-CGI shore: "$8/$1"\n-""'$region_name'"": "$9/$1}' 
echo -e "Name\tNCpG\tIntergenic\tIntron\tExon\tGene\tPromoter\tCGI_promoter\tCGI_genebody\tCGI_intergenic\tCGI_shore\t$region_name" > $dirOut/genomic.breakdown.summary
for dmr in *.bed; do
    name=$(echo $dmr | sed 's/.bed//g')
    echo "Processing $name" 
    less $dmr | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dmr.tmp
    $BEDTOOLS/intersectBed -a $dmr.tmp -b $gene -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/$name.gene.bed
    $BEDTOOLS/intersectBed -a $dmr.tmp -b $exon -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/$name.exon.bed
    $BEDTOOLS/intersectBed -a $dmr.tmp -b $promoter -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/$name.promoter.bed
    $BEDTOOLS/intersectBed -a $dmr.tmp -b $CGI_promoter -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/$name.CGI_promoter.bed
    $BEDTOOLS/intersectBed -a $dmr.tmp -b $CGI_genebody -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/$name.CGI_genebody.bed
    $BEDTOOLS/intersectBed -a $dmr.tmp -b $CGI_intergenic -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/$name.CGI_intergenic.bed
    $BEDTOOLS/intersectBed -a $dmr.tmp -b $CGI_shore -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/$name.CGI_shore.bed
    less $dirOut/DMR.$name.gene.bed | awk '/protein_coding/ {gsub("_protein_coding", "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $dirOut/$name.gene_pc.bed
    dmr_total=$(less $dmr | awk '{s=s+$3-$2}END{print s}')
    dmr_gene=$(less $dirOut/$name.gene.bed | awk '{s=s+$3-$2}END{print s}')
    dmr_exon=$(less $dirOut/$name.exon.bed | awk '{s=s+$3-$2}END{print s}')
    dmr_promoter=$(less $dirOut/$name.promoter.bed | awk '{s=s+$3-$2}END{print s}')
    dmr_CGI_promoter=$(less $dirOut/$name.CGI_promoter.bed | awk '{s=s+$3-$2}END{print s}')
    dmr_CGI_genebody=$(less $dirOut/$name.CGI_genebody.bed | awk '{s=s+$3-$2}END{print s}')
    dmr_CGI_intergenic=$(less $dirOut/$name.CGI_intergenic.bed | awk '{s=s+$3-$2}END{print s}')
    dmr_CGI_shore=$(less $dirOut/$name.CGI_shore.bed | awk '{s=s+$3-$2}END{print s}')
    if [ -f "$region" ]; then
        less $region | sed 's/chr//g'| $BEDTOOLS/intersectBed -a $dmr.tmp -b stdin -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5":"$6"-"$7}' > $dirOut/$name.$region_name.bed
        DMR_region=$(less $dirOut/$name.$region_name.bed | awk '{s=s+$3-$2}END{print s}')
        echo -e "$name\t$dmr_total\t$dmr_gene\t$dmr_exon\t$dmr_promoter\t$dmr_CGI_promoter\t$dmr_CGI_genebody\t$dmr_CGI_intergenic\t$dmr_CGI_shore\t$DMR_region\t$all_total\t$all_gene\t$all_exon\t$all_promoter\t$all_CGI_promoter\t$all_CGI_genebody\t$all_CGI_intergenic\t$all_CGI_shore\t$all_region" | awk '{print $1"\t"$2"\t"(($2-$3)/$2)/(($11-$12)/$11)"\t"(($3-$4)/$2)/(($12-$13)/$11)"\t"($4/$2)/($13/$11)"\t"($3/$2)/($12/$11)"\t"($5/$2)/($14/$11)"\t"($6/$2)/($15/$11)"\t"($7/$2)/($16/$11)"\t"($8/$2)/($17/$11)"\t"($9/$2)/($18/$11)"\t"($10/$2)/($19/$11)}' >> $dirOut/genomic.breakdown.summary
    else
        echo -e "$name\t$dmr_total\t$dmr_gene\t$dmr_exon\t$dmr_promoter\t$dmr_CGI_promoter\t$dmr_CGI_genebody\t$dmr_CGI_intergenic\t$dmr_CGI_shore\t$all_total\t$all_gene\t$all_exon\t$all_promoter\t$all_CGI_promoter\t$all_CGI_genebody\t$all_CGI_intergenic\t$all_CGI_shore" | awk '{print $1"\t"$2"\t"(($2-$3)/$2)/(($10-$11)/$10)"\t"(($3-$4)/$2)/(($11-$12)/$10)"\t"($4/$2)/($12/$10)"\t"($3/$2)/($11/$10)"\t"($5/$2)/($13/$10)"\t"($6/$2)/($14/$10)"\t"($7/$2)/($15/$10)"\t"($8/$2)/($16/$10)"\t"($9/$2)/($17/$10)}' >> $dirOut/genomic.breakdown.summary        
    fi
    rm $dmr.tmp
done
