#!/bin/sh

# IS DMRs on chrX
dirIn='/projects/mbilenky/REMC/breast/hg19/MeDIP/analysis/CpG_coverage_10bp/individual/lum/0.5'
dirOut='/home/lli/REMC/IS.DMR/chrX/lum'
#dirIn='/projects/mbilenky/REMC/breast/hg19/MeDIP/analysis/CpG_coverage_10bp/individual/myo/0.5'
#dirOut='/home/lli/REMC/IS.DMR/chrX/myo'

cd $dirIn
for name in {1000,0111,0100,1011,0010,1101,0001,1110}
do
    file=$name.DMR.BED
    echo "Processing $file"
    less $dirIn/$file | awk -v name=$name '/X/ { print $0"\t"name }' > $dirOut/$name.DMR.X.BED
done
cat $dirOut/*.DMR.X.BED > $dirOut/IS.DMR.X.bed
# intersect with genomic regions
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b $dirOut/IS.DMR.X.bed -wa > $dirOut/IS.DMR.X.CpG.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/IS.DMR.X.CpG.bed -b /home/lli/hg19/hg19v65_genes.bed -wa > $dirOut/IS.DMR.X.CpG_gene.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/IS.DMR.X.CpG.bed -b /home/lli/hg19/hg19v65_exons.bed -wa > $dirOut/IS.DMR.X.CpG_exon.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/IS.DMR.X.CpG.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa > $dirOut/IS.DMR.X.CpG_promoter.bed
echo $dirOut
wc -l $dirOut/IS.DMR.X.CpG.bed
wc -l $dirOut/IS.DMR.X.CpG_gene.bed
wc -l $dirOut/IS.DMR.X.CpG_exon.bed
wc -l $dirOut/IS.DMR.X.CpG_promoter.bed


# all IS DMRs
#dirIn='/projects/mbilenky/REMC/breast/hg19/MeDIP/analysis/CpG_coverage_10bp/individual/lum/0.5'
#dirOut='/home/lli/REMC/IS.DMR/lum'
dirIn='/projects/mbilenky/REMC/breast/hg19/MeDIP/analysis/CpG_coverage_10bp/individual/myo/0.5'
dirOut='/home/lli/REMC/IS.DMR/myo'

cd $dirIn
for name in {1000,0111,0100,1011,0010,1101,0001,1110}
do
    file=$name.DMR.BED
    echo "Processing $file"
    less $dirIn/$file | awk -v name=$name '{ print $0"\t"name }' > $dirOut/$name.DMR.BED
done
cat $dirOut/*.DMR.BED > $dirOut/IS.DMR.bed
# intersect with genomic regions
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b $dirOut/IS.DMR.bed -wa > $dirOut/IS.DMR.CpG.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/IS.DMR.CpG.bed -b /home/lli/hg19/hg19v65_genes.bed -wa > $dirOut/IS.DMR.CpG_gene.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/IS.DMR.CpG.bed -b /home/lli/hg19/hg19v65_exons.bed -wa > $dirOut/IS.DMR.CpG_exon.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/IS.DMR.CpG.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa > $dirOut/IS.DMR.CpG_promoter.bed
echo $dirOut
wc -l $dirOut/IS.DMR.CpG.bed
wc -l $dirOut/IS.DMR.CpG_gene.bed
wc -l $dirOut/IS.DMR.CpG_exon.bed
wc -l $dirOut/IS.DMR.CpG_promoter.bed

