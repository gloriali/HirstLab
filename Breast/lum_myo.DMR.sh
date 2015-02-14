#!/bin/sh

dirOut='/home/lli/REMC/UMR'
# intersect with genomic regions
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b /projects/mbilenky/REMC/breast/WGBS/pash/analysis/newDMRs/lum_myo/DMRs.p0.0005.s200.c3.-1 -wa > $dirOut/lum.UMR.CpG.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/lum.UMR.CpG.bed -b /home/lli/hg19/hg19v65_genes.bed -wa > $dirOut/lum.UMR.CpG_gene.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/lum.UMR.CpG.bed -b /home/lli/hg19/hg19v65_exons.bed -wa > $dirOut/lum.UMR.CpG_exon.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/lum.UMR.CpG.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa > $dirOut/lum.UMR.CpG_promoter.bed
echo $dirOut
wc -l $dirOut/lum.UMR.CpG.bed
wc -l $dirOut/lum.UMR.CpG_gene.bed
wc -l $dirOut/lum.UMR.CpG_exon.bed
wc -l $dirOut/lum.UMR.CpG_promoter.bed


# intersect with genomic regions
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/CG.BED -b /projects/mbilenky/REMC/breast/WGBS/pash/analysis/newDMRs/lum_myo/DMRs.p0.0005.s200.c3.1 -wa > $dirOut/myo.UMR.CpG.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/myo.UMR.CpG.bed -b /home/lli/hg19/hg19v65_genes.bed -wa > $dirOut/myo.UMR.CpG_gene.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/myo.UMR.CpG.bed -b /home/lli/hg19/hg19v65_exons.bed -wa > $dirOut/myo.UMR.CpG_exon.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/myo.UMR.CpG.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa > $dirOut/myo.UMR.CpG_promoter.bed
echo $dirOut
wc -l $dirOut/myo.UMR.CpG.bed
wc -l $dirOut/myo.UMR.CpG_gene.bed
wc -l $dirOut/myo.UMR.CpG_exon.bed
wc -l $dirOut/myo.UMR.CpG_promoter.bed

