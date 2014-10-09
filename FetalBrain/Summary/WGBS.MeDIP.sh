#!/bin/sh

# Compare WGBS and MeDIP DMRs - HuFNSC02
dirOut='/projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/'
dirWGBS='/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/'
dirMeDIP='/projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/'
mkdir -p $dirOut
cd $dirOut
# intersect DM CpGs
less $dirWGBS/DM.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.bed | awk '{if($4==1){print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3 >> "'$dirOut'""/WGBS.Cortex_GE.hyper.bed"} else{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3 >> "'$dirOut'""/WGBS.Cortex_GE.hypo.bed"}}'
less $dirMeDIP/DM.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.bed | awk '{if($4==1){print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3 >> "'$dirOut'""/MeDIP.Cortex_GE.hyper.bed"} else{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3 >> "'$dirOut'""/MeDIP.Cortex_GE.hypo.bed"}}'
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/WGBS.Cortex_GE.hyper.bed -b $dirOut/MeDIP.Cortex_GE.hyper.bed -wa > $dirOut/DM.hyper.intersect.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/WGBS.Cortex_GE.hyper.bed -b $dirOut/MeDIP.Cortex_GE.hyper.bed -v > $dirOut/DM.hyper.WGBS.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/MeDIP.Cortex_GE.hyper.bed -b $dirOut/WGBS.Cortex_GE.hyper.bed -v > $dirOut/DM.hyper.MeDIP.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/WGBS.Cortex_GE.hypo.bed -b $dirOut/MeDIP.Cortex_GE.hypo.bed -wa > $dirOut/DM.hypo.intersect.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/WGBS.Cortex_GE.hypo.bed -b $dirOut/MeDIP.Cortex_GE.hypo.bed -v > $dirOut/DM.hypo.WGBS.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/MeDIP.Cortex_GE.hypo.bed -b $dirOut/WGBS.Cortex_GE.hypo.bed -v > $dirOut/DM.hypo.MeDIP.bed
## genomic breakdown
mkdir -p $dirOut/CpG/
> $dirOut/CpG/genomic.breakdown.summary
for file in DM.*.bed
do
    name=$(echo $file | sed -e s/'.bed'//g)
    echo "Processing $name" 
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $file -b /home/lli/hg19/hg19v65_genes.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/CpG/$name.CpG_gene.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $file -b /home/lli/hg19/hg19v65_exons.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5":"$6"-"$7}' > $dirOut/CpG/$name.CpG_exon.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $file -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $dirOut/CpG/$name.CpG_promoter.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $file -b /projects/epigenomics/resources/UCSC_hg19/CGI/CGI.forProfiles.BED -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$9}' > $dirOut/CpG/$name.CpG_CGI.bed
    less $dirOut/CpG/$name.CpG_gene.bed | awk '/protein_coding/ {gsub("_protein_coding", "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $dirOut/CpG/$name.CpG_gene_pc.bed
    less $dirOut/CpG/$name.CpG_promoter.bed | awk '/protein_coding/ {gsub("_protein_coding", "", $5); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $dirOut/CpG/$name.CpG_promoter_pc.bed
    total=`wc -l $dirOut/$file | cut -d' ' -f 1`
    gene=`less $dirOut/CpG/$name.CpG_gene.bed | awk '{if(!($4 in h)){s=s+1; h[$4]=1}} END{print s}'`
    exon=`less $dirOut/CpG/$name.CpG_exon.bed | awk '{if(!($4 in h)){s=s+1; h[$4]=1}} END{print s}'`
    promoter=`less $dirOut/CpG/$name.CpG_promoter.bed | awk '{if(!($4 in h)){s=s+1; h[$4]=1}} END{print s}'`
    CGI=`wc -l $dirOut/CpG/$name.CpG_CGI.bed | cut -d' ' -f 1`
    echo -e "$name\t$total\t$gene\t$exon\t$promoter\t$CGI" | awk '{print $1"\t"$2"\t"($2-$3)/$2"\t"($3-$4)/$2"\t"$4/$2"\t"$3/$2"\t"$5/$2"\t"$6/$2}' >> $dirOut/CpG/genomic.breakdown.summary
done

# intersect DMRs
WGBS_Cortex_UMR=$dirWGBS'/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed'
WGBS_GE_UMR=$dirWGBS'/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed'
MeDIP_Cortex_UMR=$dirMeDIP'/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.s300.c4.hypo.bed'
MeDIP_GE_UMR=$dirMeDIP'/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.s300.c4.hyper.bed'
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $WGBS_Cortex_UMR -b $MeDIP_Cortex_UMR -wo > $dirOut/Cortex_UMR.intersect
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $WGBS_Cortex_UMR -b $MeDIP_Cortex_UMR -v > $dirOut/Cortex_UMR.WGBS
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $MeDIP_Cortex_UMR -b $WGBS_Cortex_UMR -v > $dirOut/Cortex_UMR.MeDIP
less $dirOut/Cortex_UMR.intersect | awk '{if(!($4 in h)){c=c+1; h[$4]=1}}END{print "Cortex_UMR intersect", c}'
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $WGBS_GE_UMR -b $MeDIP_GE_UMR -wo > $dirOut/GE_UMR.intersect
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $WGBS_GE_UMR -b $MeDIP_GE_UMR -v > $dirOut/GE_UMR.WGBS
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $MeDIP_GE_UMR -b $WGBS_GE_UMR -v > $dirOut/GE_UMR.MeDIP
less $dirOut/GE_UMR.intersect | awk '{if(!($4 in h)){c=c+1; h[$4]=1}}END{print "GE_UMR intersect", c}'
## GC density
for file in *UMR*
do
    echo "Processing $file" 
    less $file | awk '{gsub("chr", "", $1); print $1"\t"$2"\t"$3"\t"$4}' > $file.bed
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/bedtools nuc -fi /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa -bed $file.bed > GC.$file.txt
done

# intersect closest genes
Gene_WGBS_Cortex_UMR=$dirWGBS'/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.closest.gene'
Gene_WGBS_GE_UMR=$dirWGBS'/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.closest.gene'
Gene_MeDIP_Cortex_UMR=$dirMeDIP'/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.s300.c4.hypo.closest.gene'
Gene_MeDIP_GE_UMR=$dirMeDIP'/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.s300.c4.hyper.closest.gene'
awk 'NR==FNR {h[$8]=$4'_'$9; next} {if($8 in h){print $8"\t"h[$8]"\t"$4'_'$9}}' $Gene_MeDIP_Cortex_UMR $Gene_WGBS_Cortex_UMR > $dirOut/Cortex_UMR_closest_gene.intersect
awk 'NR==FNR {h[$8]=$4'_'$9; next} {if(!($8 in h)){print $8"\t"$4'_'$9}}' $Gene_MeDIP_Cortex_UMR $Gene_WGBS_Cortex_UMR > $dirOut/Cortex_UMR_closest_gene.WGBS
awk 'NR==FNR {h[$8]=$4'_'$9; next} {if(!($8 in h)){print $8"\t"$4'_'$9}}' $Gene_WGBS_Cortex_UMR $Gene_MeDIP_Cortex_UMR > $dirOut/Cortex_UMR_closest_gene.MeDIP
awk 'NR==FNR {h[$8]=$4'_'$9; next} {if($8 in h){print $8"\t"h[$8]"\t"$4'_'$9}}' $Gene_MeDIP_GE_UMR $Gene_WGBS_GE_UMR > $dirOut/GE_UMR_closest_gene.intersect
awk 'NR==FNR {h[$8]=$4'_'$9; next} {if(!($8 in h)){print $8"\t"$4'_'$9}}' $Gene_MeDIP_GE_UMR $Gene_WGBS_GE_UMR > $dirOut/GE_UMR_closest_gene.WGBS
awk 'NR==FNR {h[$8]=$4'_'$9; next} {if(!($8 in h)){print $8"\t"$4'_'$9}}' $Gene_WGBS_GE_UMR $Gene_MeDIP_GE_UMR > $dirOut/GE_UMR_closest_gene.MeDIP



