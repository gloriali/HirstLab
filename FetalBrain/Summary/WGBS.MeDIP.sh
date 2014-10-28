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

# hydroxymethylation: high mC in WGBS & low in MeDIP - HuFNSC02
dirOut='/projects/epigenomics/users/lli/FetalBrain/WGBS_MeDIP/'
cutoff=0.8
## potential CpGs
### Cortex
cd $dirOut
Cortex02_WGBS='/projects/epigenomics/users/lli/FetalBrain/WGBS/A22475.WGBS.NeurospheresCortex02.sam.bedGraph.combine'
Cortex02_MeDIP='/projects/epigenomics/users/lli/FetalBrain/MeDIP/HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174.dip'
awk 'NR==FNR {h[$1]=$4; next} {chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; id=chr":"start"-"end; if(id in h){diff=h[id]-$2; print chr"\t"start"\t"end"\t"h[id]"\t"$2"\t"diff;}}' $Cortex02_WGBS $Cortex02_MeDIP > $dirOut/Cortex02_WGBS_MeDIP.fractional
less $dirOut/Cortex02_WGBS_MeDIP.fractional | awk '{s=s+1; diff[int($6*10)]=diff[int($6*10)]+1; } END{for(key in diff){print key"\t"diff[key]; sum=sum+diff[key]}; if(s!=sum){print "ERROR!";}}' | sort -k1,1n > $dirOut/Cortex02_WGBS_MeDIP.diff.summary
less $dirOut/Cortex02_WGBS_MeDIP.fractional | awk '{if($6>="'$cutoff'"){print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}}' | sort -k1,1 -k2,2n > $dirOut/Cortex02_WGBS_MeDIP.diff.bed
### GE
GE02_WGBS='/projects/epigenomics/users/lli/FetalBrain/WGBS/A17784-A13819.WGBS.NeurospheresGE02.sam.bedGraph.combine'
GE02_MeDIP='/projects/epigenomics/users/lli/FetalBrain/MeDIP/HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166.dip'
awk 'NR==FNR {h[$1]=$4; next} {chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; id=chr":"start"-"end; if(id in h){diff=h[id]-$2; print chr"\t"start"\t"end"\t"h[id]"\t"$2"\t"diff;}}' $GE02_WGBS $GE02_MeDIP > $dirOut/GE02_WGBS_MeDIP.fractional
less $dirOut/GE02_WGBS_MeDIP.fractional | awk '{s=s+1; diff[int($6*10)]=diff[int($6*10)]+1; } END{for(key in diff){print key"\t"diff[key]; sum=sum+diff[key]}; if(s!=sum){print "ERROR!";}}' | sort -k1,1n > $dirOut/GE02_WGBS_MeDIP.diff.summary
less $dirOut/GE02_WGBS_MeDIP.fractional | awk '{if($6>="'$cutoff'"){print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}}' | sort -k1,1 -k2,2n > $dirOut/GE02_WGBS_MeDIP.diff.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/Cortex02_WGBS_MeDIP.diff.bed -b $dirOut/GE02_WGBS_MeDIP.diff.bed | wc -l 
## intersect with Cortex02 vs GE02 DM CpGs and UMRs
### WGBS DM CpGs 
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/Cortex02_WGBS_MeDIP.diff.bed -b /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DM.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.bed | wc -l
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/GE02_WGBS_MeDIP.diff.bed -b /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DM.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.bed | wc -l
### WGBS UMRs
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/Cortex02_WGBS_MeDIP.diff.bed -b /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed | wc -l
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/GE02_WGBS_MeDIP.diff.bed -b /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed | wc -l
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed -b $dirOut/Cortex02_WGBS_MeDIP.diff.bed -u | wc -l
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed -b $dirOut/GE02_WGBS_MeDIP.diff.bed -u | wc -l
### MeDIP DM CpGs 
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/Cortex02_WGBS_MeDIP.diff.bed -b /projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/DM.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.bed | wc -l
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/GE02_WGBS_MeDIP.diff.bed -b /projects/epigenomics/users/lli/FetalBrain/MeDIP/DMR/DM.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.d0.6.bed | wc -l
## genomic breakdown
> $dirOut/hydroxy.breakdown.summary
for file in Cortex02_WGBS_MeDIP.diff.bed GE02_WGBS_MeDIP.diff.bed
do
    name=$(echo $file | sed -e s/'_.*'//g)
    echo 'Processing '$name
    total=`wc -l $dirOut/$file | cut -d' ' -f 1`
    gene=`/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/$file -b /home/lli/hg19/hg19v65_genes.bed -u | wc -l | cut -d' ' -f 1`
    exon=`/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/$file -b /home/lli/hg19/hg19v65_exons.bed -u | wc -l | cut -d' ' -f 1`
    promoter=`/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/$file -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -u | wc -l | cut -d' ' -f 1`
    CGI=`/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/$file -b /projects/epigenomics/resources/UCSC_hg19/CGI/CGI.forProfiles.BED -u | wc -l | cut -d' ' -f 1`
    echo -e "$name\t$total\t$gene\t$exon\t$promoter\t$CGI" | awk '{print $1"\t"$2"\t"($2-$3)/$2"\t"($3-$4)/$2"\t"$4/$2"\t"$3/$2"\t"$5/$2"\t"$6/$2}' >> $dirOut/hydroxy.breakdown.summary
done
## Collapse potential hydroxy CpGs into regions
less $dirOut/Cortex02_WGBS_MeDIP.diff.bed | awk 'BEGIN{size=300; cut=3} {if($2<end+size && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"c"\t"l}}' > $dirOut/Cortex02_WGBS_MeDIP.region
less $dirOut/GE02_WGBS_MeDIP.diff.bed | awk 'BEGIN{size=300; cut=3} {if($2<end+size && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"c"\t"l}}' > $dirOut/GE02_WGBS_MeDIP.region
less $dirOut/Cortex02_WGBS_MeDIP.region | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $dirOut/Cortex02_WGBS_MeDIP.region.bed
less $dirOut/GE02_WGBS_MeDIP.region | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $dirOut/GE02_WGBS_MeDIP.region.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirOut/Cortex02_WGBS_MeDIP.region.bed -b $dirOut/GE02_WGBS_MeDIP.region.bed | wc -l
## Closest gene for potential hydroxy regions
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/closestBed -a $dirOut/Cortex02_WGBS_MeDIP.region.bed -b /home/lli/hg19/hg19v65_genes.bed -d > $dirOut/Cortex02_WGBS_MeDIP.region.closest.gene
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/closestBed -a $dirOut/GE02_WGBS_MeDIP.region.bed -b /home/lli/hg19/hg19v65_genes.bed -d > $dirOut/GE02_WGBS_MeDIP.region.closest.gene

