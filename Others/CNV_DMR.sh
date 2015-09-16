#!/bin/sh

cd /projects/edcc_prj2/bs-seq/a34002
paste <(zcat A34002_3_lanes_dupsFlagged.q5.f0.5mC.CpG.coverage.bedGraph.gz) <(zcat A34002_3_lanes_dupsFlagged.q5.f0.5mC.CpG.fractional.bedGraph.gz) | awk '{printf "%s\t%d\t%d\t%d\t%0.0f\t%0.0f\n", $1, $2, $3, $8, $4*$8/10, $4-$4*$8/10}' > /home/ysun/CEMT/DMR/A34002.CEMT_7.sam.bedGraph
/home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i /home/ysun/CEMT/DMR/ -o /home/ysun/CEMT/DMR/ -f A34002.CEMT_7.sam.bedGraph

cd /projects/edcc_prj2/bs-seq/a34004
paste <(zcat A34004_3_lanes_dupsFlagged.q5.f0.5mC.CpG.coverage.bedGraph.gz) <(zcat A34004_3_lanes_dupsFlagged.q5.f0.5mC.CpG.fractional.bedGraph.gz) | awk '{printf "%s\t%d\t%d\t%d\t%0.0f\t%0.0f\n", $1, $2, $3, $8, $4*$8/10, $4-$4*$8/10}' > /home/ysun/CEMT/DMR/A34004.CEMT_9.sam.bedGraph
/home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i /home/ysun/CEMT/DMR/ -o /home/ysun/CEMT/DMR/ -f A34004.CEMT_9.sam.bedGraph

cd /projects/edcc_prj2/bs-seq/a36002
paste <(zcat A36002_3_lanes_dupsFlagged.q5.f0.5mC.CpG.coverage.bedGraph.gz) <(zcat A36002_3_lanes_dupsFlagged.q5.f0.5mC.CpG.fractional.bedGraph.gz) | awk '{printf "%s\t%d\t%d\t%d\t%0.0f\t%0.0f\n", $1, $2, $3, $8, $4*$8/10, $4-$4*$8/10}' > /home/ysun/CEMT/DMR/A36002.CEMT_19.sam.bedGraph
/home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i /home/ysun/CEMT/DMR/ -o /home/ysun/CEMT/DMR/ -f A36002.CEMT_19.sam.bedGraph

cd /home/ysun/CEMT/DMR/
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i /home/ysun/CEMT/DMR/ -o /home/ysun/CEMT/DMR/ -f1 A34004.CEMT_9.sam.bedGraph.combine -f2 A36002.CEMT_19.sam.bedGraph.combine -n A34004-CEMT9_A36002-CEMT19
/home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i /home/ysun/CEMT/DMR/ -o /home/ysun/CEMT/DMR/ -f DM.A34004-CEMT9_A36002-CEMT19.m0.75.p0.0005.d0.6.bed -n A34004-CEMT9_A36002-CEMT19 -s 300 -c 3
for file in DMR.A34004-CEMT9_A36002-CEMT19*
do
    less $file | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $file.bed
done
wc -l DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.*bed
less DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.bed | awk '{s=s+$3-$2}END{print s}' # No. of bp
less DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.hyper.bed | awk '{s=s+$3-$2}END{print s}' # No. of bp
less DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.hypo.bed | awk '{s=s+$3-$2}END{print s}' # No. of bp

cd /home/ysun/CEMT/CNV/
less /projects/wtsspipeline/jira_tickets/GNRG-2612/CEMT_test/output/A24724/A24724.meta_bwa.hg19a.A24724_3_lanes_dupsFlagged.bam_ratio.BedGraph | awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4}' > /home/ysun/CEMT/CNV/A24724-CEMT_9.CNV.bedgraph
less /projects/wtsspipeline/jira_tickets/GNRG-2612/CEMT_test/output/A29507/A29507.A29507_2_lanes_dupsFlagged.bam_ratio.BedGraph | awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4}' > /home/ysun/CEMT/CNV/A29507-CEMT_19.CNV.bedgraph
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/unionBedGraph s -i /home/ysun/CEMT/CNV/A24724-CEMT_9.CNV.bedgraph /home/ysun/CEMT/CNV/A29507-CEMT_19.CNV.bedgraph -filler 2 | awk '{if(sqrt(($4-$5)^2)>1){print $1"\t"$2"\t"$3"\t"$4-$5}}' > A24724-CEMT_9-A29507-CEMT_19.CNV
less A24724-CEMT_9-A29507-CEMT_19.CNV | awk '{if($4<0)print $0}' > A24724-CEMT_9-A29507-CEMT_19.CNV.loss
less A24724-CEMT_9-A29507-CEMT_19.CNV | awk '{if($4>0)print $0}' > A24724-CEMT_9-A29507-CEMT_19.CNV.gain
wc -l A24724-CEMT_9-A29507-CEMT_19.CNV* # No. of regions
less A24724-CEMT_9-A29507-CEMT_19.CNV | awk '{s=s+$3-$2}END{print s}' # No. of bp
less A24724-CEMT_9-A29507-CEMT_19.CNV.gain | awk '{s=s+$3-$2}END{print s}' # No. of bp
less A24724-CEMT_9-A29507-CEMT_19.CNV.loss | awk '{s=s+$3-$2}END{print s}' # No. of bp

/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/ysun/CEMT/CNV/A24724-CEMT_9-A29507-CEMT_19.CNV -b /home/ysun/CEMT/DMR/DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.bed > /home/ysun/CEMT/CNV_DMR/CNV_DMR.CEMT9_CEMT19.bed
wc -l CNV_DMR.CEMT9_CEMT19.bed
less CNV_DMR.CEMT9_CEMT19.bed | awk '{s=s+$3-$2}END{print s}' # No. of bp
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/ysun/CEMT/CNV/A24724-CEMT_9-A29507-CEMT_19.CNV.gain -b /home/ysun/CEMT/DMR/DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.hyper.bed > /home/ysun/CEMT/CNV_DMR/CNV-gain_DMR-hyper.CEMT9_CEMT19.bed
wc -l CNV-gain_DMR-hyper.CEMT9_CEMT19.bed
less CNV-gain_DMR-hyper.CEMT9_CEMT19.bed | awk '{s=s+$3-$2}END{print s}' # No. of bp
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/ysun/CEMT/CNV/A24724-CEMT_9-A29507-CEMT_19.CNV.gain -b /home/ysun/CEMT/DMR/DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.hypo.bed > /home/ysun/CEMT/CNV_DMR/CNV-gain_DMR-hypo.CEMT9_CEMT19.bed
wc -l CNV-gain_DMR-hypo.CEMT9_CEMT19.bed
less CNV-gain_DMR-hypo.CEMT9_CEMT19.bed | awk '{s=s+$3-$2}END{print s}' # No. of bp
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/ysun/CEMT/CNV/A24724-CEMT_9-A29507-CEMT_19.CNV.gain -b /home/ysun/CEMT/DMR/DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.bed > /home/ysun/CEMT/CNV_DMR/CNV-gain_DMR.CEMT9_CEMT19.bed
wc -l CNV-gain_DMR.CEMT9_CEMT19.bed
less CNV-gain_DMR.CEMT9_CEMT19.bed | awk '{s=s+$3-$2}END{print s}' # No. of bp
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/ysun/CEMT/CNV/A24724-CEMT_9-A29507-CEMT_19.CNV.loss -b /home/ysun/CEMT/DMR/DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.hyper.bed > /home/ysun/CEMT/CNV_DMR/CNV-loss_DMR-hyper.CEMT9_CEMT19.bed
wc -l CNV-loss_DMR-hyper.CEMT9_CEMT19.bed
less CNV-loss_DMR-hyper.CEMT9_CEMT19.bed | awk '{s=s+$3-$2}END{print s}' # No. of bp
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/ysun/CEMT/CNV/A24724-CEMT_9-A29507-CEMT_19.CNV.loss -b /home/ysun/CEMT/DMR/DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.hypo.bed > /home/ysun/CEMT/CNV_DMR/CNV-loss_DMR-hypo.CEMT9_CEMT19.bed
wc -l CNV-loss_DMR-hypo.CEMT9_CEMT19.bed
less CNV-loss_DMR-hypo.CEMT9_CEMT19.bed | awk '{s=s+$3-$2}END{print s}' # No. of bp
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/ysun/CEMT/CNV/A24724-CEMT_9-A29507-CEMT_19.CNV.loss -b /home/ysun/CEMT/DMR/DMR.A34004-CEMT9_A36002-CEMT19.s300.c3.bed > /home/ysun/CEMT/CNV_DMR/CNV-loss_DMR.CEMT9_CEMT19.bed
wc -l CNV-loss_DMR.CEMT9_CEMT19.bed
less CNV-loss_DMR.CEMT9_CEMT19.bed | awk '{s=s+$3-$2}END{print s}' # No. of bp

# genome length excluding chrY
less /home/lli/hg19/hg19.chrom.sizes | awk 'NR<25&&NR!=21{s=s+$2}END{print s}'
# hypergeometric test R code: phyper(intersect, sample1, total-sample1, sample2, lower.tail = F, log.p = T)


# 2015. 9. 14
cd /home/ysun/CEMT/DMR/
less A36002.CEMT_19.sam.bedGraph.combine | awk '{gsub(":", "\t"); gsub("-","\t"); print $1"\t"$2"\t"$3"\t"$6"\t"$4+$5}' | sort -k1,1 -k2,2n -T /projects/epigenomics/temp/ > A36002.CEMT_19.cov.sam.bedGraph.combine
cd /home/ysun/CEMT/CNV/
less A29507-CEMT_19.CNV.bedgraph | awk '{if($4>1&&$1!="chrY"&&$1!="chrMT"&&(!($1~/GL/))){print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}}' > A29507-CEMT_19.CNV.gain
less A29507-CEMT_19.CNV.gain | awk '{len=$3-$2; if($2-len>0){print $1"\t"$2-len"\t"$2"\t"$1":"$2"-"$3}}' > A29507-CEMT_19.CNV.gain.up
less A29507-CEMT_19.CNV.gain | awk '{len=$3-$2; if($2-len>0){print $1"\t"$3"\t"$3+len"\t"$1":"$2"-"$3}}' > A29507-CEMT_19.CNV.gain.dn

/home/lli/HirstLab/Pipeline/shell/region.mean.sh -i /home/ysun/CEMT/DMR/ -o /home/ysun/CEMT/CNV/ -f A36002.CEMT_19.cov.sam.bedGraph.combine -r /home/ysun/CEMT/CNV/A29507-CEMT_19.CNV.gain -n A29507-CEMT_19.CNV.gain.5mC -w
/home/lli/HirstLab/Pipeline/shell/region.mean.sh -i /home/ysun/CEMT/DMR/ -o /home/ysun/CEMT/CNV/ -f A36002.CEMT_19.cov.sam.bedGraph.combine -r /home/ysun/CEMT/CNV/A29507-CEMT_19.CNV.gain.up -n A29507-CEMT_19.CNV.gain.up.5mC -w
/home/lli/HirstLab/Pipeline/shell/region.mean.sh -i /home/ysun/CEMT/DMR/ -o /home/ysun/CEMT/CNV/ -f A36002.CEMT_19.cov.sam.bedGraph.combine -r /home/ysun/CEMT/CNV/A29507-CEMT_19.CNV.gain.dn -n A29507-CEMT_19.CNV.gain.dn.5mC -w

cd /home/ysun/CEMT/CNV/
awk 'NR==FNR{h[$4]=$5; next} {if($4 in h){print $0"\t"h[$4]"\t"h[$4]-$5}}' A29507-CEMT_19.CNV.gain.up.5mC.weighted.mean.bed A29507-CEMT_19.CNV.gain.5mC.weighted.mean.bed > A29507-CEMT_19.CNV.up-gain.dm
awk 'NR==FNR{id=$1":"$2"-"$3; h[id]=$4} {if($4 in h){print $0"\t"h[$4]}}' A29507-CEMT_19.CNV.bedgraph A29507-CEMT_19.CNV.up-gain.dm > A29507-CEMT_19.CNV.up-gain.dm.CNV
awk 'NR==FNR{h[$4]=$5; next} {if($4 in h){print $0"\t"h[$4]"\t"h[$4]-$5}}' A29507-CEMT_19.CNV.gain.dn.5mC.weighted.mean.bed A29507-CEMT_19.CNV.gain.5mC.weighted.mean.bed > A29507-CEMT_19.CNV.dn-gain.dm
awk 'NR==FNR{id=$1":"$2"-"$3; h[id]=$4} {if($4 in h){print $0"\t"h[$4]}}' A29507-CEMT_19.CNV.bedgraph A29507-CEMT_19.CNV.dn-gain.dm > A29507-CEMT_19.CNV.dn-gain.dm.CNV


