#!/bin/sh

export PATH=/home/lli/anaconda2/bin/:$PATH
export PYTHONPATH=/home/lli/anaconda2/lib/python2.7/site-packages
java=/gsc/software/linux-x86_64-centos6/jdk1.8.0_162/jre/bin/java
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirNPC=/projects/epigenomics2/users/lli/glioma/
dirRNA=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/
dirHM=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/ChIPseq/
dir5mC=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/PBAL/
dir5hmC=/projects/epigenomics3/epigenomics3_results/users/lli/NHA/hMeDIP/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/
mkdir -p $dirOut

## ChIPseq
dirIn=$dirOut/ChIPseq/bam/
for file in $dirNPC/ChIPseq/bam/*/*NPC*.bam $dirHM/bam/*/*MGG*.bam; do
    mark=$(echo $file | sed 's/.*bam\///' | sed 's/\/.*//')
    sample=$(basename $file | sed 's/.bam//' | sed 's/\./_/' | cut -d'_' -f2,3 | sed 's/119//')
    echo $mark $sample
    mkdir -p $dirIn/$mark/
    ln -s $file $dirIn/$mark/$mark"."$sample.bam
    ln -s $file.bai $dirIn/$mark/$mark"."$sample.bam.bai
done
for file in $dirNPC/ChIPseq/wig/*/*NPC*.PET.wig.gz $dirHM/wig/*/*MGG*.wig.gz; do
    mark=$(echo $file | sed 's/.*wig\///' | sed 's/\/.*//')
    sample=$(basename $file | sed 's/.q5.F1028.PET.wig.gz//' | sed 's/H.*_MGG/MGG/' | sed 's/119//' | sed 's/input_//')
    echo $mark $sample
    mkdir -p $dirOut/ChIPseq/wig/$mark/
    ln -s $file $dirOut/ChIPseq/wig/$mark/$mark"."$sample.wig.gz
done
### FindER
FindER=/home/mbilenky/bin/Solexa_Java/FindER.1.0.1e.jar
mkdir -p $dirOut/ChIPseq/FindER/
echo -e "Mark\tSample\tN_region\tTotal_length\tAverage_length" > $dirOut/ChIPseq/FindER/ER_summary.txt
for file in $dirIn/H*/*.bam; do
    mark=$(basename $file | cut -d'.' -f1)
    sample=$(basename $file | cut -d'.' -f2)
    echo $mark $sample
    mkdir -p $dirOut/ChIPseq/FindER/$mark/
    $java -jar -Xmx12G $FindER -signalBam $file -inputBam $dirIn/input/input.$sample.bam -out $dirOut/ChIPseq/FindER/$mark/ > $dirOut/ChIPseq/FindER/$mark/$mark.$sample.log
    echo -e $mark"\t"$sample"\t"$(less $dirOut/ChIPseq/FindER/$mark/$mark.$sample.vs.input.$sample.FDR_0.05.FindER.bed.gz | wc -l)"\t"$(less $dirOut/ChIPseq/FindER/$mark/$mark.$sample.vs.input.$sample.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"$4/$3}' >> $dirOut/ChIPseq/FindER/ER_summary.txt
done
### FindER2
finder2=/home/mbilenky/bin/FindER2/finder2.jar
mkdir -p $dirOut/ChIPseq/FindER2/
for inp in $dirIn/input/*.bam; do
    sample=$(basename $inp | cut -d'.' -f2);
    echo $sample
    sig1=$(echo $dirIn/H3K27me3/H3K27me3.$sample.bam)
    sig2=$(echo $dirIn/H3K4me3/H3K4me3.$sample.bam)
    sig3=$(echo $dirIn/H3K4me1/H3K4me1.$sample.bam)
    sig4=$(echo $dirIn/H3K27ac/H3K27ac.$sample.bam)
    sig5=$(echo $dirIn/H3K36me3/H3K36me3.$sample.bam)
    sig=${sig1},${sig2},${sig3},${sig4},${sig5}
    $java -jar -Xmx25G $finder2 inputBam:$inp signalBam:$sig outDir:$dirOut/ChIPseq/FindER2/ acgtDir:/projects/epigenomics2/users/mbilenky/resources/hg19/ACGT 
done
echo -e "Mark\tSample\tN_region\tTotal_length\tAverage_length" > $dirOut/ChIPseq/FindER2/ER_summary.txt
echo -e "chr\tstart\tend\tdis\tlen\tmark\tsample\tcategory" > $dirOut/ChIPseq/FindER2/ER.nearest.dis.bed
for file in $dirOut/ChIPseq/FindER2/*FindER2.bed; do
    mark=$(basename $file | cut -d'.' -f1)
    sample=$(basename $file | cut -d'.' -f2)
    echo $mark $sample
    $BEDTOOLS/closestBed -a $file -b $file -io -d | awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$3-$2"\t""'$mark'""\t""'$sample'""\toriginal"}' >> $dirOut/ChIPseq/FindER2/ER.nearest.dis.bed
    $BEDTOOLS/mergeBed -i $file -d 500 | awk '$1 !~ /GL/ {if($3-$2>=200){print}}' > $dirOut/ChIPseq/FindER2/$mark.$sample.FindER2.merge.bed
    $BEDTOOLS/closestBed -a $dirOut/ChIPseq/FindER2/$mark.$sample.FindER2.merge.bed -b $dirOut/ChIPseq/FindER2/$mark.$sample.FindER2.merge.bed -io -d | awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$3-$2"\t""'$mark'""\t""'$sample'""\tmerge"}' >> $dirOut/ChIPseq/FindER2/ER.nearest.dis.bed
    echo -e $mark"\t"$sample"\t"$(less $dirOut/ChIPseq/FindER2/$mark.$sample.FindER2.merge.bed | wc -l)"\t"$(less $dirOut/ChIPseq/FindER2/$mark.$sample.FindER2.merge.bed | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"$4/$3}' >> $dirOut/ChIPseq/FindER2/ER_summary.txt
done
for mark in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    intervene upset -i $dirOut/ChIPseq/FindER2/$mark.*.bed --project $mark -o $dirOut/ChIPseq/FindER2/
done
cat $dirOut/ChIPseq/FindER2/H3K27ac.MGG*.FindER2.merge.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin > $dirOut/ChIPseq/FindER2/H3K27ac.MGG.union.bed

## 5mC
mkdir -p $dirOut/WGBS/
for file in $dir5mC/MGG*.combine.5mC.CpG $dirNPC/WGBS/NPC*.combine.5mC.CpG $dirNPC/WGBS/IDHwt_CEMT_23.combine.5mC.CpG; do
    ln -s $file $dirOut/WGBS/
done
echo -e "sample\tcoverage\tN" > $dirOut/WGBS/qc_5mC_coverage.txt
echo -e "sample\ttype\tfractional\tN" > $dirOut/WGBS/qc_5mC_profile.txt 
echo -e "sample\ttype\tmin\tymin\tlower\tmedian\tupper\tymax\tmax" > $dirOut/WGBS/qc_5mC_quantile.txt #ymin: 10% quantile; ymax: 90% quantile
cd $dirOut/WGBS/
for file in *.combine.5mC.CpG; do
    lib=$(echo $file | sed -e 's/.combine.5mC.CpG//g')
    echo "Processing" $lib
    less $file | awk '{c = $4 + $5; if(c >= 5000){s[5001]++} else {s[c]++}} END{for(i = 3; i <= 5001; i++){print "'$lib'""\t"i"\t"s[i]}}' >> $dirOut/WGBS/qc_5mC_coverage.txt
    less $file | awk '{s[int($6*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tgenome\t"i/100"\t"s[i]}}' >> $dirOut/WGBS/qc_5mC_profile.txt 
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | awk '{s[int($1*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tCGI\t"i/100"\t"s[i]}}' >> $dirOut/WGBS/qc_5mC_profile.txt 
    less $file | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirOut/WGBS/qc_5mC_quantile.txt
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirOut/WGBS/qc_5mC_quantile.txt
done
### promoter enhancer 5mC
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
gene=/home/lli/hg19/hg19v69_genes.bed
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
mkdir -p $dirOut/WGBS/enhancer/
for file in $dirOut/ChIPseq/FindER/H3K27ac/*.bed.gz; do
    sample=$(basename $file | cut -d'.' -f2)
    echo $sample
    $BEDTOOLS/intersectBed -a $file -b $promoter -u -f 0.5 | awk '{print $0"\t"$1":"$2"-"$3"_promoter"}' > $dirOut/WGBS/enhancer/$sample.promoter.enhancer
    $java -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $dirOut/ChIPseq/wig/H3K27ac/H3K27ac.$sample.wig.gz -r $dirOut/WGBS/enhancer/$sample.promoter.enhancer -o $dirOut/WGBS/enhancer/ -c $chr -n $sample > $dirOut/WGBS/enhancer/$sample.coverage.log
    less $dirOut/WGBS/enhancer/$sample.promoter.enhancer.$sample.coverage | sed 's/chr//g' | $BEDTOOLS/intersectBed -a stdin -b $dirOut/WGBS/$sample.combine.5mC.CpG -wa -wb | awk '{if($11+$12 >= 3){t[$5]=t[$5]+$11; c[$5]=c[$5]+$12; chr[$5]=$1; start[$5]=$2; end[$5]=$3; signal[$5]=$6}} END{for(i in chr){if(t[i]+c[i]>0){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"signal[i]"\t"c[i]/(c[i]+t[i])"\t""'$sample'"}}}' | sort -k1,1 -k2,2n > $dirOut/WGBS/enhancer/$sample.promoter.enhancer.5mC
done
cat $dirOut/WGBS/enhancer/*.promoter.enhancer.5mC > $dirOut/WGBS/enhancer/promoter.enhancer.5mC
$BEDTOOLS/intersectBed -a $dirOut/WGBS/enhancer/promoter.enhancer.5mC -b /home/lli/hg19/CG.BED -c | awk '{print $0"\t"$8/($3-$2)*1000}' > $dirOut/WGBS/enhancer/promoter.enhancer.5mC.CpG
### DMR against NPC
mkdir -p $dirOut/WGBS/DMR/
echo -e "sample\tp-value\tdelta\tm\ttotal\thyper\thypo" > $dirOut/WGBS/DMR/DM.summary.stats
echo -e "sample\tsize\tcut\tmedian_length\tmedian_N_CpG\ttotal\thyper\thypo" > $dirOut/WGBS/DMR/DMR.summary.stats
pth=0.0005; delta=0.6; m=0.75; cov=3; size=500; cut=3; lib1=MGG_control
for lib2 in MGG_vitc NPC_Cortex02 NPC_Cortex04 NPC_GE02 NPC_GE04 IDHwt_CEMT_23; do
    name=$lib1'_'$lib2; echo $name
    /home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirOut/WGBS/ -o $dirOut/WGBS/DMR/ -f1 $lib1.combine.5mC.CpG -f2 $lib2.combine.5mC.CpG -n $name -p $pth -d $delta -m $m -c $cov
    /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut/WGBS/DMR/ -o $dirOut/WGBS/DMR/ -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
done
$BEDTOOLS/intersectBed -a $dirOut/WGBS/DMR/DMR.MGG_control_NPC_Cortex02.s500.c3.hyper.bed -b $dirOut/WGBS/DMR/DMR.MGG_control_NPC_Cortex04.s500.c3.hyper.bed | $BEDTOOLS/intersectBed -a stdin -b $dirOut/WGBS/DMR/DMR.MGG_control_NPC_GE02.s500.c3.hyper.bed | $BEDTOOLS/intersectBed -a stdin -b $dirOut/WGBS/DMR/DMR.MGG_control_NPC_GE04.s500.c3.hyper.bed > $dirOut/WGBS/DMR/DMR.MGG_control_NPC.s500.c3.hyper.bed
less $dirOut/WGBS/DMR/DMR.MGG_control_NPC.s500.c3.hyper.bed | sed 's/chr//g' > $dirOut/WGBS/DMR/DMR.MGG_control_NPC.hyper.bed
less $dirOut/WGBS/DMR/DMR.MGG_control_NPC.s500.c3.hypo.bed | sed 's/chr//g' > $dirOut/WGBS/DMR/DMR.MGG_control_NPC.hypo.bed
$BEDTOOLS/intersectBed -a $dirOut/WGBS/DMR/DMR.MGG_control_NPC_Cortex02.s500.c3.hypo.bed -b $dirOut/WGBS/DMR/DMR.MGG_control_NPC_Cortex04.s500.c3.hypo.bed | $BEDTOOLS/intersectBed -a stdin -b $dirOut/WGBS/DMR/DMR.MGG_control_NPC_GE02.s500.c3.hypo.bed | $BEDTOOLS/intersectBed -a stdin -b $dirOut/WGBS/DMR/DMR.MGG_control_NPC_GE04.s500.c3.hypo.bed > $dirOut/WGBS/DMR/DMR.MGG_control_NPC.s500.c3.hypo.bed
enhancer=$dirOut/ChIPseq/FindER2/H3K27ac.MGG.union.bed
/home/lli/HirstLab/Pipeline/shell/region.intersect.sh -d $dirOut/WGBS/DMR/ -r $enhancer -n "enhancer"
for file in $dirNPC/WGBS/DMR/DMR.IDHmut_CEMT_*_NPC.hyper.bed $dirNPC/WGBS/DMR/DMR.IDHmut_CEMT_*_NPC.hypo.bed; do
    ln -s $file $dirOut/WGBS/DMR/
done
rm $dirOut/WGBS/DMR/*CEMT_21*
intervene upset -i $dirOut/WGBS/DMR/DMR.MGG_control_NPC.s500.c3.hyper.bed $dirOut/WGBS/DMR/DMR.IDHmut_CEMT_*_NPC.hyper.bed --project hyper -o $dirOut/WGBS/DMR/
intervene upset -i $dirOut/WGBS/DMR/DMR.MGG_control_NPC.s500.c3.hypo.bed $dirOut/WGBS/DMR/DMR.IDHmut_CEMT_*_NPC.hypo.bed --project hypo -o $dirOut/WGBS/DMR/
### DMR against NB141
mkdir -p $dirOut/WGBS/NB141/
ln -s /projects/epigenomics3/epigenomics3_results/users/lli/glioma/WGBS/Normal.NB141.combine.5mC.CpG $dirOut/WGBS/Normal.NB141.combine.5mC.CpG
echo -e "sample\tp-value\tdelta\tm\ttotal\thyper\thypo" > $dirOut/WGBS/NB141/DM.summary.stats
echo -e "sample\tsize\tcut\tmedian_length\tmedian_N_CpG\ttotal\thyper\thypo" > $dirOut/WGBS/NB141/DMR.summary.stats
pth=0.0005; delta=0.6; m=0.75; cov=3; size=500; cut=3
lib1=MGG_control; lib2=Normal.NB141; name=$lib1'_'$lib2; echo $name
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirOut/WGBS/ -o $dirOut/WGBS/NB141/ -f1 $lib1.combine.5mC.CpG -f2 $lib2.combine.5mC.CpG -n $name -p $pth -d $delta -m $m -c $cov
/home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut/WGBS/NB141/ -o $dirOut/WGBS/NB141/ -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
enhancer=$dirOut/ChIPseq/FindER2/H3K27ac.MGG.union.bed
/home/lli/HirstLab/Pipeline/shell/region.intersect.sh -d $dirOut/WGBS/NB141/ -r $enhancer -n "enhancer"

## hMeDIP
mkdir -p $dirOut/hMeDIP/bam/
for file in $dir5hmC/bam/MGG*realign*.bam; do
    ln -s $file $dirOut/hMeDIP/bam/
    ln -s $file.bai $dirOut/hMeDIP/bam/
done
mkdir -p $dirOut/hMeDIP/bw/
for file in $dir5hmC/bw/MGG*.bw; do
    ln -s $file $dirOut/hMeDIP/bw/
done
### FindER2
finder2=/home/mbilenky/bin/FindER2/finder2.jar
mkdir -p $dirOut/hMeDIP/FindER2/
echo -e "Sample\tN_region\tTotal_length\tAverage_length" > $dirOut/ChIPseq/FindER2/ER_summary.txt
for bam in $dirOut/hMeDIP/bam/*.bam; do
    sample=$(basename $bam | cut -d'.' -f1);
    echo $sample
    $java -jar -Xmx25G $finder2 inputBam:$dirOut/ChIPseq/bam/input/input.$sample.bam signalBam:$bam outDir:$dirOut/hMeDIP/FindER2/ acgtDir:/projects/epigenomics2/users/mbilenky/resources/hg19/ACGT 
    echo -e $sample"\t"$(less $dirOut/hMeDIP/FindER2/$sample.realign.FindER2.bed | wc -l)"\t"$(less $dirOut/hMeDIP/FindER2/$sample.realign.FindER2.bed | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"$3/$2}' >> $dirOut/hMeDIP/FindER2/ER_summary.txt
done
### unique: FC >= 2 & RPKM >= 5
CG=/home/lli/hg19/CG.BED
cat $dirOut/hMeDIP/FindER2/*.realign.FindER2.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin > $dirOut/hMeDIP/FindER2/ER.union.bed
multiBigwigSummary BED-file -b $dirOut/hMeDIP/bw/MGG_control.realign.bw $dirOut/hMeDIP/bw/MGG_vitc.realign.bw --BED $dirOut/hMeDIP/FindER2/ER.union.bed --labels MGG_control MGG_vitc -out $dirOut/hMeDIP/FindER2/ER.union.matrix.npz --outRawCounts $dirOut/hMeDIP/FindER2/ER.union.matrix.RPKM
less $dirOut/hMeDIP/FindER2/ER.union.matrix.RPKM | awk 'NR>1{fc=($4+0.0001)/($5+0.0001); if($4>=5&&fc>=2){print $0"\t"fc}}' | $BEDTOOLS/intersectBed -a stdin -b $CG -c > $dirOut/hMeDIP/FindER2/MGG_control.unique.bed
less $dirOut/hMeDIP/FindER2/ER.union.matrix.RPKM | awk 'NR>1{fc=($5+0.0001)/($4+0.0001); if($5>=5&&fc>=2){print $0"\t"fc}}' | $BEDTOOLS/intersectBed -a stdin -b $CG -c > $dirOut/hMeDIP/FindER2/MGG_vitc.unique.bed
for file in $dirOut/hMeDIP/FindER2/*unique.bed; do
    less $file | awk '{print "chr"$0}' > $file.bed
    less $file | awk '{if($7>=5){print "chr"$0}}' > $file.c5.bed
done
enhancer=$dirOut/ChIPseq/FindER2/H3K27ac.MGG.union.bed
/home/lli/HirstLab/Pipeline/shell/region.intersect.sh -d $dirOut/hMeDIP/FindER2/ -r $enhancer -n "enhancer"
intervene upset -i $dirOut/WGBS/DMR/DMR.MGG_control_NPC.hyp*.bed $dirOut/hMeDIP/FindER2/*unique.bed --project DMR.DhMR -o $dirOut/hMeDIP/FindER2/
intervene upset -i $dirOut/WGBS/DMR/DMR.IDHmut_CEMT_*hyper.bed $dirOut/hMeDIP/FindER2/*unique.bed.bed --project CEMT_hyper.DhMR -o $dirOut/hMeDIP/FindER2/
intervene upset -i $dirOut/WGBS/DMR/DMR.IDHmut_CEMT_*hypo.bed $dirOut/hMeDIP/FindER2/*unique.bed.bed --project CEMT_hypo.DhMR -o $dirOut/hMeDIP/FindER2/
### Homer
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
for file in $dirOut/hMeDIP/FindER2/*unique.bed; do
     name=$(basename $file | sed 's/.bed//g')
     echo $name
     mkdir -p $dirOut/hMeDIP/FindER2/homer/$name/
     /home/lli/bin/homer/bin/findMotifsGenome.pl <(less $file | awk '{print "chr"$0}') hg19 $dirOut/hMeDIP/FindER2/homer/$name/ -size 200 -len 8 
done
mkdir -p $dirOut/hMeDIP/FindER2/homer/MGG_control.MGG_vitc/
/home/lli/bin/homer/bin/findMotifsGenome.pl <(less $dirOut/hMeDIP/FindER2/MGG_control.unique.bed | awk '{print "chr"$0}') hg19 $dirOut/hMeDIP/FindER2/homer/MGG_control.MGG_vitc/ -size 200 -len 8 -bg <(less $dirOut/hMeDIP/FindER2/MGG_vitc.unique.bed | awk '{print "chr"$0}')
mkdir -p $dirOut/hMeDIP/FindER2/homer/MGG_vitc.MGG_control/
/home/lli/bin/homer/bin/findMotifsGenome.pl <(less $dirOut/hMeDIP/FindER2/MGG_vitc.unique.bed | awk '{print "chr"$0}') hg19 $dirOut/hMeDIP/FindER2/homer/MGG_vitc.MGG_control/ -size 200 -len 8 -bg <(less $dirOut/hMeDIP/FindER2/MGG_control.unique.bed | awk '{print "chr"$0}')
### intersect with enhancer
mkdir -p $dirOut/hMeDIP/FindER2/enhancer/
echo -e "Name\tN_total\tlength_total\tN_enhancer\tlength_enhancer\tpercent" > $dirOut/hMeDIP/FindER2/enhancer/ER_enhancer_summary.txt
for file in $dirOut/hMeDIP/FindER2/*unique.bed; do
    name=$(basename $file | sed 's/.bed//g')
    sample=$(basename $file | cut -d'.' -f1)
    echo $name $sample
    $BEDTOOLS/intersectBed -a $file -b $enhancer -u > $dirOut/hMeDIP/FindER2/enhancer/$name.enhancer.bed
    echo -e $name"\t"$(less $file | wc -l)"\t"$(less $file | awk '{s=s+$3-$2}END{print s}')"\t"$(less $dirOut/hMeDIP/FindER2/enhancer/$name.enhancer.bed | wc -l)"\t"$(less $dirOut/hMeDIP/FindER2/enhancer/$name.enhancer.bed | awk '{s=s+$3-$2}END{print s}') | awk '{print $0"\t"$5/$3}' >> $dirOut/hMeDIP/FindER2/enhancer/ER_enhancer_summary.txt
done
intervene upset -i $dirOut/WGBS/DMR/DMR.MGG_control_NPC.hyper.bed $dirOut/hMeDIP/FindER2/enhancer/*.bed --project DMR.DhMR_enhancer -o $dirOut/hMeDIP/FindER2/
#### Homer
for file in $dirOut/hMeDIP/FindER2/enhancer/*.bed; do
     name=$(basename $file | sed 's/.enhancer.bed//g')
     echo $name
     mkdir -p $dirOut/hMeDIP/FindER2/enhancer/$name/
     /home/lli/bin/homer/bin/findMotifsGenome.pl <(less $file | awk '{print "chr"$0}') hg19 $dirOut/hMeDIP/FindER2/enhancer/$name/ -size 200 -len 8 
done
mkdir -p $dirOut/hMeDIP/FindER2/enhancer/homer/MGG_control.MGG_vitc/
/home/lli/bin/homer/bin/findMotifsGenome.pl <(less $dirOut/hMeDIP/FindER2/enhancer/MGG_control.unique.enhancer.bed | awk '{print "chr"$0}') hg19 $dirOut/hMeDIP/FindER2/enhancer/homer/MGG_control.MGG_vitc/ -size 200 -len 8 -bg <(less $dirOut/hMeDIP/FindER2/enhancer/MGG_vitc.unique.enhancer.bed | awk '{print "chr"$0}')
mkdir -p $dirOut/hMeDIP/FindER2/enhancer/homer/MGG_vitc.MGG_control/
/home/lli/bin/homer/bin/findMotifsGenome.pl <(less $dirOut/hMeDIP/FindER2/enhancer/MGG_vitc.unique.enhancer.bed | awk '{print "chr"$0}') hg19 $dirOut/hMeDIP/FindER2/enhancer/homer/MGG_vitc.MGG_control/ -size 200 -len 8 -bg <(less $dirOut/hMeDIP/FindER2/enhancer/MGG_control.unique.enhancer.bed | awk '{print "chr"$0}')

## RNAseq
mkdir -p $dirOut/RNAseq/DEfine/
for file in $dirRNA/RPKM/MGG*/coverage/*.G.A.rpkm.pc $dirNPC/RNAseq/NPC_RPKM/*GE04/coverage/*.G.A.rpkm.pc; do
    name=$(basename $file | sed 's/119//' | sed 's/A15299.GE04/NPC_GE04/')
    ln -s $file $dirOut/RNAseq/$name
done
echo -e "%/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab
addpath /home/mbilenky/matlab -end" > $dirOut/RNAseq/DEfine/DEfine.MGG.m
sample1=MGG_control
for sample2 in MGG_vitc NPC_GE04; do
    name1=$sample1; name2=$sample2
    echo $sample1 $sample2
    echo -e "
%%%%%%%%%%%%%%%%%%%%%%%%%% $name1 vs $name2 %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
out=true; corr1=true; corr2=true; rpkm=true; figs=true; RPKMmin=0.005; Nmin=25; eps=0.0001; maxLim=3.5; fdr=0.05; 
[idl,gl]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.length','%s %f');
[id,gc]=textread('/projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.pc.EnsID.GC','%s %f');
dirIn='$dirOut/RNAseq/';
dirOut='$dirOut/RNAseq/DEfine/';
sample1='$sample1'; name1='$name1';
sample2='$sample2'; name2='$name2';
[idl,n1,r1,rmi,ra,rma]=textread(strcat(dirIn, sample1, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[id,n2,r2,rmi,ra,rma]=textread(strcat(dirIn, sample2, '.G.A.rpkm.pc'),'%s %f %f %f %f %f');
[C,ix,ixl]=intersect(id,idl);
[cc,nfup,nfdn]=DEfine(idl(ixl), r1(ixl), r2(ix), n1(ixl), n2(ix), [gl(ixl), gc(ix)], dirOut, name1, name2, out, figs, fdr, corr1, corr2, rpkm, RPKMmin, Nmin, eps, maxLim);
" >> $dirOut/RNAseq/DEfine/DEfine.MGG.m
done
### trackHub
chr=/projects/epigenomics2/resources/UCSC_chr/hg19.bwa2ucsc.names
chrsize=/home/lli/hg19/hg19.chrom.sizes
dirWig=$dirOut/RNAseq/wig/
dirHub=/gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/RNAseqHub/hg19/
mkdir -p $dirWig; mkdir -p $dirHub
cp /gsc/www/bcgsc.ca/downloads/mb/BrainHubs/HistoneHub/genomes.txt /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/RNAseqHub/
echo -e "hub VitC_gliomaHub_RNAseq
shortLabel VitC_glioma Hub (RNAseq)
longLabel Hub to display VitC glioma data at UCSC (RNAseq)
genomesFile genomes.txt
email lli@bcgsc.ca" > /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/RNAseqHub/hub.txt
> $dirHub/trackDb.txt
for bam in /projects/epigenomics3/epigenomics3_results/users/lli/NHA/RNAseq/bam/MGG*.bam; do
    name=$(basename $bam | cut -d'.' -f1); echo $name
    /home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirWig -F:1028,-q:5,-n:$name,-chr:$chr,-cp
    /home/lli/HirstLab/Pipeline/UCSC/wigToBigWig $dirWig/$name.q5.F1028.PET.wig.gz $chrsize $dirHub/$name.q5.F1028.PET.bw
    if [[ "$name" =~ "control" ]]; then
        color="255,0,0"
    else
        color="0,0,255"
    fi
    echo -e "
track $name
shortLabel $name
longLabel RNAseq $name
type bigWig
visibility full
maxHeightPixels 70:70:32
configurable on
autoScale on
alwaysZero on
priority 0.1
bigDataUrl $name.q5.F1028.PET.bw
color $color
" >> $dirHub/trackDb.txt
done
