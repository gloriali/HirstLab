#!/bin/sh

dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/bam/
mkdir -p $dirOut
cd /projects/glioma_dart/glimoa_ctcf/
for bam in *.bam; do
	ln -s /projects/glioma_dart/glimoa_ctcf/$bam $dirOut/$bam
done
cat mgh59_idh1m_ctcf_10/subfile1.bin mgh59_idh1m_ctcf_10/subfile2.bin > $dirOut/mgh59_idh1m_ctcf_10.bam
cat mgh76_ctcf_19/subfile1.bin mgh76_ctcf_19/subfile2.bin > $dirOut/mgh76_ctcf_19.bam
cat mgh80_ctcf_21/subfile1.bin mgh80_ctcf_21/subfile2.bin > $dirOut/mgh80_ctcf_21.bam
cat mgh81_idh1m_ctcf_22/subfile1\ \(1\).bin mgh81_idh1m_ctcf_22/subfile2.bin > $dirOut/mgh81_idh1m_ctcf_22.bam

## QC
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
chr=/projects/epigenomics/resources/UCSC_chr/hg19.bwa2ucsc.names
chrsize=/home/lli/hg19/hg19.chrom.sizes
dirIn=/projects/epigenomics2/users/lli/glioma/CTCF/bam/
dirWig=/projects/epigenomics2/users/lli/glioma/CTCF/wig/
dirBW=/gsc/www/bcgsc.ca/downloads/mb/CTCF_glioma/hg19/
mkdir -p $dirWig
mkdir -p $dirBW
cp /gsc/www/bcgsc.ca/downloads/mb/BrainHubs/HistoneHub/genomes.txt /gsc/www/bcgsc.ca/downloads/mb/CTCF_glioma/
echo -e "hub CTCF_glioma
shortLabel CTCF_glioma 
longLabel Hub to display CTCF glioma data at UCSC
genomesFile genomes.txt
email lli@bcgsc.ca" > /gsc/www/bcgsc.ca/downloads/mb/CTCF_glioma/hub.txt
> $dirBW/trackDb.txt
for bam in $dirIn/*.bam; do
    name=$(basename $bam | sed 's/.bam//g')
    sample=$(echo $name | cut -d'_' -f1)
    idh=$(less $dirIn/../sample.info | awk '{if($2=="'$sample'")print $1}' | awk '{if($1~/GBM/){print "IDHwt"}else{print "IDHmut"}}')
    if [ $idh == "IDHmut" ]; then
        color=255,0,0
    else
        color=0,0,255
    fi
    echo $name $sample $idh $color
    $samtools index $bam
    $bamstats -g 2864785220 -q 10 -b $bam  > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn/ $name $dirIn/$name.bamstats
    /home/mbilenky/bin/PETLengthDist.sh $bam 5 $dirIn 10
    /home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirWig -F:1028,-q:5,-n:$name,-chr:$chr,-cp
    /home/lli/HirstLab/Pipeline/UCSC/wigToBigWig $dirWig/$name.q5.F1028.PET.wig.gz $chrsize $dirBW/$name.q5.F1028.PET.bw
    echo -e "
track $name
shortLabel $name
longLabel $idh  $name
type bigWig
visibility full
maxHeightPixels 70:70:32
configurable on
autoScale on
alwaysZero on
priority 0.1
bigDataUrl $name.q5.F1028.PET.bw
color $color
" >> $dirBW/trackDb.txt
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn

## FindER
JAVA=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/PET_200/map_old/
dirIn=/projects/epigenomics2/users/lli/glioma/CTCF/bam/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/
mkdir -p $dirOut
cd $dirIn
for bam in *ctcf*.bam; do
	echo $bam
	$JAVA -jar -Xmx10G /home/mbilenky/bin/Solexa_Java/FindER.0.9.3b.jar -i $bam -r $reg -o $dirOut -v -m $map -info -TF > $dirOut/FindER.$bam.log
done
for bam in *h3k27ac*.bam; do
	echo $bam
	$JAVA -jar -Xmx10G /home/mbilenky/bin/Solexa_Java/FindER.0.9.3b.jar -i $bam -r $reg -o $dirOut -v -m $map -info -minER 300 > $dirOut/FindER.$bam.log
done
### CTCF changes
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/
cd $dirOut
for f1 in FindER_scan.mgh7530_ctcf_13.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh6971_ctcf_11.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh7770_ctcf_15.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh80_ctcf_21.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh76_ctcf_19.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh79_ctcf_20.pctl_0.1.FDR_0.05.bed.gz; do
    s1=$(echo $f1 | sed 's/FindER_scan.//g' | cut -d'_' -f1)
    > $dirOut/CTCF_loss_all.$s1.bed
    i=0
    for f2 in FindER_scan.mgh7478_idhm_ctcf_12.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh59_idh1m_ctcf_10.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh17m_ctcf_17.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh18m_ctcf_18.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh81_idh1m_ctcf_22.pctl_0.1.FDR_0.05.bed.gz; do
        s2=$(echo $f2 | sed 's/FindER_scan.//g' | cut -d'_' -f1)
        echo $s1 $s2
        i=$((i+1))
        $BEDTOOLS/intersectBed -a $f1 -b $f2 -v | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t""'$s2'"}' >> $dirOut/CTCF_loss_all.$s1.bed
    done
    less $dirOut/CTCF_loss_all.$s1.bed | awk '{print $4}' | sort | uniq -c | awk '{if($1=="'$i'"){gsub(":", "\t"); gsub("-", "\t"); print $2"\t"$3"\t"$4"\t"$2":"$3"-"$4}}' > $dirOut/CTCF_loss.$s1.bed 
done
> $dirOut/CTCF_loss.intersect.bed
for file in CTCF_loss.*.bed; do
    s=$(echo $file | cut -d'.' -f2)
    $BEDTOOLS/intersectBed -a CTCF_loss.mgh6971.bed -b $file -u | awk '{print $0"\t""'$s'"}' >> $dirOut/CTCF_loss.intersect.bed
done
less $dirOut/CTCF_loss.intersect.bed | awk '{print $4}' | sort | uniq -c | awk '{if($1>=4){gsub(":", "\t"); gsub("-", "\t"); print $2"\t"$3"\t"$4"\t"$2":"$3"-"$4}}' > $dirOut/CTCF_IDHwt_unique.bed
for f1 in FindER_scan.mgh7478_idhm_ctcf_12.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh59_idh1m_ctcf_10.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh17m_ctcf_17.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh18m_ctcf_18.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh81_idh1m_ctcf_22.pctl_0.1.FDR_0.05.bed.gz; do
    s1=$(echo $f1 | sed 's/FindER_scan.//g' | cut -d'_' -f1)
    > $dirOut/CTCF_gain_all.$s1.bed
    > $dirOut/CTCF_retained_all.$s1.bed
    i=0
    for f2 in FindER_scan.mgh7530_ctcf_13.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh6971_ctcf_11.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh7770_ctcf_15.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh80_ctcf_21.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh76_ctcf_19.pctl_0.1.FDR_0.05.bed.gz FindER_scan.mgh79_ctcf_20.pctl_0.1.FDR_0.05.bed.gz; do
        s2=$(echo $f2 | sed 's/FindER_scan.//g' | cut -d'_' -f1)
        echo $s1 $s2
        i=$((i+1))
        $BEDTOOLS/intersectBed -a $f1 -b $f2 -v | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t""'$s2'"}' >> $dirOut/CTCF_gain_all.$s1.bed
        $BEDTOOLS/intersectBed -a $f1 -b $f2 -u | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t""'$s2'"}' >> $dirOut/CTCF_retained_all.$s1.bed
    done
    less $dirOut/CTCF_gain_all.$s1.bed | awk '{print $4}' | sort | uniq -c | awk '{if($1=="'$i'"){gsub(":", "\t"); gsub("-", "\t"); print $2"\t"$3"\t"$4"\t"$2":"$3"-"$4}}' > $dirOut/CTCF_gain.$s1.bed 
    less $dirOut/CTCF_retained_all.$s1.bed | awk '{print $4}' | sort | uniq -c | awk '{if($1=="'$i'"){gsub(":", "\t"); gsub("-", "\t"); print $2"\t"$3"\t"$4"\t"$2":"$3"-"$4}}' > $dirOut/CTCF_retained.$s1.bed 
done
> $dirOut/CTCF_gain.intersect.bed
for file in CTCF_gain.*.bed; do
    s=$(echo $file | cut -d'.' -f2)
    $BEDTOOLS/intersectBed -a CTCF_gain.mgh17m.bed -b $file -u | awk '{print $0"\t""'$s'"}' >> $dirOut/CTCF_gain.intersect.bed
done
less $dirOut/CTCF_gain.intersect.bed | awk '{print $4}' | sort | uniq -c | awk '{if($1>=4){gsub(":", "\t"); gsub("-", "\t"); print $2"\t"$3"\t"$4"\t"$2":"$3"-"$4}}' > $dirOut/CTCF_IDHmut_unique.bed
> $dirOut/CTCF_retained.intersect.bed
for file in CTCF_retained.*.bed; do
    s=$(echo $file | cut -d'.' -f2)
    $BEDTOOLS/intersectBed -a CTCF_retained.mgh17m.bed -b $file -u | awk '{print $0"\t""'$s'"}' >> $dirOut/CTCF_retained.intersect.bed
done
less $dirOut/CTCF_retained.intersect.bed | awk '{print $4}' | sort | uniq -c | awk '{if($1>=4){gsub(":", "\t"); gsub("-", "\t"); print $2"\t"$3"\t"$4"\t"$2":"$3"-"$4}}' > $dirOut/CTCF_IDHmut_retained.bed

## 5mC
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/WGBS/
mkdir -p $dirOut
cd $dir5mC
for CTCF in /projects/epigenomics2/users/lli/glioma/CTCF/FindER/CTCF_IDH*.bed; do
    echo -e "chr\tstart\tend\tID\tfractional\tsample" > $dirOut/$(basename $CTCF).5mC
    for file in *.combine.5mC.CpG; do
        sample=$(echo $file | sed 's/.5mC.CpG.combine.5mC.CpG//g');
        echo $(basename $CTCF) $sample;
        $BEDTOOLS/intersectBed -a $file -b <(less $CTCF | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}') -wa -wb | awk '{t[$10]=t[$10]+$4; c[$10]=c[$10]+$5; chr[$10]=$7; start[$10]=$8; end[$10]=$9} END{for(i in chr){f=c[i]/(c[i]+t[i]); print "chr"chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"f"\t""'$sample'"}}' | sort -k1,1 -k 2,2n >> $dirOut/$(basename $CTCF).5mC
    done
done

## H3K36me3
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
responder=/projects/epigenomics2/users/lli/glioma/CTCF/WGBS/responder.bed
nonresponder=/projects/epigenomics2/users/lli/glioma/CTCF/WGBS/nonresponder.bed
gain=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/CTCF_IDHmut_unique.bed
loss=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/CTCF_IDHwt_unique.bed
retained=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/CTCF_IDHmut_retained.bed
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/H3K36me3/
dirK36=/projects/epigenomics2/users/lli/glioma/ChIPseq/bam/H3K36me3/
echo -e "chr\tstart\tend\tID\tmut\twt\tdiff\tN\tsample\tRPKM" > $dirOut/resonpder.H3K36me3.bed
echo -e "chr\tstart\tend\tID\tmut\twt\tdiff\tN\tsample\tRPKM" > $dirOut/nonresonpder.H3K36me3.bed
echo -e "chr\tstart\tend\tID\tN\tsample\tRPKM" > $dirOut/gain.H3K36me3.bed
echo -e "chr\tstart\tend\tID\tN\tsample\tRPKM" > $dirOut/loss.H3K36me3.bed
echo -e "chr\tstart\tend\tID\tsample\tRPKM" > $dirOut/retained.H3K36me3.bed
mkdir -p $dirOut
cd $dirK36
for bam in *.bam; do
    sample=$(echo $bam | cut -d'.' -f2)
    echo $sample
    depth=$($samtools view -q 5 -F 1028 $bam | wc -l)
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a $responder -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$8*10^9/($3-$2)/"'$depth'"}' >> $dirOut/resonpder.H3K36me3.bed
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a $nonresponder -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$8*10^9/($3-$2)/"'$depth'"}' >> $dirOut/nonresonpder.H3K36me3.bed    
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a <(less $gain | sed 's/chr//g') -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$5*10^9/($3-$2)/"'$depth'"}' >> $dirOut/gain.H3K36me3.bed    
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a <(less $loss | sed 's/chr//g') -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$5*10^9/($3-$2)/"'$depth'"}' >> $dirOut/loss.H3K36me3.bed    
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a <(less $retained | sed 's/chr//g') -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$5*10^9/($3-$2)/"'$depth'"}' >> $dirOut/retained.H3K36me3.bed    
done

## H3K27me3
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
responder=/projects/epigenomics2/users/lli/glioma/CTCF/WGBS/responder.bed
nonresponder=/projects/epigenomics2/users/lli/glioma/CTCF/WGBS/nonresponder.bed
gain=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/CTCF_IDHmut_unique.bed
loss=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/CTCF_IDHwt_unique.bed
retained=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/CTCF_IDHmut_retained.bed
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/H3K27me3/
dirK27=/projects/epigenomics2/users/lli/glioma/ChIPseq/bam/H3K27me3/
echo -e "chr\tstart\tend\tID\tmut\twt\tdiff\tN\tsample\tRPKM" > $dirOut/resonpder.H3K27me3.bed
echo -e "chr\tstart\tend\tID\tmut\twt\tdiff\tN\tsample\tRPKM" > $dirOut/nonresonpder.H3K27me3.bed
echo -e "chr\tstart\tend\tID\tN\tsample\tRPKM" > $dirOut/gain.H3K27me3.bed
echo -e "chr\tstart\tend\tID\tN\tsample\tRPKM" > $dirOut/loss.H3K27me3.bed
echo -e "chr\tstart\tend\tID\tsample\tRPKM" > $dirOut/retained.H3K27me3.bed
mkdir -p $dirOut
cd $dirK27
for bam in *.bam; do
    sample=$(echo $bam | cut -d'.' -f2)
    echo $sample
    depth=$($samtools view -q 5 -F 1028 $bam | wc -l)
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a $responder -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$8*10^9/($3-$2)/"'$depth'"}' >> $dirOut/resonpder.H3K27me3.bed
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a $nonresponder -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$8*10^9/($3-$2)/"'$depth'"}' >> $dirOut/nonresonpder.H3K27me3.bed    
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a <(less $gain | sed 's/chr//g') -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$5*10^9/($3-$2)/"'$depth'"}' >> $dirOut/gain.H3K27me3.bed    
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a <(less $loss | sed 's/chr//g') -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$5*10^9/($3-$2)/"'$depth'"}' >> $dirOut/loss.H3K27me3.bed    
    $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a <(less $retained | sed 's/chr//g') -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$5*10^9/($3-$2)/"'$depth'"}' >> $dirOut/retained.H3K27me3.bed    
done

## H3K27ac
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
gain=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/CTCF_IDHmut_unique.bed
loss=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/CTCF_IDHwt_unique.bed
RPKM=/projects/epigenomics2/users/lli/glioma/RNAseq/RPKM/glioma.RPKM
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/H3K27ac/
dirK27=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/
mkdir -p $dirOut
cd $dirK27
less /projects/epigenomics/resources/Ensembl/hg19v69/hg19v69_genes.TSS_2000.pc.bed | awk '{print $1"\t"$2+2000"\t"$2+2001"\t"$4}' > /home/lli/hg19/hg19v69_genes.TSS.pc.bed
for file in $gain $loss; do
    $BEDTOOLS/closestBed -a <(less /home/lli/hg19/hg19v69_genes.TSS.pc.bed | sort -k1,1 -k2,2n) -b <(less $file | sort -k1,1 -k2,2n) -iu -D b | awk '{if($9<=1000000&&$6!=-1){print $1"\t"$2"\t"$3"\t"$4"\t"$8}}' > $dirOut/TSS_iu.$(basename $file)
    $BEDTOOLS/closestBed -a <(less /home/lli/hg19/hg19v69_genes.TSS.pc.bed | sort -k1,1 -k2,2n) -b <(less $file | sort -k1,1 -k2,2n) -id -D b | awk '{if($9>=-1000000&&$6!=-1){print $1"\t"$2"\t"$3"\t"$4"\t"$8}}' > $dirOut/TSS_id.$(basename $file)
    echo -e "ENSG\tCTCF\tenhancer\tsample\td" > $dirOut/H3K27ac.TSS.$(basename $file)
    for enhancer in CEMT*.bed.gz; do
        sample=$(echo $enhancer | cut -d'.' -f1)
        echo $sample
        $BEDTOOLS/closestBed -a <(less $enhancer | sort -k1,1 -k2,2n | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}') -b <(less $file | sort -k1,1 -k2,2n) -iu -D b | awk '{if($9<=1000000&&$6!=-1){print $1"\t"$2"\t"$3"\t"$4"\t"$8}}' > $dirOut/H3K27ac_iu.$sample.$(basename $file)
        $BEDTOOLS/closestBed -a <(less $enhancer | sort -k1,1 -k2,2n | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}') -b <(less $file | sort -k1,1 -k2,2n) -id -D b | awk '{if($9>=-1000000&&$6!=-1){print $1"\t"$2"\t"$3"\t"$4"\t"$8}}' > $dirOut/H3K27ac_id.$sample.$(basename $file)
        join <(less $dirOut/H3K27ac_iu.$sample.$(basename $file) | sort -k5,5) <(less $dirOut/TSS_id.$(basename $file) | sort -k5,5) -1 5 -2 5 | awk '{d=sqrt(($7-($3+$4)/2)^2); if(d<=1000000){print $9"\t"$1"\t"$5"\t""'$sample'""\t"d}}' >> $dirOut/H3K27ac.TSS.$(basename $file)
        join <(less $dirOut/H3K27ac_id.$sample.$(basename $file) | sort -k5,5) <(less $dirOut/TSS_iu.$(basename $file) | sort -k5,5) -1 5 -2 5 | awk '{d=sqrt(($7-($3+$4)/2)^2); if(d<=1000000){print $9"\t"$1"\t"$5"\t""'$sample'""\t"d}}' >> $dirOut/H3K27ac.TSS.$(basename $file)
    done
    join <(less $dirOut/H3K27ac.TSS.$(basename $file) | sort -k1,1) <(less $RPKM | sort -k1,1) | sed 's/ /\t/g' > $dirOut/H3K27ac.TSS.$(basename $file).RPKM
done

### Homer for CTCF BS
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics2/users/lli/glioma/CTCF/FindER/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/Homer/
mkdir -p $dirOut
cd $dirIn
gunzip FindER_scan.*ctcf*.bed.gz 
for file in CTCF_IDH*.bed FindER_scan.*ctcf*.bed CTCF_gain.m*.bed  CTCF_loss.m*.bed  CTCF_retained.m*.bed; do
	name=$(echo $file | sed 's/.bed//g' | sed 's/_ctcf.*//g' | sed 's/_idh.*//g')
	echo "Processing "$name
	mkdir -p $dirOut/$name/
	/home/lli/bin/homer/bin/findMotifsGenome.pl $file hg19 $dirOut/$name/ -size 200 -len 6,10,15,20 
done
/home/lli/bin/homer/bin/annotatePeaks.pl $dirIn/CTCF_IDHmut_retained.bed hg19 -m $dirOut/CTCF_IDHmut_retained/knownResults/known1.motif -nmotifs -mbed $dirOut/CTCF_IDHmut_retained_long.BS.bed > $dirOut/CTCF_IDHmut_retained_long.annotate
/home/lli/bin/homer/bin/annotatePeaks.pl $dirIn/CTCF_IDHmut_retained.bed hg19 -m $dirOut/CTCF_IDHmut_retained/homerResults/motif1.motif -nmotifs -mbed $dirOut/CTCF_IDHmut_retained_short.BS.bed > $dirOut/CTCF_IDHmut_retained_short.annotate
/home/lli/bin/homer/bin/annotatePeaks.pl $dirIn/CTCF_IDHmut_unique.bed hg19 -m $dirOut/CTCF_IDHmut_unique/homerResults/motif9.motif -nmotifs -mbed $dirOut/CTCF_IDHmut_unique_short.BS.bed > $dirOut/CTCF_IDHmut_unique_short.annotate
/home/lli/bin/homer/bin/annotatePeaks.pl $dirIn/CTCF_IDHwt_unique.bed hg19 -m $dirOut/CTCF_IDHwt_unique/homerResults/motif1.motif -nmotifs -mbed $dirOut/CTCF_IDHwt_unique_short.BS.bed > $dirOut/CTCF_IDHwt_unique_short.annotate
#### 5mC
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/WGBS/
mkdir -p $dirOut
cd $dir5mC
for CTCF in /projects/epigenomics2/users/lli/glioma/CTCF/Homer/*.bed; do
    name=$(basename $CTCF | cut -d'.' -f1)
    echo -e "chr\tstart\tend\tID\tfractional\tsample" > $dirOut/$name.5mC
    for file in *.combine.5mC.CpG; do
        sample=$(echo $file | sed 's/.5mC.CpG.combine.5mC.CpG//g');
        echo $name $sample;
        $BEDTOOLS/intersectBed -a $file -b <(less $CTCF | awk 'NR>1{gsub("chr", ""); print $1"\t"$2-20"\t"$3+20"\t"$1":"$2"-"$3}') -wa -wb | awk '{t[$10]=t[$10]+$4; c[$10]=c[$10]+$5; chr[$10]=$7; start[$10]=$8; end[$10]=$9} END{for(i in chr){if(c[i]+t[i]>0){f=c[i]/(c[i]+t[i]); print "chr"chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"f"\t""'$sample'"}}}' | sort -k1,1 -k 2,2n >> $dirOut/$name.5mC
    done
done
#### H3K36me3
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/H3K36me3/
dirK36=/projects/epigenomics2/users/lli/glioma/ChIPseq/bam/H3K36me3/
mkdir -p $dirOut
cd $dirK36
for CTCF in /projects/epigenomics2/users/lli/glioma/CTCF/Homer/*.bed; do
    name=$(basename $CTCF | cut -d'.' -f1)
    echo -e "chr\tstart\tend\tID\tN\tsample\tRPKM" > $dirOut/$name.H3K36me3.bed
    for bam in *.bam; do
        sample=$(echo $bam | cut -d'.' -f2)
        echo $name $sample
        depth=$($samtools view -q 5 -F 1028 $bam | wc -l)
        $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a <(less $CTCF | awk 'NR>1{gsub("chr", ""); print $1"\t"$2-20"\t"$3+20"\t"$1":"$2"-"$3}') -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$5*10^9/($3-$2)/"'$depth'"}' >> $dirOut/$name.H3K36me3.bed    
    done
done
#### RNA-seq
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/RNAseq/
mkdir -p $dirOut
for CTCF in /projects/epigenomics2/users/lli/glioma/CTCF/Homer/*.bed; do
    name=$(basename $CTCF | cut -d'.' -f1)
    echo -e "chr\tstart\tend\tID\tN\tsample\tRPKM" > $dirOut/$name.RNAseq.bed
    for dirRNA in /projects/epigenomics2/users/lli/glioma/RNAseq/bam/ /projects/epigenomics/users/lli/FetalBrain/RNAseq/bam/; do
        cd $dirRNA
        for bam in *.bam; do
            sample=$(echo $bam | cut -d'.' -f2)
            echo $name $sample
            depth=$($samtools view -q 5 -F 1028 $bam | wc -l)
            $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a <(less $CTCF | awk 'NR>1{gsub("chr", ""); print $1"\t"$2-20"\t"$3+20"\t"$1":"$2"-"$3}') -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$5*10^9/($3-$2)/"'$depth'"}' >> $dirOut/$name.RNAseq.bed    
        done
    done
done

## ===================================================================
## MACS2
export PATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/bin:$PATH
export PYTHONPATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/lib/python2.7/site-packages:$PYTHONPATH
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
dirIn=/projects/epigenomics2/users/lli/glioma/CTCF/bam/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/MACS2/
mkdir -p $dirOut/
echo -e "Sample\tIDH\tN_peak\tTotal_length" > $dirOut/peak_summary.txt
cd $dirIn
for file in $dirIn/*.bam; do
    name=$(basename $file | sed 's/.bam//g')
    sample=$(echo $name | cut -d'_' -f1)
    idh=$(less $dirIn/../sample.info | awk '{if($2=="'$sample'")print $1}' | awk '{if($1~/GBM/){print "IDHwt"}else{print "IDHmut"}}')
    echo $name $sample $idh
    $samtools view -F 1028 -q 5 -b $file > $dirIn/$name.q5.F1028.bam
    macs2 callpeak -f BAMPE -g hs -t $dirIn/$name.q5.F1028.bam -q 0.01 -n $idh"_"$name --outdir $dirOut --call-summits
    echo -e $name"\t"$idh"\t"$(less $dirOut/"$idh"_"$name"_peaks.narrowPeak | wc -l)"\t"$(less $dirOut/"$idh"_"$name"_peaks.narrowPeak | awk '{s=s+$3-$2}END{print s}') >> $dirOut/peak_summary.txt
done

### CTCF changes 
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics2/users/lli/glioma/CTCF/MACS2/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/MACS2/unique/
mkdir -p $dirOut
cd $dirIn
cat IDHmut*ctcf*summits.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -d 50 -c 5 -o mean,count > $dirIn/IDHmut_ctcf.summits.bed
cat IDHwt*ctcf*summits.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -d 50 -c 5 -o mean,count > $dirIn/IDHwt_ctcf.summits.bed
echo -e "Sample\tIDH\tN_peak\tTotal_length" > $dirOut/unique_peak_summary.txt
for file in $dirIn/*mgh*ctcf*summits.bed; do
    name=$(basename $file | sed 's/_summits.bed//g' | sed 's/IDH.*t_//g')
    idh=$(basename $file | cut -d'_' -f1)
    ct=$(echo $idh | awk '{if($1=="IDHmut"){print "IDHwt_ctcf.summits.bed"} else{print "IDHmut_ctcf.summits.bed"}}')
    echo $name $idh $ct
    $BEDTOOLS/windowBed -a $file -b $dirIn/$ct -w 1000 -v > $dirOut/"$idh"_"$name"_unique.summits.bed
    echo -e $name"\t"$idh"\t"$(less $dirOut/"$idh"_"$name"_unique.summits.bed | wc -l)"\t"$(less $dirOut/"$idh"_"$name"_unique.summits.bed | awk '{s=s+$3-$2}END{print s}') >> $dirOut/unique_peak_summary.txt
done
$BEDTOOLS/windowBed -a $dirIn/IDHmut_ctcf.summits.bed -b $dirIn/IDHwt_ctcf.summits.bed -w 50 -u | awk '$5>1 {mid=int(($2+$3)/2); print $1"\t"mid-50"\t"mid+50"\tchr"$1":"$2"-"$3"\t"$4"\t"$5}' > $dirIn/CTCF_IDHmut_retained.bed
cat $dirOut/IDHmut*summits.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -d 50 -c 5 -o mean,count | awk '$5>1 {mid=int(($2+$3)/2); print $1"\t"mid-50"\t"mid+50"\tchr"$1":"$2"-"$3"\t"$4"\t"$5}' > $dirIn/CTCF_IDHmut_gain.bed
cat $dirOut/IDHwt*summits.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -d 50 -c 5 -o mean,count | awk '$5>1 {mid=int(($2+$3)/2); print $1"\t"mid-50"\t"mid+50"\tchr"$1":"$2"-"$3"\t"$4"\t"$5}' > $dirIn/CTCF_IDHmut_loss.bed

### 5mC
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/MACS2/WGBS/
mkdir -p $dirOut
cd $dir5mC
for CTCF in /projects/epigenomics2/users/lli/glioma/CTCF/MACS2/CTCF*.bed; do
    name=$(basename $CTCF | sed 's/.bed//g')
    echo -e "chr\tstart\tend\tID\tsummit\tN\tfractional\tsample" > $dirOut/$name.5mC
    for file in *.combine.5mC.CpG; do
        sample=$(echo $file | sed 's/.combine.5mC.CpG//g');
        echo $name $sample;
        $BEDTOOLS/intersectBed -a $file -b $CTCF -wa -wb | awk '{t[$10]=t[$10]+$4; c[$10]=c[$10]+$5; chr[$10]=$7; start[$10]=$8; end[$10]=$9; signal[$10]=$11; n[$10]=$12} END{for(i in chr){if(c[i]+t[i]>0){f=c[i]/(c[i]+t[i]); print "chr"chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"signal[i]"\t"n[i]"\t"f"\t""'$sample'"}}}' | sort -k1,1 -k 2,2n >> $dirOut/$name.5mC
    done
done

### H3K36me3
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
dirK36=/projects/epigenomics2/users/lli/glioma/ChIPseq/bam/H3K36me3/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/MACS2/H3K36me3/
mkdir -p $dirOut
cd $dirK36
for CTCF in /projects/epigenomics2/users/lli/glioma/CTCF/MACS2/CTCF*.bed; do
    name=$(basename $CTCF | sed 's/.bed//g')
    echo -e "chr\tstart\tend\tID\tsummit\tN\tN_reads\tsample\tRPKM" > $dirOut/$name.H3K36me3.bed
    for bam in *.bam; do
        sample=$(echo $bam | cut -d'.' -f2)
        echo $name $sample
        depth=$($samtools view -q 5 -F 1028 $bam | wc -l)
        $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a $CTCF -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$7*10^9/($3-$2)/"'$depth'"}' >> $dirOut/$name.H3K36me3.bed    
    done
done

### RNA-seq
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/MACS2/RNAseq/
mkdir -p $dirOut
for CTCF in /projects/epigenomics2/users/lli/glioma/CTCF/MACS2/CTCF*.bed; do
    name=$(basename $CTCF | sed 's/.bed//g')
    echo -e "chr\tstart\tend\tID\tsummit\tN\tN_reads\tsample\tRPKM" > $dirOut/$name.RNAseq.bed
    for dirRNA in /projects/epigenomics2/users/lli/glioma/RNAseq/bam/ /projects/epigenomics/users/lli/FetalBrain/RNAseq/bam/; do
        cd $dirRNA
        for bam in *.bam; do
            sample=$(echo $bam | cut -d'.' -f2)
            echo $name $sample
            depth=$($samtools view -q 5 -F 1028 $bam | wc -l)
            $samtools view -q 5 -F 1028 -b $bam | $BEDTOOLS/coverageBed -a $CTCF -b stdin -counts | awk '{print $0"\t""'$sample'""\t"$7*10^9/($3-$2)/"'$depth'"}' >> $dirOut/$name.RNAseq.bed    
        done
    done
done

### Homer 
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics2/users/lli/glioma/CTCF/MACS2/
dirOut=/projects/epigenomics2/users/lli/glioma/CTCF/MACS2/Homer/
mkdir -p $dirOut
cd $dirIn
for file in CTCF*.bed; do
    name=$(echo $file | sed 's/.bed//g')
    echo $name
    mkdir -p $dirOut/$name/
    /home/lli/bin/homer/bin/findMotifsGenome.pl <(less $file | awk '{print "chr"$0}') hg19 $dirOut/$name/ -size 200 -len 6,10,15,20 
done
