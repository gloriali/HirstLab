#!/bin/sh

# ============= DNAme =============
## check coverage profile and 5mC profile
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/AML/DNAme/
echo -e "sample\tcoverage\tN" > $dirIn/qc_5mC_coverage.txt
echo -e "sample\ttype\tfractional\tN" > $dirIn/qc_5mC_profile.txt 
echo -e "sample\ttype\tmin\tymin\tlower\tmedian\tupper\tymax\tmax" > $dirIn/qc_5mC_quantile.txt #ymin: 10% quantile; ymax: 90% quantile
cd $dirIn
for file in *.combine.5mC.CpG; do
    lib=$(echo $file | sed -e 's/.combine.5mC.CpG//g')
    echo "Processing" $lib
    less $file | awk '{c = $4 + $5; if(c >= 5000){s[5001]++} else {s[c]++}} END{for(i = 1; i <= 5001; i++){print "'$lib'""\t"i"\t"s[i]}}' >> $dirIn/qc_5mC_coverage.txt
    less $file | awk '{s[int($6*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tgenome\t"i/100"\t"s[i]}}' >> $dirIn/qc_5mC_profile.txt 
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | awk '{s[int($1*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tCGI\t"i/100"\t"s[i]}}' >> $dirIn/qc_5mC_profile.txt 
    less $file | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc_5mC_quantile.txt
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc_5mC_quantile.txt
done

## 5mC matrix
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/AML/DNAme/
cd $dirIn
less /home/lli/hg19/CG.strand | awk '{if(NR%2){print $2}}' | sort > x
less /home/lli/hg19/CGI.forProfiles.BED | awk '{print $4}' | sort > a
header="ID"
for file in *combine.5mC.CpG; do
    sample=$(echo $file | sed -e 's/.combine.5mC.CpG//g')
    header=$header" "$sample
    echo $sample
    less $file | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1 | join x - > y
    mv y x
    less $file | $BEDTOOLS/intersectBed -a /home/lli/hg19/CGI.forProfiles.BED -b stdin -wa -wb | awk '{t[$4]=t[$4]+$8; c[$4]=c[$4]+$9} END{for(i in t){if(t[i]+c[i]>0){print i"\t"c[i]/(c[i]+t[i])}}}' | sort -k1,1 | join a - > b
    mv b a
done
echo -e $header | cat - x > matrix_genome.5mC
echo -e $header | cat - a > matrix_CGI.5mC
rm x a

## DMR between AML and Erythroid
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/AML/DNAme/
dirOut=$dirIn/DMR/; mkdir -p $dirOut
echo -e "sample\tp-value\tdelta\tm\ttotal\thyper\thypo" > $dirOut/DM.summary.stats
echo -e "sample\tsize\tcut\tmedian_length\tmedian_N_CpG\ttotal\thyper\thypo" > $dirOut/DMR.summary.stats
pth=0.0005; delta=0.6; m=0.75; cov=3; size=500; cut=3
file2=PX0681_TATAAT.Normal_Erythroid.combine.5mC.CpG; lib2=Normal_Erythroid
cd $dirIn
for file1 in *AML*.combine.5mC.CpG; do
    lib1=$(echo $file1 | sed -e 's/.combine.5mC.CpG//g')
    name=$lib1'-'$lib2
    echo $lib1 $name
    /home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut -f1 $file1 -f2 $file2 -n $name -p $pth -d $delta -m $m -c $cov
    /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
done
cat $dirOut/DMR.*IDH*R*.s500.c3.hyper.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -c 1 -o count | awk '{if($4>2)print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/DMR.AML_IDHmut.Normal_Erythroid.hyper.bed
cat $dirOut/DMR.*IDH*R*.s500.c3.hypo.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -c 1 -o count | awk '{if($4>2)print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/DMR.AML_IDHmut.Normal_Erythroid.hypo.bed
cat $dirOut/DMR.*IDHwt*.s500.c3.hyper.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -c 1 -o count | awk '{if($4>1)print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/DMR.AML_IDHwt.Normal_Erythroid.hyper.bed
cat $dirOut/DMR.*IDHwt*.s500.c3.hypo.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -c 1 -o count | awk '{if($4>1)print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/DMR.AML_IDHwt.Normal_Erythroid.hypo.bed
### genomic enrichment
enhancer=/projects/epigenomics3/epigenomics3_results/users/lli/AML/ChIPseq/FindER2/H3K27ac.IDHmut.bed
/home/lli/HirstLab/Pipeline/shell/DMR.intersect.sh -d $dirOut -r $enhancer -n enhancer
### DMR upSet
export PATH=/home/lli/anaconda2/bin/:$PATH
export PYTHONPATH=/home/lli/anaconda2/lib/python2.7/site-packages
intervene upset -i $dirOut/DMR.*.s500.c3.hyper.bed --project DMR.AML_Erythroid.hyper -o $dirOut
intervene upset -i $dirOut/DMR.*.s500.c3.hypo.bed --project DMR.AML_Erythroid.hypo -o $dirOut


# ============= ChIPseq =============

## FindER2
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/AML/ChIPseq/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/AML/ChIPseq/FindER2/
dirAli=/projects/epigenomics3/epigenomics3_results/users/alorzadeh/Finder2.0/Finder2_2019/
ID=/projects/epigenomics3/epigenomics3_results/users/lli/AML/ChIPseq/ID.tsv
for file in $dirAli/*[0-9]-[0-9]*.FindER2.bed $dirAli/*[0-9]_[0-9]*_[AB]*.FindER2.bed; do
    sample=$(basename $file | sed 's/_H.*.bed//' | sed 's/B0[1-6]_//' | sed 's/_A.*//' | sed 's/_B.*//' | sed 's/-AML//')
    mark=$(basename $file | sed 's/.*_H/H/' | sed 's/.FindER2.bed//')
    name=$(less $ID | awk '$3 ~ "'$sample'" {print $1}')
    echo $sample $mark $name
    ln -s $file $dirOut/$mark.$name.FindER2.bed
done
cat $dirOut/H3K27ac.*AML_IDH*R*.FindER2.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -c 1 -o count | awk '{if($4>2)print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/H3K27ac.IDHmut.bed

