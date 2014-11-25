#!/bin/bash

# bamstats path: /projects/analysis/analysis*/$lib/*/hg19a/bwa/*.bamstats
bamstats=(
"A02879"
"A02880"
"A02881"
"A03484"
"A22475"
)
for bamstats in ${bamstats[*]}
do
    echo $bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $bamstats
done

# bamstats path for inx libraries
less /projects/epigenomics/users/lli/FetalBrain/FetalBrainID_inx.txt | awk '{cmd="/home/lli/HirstLab/Pipeline/shell/bamstats2report.sh " $1 " " $2; print cmd; system(cmd)}'
/home/acarles/bin/bamstats2report.sh A03481 /projects/analysis/analysis8/INX216/702RHABXX_7/hg19a/A03481/bwa/702RHABXX_7_GATCAG.bamstats
/home/acarles/bin/bamstats2report.sh A03483 /projects/analysis8/INX216/702RHABXX_7/hg19a/A03483/bwa/702RHABXX_7_TAGCTT.bamstats

# bamstats path for merged bams: /projects/analysis6/$lib/merge_bwa/hg19a/*.bamstats
bamstats=(
"HS2774"
"HS2775"
"HS2776"
"HS2777"
"HS2778"
"HS2779"
"HS2780"
"HS2781"
"HS2787"
"HS2788"
"HS2789"
"HS2790"
)
for bamstats in ${bamstats[*]}
do
    echo $bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $bamstats
done

# For libraries without bamstats
less /projects/epigenomics/users/lli/FetalBrain/FetalBrainID_other.txt | awk '{cmd="/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py -g 2864785220 -q 10 -b " $2 " > /projects/epigenomics/users/lli/tmp/" $1 ".bamstats"; print cmd; system(cmd); cmd2="/home/acarles/bin/bamstats2report.sh " $1 " /projects/epigenomics/users/lli/tmp/" $1 ".bamstats"; print cmd2; system(cmd2)}'

# combine reports together
dir=/projects/epigenomics/users/acarles/qc
cd $dir
cp /projects/epigenomics/users/acarles/qc/summaryTemplate.11 FetalBrain_QCsummary.xls
for file in report_*.20141124 report_*.20141125; do
less FetalBrain_QCsummary.xls | sort -k1,1 > x;
echo "$file"
rfile=$dir/$file
c=$(less $rfile | wc -l);
if [ "$c" -eq "11" ]; then
less $rfile | awk '{c=c+1; print c"\t"$3}' | sort -k1,1 > y;
join x y | sed 's/ /\t/g' | sort -k1,1n > FetalBrain_QCsummary.xls ;
rm -rf x y;
else
  echo "$file : report has a wrong format. Skipping."
fi
done
rm report_*.20141124 report_*.20141125
