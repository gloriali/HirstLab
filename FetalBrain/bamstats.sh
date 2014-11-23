#!/bin/bash

# bamstats report for FetalBrain libraries

bamstats=(
"A02879 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A02880 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A02881 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03269 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03271 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03272 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03273 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03275 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03277 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03278 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03279 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03281 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03282 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03283 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03284 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03285 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03287 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03288 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03289 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03473 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A03474 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A03475 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A03476 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A03477 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03478 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03479 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03480 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03481 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03483 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03484 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A03485 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03486 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03487 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03488 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03489 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03491 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03493 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03494 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03495 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03496 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03497 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A03499 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A04599 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A07825 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A15295 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A15298 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A15299 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A16057 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A17784 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A13819 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A19303 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A19304 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A19305 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A19306 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A19307 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A19308 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A19309 /projects/analysis/analysis*/*/hg19a/$lib/bwa*/*.bamstats"
"A22475 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A22476 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"A22477 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2774 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2775 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2776 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2777 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2778 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2779 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2780 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2781 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2787 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2788 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2789 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"HS2790 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"M01577 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"M01578 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"M01579 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"M01580 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"M01581 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"M01582 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"M01587 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"M05800 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"M05801 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
"M05802 /projects/analysis/analysis*/$lib/*/hg19a/bwa*/*.bamstats"
)

for bamstats in ${bamstats[*]}
do
    echo $bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $bamstats
done

# combine reports together
dir=/projects/epigenomics/users/acarles/qc

cd $dir
cp /projects/epigenomics/users/acarles/qc/summaryTemplate.11 summary.xls
 
for file in report_*.20141122; do
less summary.xls | sort -k1,1 > x;
echo "$file"
rfile=$dir/$file
c=$(less $rfile | wc -l);
if [ "$c" -eq "11" ]; then
less $rfile | awk '{c=c+1; print c"\t"$3}' | sort -k1,1 > y;
join x y | sed 's/ /\t/g' | sort -k1,1n > summary.xls ;
rm -rf x y;
else
  echo "$file : report has a wrong format. Skipping."
fi
done

