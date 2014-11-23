#!/bin/bash

# bamstats report for FetalBrain libraries

bamstats=(
"A02879 ChIPseq"
"A02880 ChIPseq"
"A02881 ChIPseq"
"A03269 ChIPseq"
"A03271 ChIPseq"
"A03272 ChIPseq"
"A03273 ChIPseq"
"A03275 ChIPseq"
"A03277 ChIPseq"
"A03278 ChIPseq"
"A03279 ChIPseq"
"A03281 ChIPseq"
"A03282 ChIPseq"
"A03283 ChIPseq"
"A03284 ChIPseq"
"A03285 ChIPseq"
"A03287 ChIPseq"
"A03288 ChIPseq"
"A03289 ChIPseq"
"A03473 RNAseq"
"A03474 RNAseq"
"A03475 RNAseq"
"A03476 RNAseq"
"A03477 ChIPseq"
"A03478 ChIPseq"
"A03479 ChIPseq"
"A03480 ChIPseq"
"A03481 ChIPseq"
"A03483 ChIPseq"
"A03484 RNAseq"
"A03485 ChIPseq"
"A03486 ChIPseq"
"A03487 ChIPseq"
"A03488 ChIPseq"
"A03489 ChIPseq"
"A03491 ChIPseq"
"A03493 ChIPseq"
"A03494 ChIPseq"
"A03495 ChIPseq"
"A03496 ChIPseq"
"A03497 ChIPseq"
"A03499 ChIPseq"
"A04599 RNAseq"
"A07825 RNAseq"
"A15295 RNAseq"
"A15298 RNAseq"
"A15299 RNAseq"
"A16057 ChIPseq"
"A17784 WGBS"
"A13819 WGBS"
"A19303 ChIPseq"
"A19304 ChIPseq"
"A19305 ChIPseq"
"A19306 ChIPseq"
"A19307 ChIPseq"
"A19308 ChIPseq"
"A19309 ChIPseq"
"A22475 WGBS"
"A22476 WGBS"
"A22477 WGBS"
"HS2774 MRE"
"HS2775 MeDIP"
"HS2776 MRE"
"HS2777 MeDIP"
"HS2778 MRE"
"HS2779 MeDIP"
"HS2780 MRE"
"HS2781 MeDIP"
"HS2787 MRE"
"HS2788 MeDIP"
"HS2789 MRE"
"HS2790 MeDIP"
"M01577 microRNAseq"
"M01578 microRNAseq"
"M01579 microRNAseq"
"M01580 microRNAseq"
"M01581 microRNAseq"
"M01582 microRNAseq"
"M01587 microRNAseq"
"M05800 microRNAseq"
"M05801 microRNAseq"
"M05802 microRNAseq"
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

