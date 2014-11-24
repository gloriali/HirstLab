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

# bamstats path: /projects/analysis/analysis*/*/hg19a/$lib/bwa/*.bamstats
bamstats=(
"A03269"
"A03271"
"A03272"
"A03273"
"A03275"
"A03277"
"A03278"
"A03279"
"A03281"
"A03282"
"A03283"
"A03284"
"A03285"
"A03287"
"A03288"
"A03289"
"A03473"
"A03474"
"A03475"
"A03476"
"A03477"
"A03478"
"A03479"
"A03480"
"A03481"
"A03483"
"A03485"
"A03486"
"A03487"
"A03488"
"A03489"
"A03491"
"A03493"
"A03494"
"A03495"
"A03496"
"A03497"
"A03499"
"A04599"
"A07825"
"A15295"
"A15298"
"A15299"
"A16057"
"A17784"
"A13819"
"A19303"
"A19304"
"A19305"
"A19306"
"A19307"
"A19308"
"A19309"
"A22476"
"A22477"
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
"M01577"
"M01578"
"M01579"
"M01580"
"M01581"
"M01582"
"M01587"
"M05800"
"M05801"
"M05802"
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
