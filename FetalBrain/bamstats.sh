#!/bin/bash

# bamstats report for FetalBrain libraries

less /projects/epigenomics/users/lli/FetalBrain/FetalBrainID.txt | awk '{print $1; system("/home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $1 $2")}'

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

