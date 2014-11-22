#!/bin/bash

#=====================
#Generate a QC report:
#=====================

dir=/projects/epigenomics/users/acarles/qc

cd $dir
cp /projects/epigenomics/users/acarles/qc/summaryTemplate.11 summary.xls
 
unset files
files=(
report_A27333_C34N5ACXX_1.20131224
report_A27333_C34N5ACXX_2.20131224
report_A27334_C34N5ACXX_3.20131224
report_A27334_C34N5ACXX_4.20131224
report_A27713_C2MEJACXX_5.20131224
report_A27713_C2MEJACXX_6.20131224
report_A27714_C2MPGACXX_1.20131224
report_A27714_C2MPGACXX_2.20131224
report_A27715_C2MPGACXX_3.20131224
report_A27715_C2MPGACXX_4.20131224
report_A27716_C2MPGACXX_5.20131224
report_A27716_C2MPGACXX_6.20131224
report_A27721_C2R88ACXX_5.20131224
report_A27721_C2R88ACXX_6.20131224
)


for file in ${files[*]}; do
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
