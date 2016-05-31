#!/bin/bash

#=====================
#Generate a QC report:
#bamstats2report.combine.sh <input directory> <output directory> 
#=====================

dirIn=$1
dirOut=$2

cp /projects/epigenomics/users/acarles/qc/summaryTemplate.11 $dirOut/summary.xls

cd $dirIn
for file in report*; do
    less $dirOut/summary.xls | sort -k1,1 > x;
    echo "$file"
    rfile=$dirIn/$file
    c=$(less $rfile | wc -l);
    if [ "$c" -eq "11" ]; then
        less $rfile | awk '{c=c+1; print c"\t"$3}' | sort -k1,1 > y;
        join x y | sed 's/ /\t/g' | sort -k1,1n > $dirOut/summary.xls ;
        rm -rf x y;
    else
        echo "$file : report has a wrong format. Skipping."
    fi
done
