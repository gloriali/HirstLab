#!/bin/bash

# Running bamstats2report.sh for a list of bamstats files

# Annaick Carles
# November 2014

# For libraries successfully found bamstats file automatically 
bamstats=(
"A25282" 
"A25281"
)
for bamstats in ${bamstats[*]}
do
    echo $bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $bamstats
done

# For libraries fail, specify location of bamstats manually
bamstats=(
"A25282 ../bams/A25282_D2F01ACXX_2.bamstats" 
"A25281 ../bams/A25281_C2B1NACXX_3.bamstats"
)
for bamstats in ${bamstats[*]}; do
echo $bamstats
/home/acarles/bin/bamstats2report.sh $bamstats
done

