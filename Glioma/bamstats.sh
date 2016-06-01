#!/bin/sh

# generate QC reports from bamstats files
dirIn=/projects/edcc_new/reference_epigenomes/
dirOut=/projects/epigenomics2/users/lli/glioma/qc/
mkdir -p $dirOut
less /projects/epigenomics2/users/lli/glioma/glioma_libraries.tsv | awk 'NR > 1 {cmd="readlink -f ""'$dirIn'"$1"/bams/"$2"/"$3"*.bam"; cmd | getline bam; print $0"\t"bam"stats"}' > $dirOut/glioma.bamstats
# manually curate multiple bam files cases etc. 
less $dirOut/glioma.bamstats | awk '{cmd="/home/lli/HirstLab/Pipeline/shell/bamstats2report.sh ""'$dirOut'"" "$3" "$4; system(cmd)}' > $dirOut/bamstas2report.log
less $dirOut/bamstas2report.log | grep 'ERROR' > $dirOut/bamstas2report.error
# generate bamstats files for libraries without
less $dirOut/bamstas2report.error | awk '{print $6; gsub("bamstats", "bam"); cmd="/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py -g 2864785220 -q 10 -b " $7 " > ""'$dirOut'" $6 ".bamstats"; system(cmd); cmd2="/home/lli/HirstLab/Pipeline/shell/bamstats2report.sh ""'$dirOut'"" "$6" ""'$dirOut'" $6 ".bamstats"; system(cmd2)}'
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut
# add metadata info
awk 'NR==FNR {donor[$3]=$1; assay[$3]=$2; next} {print $0"\t"donor[$1]"\t"assay[$1]}' /projects/epigenomics2/users/lli/glioma/glioma_libraries.tsv $dirOut/summary.txt > $dirOut/summary.tsv

