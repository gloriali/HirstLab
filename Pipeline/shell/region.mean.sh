#!/bin/sh

# Compute regional mean / sum / weighted mean 
if [ "$1" == "-h" ] ; then
    echo -e "Compute regional mean / sum / weighted mean 
Usage: `basename $0` -i <dirIn> -o <dirOut> -f <file> -r <region> -n <name> -m [-s] [-w]
    <dirIn>: input directory for the source file, optional
    <dirOut>: output directory
    <file>: input source file, format chr\tstart\tend\tvalue\t[weight], weight column is only needed for -w
    <region>: region bed file, format chr\tstart\tend\tID
    <name>: output file name
    -m: compute mean
    -s: compute sum
    -w: compute weighted mean"
    exit 0
fi

while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f) file="$2"; shift;;
        -r) region="$2"; shift;;
        -n) name="$2"; shift;;
        -m) MEAN=YES; shift;;
        -s) SUM=YES; shift;;
        -w) WEIGHT=YES; shift;;
    esac
    shift
done

/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $dirIn/$file -b $region -wa -wb > $dirOut/$name.intersect

if [ "$MEAN" = YES ]; then
    echo "Processing "$file ": compute regional mean for "$region ", output to "$dirOut/$name.mean.bed
    less $dirOut/$name.intersect | awk '{sum[$8]=sum[$8]+$4; count[$8]=count[$8]+1; chr[$8]=$5; start[$8]=$6; end[$8]=$7;} END {for(i in chr){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"sum[i]/count[i]}}' | sort -k1,1 -k 2,2n -T /projects/epigenomics/users/lli/tmp/ > $dirOut/$name.mean.bed
fi

if [ "$SUM" = YES ]; then
    echo "Processing "$file ": compute regional sum for "$region ", output to "$dirOut/$name.sum.bed
    less $dirOut/$name.intersect | awk '{sum[$8]=sum[$8]+$4; chr[$8]=$5; start[$8]=$6; end[$8]=$7;} END {for(i in chr){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"sum[i]}}' | sort -k1,1 -k 2,2n -T /projects/epigenomics/users/lli/tmp/ > $dirOut/$name.sum.bed
fi

if [ "$WEIGHT" = YES ]; then
    echo "Processing "$file ": compute regional weighted mean for "$region ", output to "$dirOut/$name.weighted.mean.bed
    less $dirOut/$name.intersect | awk '{sum[$9]=sum[$9]+$4*$5; weight[$9]=weight[$9]+$5; chr[$9]=$6; start[$9]=$7; end[$9]=$8;} END {for(i in chr){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"sum[i]/weight[i]}}' | sort -k1,1 -k 2,2n -T /projects/epigenomics/users/lli/tmp/ > $dirOut/$name.sum.bed
fi
