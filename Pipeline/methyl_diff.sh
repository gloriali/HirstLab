#!/bin/sh

# identify DM CpGs from WGBS coverage with methyl_diff
pth=0.0005
delta=0.6
dirIn=''
while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -f1) file1="$2"; shift;;
        -f2) file2="$2"; shift;;
        -n) name="$2"; shift;;
        -p) pth="$2"; shift;;
        -d) delta="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut

echo "Running methyl_diff; output to "$dirOut/DM.$name.p$pth.d$delta.bed

awk 'NR==FNR {h[$1]=$2"\t"$3; next} {if($1 in h){print $1"\t"$2"\t"$3"\t"h[$1]}}' $dirIn/$file2 $dirIn/$file1 > $dirOut/$name.join
less $dirOut/$name.join | awk '{if($2+$3+$4+$5>100000){print "0\t0\t0\t0"} else{print $2"\t"$3"\t"$4"\t"$5}}' > $dirOut/$name.input
less $dirOut/$name.input | /home/mbilenky/methyl_diff-methyl_diff/methyl_diff > $dirOut/$name.output
paste $dirOut/$name.join $dirOut/$name.output | awk '{cov1=$2+$3; cov2=$4+$5; m1=$2/cov1; m2=$4/cov2; print $1"\t"cov1"\t"m1"\t"cov2"\t"m2"\t"$6}' > $dirOut/$name.diff
less $dirOut/$name.diff | awk 'BEGIN{pth="'$pth'"+0; delta="'$delta'"+0} {d=$3-$5; chr=gensub(":.*", "", "g", $1); end=gensub(".*-", "", "g", $1); start=end-2; if(1-$6<pth && d>delta){print chr"\t"start"\t"end"\t1\t"$3"\t"$5} else if($6<pth && d<-delta){print chr"\t"start"\t"end"\t-1\t"$3"\t"$5}}' | sort -k1,1 -k2,2n > $dirOut/DM.$name.p$pth.d$delta.bed
nC1=($(wc -l $dirIn/$file1)); nC2=($(wc -l $dirIn/$file2)); nC=($(wc -l $dirOut/$name.join));
dm=($(wc -l $dirOut/DM.$name.p$pth.d$delta.bed));
hyper=($(less $dirOut/DM.$name.p$pth.d$delta.bed | awk '{if($4==1){c=c+1}} END{print c}'))
hypo=($(less $dirOut/DM.$name.p$pth.d$delta.bed | awk '{if($4==-1){c=c+1}} END{print c}'))
echo "No. of CpGs in "$file1": "$nC1
echo "No. of CpGs in "$file2": "$nC2
echo "No. of CpGs shared: "$nC
echo -e $name"\t"$pth"\t"$delta"\t"$dm"\t"$hyper"\t"$hypo >> $dirOut/DM.summary.stats
rm $dirOut/$name.input $dirOut/$name.output 
