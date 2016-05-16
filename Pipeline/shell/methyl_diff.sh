#!/bin/sh

if [ "$1" == "-h" ] ; then
    echo -e "Identify DM CpGs from WGBS coverage with methyl_diff
Usage: `basename $0` -i <dirIn> -o <dirOut> -f1 <file1> -f2 <file2> -n <name> -p <p-value> -d <delta> -m <methylation> -c <coverage>
    <dirIn>: input directory
    <dirOut>: ourput directory
    <file1>: input for sample1, format: chr\tposition\tstrand\tconverted\tunconverted (novoalign output .5mC.CpG file)
    <file2>: input for sample2, format: chr\tposition\tstrand\tconverted\tunconverted (novoalign output .5mC.CpG file)
    <name>: output file name
    <p-value>: p-value cutoff for methyl_diff, default to 0.0005
    <delta>: min difference in fractional methylation, default to 0.6
    <methylation>: min fractional methylation for the hyper, default to 0.75
    <coverage>: min coverage for the CpG, default to 3"
    exit 0
fi

pth=0.0005
delta=0.6
m=0.75
cov=3
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
        -m) m="$2"; shift;;
        -c) cov="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut

awk 'NR==FNR {gsub("chr", ""); id=$1":"$2; if(($4+$5 >= "'$cov'"+0)&&($4+$5 <= 5000)){h[id]=$5"\t"$4}; next} {gsub("chr", ""); id=$1":"$2; if((id in h)&&($4+$5 >= "'$cov'"+0)&&($4+$5 <= 5000)){print id"\t"$5"\t"$4"\t"h[id]}}' $dirIn/$file2 $dirIn/$file1 > $dirOut/$name.join
nC1=($(wc -l $dirIn/$file1)); nC2=($(wc -l $dirIn/$file2)); nC=($(wc -l $dirOut/$name.join));
echo "No. of CpGs in file1: "$nC1
echo "No. of CpGs in file2: "$nC2
echo "No. of CpGs shared, and with enough coverage: "$nC

echo "Running methyl_diff; output to "$dirOut/DM.$name.p$pth.d$delta.bed
less $dirOut/$name.join | awk '{print $2"\t"$3"\t"$4"\t"$5}' > $dirOut/$name.input
less $dirOut/$name.input | /home/mbilenky/methyl_diff-methyl_diff/methyl_diff > $dirOut/$name.output
l1=($(wc -l $dirOut/$name.input)); l2=($(wc -l $dirOut/$name.output));
if [ "$l1" != "$l2" ]; then
    echo "methyl_diff input has" $l1 "lines, output has" $l2 "lines."; 
    exit 1
fi
paste $dirOut/$name.join $dirOut/$name.output | awk '{cov1=$2+$3; cov2=$4+$5; m1=$2/cov1; m2=$4/cov2; print $1"\t"cov1"\t"m1"\t"cov2"\t"m2"\t"$6}' > $dirOut/$name.diff
less $dirOut/$name.diff | awk 'BEGIN{pth="'$pth'"+0; delta="'$delta'"+0; m="'$m'"+0} {d=$3-$5; chr=gensub(":.*", "", "g", $1); start=gensub(".*:", "", "g", $1); end=start+1; if(1-$6<pth && d>delta && $3>m){print chr"\t"start"\t"end"\t1\t"$3"\t"$5} else if($6<pth && d<-delta && $5>m){print chr"\t"start"\t"end"\t-1\t"$3"\t"$5}}' | sort -k1,1 -k2,2n > $dirOut/DM.$name.m$m.p$pth.d$delta.bed
dm=($(wc -l $dirOut/DM.$name.m$m.p$pth.d$delta.bed)); hyper=($(less $dirOut/DM.$name.m$m.p$pth.d$delta.bed | awk '{if($4==1){c=c+1}} END{print c}')); hypo=($(less $dirOut/DM.$name.m$m.p$pth.d$delta.bed | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $name"\t"$pth"\t"$delta"\t"$m"\t"$dm"\t"$hyper"\t"$hypo >> $dirOut/DM.summary.stats
rm $dirOut/$name.join $dirOut/$name.input $dirOut/$name.output $dirOut/$name.diff
