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

# combine coverage from opposite strands  
for file in $file1 $file2
do
    echo "Combining strands for "$file
    less $dirIn/$file | awk '{print "chr"$1":"$2"\t"$0}' > $dirOut/$file.temp
    join $dirOut/$file.temp /home/lli/hg19/CG.strand > $dirOut/$file.temp.join
    less $dirOut/$file.temp.join | awk '{if($8==id){c=c+$6; t=t+$7} else{id=$8; c=$6; t=$7; print id"\t"c"\t"t}} END{print id"\t"c"\t"t}' > $dirOut/$file.combine
done

echo "Running methyl_diff; output to "$dirOut/DM.$name.p$pth.d$delta.bed

join $dirOut/$file1.combine $dirOut/$file2.combine > $dirOut/$name.join
less $dirOut/$name.join | awk '{print $2"\t"$3"\t"$4"\t"$5}' > $dirOut/$name.input
less $dirOut/$name.input | /home/mbilenky/methyl_diff-methyl_diff/methyl_diff > $dirOut/$name.output
paste $dirOut/$name.join $dirOut/$name.output | awk '{cov1=$2+$3; cov2=$4+$5; m1=$2/cov1; m2=$4/cov2; print $1"\t"cov1"\t"m1"\t"cov2"\t"m2"\t"$6}' > $dirOut/$name.diff
less $dirOut/$name.diff | awk 'BEGIN{pth="'$pth'"+0; delta="'$delta'"+0} {d=$3-$5; chr=gensub(":.*", "", "g", $1); end=gensub(".*-", "", "g", $1); start=end-2; if(1-$6<pth && d>delta){print chr"\t"start"\t"end"\t1\t"$3"\t"$5} else if($6<pth && d<-delta){print chr"\t"start"\t"end"\t-1\t"$3"\t"$5}}' > $dirOut/DM.$name.p$pth.d$delta.bed
nC1=($(wc -l $dirOut/$file1.combine)); nC2=($(wc -l $dirOut/$file2.combine)); nC=($(wc -l $dirOut/$name.join));
dm=($(wc -l $dirOut/DM.$name.p$pth.d$delta.bed));
hyper=($(less $dirOut/DM.$name.p$pth.d$delta.bed | awk '{if($4==1){c=c+1}} END{print c}'))
hypo=($(less $dirOut/DM.$name.p$pth.d$delta.bed | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $name"\t"$pth"\t"$delta"\t"$nC1"\t"$nC2"\t"$nC"\t"$dm"\t"$hyper"\t"$hypo >> $dirOut/DM.summary.stats
rm $dirOut/$file1.* $dirOut/$file2.* $dirOut/$name.* 
