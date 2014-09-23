#!/bin/sh

# DM CpGs between Cortex vs GE WGBS
dirIn=/projects/epigenomics/users/lli/FetalBrain/WGBS/
dirOut=/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/
m=0.75 # fractional methylation in at least one sample need to > m  
cd $dirIn
for file in *.sam.bedGraph
do
    /home/lli/bin/shell/WGBS.combine.sh -i $dirIn -o $dirIn -f $file 
done

/home/lli/bin/shell/methyl_diff.sh -i $dirIn -o $dirOut -f1 A22475.WGBS.NeurospheresCortex02.sam.bedGraph.combine -f2 A17784-A13819.WGBS.NeurospheresGE02.sam.bedGraph.combine -n Cortex-HuFNSC02_GE-HuFNSC02
/home/lli/bin/shell/methyl_diff.sh -i $dirIn -o $dirOut -f1 A22477.WGBS.NeurospheresCortex04.sam.bedGraph.combine -f2 A22476.WGBS.NeurospheresGE04.sam.bedGraph.combine -n Cortex-HuFNSC04_GE-HuFNSC04

# Tuning parameters
cd $dirOut
> $dirOut/DM.summary.stats  # sample, p, delta, No.of DM CpGs, No.of hypermethylated DM CpGs, No.of hypomethylated DM CpGs    
> $dirOut/DMR.summary.stats # sample, size, cut, Median length of DMRs, Median No.of CpGs per DMR, No.of DMRs, No.of hypermethylated DMRs, No.of hypomethylated DMRs    
for file in *.diff
do
    for pth in 0.05 0.01 0.005 0.001 0.0005
    do
        for delta in 0.5 0.6
        do
            name=$(echo $file | sed -e s/'.diff'//g)
            echo "Processing "$name", p="$pth "delta="$delta 
            less $dirOut/$name.diff | awk 'BEGIN{pth="'$pth'"+0; delta="'$delta'"+0; m="'$m'"+0} {d=$3-$5; chr=gensub(":.*", "", "g", $1); end=gensub(".*-", "", "g", $1); start=end-2; if(1-$6<pth && d>delta && $3>m){print chr"\t"start"\t"end"\t1\t"$3"\t"$5} else if($6<pth && d<-delta && $5>m){print chr"\t"start"\t"end"\t-1\t"$3"\t"$5}}' | sort -k1,1 -k2,2n > $dirOut/DM.$name.m$m.p$pth.d$delta.bed
            dm=($(wc -l $dirOut/DM.$name.m$m.p$pth.d$delta.bed)); hyper=($(less $dirOut/DM.$name.m$m.p$pth.d$delta.bed | awk '{if($4==1){c=c+1}} END{print c}')); hypo=($(less $dirOut/DM.$name.m$m.p$pth.d$delta.bed | awk '{if($4==-1){c=c+1}} END{print c}'))
            echo -e $name"\t"$pth"\t"$delta"\t"$dm"\t"$hyper"\t"$hypo >> $dirOut/DM.summary.stats
        done
    done
done

for file in DM.*.bed
do
    name=$(echo $file | sed -e s/'DM.'//g)
    name=$(echo $name | sed -e s/'.bed'//g)
    for size in 100 200 300 500 
    do
        echo "Processing "$name", size = "$size
        /home/lli/bin/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f $file -n $name -s $size
    done
done



