#!/bin/sh

# DM CpGs between Cortex vs GE WGBS
/home/lli/bin/shell/methyl_diff.sh -i /projects/epigenomics/users/lli/FetalBrain/WGBS -o /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/ -f1 A22475.WGBS.NeurospheresCortex02.sam.bedGraph -f2 A17784-A13819.WGBS.NeurospheresGE02.sam.bedGraph -n Cortex-HuFNSC02_GE-HuFNSC02
/home/lli/bin/shell/methyl_diff.sh -i /projects/epigenomics/users/lli/FetalBrain/WGBS -o /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/ -f1 A22477.WGBS.NeurospheresCortex04.sam.bedGraph -f2 A22476.WGBS.NeurospheresGE04.sam.bedGraph -n Cortex-HuFNSC04_GE-HuFNSC04

# Tuning parameters
dirOut=/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR
cd $dirOut
delta=0.6
for pth in 0.05 0.01 0.005 0.001
do
    for file in *.diff
    do
        name=$(echo $file | sed -e s/'.diff'//g)
        echo "Processing "$name", p="$pth "delta="$delta 
        less $dirOut/$name.diff | awk 'BEGIN{pth="'$pth'"+0; delta="'$delta'"+0} {d=$3-$5; chr=gensub(":.*", "", "g", $1); end=gensub(".*-", "", "g", $1); start=end-2; if(1-$6<pth && d>delta){print chr"\t"start"\t"end"\t1\t"$3"\t"$5} else if($6<pth && d<-delta){print chr"\t"start"\t"end"\t-1\t"$3"\t"$5}}' | sort -k1,1 -k2,2n > $dirOut/DM.$name.p$pth.d$delta.bed
        nC1=($(wc -l $dirOut/$file1.combine)); nC2=($(wc -l $dirOut/$file2.combine)); nC=($(wc -l $dirOut/$name.join));
        dm=($(wc -l $dirOut/DM.$name.p$pth.d$delta.bed));
        hyper=($(less $dirOut/DM.$name.p$pth.d$delta.bed | awk '{if($4==1){c=c+1}} END{print c}'))
        hypo=($(less $dirOut/DM.$name.p$pth.d$delta.bed | awk '{if($4==-1){c=c+1}} END{print c}'))
        echo -e $name"\t"$pth"\t"$delta"\t"$nC1"\t"$nC2"\t"$nC"\t"$dm"\t"$hyper"\t"$hypo >> $dirOut/DM.summary.stats
    done
done

