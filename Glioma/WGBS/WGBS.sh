#!/bin/sh

# DMR between glioma and NPCs: pairwise between each glioma and all 4 NPCs
dirIn='/projects/epigenomics2/users/lli/glioma/WGBS/'
dirOut='/home/lli/glioma/WGBS/DMR/'
mkdir -p $dirOut/intermediate/
less /projects/edcc_prj2/bs-seq/a54762/A54762_3_lanes_dupsFlagged.q5.5mC.CpG | awk '{gsub(/chr/, ""); print $0}' > $dirIn/CEMT_47.5mC.CpG # CEMT_47 had a different format
cd $dirIn
> $dirOut/intermediate/DM.summary.stats
> $dirOut/intermediate/DMR.summary.stats
pth=0.0005
delta=0.6
m=0.75
cov=3
size=500  
cut=3
for file1 in CEMT*.5mC.CpG; do
    lib1=$(echo $file1 | sed -e 's/.5mC.CpG//g')
    echo -e "\n\n"$lib1
    for file2 in NPC*.5mC.CpG; do
        lib2=$(echo $file2 | sed -e 's/.5mC.CpG//g')
        name=$lib1'_'$lib2
        echo -e "\n"$name
        /home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut/intermediate/ -f1 $file1 -f2 $file2 -n $name -p $pth -d $delta -m $m -c $cov
        /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut/intermediate/ -o $dirOut/intermediate/ -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
        less $dirOut/intermediate/DMR.$name.s$size.c$cut.hyper | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $dirOut/intermediate/DMR.$name.s$size.c$cut.hyper.bed 
        less $dirOut/intermediate/DMR.$name.s$size.c$cut.hypo | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $dirOut/intermediate/DMR.$name.s$size.c$cut.hypo.bed 
    done
done

# For each glioma sample, take the intersect of DMRs compared to all 4 NPCs
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirIn='/home/lli/glioma/WGBS/DMR/'
dirOut='/projects/epigenomics2/users/lli/glioma/WGBS/DMR/'
size=500  
cut=3
> $dirOut/DMR.summary.stats
cd /projects/epigenomics2/users/lli/glioma/WGBS/
for file1 in CEMT*.5mC.CpG; do
    lib1=$(echo $file1 | sed -e 's/.5mC.CpG//g')
    echo -e $lib1
    $BEDTOOLS/intersectBed -a $dirIn/intermediate/DMR.$lib1'_'NPC.Cortex02.s$size.c$cut.hyper.bed -b $dirIn/intermediate/DMR.$lib1'_'NPC.Cortex04.s$size.c$cut.hyper.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirIn/intermediate/DMR.$lib1'_'NPC.Cortex.hyper.bed
    $BEDTOOLS/intersectBed -a $dirIn/intermediate/DMR.$lib1'_'NPC.GE02.s$size.c$cut.hyper.bed -b $dirIn/intermediate/DMR.$lib1'_'NPC.GE04.s$size.c$cut.hyper.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirIn/intermediate/DMR.$lib1'_'NPC.GE.hyper.bed
    $BEDTOOLS/intersectBed -a $dirIn/intermediate/DMR.$lib1'_'NPC.Cortex.hyper.bed -b $dirIn/intermediate/DMR.$lib1'_'NPC.GE.hyper.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$3-$2}' > $dirOut/DMR.$lib1'_'NPC.hyper
    less $dirOut/DMR.$lib1'_'NPC.hyper | awk '$1 !~ /GL/ {print "chr"$1"\t"$2"\t"$3"\t"$4}' > $dirOut/DMR.$lib1'_'NPC.hyper.bed
    $BEDTOOLS/intersectBed -a $dirIn/intermediate/DMR.$lib1'_'NPC.Cortex02.s$size.c$cut.hypo.bed -b $dirIn/intermediate/DMR.$lib1'_'NPC.Cortex04.s$size.c$cut.hypo.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirIn/intermediate/DMR.$lib1'_'NPC.Cortex.hypo.bed
    $BEDTOOLS/intersectBed -a $dirIn/intermediate/DMR.$lib1'_'NPC.GE02.s$size.c$cut.hypo.bed -b $dirIn/intermediate/DMR.$lib1'_'NPC.GE04.s$size.c$cut.hypo.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirIn/intermediate/DMR.$lib1'_'NPC.GE.hypo.bed
    $BEDTOOLS/intersectBed -a $dirIn/intermediate/DMR.$lib1'_'NPC.Cortex.hypo.bed -b $dirIn/intermediate/DMR.$lib1'_'NPC.GE.hypo.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$3-$2}' > $dirOut/DMR.$lib1'_'NPC.hypo
    less $dirOut/DMR.$lib1'_'NPC.hypo | awk '$1 !~ /GL/ {print "chr"$1"\t"$2"\t"$3"\t"$4}' > $dirOut/DMR.$lib1'_'NPC.hypo.bed
    hyper=($(wc -l $dirOut/DMR.$lib1'_'NPC.hyper))
    hypo=($(wc -l $dirOut/DMR.$lib1'_'NPC.hypo))
    length_hyper=($(less $dirOut/DMR.$lib1'_'NPC.hyper | awk '{len=len+$5}END{print len}'))
    length_hypo=($(less $dirOut/DMR.$lib1'_'NPC.hypo | awk '{len=len+$5}END{print len}'))
    echo -e $lib1"_NPCs\t"$size"\t"$cut"\t"$hyper"\t"$hypo"\t"$length_hyper"\t"$length_hypo >> $dirOut/DMR.summary.stats
done

