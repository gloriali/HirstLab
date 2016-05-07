#!/bin/sh

# DMR between glioma and NPCs: pairwise between each glioma and all 4 NPCs and intersect
BEDTOOLS='/projects/epigenomics/software/bedtools-2.23.0/bin/'
dirIn='/projects/epigenomics2/users/lli/glioma/WGBS/'
dirOut=$dirIn/DMR/
mkdir -p $dirOut/intermediate/
cd $dirIn
> $dirOut/intermediate/DM.summary.stats
> $dirOut/intermediate/DMR.summary.stats
> $dirOut/DMR.summary.stats
pth=0.0005
delta=0.6
m=0.75
cov=3
size=500  
cut=3     
for file1 in CEMT*.5mC.CpG; do
    lib1=$(echo $file1 | sed -e 's/.5mC.CpG//g')
    echo $lib1
    for file2 in NPC*.5mC.CpG; do
        lib2=$(echo $file2 | sed -e 's/.5mC.CpG//g')
        echo $lib2
        /home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut/intermediate/ -f1 $file1 -f2 $file2 -n $lib1'_'$lib2 -p $pth -d $delta -m $m -c $cov
        /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirIn -o $dirOut/intermediate/ -f DM.$name.*.bed -n $lib1_$lib2 -s $size -c $cut
        less $dirOut/intermediate/DMR.$name.s$size.c$cut.hyper | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $dirOut/intermediate/DMR.$name.s$size.c$cut.hyper.bed 
        less $dirOut/intermediate/DMR.$name.s$size.c$cut.hypo | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $dirOut/intermediate/DMR.$name.s$size.c$cut.hypo.bed 
    done
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex02.s$size.c$cut.hyper.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex04.s$size.c$cut.hyper.bed | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex.hyper.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.GE02.s$size.c$cut.hyper.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.GE04.s$size.c$cut.hyper.bed | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC.GE.hyper.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex.hyper.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.hyper.GE.bed | awk '{print $0"\t"$1":"$2"-"$3"\t"$3-$2}' > $dirOut/DMR.$lib1'_'NPC.hyper
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex02.s$size.c$cut.hypo.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex04.s$size.c$cut.hypo.bed | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex.hypo.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.GE02.s$size.c$cut.hypo.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.GE04.s$size.c$cut.hypo.bed | awk '{print $0"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC.GE.hypo.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex.hypo.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.hypo.GE.bed | awk '{print $0"\t"$1":"$2"-"$3"\t"$3-$2}' > $dirOut/DMR.$lib1'_'NPC.hypo
    hyper=($(wc -l $dirOut/DMR.$lib1'_'NPC.hyper))
    hypo=($(wc -l $dirOut/DMR.$lib1'_'NPC.hypo))
    length_hyper=($(less $dirOut/DMR.$lib1'_'NPC.hyper | awk '{len=len+$5}END{print$5}'))
    length_hypo=($(less $dirOut/DMR.$lib1'_'NPC.hypo | awk '{len=len+$5}END{print$5}'))
    echo -e $name"\t"$size"\t"$cut"\t"$hyper"\t"$hypo"\t"$lenght_hyper"\t"$length_hypo >> $dirOut/DMR.summary.stats
done

