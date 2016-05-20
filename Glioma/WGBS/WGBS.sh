#!/bin/sh

# check coverage profile
dirIn='/projects/epigenomics2/users/lli/glioma/WGBS/'
cd $dirIn
for file in *.5mC.CpG; do
    echo "Coverage profile for" $file
    less $file | awk 'BEGIN{for(i=1;i<=5001;i++){s[i]=0}} {c=$4+$5; if(c>=5000){s[5001]++} else {s[c]++}} END{for(i=1;i<=5001;i++){print i"\t"s[i]}}' > $file.coverage.txt
done

# DMR between glioma and NPCs: pairwise between each glioma and all 4 NPCs
dirIn='/projects/epigenomics2/users/lli/glioma/WGBS/'
dirOut='/projects/epigenomics2/users/lli/glioma/WGBS/DMR/'
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
mkdir -p $dirOut/intermediate/
echo -e "sample\tmin\t10%\t25%\tmedian\t75%\t90%\tmax" > $dirIn/qc.5mC.quantile # QC: genome-wide and CGI methylation level summary
echo -e "sample\tp-value\tdelta\tm\ttotal\thyper\thypo" > $dirOut/intermediate/DM.summary.stats
echo -e "sample\tsize\tcut\tmedian_length\tmedian_N_CpG\ttotal\thyper\thypo" > $dirOut/intermediate/DMR.summary.stats
pth=0.0005
delta=0.6
m=0.75
cov=3
size=500  
cut=3
cd $dirIn
for file in *.5mC.CpG; do
    lib=$(echo $file | sed -e 's/.5mC.CpG//g')
    echo -e "Combining strand for" $lib
    /home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i $dirIn -o $dirIn -f $file 
    less $dirIn/$file.combine.5mC.CpG | awk '{if($4+$5 > 0){print $5/($4+$5)}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""_genome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc.5mC.quantile
    less $dirIn/$file.combine.5mC.CpG | awk '{gsub("chr", ""); if($4+$5 > 0){print $1"\t"$2"\t"$2+2"\t"$1":"$2"\t"$4"\t"$5}}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""_CGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc.5mC.quantile
done    
for file1 in CEMT*.combine.5mC.CpG; do
    lib1=$(echo $file1 | sed -e 's/.5mC.CpG.combine.5mC.CpG//g')
    echo -e "\n\n"$lib1
    for file2 in NPC*.combine.5mC.CpG; do
        lib2=$(echo $file2 | sed -e 's/.5mC.CpG.combine.5mC.CpG//g')
        name=$lib1'_'$lib2
        echo -e "\n"$name
        /home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut/intermediate/ -f1 $file1 -f2 $file2 -n $name -p $pth -d $delta -m $m -c $cov
        /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut/intermediate/ -o $dirOut/intermediate/ -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
    done
done

# For each glioma sample, take the intersect of DMRs compared to all 4 NPCs
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirOut='/projects/epigenomics2/users/lli/glioma/WGBS/DMR/'
size=500  
cut=3
echo -e "sample\tsize\tcut\thyper\thypo\tlength_hyper\tlength_hypo" > $dirOut/DMR.summary.stats
cd /projects/epigenomics2/users/lli/glioma/WGBS/
for file1 in CEMT*.combine.5mC.CpG; do
    lib1=$(echo $file1 | sed -e 's/.5mC.CpG.combine.5mC.CpG//g')
    echo -e $lib1
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex02.s$size.c$cut.hyper.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex04.s$size.c$cut.hyper.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex.hyper.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.GE02.s$size.c$cut.hyper.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.GE04.s$size.c$cut.hyper.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC.GE.hyper.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex.hyper.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.GE.hyper.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$3-$2}' > $dirOut/DMR.$lib1'_'NPC.hyper
    less $dirOut/DMR.$lib1'_'NPC.hyper | awk '!/GL/ {print "chr"$1"\t"$2"\t"$3"\t"$4}' > $dirOut/DMR.$lib1'_'NPC.hyper.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex02.s$size.c$cut.hypo.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex04.s$size.c$cut.hypo.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex.hypo.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.GE02.s$size.c$cut.hypo.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.GE04.s$size.c$cut.hypo.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC.GE.hypo.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC.Cortex.hypo.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC.GE.hypo.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$3-$2}' > $dirOut/DMR.$lib1'_'NPC.hypo
    less $dirOut/DMR.$lib1'_'NPC.hypo | awk '!/GL/ {print "chr"$1"\t"$2"\t"$3"\t"$4}' > $dirOut/DMR.$lib1'_'NPC.hypo.bed
    hyper=($(wc -l $dirOut/DMR.$lib1'_'NPC.hyper))
    hypo=($(wc -l $dirOut/DMR.$lib1'_'NPC.hypo))
    length_hyper=($(less $dirOut/DMR.$lib1'_'NPC.hyper | awk '{len=len+$6}END{print len}'))
    length_hypo=($(less $dirOut/DMR.$lib1'_'NPC.hypo | awk '{len=len+$6}END{print len}'))
    echo -e $lib1"_NPCs\t"$size"\t"$cut"\t"$hyper"\t"$hypo"\t"$length_hyper"\t"$length_hypo >> $dirOut/DMR.summary.stats
done

