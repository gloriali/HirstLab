#!/bin/sh

# combine MeDIP calls for all chrs
#cd /projects/epigenomics/users/lli/FetalBrain/MeDIP/CG_25_around_chr
#out='/projects/epigenomics/users/lli/FetalBrain/MeDIP'
#for name in "HS2788.MeDIP.Brain01.q5.F1028.SET_174" "HS2790.MeDIP.Brain02.q5.F1028.SET_174" "HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174" "HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174" "HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157" "HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166"
#do
#    echo "$name"
#    cat ./$name/*/*.dip > $out/$name.dip
#done
#
## test parameters
#cd /projects/epigenomics/users/lli/FetalBrain/MeDIP/
#dirIn='/projects/epigenomics/users/lli/FetalBrain/MeDIP/'
#dirOut=$dirIn/DMR/
#mkdir -p $dirOut
#> $dirOut/DM.summary.stats.test # samples, delta, No.of DM CpGs, No.of hypermethylated DM CpGs, No.of hypomethylated DM CpGs    
#> $dirOut/DMR.summary.stats.test # samples, delta, size, cut, Median length of DMRs, Median No.of CpGs per DMR, No.of DMRs, No.of hypermethylated DMRs, No.of hypomethylated DMRs    
#cut=3     # minimum number of CpGs
#lib1="HS2779"; cell1="Cortex"; donor1="HuFNSC02"; name1="HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174";
#lib2="HS2781"; cell2="GE"; donor2="HuFNSC02"; name2="HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166";
#name=$cell1"-"$donor1"_"$cell2"-"$donor2
#for m in 0.75 0.8 0.9
#do
#    for delta in 0.5 0.6 0.7 0.8 0.9
#    do
#        dm=DM.$name.m$m.d$delta.bed
#        paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta -v m=$m '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta && $2>m)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta && $4>m)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirOut/$dm
#        less $dirOut/$dm | grep 'Bad line'
#        Ndm=($(wc -l $dirOut/$dm)); Nhyper=($(less $dirOut/$dm | awk '{if($4==1){c=c+1}} END{print c}')); Nhypo=($(less $dirOut/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
#        echo -e $name"\t"$m"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirOut/DM.summary.stats.test
#        for size in 200 500 1000 1500 2000 2500 3000
#        do
#            echo "Processing "$name", m = "$m", delta = "$delta", size = "$size
#            less $dirOut/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirOut/DMR.$name.m$m.d$delta.s$size.c$cut
#            dmr=($(less $dirOut/DMR.$name.m$m.d$delta.s$size.c$cut | awk '{count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirOut'""/DMR.""'$name'"".m""'$m'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirOut'""/DMR.""'$name'"".m""'$m'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print count}'))
#            hyper=($(wc -l $dirOut/DMR.$name.m$m.d$delta.s$size.c$cut.hyper))
#            hypo=($(wc -l $dirOut/DMR.$name.m$m.d$delta.s$size.c$cut.hypo))
#            length=($(less $dirOut/DMR.$name.m$m.d$delta.s$size.c$cut | sort -k7,7n | awk '{len[NR]=$7} END{if(NR%2){print len[(NR+1)/2]} else{print (len[(NR/2)]+len[(NR/2)+1])/2}}'))
#            count=($(less $dirOut/DMR.$name.m$m.d$delta.s$size.c$cut | sort -k6,6n | awk '{count[NR]=$6} END{if(NR%2){print count[(NR+1)/2]} else{print (count[(NR/2)]+count[(NR/2)+1])/2}}')) 
#            echo -e $name"\t"$m"\t"$delta"\t"$size"\t"$cut"\t"$length"\t"$count"\t"$dmr"\t"$hyper"\t"$hypo >> $dirOut/DMR.summary.stats.test
#        done
#    done
#done

# DMRs
m=0.75    # fractional methylation of one sample need to > m 
delta=0.6 # minimum difference in fractional calls to call DM CpG
size=1200  # max distance between two consecutive CpGs
cut=3     # minimum number of CpGs
cd /projects/epigenomics/users/lli/FetalBrain/MeDIP/
dirIn='/projects/epigenomics/users/lli/FetalBrain/MeDIP'
dirOut=$dirIn/DMR
mkdir -p $dirOut
#> $dirOut/DMR.summary.stats # samples, size, cut, Average length of DMRs, Average No.of CpGs per DMR, No.of DMRs, No.of hypermethylated DMRs, No.of hypomethylated DMRs    
#> $dirOut/DM.summary.stats  # samples, m, delta, No.of DM CpGs, No.of hypermethylated DM CpGs, No.of hypomethylated DM CpGs    

## Brain01 vs Brain02
lib1="HS2788"; cell1="Brain"; donor1="HuFNSC01"; name1="HS2788.MeDIP.Brain01.q5.F1028.SET_174";
lib2="HS2790"; cell2="Brain"; donor2="HuFNSC02"; name2="HS2790.MeDIP.Brain02.q5.F1028.SET_174";
name=$cell1"-"$donor1"_"$cell2"-"$donor2
dm=DM.$name.m$m.d$delta.bed
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: m=$m, delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta -v m=$m '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta && $2>m)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta && $4>m)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirOut/$dm
less $dirOut/$dm | grep 'Bad line'
Ndm=($(wc -l $dirOut/$dm)); Nhyper=($(less $dirOut/$dm | awk '{if($4==1){c=c+1}} END{print c}')); Nhypo=($(less $dirOut/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $name"\t"$m"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirOut/DM.summary.stats
/home/lli/bin/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f $dm -n $name.m$m.d$delta -s $size -c $cut

## Cortex01 vs Cortex02
lib1="HS2775"; cell1="Cortex"; donor1="HuFNSC01"; name1="HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174";
lib2="HS2779"; cell2="Cortex"; donor2="HuFNSC02"; name2="HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174";
name=$cell1"-"$donor1"_"$cell2"-"$donor2
dm=DM.$name.m$m.d$delta.bed
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: m=$m, delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta -v m=$m '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta && $2>m)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta && $4>m)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirOut/$dm
less $dirOut/$dm | grep 'Bad line'
Ndm=($(wc -l $dirOut/$dm)); Nhyper=($(less $dirOut/$dm | awk '{if($4==1){c=c+1}} END{print c}')); Nhypo=($(less $dirOut/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $name"\t"$m"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirOut/DM.summary.stats
/home/lli/bin/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f $dm -n $name.m$m.d$delta -s $size -c $cut

## GE01 vs GE02
lib1="HS2777"; cell1="GE"; donor1="HuFNSC01"; name1="HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157";
lib2="HS2781"; cell2="GE"; donor2="HuFNSC02"; name2="HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166";
name=$cell1"-"$donor1"_"$cell2"-"$donor2
dm=DM.$name.m$m.d$delta.bed
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: m=$m, delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta -v m=$m '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta && $2>m)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta && $4>m)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirOut/$dm
less $dirOut/$dm | grep 'Bad line'
Ndm=($(wc -l $dirOut/$dm)); Nhyper=($(less $dirOut/$dm | awk '{if($4==1){c=c+1}} END{print c}')); Nhypo=($(less $dirOut/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $name"\t"$m"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirOut/DM.summary.stats
/home/lli/bin/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f $dm -n $name.m$m.d$delta -s $size -c $cut

## Cortex01 vs GE01
lib1="HS2775"; cell1="Cortex"; donor1="HuFNSC01"; name1="HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174";
lib2="HS2777"; cell2="GE"; donor2="HuFNSC01"; name2="HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157";
name=$cell1"-"$donor1"_"$cell2"-"$donor2
dm=DM.$name.m$m.d$delta.bed
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: m=$m, delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta -v m=$m '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta && $2>m)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta && $4>m)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirOut/$dm
less $dirOut/$dm | grep 'Bad line'
Ndm=($(wc -l $dirOut/$dm)); Nhyper=($(less $dirOut/$dm | awk '{if($4==1){c=c+1}} END{print c}')); Nhypo=($(less $dirOut/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $name"\t"$m"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirOut/DM.summary.stats
/home/lli/bin/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f $dm -n $name.m$m.d$delta -s $size -c $cut

## Cortex02 vs GE02
lib1="HS2779"; cell1="Cortex"; donor1="HuFNSC02"; name1="HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174";
lib2="HS2781"; cell2="GE"; donor2="HuFNSC02"; name2="HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166";
name=$cell1"-"$donor1"_"$cell2"-"$donor2
dm=DM.$name.m$m.d$delta.bed
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: m=$m, delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta -v m=$m '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta && $2>m)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta && $4>m)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirOut/$dm
less $dirOut/$dm | grep 'Bad line'
Ndm=($(wc -l $dirOut/$dm)); Nhyper=($(less $dirOut/$dm | awk '{if($4==1){c=c+1}} END{print c}')); Nhypo=($(less $dirOut/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $name"\t"$m"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirOut/DM.summary.stats
/home/lli/bin/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f $dm -n $name.m$m.d$delta -s $size -c $cut
