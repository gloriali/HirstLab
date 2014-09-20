#!/bin/sh

# combine MeDIP calls for all chrs
cd /projects/epigenomics/users/lli/FetalBrain/MeDIP/CG_25_around_chr
out='/projects/epigenomics/users/lli/FetalBrain/MeDIP'
for name in "HS2788.MeDIP.Brain01.q5.F1028.SET_174" "HS2790.MeDIP.Brain02.q5.F1028.SET_174" "HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174" "HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174" "HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157" "HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166"
do
    echo "$name"
    cat ./$name/*/*.dip > $out/$name.dip
done

# test parameters
cd /projects/epigenomics/users/lli/FetalBrain/MeDIP/
dirIn='/projects/epigenomics/users/lli/FetalBrain/MeDIP'
dirDM=$dirIn/DMR
mkdir -p $dirDM
>$dirDM/DMR.summary.stats.test # samples, delta, size, cut, Average length of DMRs, Average No.of CpGs per DMR, No.of DMRs, No.of hypermethylated DMRs, No.of hypomethylated DMRs    
delta=0.6 # minimum difference in fractional calls to call DM CpG
cut=3     # minimum number of CpGs
lib1="HS2777"; cell1="GE"; donor1="HuFNSC01"; name1="HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157";
lib2="HS2781"; cell2="GE"; donor2="HuFNSC02"; name2="HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
for size in 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500
do
    less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
    dmr=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print count}'))
    hyper=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hyper))
    hypo=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hypo))
    length=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k7,7n | awk '{len[NR]=$7} END{if(NR%2){print len[(NR+1)/2]} else{print (len[(NR/2)]+len[(NR/2)+1])/2}}'))
    count=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k6,6n | awk '{count[NR]=$6} END{if(NR%2){print count[(NR+1)/2]} else{print (count[(NR/2)]+count[(NR/2)+1])/2}}')) 
    echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$size"\t"$cut"\t"$length"\t"$count"\t"$dmr"\t"$hyper"\t"$hypo >> $dirDM/DMR.summary.stats.test
done
less $dirDM/DMR.summary.stats.test

# DMRs
delta=0.6 # minimum difference in fractional calls to call DM CpG
size=500  # max distance between two consecutive CpGs
cut=3     # minimum number of CpGs
cd /projects/epigenomics/users/lli/FetalBrain/MeDIP/
dirIn='/projects/epigenomics/users/lli/FetalBrain/MeDIP'
dirDM=$dirIn/DMR
mkdir -p $dirDM
>$dirDM/DMR.summary.stats # samples, delta, size, cut, Average length of DMRs, Average No.of CpGs per DMR, No.of DMRs, No.of hypermethylated DMRs, No.of hypomethylated DMRs    
>$dirDM/DM.summary.stats # samples, delta, No.of DM CpGs, No.of hypermethylated DM CpGs, No.of hypomethylated DM CpGs    

## Brain01 vs Brain02
lib1="HS2788"; cell1="Brain"; donor1="HuFNSC01"; name1="HS2788.MeDIP.Brain01.q5.F1028.SET_174";
lib2="HS2790"; cell2="Brain"; donor2="HuFNSC02"; name2="HS2790.MeDIP.Brain02.q5.F1028.SET_174";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
Ndm=($(wc -l $dirDM/$dm))
Nhyper=($(less $dirDM/$dm | awk '{if($4==1){c=c+1}} END{print c}'))
Nhypo=($(less $dirDM/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirDM/DM.summary.stats
less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
dmr=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print count}'))
hyper=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hyper))
hypo=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hypo))
length=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k7,7n | awk '{len[NR]=$7} END{if(NR%2){print len[(NR+1)/2]} else{print (len[(NR/2)]+len[(NR/2)+1])/2}}'))
count=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k6,6n | awk '{count[NR]=$6} END{if(NR%2){print count[(NR+1)/2]} else{print (count[(NR/2)]+count[(NR/2)+1])/2}}')) 
echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$size"\t"$cut"\t"$length"\t"$count"\t"$dmr"\t"$hyper"\t"$hypo >> $dirDM/DMR.summary.stats.test

## Cortex01 vs Cortex02
lib1="HS2775"; cell1="Cortex"; donor1="HuFNSC01"; name1="HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174";
lib2="HS2779"; cell2="Cortex"; donor2="HuFNSC02"; name2="HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
Ndm=($(wc -l $dirDM/$dm))
Nhyper=($(less $dirDM/$dm | awk '{if($4==1){c=c+1}} END{print c}'))
Nhypo=($(less $dirDM/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirDM/DM.summary.stats
less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
dmr=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print count}'))
hyper=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hyper))
hypo=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hypo))
length=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k7,7n | awk '{len[NR]=$7} END{if(NR%2){print len[(NR+1)/2]} else{print (len[(NR/2)]+len[(NR/2)+1])/2}}'))
count=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k6,6n | awk '{count[NR]=$6} END{if(NR%2){print count[(NR+1)/2]} else{print (count[(NR/2)]+count[(NR/2)+1])/2}}')) 
echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$size"\t"$cut"\t"$length"\t"$count"\t"$dmr"\t"$hyper"\t"$hypo >> $dirDM/DMR.summary.stats.test

## GE01 vs GE02
lib1="HS2777"; cell1="GE"; donor1="HuFNSC01"; name1="HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157";
lib2="HS2781"; cell2="GE"; donor2="HuFNSC02"; name2="HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
Ndm=($(wc -l $dirDM/$dm))
Nhyper=($(less $dirDM/$dm | awk '{if($4==1){c=c+1}} END{print c}'))
Nhypo=($(less $dirDM/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirDM/DM.summary.stats
less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
dmr=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print count}'))
hyper=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hyper))
hypo=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hypo))
length=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k7,7n | awk '{len[NR]=$7} END{if(NR%2){print len[(NR+1)/2]} else{print (len[(NR/2)]+len[(NR/2)+1])/2}}'))
count=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k6,6n | awk '{count[NR]=$6} END{if(NR%2){print count[(NR+1)/2]} else{print (count[(NR/2)]+count[(NR/2)+1])/2}}')) 
echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$size"\t"$cut"\t"$length"\t"$count"\t"$dmr"\t"$hyper"\t"$hypo >> $dirDM/DMR.summary.stats.test

## Cortex01 vs GE01
lib1="HS2775"; cell1="Cortex"; donor1="HuFNSC01"; name1="HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174";
lib2="HS2777"; cell2="GE"; donor2="HuFNSC01"; name2="HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
Ndm=($(wc -l $dirDM/$dm))
Nhyper=($(less $dirDM/$dm | awk '{if($4==1){c=c+1}} END{print c}'))
Nhypo=($(less $dirDM/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirDM/DM.summary.stats
less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
dmr=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print count}'))
hyper=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hyper))
hypo=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hypo))
length=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k7,7n | awk '{len[NR]=$7} END{if(NR%2){print len[(NR+1)/2]} else{print (len[(NR/2)]+len[(NR/2)+1])/2}}'))
count=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k6,6n | awk '{count[NR]=$6} END{if(NR%2){print count[(NR+1)/2]} else{print (count[(NR/2)]+count[(NR/2)+1])/2}}')) 
echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$size"\t"$cut"\t"$length"\t"$count"\t"$dmr"\t"$hyper"\t"$hypo >> $dirDM/DMR.summary.stats.test

## Cortex02 vs GE02
lib1="HS2779"; cell1="Cortex"; donor1="HuFNSC02"; name1="HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174";
lib2="HS2781"; cell2="GE"; donor2="HuFNSC02"; name2="HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
Ndm=($(wc -l $dirDM/$dm))
Nhyper=($(less $dirDM/$dm | awk '{if($4==1){c=c+1}} END{print c}'))
Nhypo=($(less $dirDM/$dm | awk '{if($4==-1){c=c+1}} END{print c}'))
echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirDM/DM.summary.stats
less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
dmr=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print count}'))
hyper=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hyper))
hypo=($(wc -l $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut.hypo))
length=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k7,7n | awk '{len[NR]=$7} END{if(NR%2){print len[(NR+1)/2]} else{print (len[(NR/2)]+len[(NR/2)+1])/2}}'))
count=($(less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | sort -k6,6n | awk '{count[NR]=$6} END{if(NR%2){print count[(NR+1)/2]} else{print (count[(NR/2)]+count[(NR/2)+1])/2}}')) 
echo -e $cell1"-"$donor1"_"$cell2"-"$donor2"\t"$delta"\t"$size"\t"$cut"\t"$length"\t"$count"\t"$dmr"\t"$hyper"\t"$hypo >> $dirDM/DMR.summary.stats.test

