#!/bin/sh

# combine MeDIP calls for all chrs
cd /projects/epigenomics/lli/FetalBrain/MeDIP/CG_25_around_chr
out='/projects/epigenomics/lli/FetalBrain/MeDIP'
for name in "HS2788.MeDIP.Brain01.q5.F1028.SET_174" "HS2790.MeDIP.Brain02.q5.F1028.SET_174" "HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174" "HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174" "HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157" "HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166"
do
    echo "$name"
    cat ./$name/*/*.dip > $out/$name.dip
done

# DMRs
delta=0.6 # minimum difference in fractional calls to call DM CpG
size=500  # max distance between two consecutive CpGs
cut=3     # minimum number of CpGs
cd /projects/epigenomics/lli/FetalBrain/MeDIP/
dirIn='/projects/epigenomics/lli/FetalBrain/MeDIP'
dirDM=$dirIn/DMR
mkdir -p $dirDM
>$dirDM/DMR.summary.stats # samples, delta, size, cut, Average length of DMRs, Average No.of CpGs per DMR, No.of DMRs, No.of hypermethylated DMRs, No.of hypomethylated DMRs    

## Brain01 vs Brain02
lib1="HS2788"; cell1="Brain"; donor1="HuFNSC01"; name1="HS2788.MeDIP.Brain01.q5.F1028.SET_174";
lib2="HS2790"; cell2="Brain"; donor2="HuFNSC02"; name2="HS2790.MeDIP.Brain02.q5.F1028.SET_174";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{l=l+$7; c=c+$6; count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print "'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'""\t""'$delta'""\t""'$size'""\t""'$cut'""\t"l/count"\t"c/count"\t"count"\t"hyper"\t"count-hyper}' >> $dirDM/DMR.summary.stats

## Cortex01 vs Cortex02
lib1="HS2775"; cell1="Cortex"; donor1="HuFNSC01"; name1="HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174";
lib2="HS2779"; cell2="Cortex"; donor2="HuFNSC02"; name2="HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{l=l+$7; c=c+$6; count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print "'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'""\t""'$delta'""\t""'$size'""\t""'$cut'""\t"l/count"\t"c/count"\t"count"\t"hyper"\t"count-hyper}' >> $dirDM/DMR.summary.stats

## GE01 vs GE02
lib1="HS2777"; cell1="GE"; donor1="HuFNSC01"; name1="HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157";
lib2="HS2781"; cell2="GE"; donor2="HuFNSC02"; name2="HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{l=l+$7; c=c+$6; count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print "'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'""\t""'$delta'""\t""'$size'""\t""'$cut'""\t"l/count"\t"c/count"\t"count"\t"hyper"\t"count-hyper}' >> $dirDM/DMR.summary.stats

## Cortex01 vs GE01
lib1="HS2775"; cell1="Cortex"; donor1="HuFNSC01"; name1="HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174";
lib2="HS2777"; cell2="GE"; donor2="HuFNSC01"; name2="HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{l=l+$7; c=c+$6; count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print "'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'""\t""'$delta'""\t""'$size'""\t""'$cut'""\t"l/count"\t"c/count"\t"count"\t"hyper"\t"count-hyper}' >> $dirDM/DMR.summary.stats

## Cortex02 vs GE02
lib1="HS2779"; cell1="Cortex"; donor1="HuFNSC02"; name1="HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174";
lib2="HS2781"; cell2="GE"; donor2="HuFNSC02"; name2="HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed 
echo -e DMRs between $cell1"-"$donor1 and $cell2"-"$donor2: delta=$delta, size=$size, cut=$cut, output: chr"\t"start"\t"end"\t"ID"\t"DM"\t"No.of CpGs"\t"length
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
less $dirDM/$dm | awk 'BEGIN{size="'$size'"+0; cut="'$cut'"+0} {if($2<end+size && $4==dm && $1==chr){end=$3;c=c+1} else {if(end!=null){if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}; chr=$1;start=$2;end=$3;c=1;dm=$4}}END{if(c>cut){l=end-start;print chr"\t"start"\t"end"\t"chr":"start"-"end"\t"dm"\t"c"\t"l}}' > $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut
less $dirDM/DMR.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.s$size.c$cut | awk '{l=l+$7; c=c+$6; count++; if($5==1){hyper++;hyperlen=hyperlen+$7; print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hyper"} else {print $0 >> "'$dirDM'""/DMR.""'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'"".d""'$delta'"".s""'$size'"".c""'$cut'"".hypo"}} END {print "'$cell1'""-""'$donor1'""_""'$cell2'""-""'$donor2'""\t""'$delta'""\t""'$size'""\t""'$cut'""\t"l/count"\t"c/count"\t"count"\t"hyper"\t"count-hyper}' >> $dirDM/DMR.summary.stats

