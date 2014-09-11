#!/bin/sh

# combine MeDIP calls for all chrs
cd /projects/epigenomics/lli/FetalBrain/MeDIP/CG_25_around_chr
out='/projects/epigenomics/lli/FetalBrain/MeDIP'
for name in "HS2788.MeDIP.Brain01.q5.F1028.SET_174" "HS2790.MeDIP.Brain02.q5.F1028.SET_174" "HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174" "HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174" "HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157" "HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166"
do
    echo "$name"
    cat ./$name/*/*.dip > $out/$name.dip
done

# DM CpGs 
cd /projects/epigenomics/lli/FetalBrain/MeDIP/
delta=0.4 # difference in fractional calls to call DM CpG
dirIn='/projects/epigenomics/lli/FetalBrain/MeDIP'
dirDM=$dirIn/DMR
mkdir -p $dirDM

## Brain01 vs Brain02
lib1="HS2788"; cell1="Brain"; donor1="HuFNSC01"; name1="HS2788.MeDIP.Brain01.q5.F1028.SET_174";
lib2="HS2790"; cell2="Brain"; donor2="HuFNSC02"; name2="HS2790.MeDIP.Brain02.q5.F1028.SET_174";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed
echo $cell1"-"$donor1"_"$cell2"-"$donor2: delta=$delta
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4; else print chr"\t"start"\t"end"\t0\t"$2"\t"$4}' > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
less $dirDM/$dm | awk 'BEGIN {print "DM CpGs: "} {if($4==1)hyper++; if($4==-1)hypo++} END {print "No.of hypermethylated CpGs: ", hyper, "\nNo.of hypomethylated CpGs:  ", hypo}'
/gsc/software/linux-x86_64-centos5/python-2.7.5/bin/python /home/lli/bin/python/MeDIP.DMR.py -i $dirDM/$dm -o $dirDM
less $dirDM/*$cell1"-"$donor1"_"$cell2"-"$donor2*s*c* | awk '{l=l+$7; c=c+$6; count++; if($5==1){hyper++;hyperlen=hyperlen+$7}} END {print "Average length of DMRs:     "l/count, "\nAverage No.of CpGs per DMR: "c/count, "\nDMRs:                       "count, "regions\t", l, "bases", "\nhypermethylated DMRs:       "hyper, "regions\t", hyperlen, "bases", "\nhypomethylated DMRs:        "count-hyper, "regions\t", l-hyperlen, "bases\n\n"}'

## Cortex01 vs Cortex02
lib1="HS2775"; cell1="Cortex"; donor1="HuFNSC01"; name1="HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174";
lib2="HS2779"; cell2="Cortex"; donor2="HuFNSC02"; name2="HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed
echo $cell1"-"$donor1"_"$cell2"-"$donor2: delta=$delta
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4; else print chr"\t"start"\t"end"\t0\t"$2"\t"$4}' > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
less $dirDM/$dm | awk 'BEGIN {print "DM CpGs: "} {if($4==1)hyper++; if($4==-1)hypo++} END {print "No.of hypermethylated CpGs: ", hyper, "\nNo.of hypomethylated CpGs:  ", hypo}'
/gsc/software/linux-x86_64-centos5/python-2.7.5/bin/python /home/lli/bin/python/MeDIP.DMR.py -i $dirDM/$dm -o $dirDM
less $dirDM/*$cell1"-"$donor1"_"$cell2"-"$donor2*s*c* | awk '{l=l+$7; c=c+$6; count++; if($5==1){hyper++;hyperlen=hyperlen+$7}} END {print "Average length of DMRs:     "l/count, "\nAverage No.of CpGs per DMR: "c/count, "\nDMRs:                       "count, "regions\t", l, "bases", "\nhypermethylated DMRs:       "hyper, "regions\t", hyperlen, "bases", "\nhypomethylated DMRs:        "count-hyper, "regions\t", l-hyperlen, "bases\n\n"}'

## GE01 vs GE02
lib1="HS2777"; cell1="GE"; donor1="HuFNSC01"; name1="HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157";
lib2="HS2781"; cell2="GE"; donor2="HuFNSC02"; name2="HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166";
dm=DM.$cell1"-"$donor1"_"$cell2"-"$donor2.d$delta.bed
echo $cell1"-"$donor1"_"$cell2"-"$donor2: delta=$delta
paste $dirIn/$name1.dip $dirIn/$name2.dip | awk -v delta=$delta '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4; else print chr"\t"start"\t"end"\t0\t"$2"\t"$4}' > $dirDM/$dm
less $dirDM/$dm | grep 'Bad line'
less $dirDM/$dm | awk 'BEGIN {print "DM CpGs: "} {if($4==1)hyper++; if($4==-1)hypo++} END {print "No.of hypermethylated CpGs: ", hyper, "\nNo.of hypomethylated CpGs:  ", hypo}'
/gsc/software/linux-x86_64-centos5/python-2.7.5/bin/python /home/lli/bin/python/MeDIP.DMR.py -i $dirDM/$dm -o $dirDM
less $dirDM/*$cell1"-"$donor1"_"$cell2"-"$donor2*s*c* | awk '{l=l+$7; c=c+$6; count++; if($5==1){hyper++;hyperlen=hyperlen+$7}} END {print "Average length of DMRs:     "l/count, "\nAverage No.of CpGs per DMR: "c/count, "\nDMRs:                       "count, "regions\t", l, "bases", "\nhypermethylated DMRs:       "hyper, "regions\t", hyperlen, "bases", "\nhypomethylated DMRs:        "count-hyper, "regions\t", l-hyperlen, "bases\n\n"}'

