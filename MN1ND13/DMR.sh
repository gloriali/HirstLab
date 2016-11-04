#!/bin/sh

BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirIn=/projects/epigenomics2/users/lli/MN1ND13/DMR/
mkdir -p $dirIn
ln -s /projects/edcc_prj2/bs-seq/PX0409/PX0409_ACCCAG_4_lanes_dupsFlagged.ClipOverlap.Cmethyl.cons.bed.CpG.txt.gz /projects/epigenomics2/users/lli/MN1ND13/DMR/A04_MNID13_DP.bed.CpG.txt.gz
ln -s /projects/edcc_prj2/bs-seq/PX0409/PX0409_AGCGCT_4_lanes_dupsFlagged.ClipOverlap.Cmethyl.cons.bed.CpG.txt.gz /projects/epigenomics2/users/lli/MN1ND13/DMR/A05_MNID13_DN.bed.CpG.txt.gz
echo -e "sample\tmin\tymin\tlower\tmedian\tupper\tymax\tmax" > $dirIn/qc.5mC.quantile # QC: genome-wide and CGI methylation level summary; ymin: 10% quantile, ymax: 90% quantile
cd $dirIn
for file in *.bed.CpG.txt.gz; do
	lib=$(echo $file | sed 's/.bed.CpG.txt.gz//g')
    echo "Coverage profile for" $lib
    less $file | awk 'BEGIN{for(i=1;i<=5001;i++){s[i]=0}} {c=$9; if(c>=5000){s[5001]++} else {s[c]++}} END{for(i=1;i<=5001;i++){print i"\t"s[i]}}' > $file.coverage.txt
	/home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i $dirIn -o $dirIn -f $file
	less $dirIn/$file.combine.5mC.CpG | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc.5mC.quantile
    less $dirIn/$file.combine.5mC.CpG | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc.5mC.quantile
done
pth=0.0005; delta=0.6; m=0.75; cov=3; size=500; cut=3;
#pth=0.0005; delta=0.5; m=0.6; cov=3; size=500; cut=3;
#pth=0.001; delta=0.6; m=0.75; cov=3; size=500; cut=3;
#pth=0.001; delta=0.5; m=0.6; cov=3; size=500; cut=3;
#pth=0.005; delta=0.6; m=0.75; cov=3; size=500; cut=3;
#pth=0.005; delta=0.5; m=0.6; cov=3; size=500; cut=3;
file1=A04_MNID13_DP.bed.CpG.txt.gz.combine.5mC.CpG;
file2=A05_MNID13_DN.bed.CpG.txt.gz.combine.5mC.CpG;
name=MNID13_DP_MNID13_DN;
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirIn -f1 $file1 -f2 $file2 -n $name -p $pth -d $delta -m $m -c $cov
/home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirIn -o $dirIn -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
less $file1 | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1 > DP.tmp
less $file2 | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1 > DN.tmp
join DP.tmp DN.tmp | sed 's/ /\t/g' | awk 'BEGIN{for(i=1;i<=21;i++){s[i]=0}} {d=int(($2-$3)*10+11); s[d]++} END{for(i=1;i<=21;i++){print (i-11)/10"\t"s[i]}}' > DP_DN.5mC.delta
