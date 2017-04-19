#!/bin/sh

BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirOut=/projects/epigenomics3/epigenomics3_results/alorzadeh/Sarcoma/WGBS/
dirIn=/projects/epigenomics2/SARCOMA/WGBS/
dirNPC=/projects/epigenomics2/users/lli/glioma/WGBS/
dirH1=/projects/epigenomics/ep50/5mC_new/SBS/
mkdir -p $dirOut
echo -e "sample\ttype\tmin\tymin\tlower\tmedian\tupper\tymax\tmax" > $dirOut/qc.5mC.quantile
## NPC
cd $dirNPC
for file in NPC.*.combine.5mC.CpG; do
	lib=$(echo $file | sed 's/.5mC.CpG.combine.5mC.CpG//g')
	echo "Coverage profile for" $lib
	ln -s $dirNPC/$file $dirOut/$file 
	less $dirOut/$file | awk 'BEGIN{for(i=1;i<=5001;i++){s[i]=0}} {c=$4+$5; if(c>=5000){s[5001]++} else {s[c]++}} END{for(i=1;i<=5001;i++){print i"\t"s[i]}}' > $dirOut/$lib.coverage.txt
	less $dirOut/$file | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirOut/qc.5mC.quantile
	less $dirOut/$file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirOut/qc.5mC.quantile
	less $dirOut/$file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.hg19v69_TSS2000.bed -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print i"\t"c[i]/(c[i]+t[i])"\t""'$lib'"}}' >> $dirOut/CGI.promoter.5mC
done
cd $dirOut
cat NPC*.combine.5mC.CpG | awk '{id=$1":"$2"-"$3; chr[id]=$1; start[id]=$2; end[id]=$3; t[id]=t[id]+$4; c[id]=c[id]+$5}END{for(i in chr){print chr[i]"\t"start[i]"\t"end[i]"\t"t[i]"\t"c[i]"\t"c[i]/(c[i]+t[i])}}' | sort -k1,1 -k2,2n -T /projects/epigenomics/temp/ > $dirOut/NPC.combine.5mC.CpG
## SS
ln -s $dirIn/PX0343_AACCCC_10_lanes_dupsFlagged.Cmethyl.cons.bed.CpG.txt.gz $dirOut/DG1344.Cmethyl.cons.bed.CpG.txt.gz
ln -s $dirIn/PX0343_ACCCAG_10_lanes_dupsFlagged.Cmethyl.cons.bed.CpG.txt.gz $dirOut/DG1348.Cmethyl.cons.bed.CpG.txt.gz
ln -s $dirIn/PX0343_AGCGCT_10_lanes_dupsFlagged.Cmethyl.cons.bed.CpG.txt.gz $dirOut/DG1349.Cmethyl.cons.bed.CpG.txt.gz
ln -s $dirIn/PX0343_ATCACG_10_lanes_dupsFlagged.Cmethyl.cons.bed.CpG.txt.gz $dirOut/DG1346.Cmethyl.cons.bed.CpG.txt.gz
ln -s $dirIn/PX0343_GATCAG_10_lanes_dupsFlagged.Cmethyl.cons.bed.CpG.txt.gz $dirOut/DG1343a.Cmethyl.cons.bed.CpG.txt.gz
cd $dirOut
for file in *.Cmethyl.cons.bed.CpG.txt.gz; do
	lib=$(echo $file | sed 's/.Cmethyl.cons.bed.CpG.txt.gz//g');
	echo $lib;
	/home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i $dirOut -o $dirOut -f $file
	echo "Coverage profile for" $lib
	less $dirOut/$file.combine.5mC.CpG | awk 'BEGIN{for(i=1;i<=5001;i++){s[i]=0}} {c=$4+$5; if(c>=5000){s[5001]++} else {s[c]++}} END{for(i=1;i<=5001;i++){print i"\t"s[i]}}' > $dirOut/$lib.coverage.txt
	less $dirOut/$file.combine.5mC.CpG | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirOut/qc.5mC.quantile
	less $dirOut/$file.combine.5mC.CpG | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirOut/qc.5mC.quantile
	less $dirOut/$file.combine.5mC.CpG | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.hg19v69_TSS2000.bed -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print i"\t"c[i]/(c[i]+t[i])"\t""'$lib'"}}' >> $dirOut/CGI.promoter.5mC
done
## H1: E003 and E006
cd $dirH1
for file in */E003_ESC.H1_combined */E006_ESDR.H1.MSC_combined; do
	chr=$(echo $file | sed 's/chr//g' | sed 's/\/.*//g')
	lib=$(basename $file | cut -d'_' -f1)
	echo $lib $chr
	less $file | awk '{printf "%s\t%d\t%s\t%.0f\t%.0f\n", "'$chr'", $1, ".", $2-$2*$3, $2*$3}' >> $dirOut/$lib.5mC.CpG
done
cd $dirOut
for file in E00*.5mC.CpG; do
	lib=$(echo $file | sed 's/.5mC.CpG//g')
	echo $lib
	/home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i $dirOut -o $dirOut -f $file -format novo5mC
	echo "Coverage profile for" $lib
	less $dirOut/$file.combine.5mC.CpG | awk 'BEGIN{for(i=1;i<=5001;i++){s[i]=0}} {c=$4+$5; if(c>=5000){s[5001]++} else {s[c]++}} END{for(i=1;i<=5001;i++){print i"\t"s[i]}}' > $dirOut/$lib.coverage.txt
	less $dirOut/$file.combine.5mC.CpG | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirOut/qc.5mC.quantile
	less $dirOut/$file.combine.5mC.CpG | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirOut/qc.5mC.quantile
	less $dirOut/$file.combine.5mC.CpG | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.hg19v69_TSS2000.bed -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print i"\t"c[i]/(c[i]+t[i])"\t""'$lib'"}}' >> $dirOut/CGI.promoter.5mC
done

