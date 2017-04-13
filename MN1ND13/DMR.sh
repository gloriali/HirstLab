#!/bin/sh

## global QC
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics2/PING/DMR/
mkdir -p $dirIn
ln -s /projects/edcc_prj2/bs-seq/PX0409/PX0409_ACCCAG_4_lanes_dupsFlagged.ClipOverlap.Cmethyl.cons.bed.CpG.txt.gz /projects/epigenomics2/PING/DMR/A04_MN1ND13_DP.bed.CpG.txt.gz
ln -s /projects/edcc_prj2/bs-seq/PX0409/PX0409_AGCGCT_4_lanes_dupsFlagged.ClipOverlap.Cmethyl.cons.bed.CpG.txt.gz /projects/epigenomics2/PING/DMR/A05_MN1ND13_DN.bed.CpG.txt.gz
echo -e "sample\ttype\tmin\tymin\tlower\tmedian\tupper\tymax\tmax" > $dirIn/qc.5mC.quantile # QC: genome-wide and CGI methylation level summary; ymin: 10% quantile, ymax: 90% quantile
cd $dirIn
for file in *.bed.CpG.txt.gz; do
    lib=$(echo $file | sed 's/.bed.CpG.txt.gz//g')
    /home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i $dirIn -o $dirIn -f $file
    echo "Coverage profile for" $lib
    less $dirIn/$file.combine.5mC.CpG | awk 'BEGIN{for(i=1;i<=5001;i++){s[i]=0}} {c=$4+$5; if(c>=5000){s[5001]++} else {s[c]++}} END{for(i=1;i<=5001;i++){print i"\t"s[i]}}' > $file.coverage.txt
    less $dirIn/$file.combine.5mC.CpG | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc.5mC.quantile
    less $dirIn/$file.combine.5mC.CpG | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc.5mC.quantile
done
join <(less A04_MN1ND13_DP.bed.CpG.txt.gz.combine.5mC.CpG | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1) <(less A05_MN1ND13_DN.bed.CpG.txt.gz.combine.5mC.CpG | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1) | sed 's/ /\t/g' | awk 'BEGIN{for(i=1;i<=21;i++){s[i]=0}} {d=int(($2-$3)*10+11); s[d]++} END{for(i=1;i<=21;i++){print (i-11)/10"\t"s[i]}}' > DP_DN.5mC.delta.summary
> $dirIn/CGI.5mC
for file in *.combine.5mC.CpG; do
    lib=$(echo $file | sed 's/.combine//g' | sed 's/.5mC//g' | sed 's/.CpG//g' | sed 's/.gz//g' | sed 's/.bed.txt//g' | sed 's/CEMT_32_//g' | sed 's/A0._MN1ND13_//g')
    echo $lib
    less $dirIn/$file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print i"\t"c[i]/(c[i]+t[i])"\t""'$lib'"}}' >> $dirIn/CGI.5mC
done

## DMRs
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics2/PING/DMR/
file1=A04_MN1ND13_DP.bed.CpG.txt.gz.combine.5mC.CpG;
file2=A05_MN1ND13_DN.bed.CpG.txt.gz.combine.5mC.CpG;
name=DP_DN;
#pth=0.0005; delta=0.6; m=0.75; cov=3; size=500; cut=3;
#pth=0.0005; delta=0.5; m=0.6; cov=3; size=500; cut=3;
#pth=0.001; delta=0.6; m=0.75; cov=3; size=500; cut=3;
#pth=0.001; delta=0.5; m=0.6; cov=3; size=500; cut=3;
#pth=0.005; delta=0.6; m=0.75; cov=3; size=500; cut=3;
pth=0.005; delta=0.5; m=0.6; cov=3; size=500; cut=3;
#pth=0.01; delta=0.5; m=0.6; cov=3; size=500; cut=3;
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirIn -f1 $file1 -f2 $file2 -n $name -p $pth -d $delta -m $m -c $cov
/home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirIn -o $dirIn -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut

## Compare to CD34+ cells
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics2/PING/DMR/
ln -s /projects/edcc_prj2/bs-seq/a34042/A34042_3_lanes_dupsFlagged.q5.f0.5mC.CpG $dirIn/CEMT_32_CD34.5mC.CpG
file=CEMT_32_CD34.5mC.CpG
lib=$(echo $file | sed 's/.5mC.CpG//g' | sed 's/CEMT_32_//g')
/home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i $dirIn -o $dirIn -f $file -format novo5mC
less $dirIn/$file.combine.5mC.CpG | awk 'BEGIN{for(i=1;i<=5001;i++){s[i]=0}} {c=$4+$5; if(c>=5000){s[5001]++} else {s[c]++}} END{for(i=1;i<=5001;i++){print i"\t"s[i]}}' > $dirIn/$file.coverage.txt
less $dirIn/$file.combine.5mC.CpG | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc.5mC.quantile
less $dirIn/$file.combine.5mC.CpG | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc.5mC.quantile
join <(less A04_MN1ND13_DP.bed.CpG.txt.gz.combine.5mC.CpG | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1) <(less CEMT_32_CD34.5mC.CpG.combine.5mC.CpG | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1) | sed 's/ /\t/g' | awk 'BEGIN{for(i=1;i<=21;i++){s[i]=0}} {d=int(($2-$3)*10+11); s[d]++} END{for(i=1;i<=21;i++){print (i-11)/10"\t"s[i]}}' > DP_CD34.5mC.delta.summary
join <(less A05_MN1ND13_DN.bed.CpG.txt.gz.combine.5mC.CpG | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1) <(less CEMT_32_CD34.5mC.CpG.combine.5mC.CpG | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1) | sed 's/ /\t/g' | awk 'BEGIN{for(i=1;i<=21;i++){s[i]=0}} {d=int(($2-$3)*10+11); s[d]++} END{for(i=1;i<=21;i++){print (i-11)/10"\t"s[i]}}' > DN_CD34.5mC.delta.summary
file1=A04_MN1ND13_DP.bed.CpG.txt.gz.combine.5mC.CpG;
file2=CEMT_32_CD34.5mC.CpG.combine.5mC.CpG;
name=DP_CD34;
pth=0.005; delta=0.5; m=0.6; cov=3; size=500; cut=3;
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirIn -f1 $file1 -f2 $file2 -n $name -p $pth -d $delta -m $m -c $cov
/home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirIn -o $dirIn -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
file1=A05_MN1ND13_DN.bed.CpG.txt.gz.combine.5mC.CpG;
file2=CEMT_32_CD34.5mC.CpG.combine.5mC.CpG;
name=DN_CD34;
pth=0.005; delta=0.5; m=0.6; cov=3; size=500; cut=3;
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirIn -f1 $file1 -f2 $file2 -n $name -p $pth -d $delta -m $m -c $cov
/home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirIn -o $dirIn -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
/home/lli/HirstLab/Pipeline/shell/DMR.intersect.sh -d $dirIn

### intersect with DE genes
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
dirDE=/projects/epigenomics2/PING/DEFINE/FDR_0.01/
dirIn=/projects/epigenomics2/PING/DMR/
dirOut=$dirIn/DE/
mkdir -p $dirOut
echo -e "Sample\tDM\tDE\tN_DM_promoter\tN_DE\tN_intersect\tp_Fisher\tPercent_intersect" > $dirOut/DMR.DE.summary
n_total=19865 
cd $dirIn
for dmr in DMR.DP_DN*.bed; do
    lib=$(echo $dmr | sed -e 's/DMR.//g' | sed -e 's/.s500.c3.*.bed//g')
    dm=$(echo $dmr | sed -e 's/.bed//g' | sed -e 's/.*c3.//g')
    echo $lib $dm
    less $dmr | awk '{gsub("chr", ""); print}' | $BEDTOOLS/intersectBed -a stdin -b $promoter -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | sort -k5,5 > $dirOut/DMR.$lib.$dm.promoter.bed
    for file in $dirDE/*FDR_0.01.rmin_0.005.Nmin_25; do
        de=$(basename $file | sed -e 's/.CB_DP_CB_DN.FDR_0.01.rmin_0.005.Nmin_25//g' | sed 's/DEFINE//g')
        echo $de
        join $dirOut/DMR.$lib.$dm.promoter.bed <(less $file | sort -k1,1) -1 5 -2 1 | sed 's/ /\t/g' > $dirOut/DMR.$lib.$dm.$de
        n_dm=$(less $dirOut/DMR.$lib.$dm.promoter.bed | wc -l)
        n_de=$(less $file | wc -l)
        n_intersect=$(less $dirOut/DMR.$lib.$dm.$de | wc -l)
        p=$(echo "phyper($n_intersect, $n_dm, $n_total - $n_dm, $n_de, lower.tail = F)" | $R - | sed -e 's/\[1\] //g')
        echo -e "$lib\t$dm\t$de\t$n_dm\t$n_de\t$n_intersect\t$p" | awk '{print $0"\t"$6/$5}' >> $dirOut/DMR.DE.summary
    done
done

### intersect with SE
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
CG=/home/lli/hg19/CG.BED
dirSE=/projects/epigenomics2/PING/Superenhancers/
dirIn=/projects/epigenomics2/PING/DMR/
dirOut=$dirIn/SE/
mkdir -p $dirOut
echo -e "Sample\tDM\tSE\tN_DMR\tN_SE\tN_intersect\tp_Fisher\tFC" > $dirOut/DMR.SE.summary
CG_total=$(less $CG | wc -l) 
cd $dirIn
for dmr in DMR.DP_DN*.bed; do
    lib=$(echo $dmr | sed -e 's/DMR.//g' | sed -e 's/.s500.c3.*.bed//g')
    dm=$(echo $dmr | sed -e 's/.bed//g' | sed -e 's/.*c3.//g')
    echo $lib $dm
    for file in $dirSE/*_SuperEnhancers.bed; do
        se=$(basename $file | sed 's/_H3K27Ac_Gateway_SuperEnhancers.bed//g' | sed 's/.*MN1ND13_//g')
        echo $se
        $BEDTOOLS/intersectBed -a $dmr -b $file -wa -wb > $dirOut/DMR.$lib.$dm.SE.$se
		n_dm=$(less $dmr | wc -l)
        n_se=$(less $file | wc -l)
        n_intersect=$(less $dirOut/DMR.$lib.$dm.SE.$se | wc -l)
		CG_se=$(less $file | awk '{gsub("chr", ""); print $0}' | $BEDTOOLS/intersectBed -a $CG -b stdin -u | wc -l)
        CG_dm=$(less $dmr | awk '{gsub("chr", ""); print $0}' | $BEDTOOLS/intersectBed -a $CG -b stdin -u | wc -l)
        CG_intersect=$(less $dirOut/DMR.$lib.$dm.SE.$se | awk '{gsub("chr", ""); print $0}' | $BEDTOOLS/intersectBed -a $CG -b stdin -u | wc -l)
        p=$(echo "phyper($CG_intersect, $CG_dm, $CG_total - $CG_dm, $CG_se, lower.tail = F)" | $R - | sed -e 's/\[1\] //g')
		FC=$(echo -e "$CG_intersect\t$CG_dm\t$CG_se\t$CG_total" | awk '{print ($1/$2)/($3/$4)}')
        echo -e "$lib\t$dm\t$se\t$n_dm\t$n_se\t$n_intersect\t$p\t$FC" >> $dirOut/DMR.SE.summary
    done
done

