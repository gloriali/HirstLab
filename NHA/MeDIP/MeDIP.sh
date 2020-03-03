#!/bin/sh

# QC
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/bam/
dirHub=/gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/MeDIP/hg19/
mkdir -p $dirHub
cp /gsc/www/bcgsc.ca/downloads/mb/BrainHubs/HistoneHub/genomes.txt /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/MeDIP/
echo -e "hub VitC_gliomaHub_MeDIP
shortLabel VitC_glioma Hub (MeDIP)
longLabel Hub to display VitC glioma data at UCSC (MeDIP)
genomesFile genomes.txt
email lli@bcgsc.ca" > /gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/MeDIP/hub.txt
> $dirHub/trackDb.txt
function qc {
    export PATH=/home/lli/anaconda2/bin/:$PATH
    export PYTHONPATH=/home/lli/anaconda2/lib/python2.7/site-packages
    sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
    bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/sambamba-bamStats
    chr=/projects/epigenomics2/resources/UCSC_chr/hg19.bwa2ucsc.names
    chrsize=/home/lli/hg19/hg19.chrom.sizes
    dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/bam/
    dirBW=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/bw/
    dirWig=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/wig/
    dirHub=/gsc/www/bcgsc.ca/downloads/mb/VitC_glioma/MeDIP/hg19/
    mkdir -p $dirBW; mkdir -p $dirWig; mkdir -p $dirHub
    bam=$1
    name=$(basename $bam | sed 's/.bam//g')
    echo $name 
    $sambamba index $bam -t 8
    $sambamba flagstat $bam -t 8 > $dirIn/$name.flagstat
    $bamstats -g 2864785220 -t 8 $bam > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirIn $name $dirIn/$name.bamstats
    /home/mbilenky/bin/PETLengthDist.sh $bam 5 $dirIn 10
    bamCoverage -b $bam -o $dirBW/$name.bw --normalizeUsing RPKM --ignoreDuplicates --samFlagExclude 1028 --minMappingQuality 5 --binSize 20 --extendReads --numberOfProcessors 8
    /home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirWig -F:1028,-q:5,-n:$name,-chr:$chr,-cp
    /home/lli/HirstLab/Pipeline/UCSC/wigToBigWig $dirWig/$name.q5.F1028.PET.wig.gz $chrsize $dirHub/$name.q5.F1028.PET.bw
    if [[ "$name" =~ "control" ]]; then
        color="255,0,0"
    else
        color="0,0,255"
    fi
    echo -e "
track $name
shortLabel $name
longLabel MeDIP $name
type bigWig
visibility full
maxHeightPixels 70:70:32
configurable on
autoScale on
alwaysZero on
priority 0.1
bigDataUrl $name.q5.F1028.PET.bw
color $color
" >> $dirHub/trackDb.txt
}
export -f qc
for bam in $dirIn/*.bam; do
    qc $bam
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirIn $dirIn

# MACS2
export PATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/bin:$PATH
export PYTHONPATH=/projects/edcc_new/reference_epigenomes/housekeeping/bin/anaconda/lib/python2.7/site-packages:$PYTHONPATH
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/MACS2/
mkdir -p $dirOut
echo -e "Sample\tN_region\tTotal_length\tTotal_reads\tReads_in_peaks\tAverage_length\tPercent_reads_in_peaks" > $dirOut/ER_summary.txt
for file in $dirIn/*.bam; do
    sample=$(basename $file | sed 's/.bam//g')
    echo $sample
    macs2 callpeak -f BAMPE -g hs -t $file -q 0.05 -n $sample --outdir $dirOut
    n_all=$($sambamba view $file -c -F "not (unmapped or duplicate) and mapping_quality >= 5")
    n_peak=$($sambamba view $file -c -F "not (unmapped or duplicate) and mapping_quality >= 5" -L $dirOut/$sample"_peaks.narrowPeak")
    echo -e $sample"\t"$(less $dirOut/$sample"_peaks.narrowPeak" | wc -l)"\t"$(less $dirOut/$sample"_peaks.narrowPeak" | awk '{s=s+$3-$2}END{print s}')"\t"$n_all"\t"$n_peak | awk '{print $0"\t"int($3/$2)"\t"$5/$4}' >> $dirOut/ER_summary.txt
done
export PATH=/home/lli/anaconda2/bin/:$PATH
export PYTHONPATH=/home/lli/anaconda2/lib/python2.7/site-packages
intervene upset -i $dirOut/*.narrowPeak --project upSet.MACS2 -o $dirOut 

# FindER
java=/home/mbilenky/jdk1.8.0_92/jre/bin/java
sambamba=/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5
FindER=/home/mbilenky/bin/Solexa_Java/FindER.0.9.3b.jar
reg=/projects/epigenomics/nci_rt/ChIPseq_RT/FindER/chr.regions
map=/projects/mbilenky/FindER/synthetic/hg19/PET_200/75bp/map/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/FindER/
mkdir -p $dirOut
echo -e "Sample\tN_region\tTotal_length\tTotal_reads\tReads_in_peaks\tAverage_length\tPercent_reads_in_peaks" > $dirOut/ER_summary.txt
for bam in $dirIn/*.bam; do
    sample=$(basename $bam | sed 's/.bam//g')
    echo $sample
    $java -jar -Xmx10G $FindER -i $bam -r $reg -o $dirOut -v -m $map -minER 200 -info -cp > $dirOut/FindER_scan.$sample.log
    n_all=$($sambamba view $bam -c -F "not (unmapped or duplicate) and mapping_quality >= 5")
    n_peak=$($sambamba view $bam -c -F "not (unmapped or duplicate) and mapping_quality >= 5" -L <(less $dirOut/FindER_scan.$sample.pctl_0.1.FDR_0.05.bed | sed 's/chr//g'))
    echo -e $sample"\t"$(less $dirOut/FindER_scan.$sample.pctl_0.1.FDR_0.05.bed | wc -l)"\t"$(less $dirOut/FindER_scan.$sample.pctl_0.1.FDR_0.05.bed | awk '{s=s+$3-$2}END{print s}')"\t"$n_all"\t"$n_peak | awk '{print $0"\t"int($3/$2)"\t"$5/$4}' >> $dirOut/ER_summary.txt
done
export PATH=/home/lli/anaconda2/bin/:$PATH
export PYTHONPATH=/home/lli/anaconda2/lib/python2.7/site-packages
intervene upset -i $dirOut/FindER_scan.*.bed --project upSet.FindER -o $dirOut 

# FindER2
java=/gsc/software/linux-x86_64-centos6/jdk1.8.0_162/jre/bin/java
finder2=/home/mbilenky/bin/FindER2/finder2.jar
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/bam/
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/FindER2/
mkdir -p $dirOut
echo -e "Sample\tN_region\tTotal_length\tTotal_reads\tReads_in_peaks\tAverage_length\tPercent_reads_in_peaks" > $dirOut/ER_summary.txt
for bam in $dirIn/*.bam; do
    sample=$(basename $bam | sed 's/.bam//g')
    inp=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/ChIPseq/bam/input/input.$(echo $sample | cut -f1 -d'.').bam
    echo $sample $(basename $inp)
    $java -jar -Xmx25G $finder2 inputBam:$inp signalBam:$bam outDir:$dirOut acgtDir:/projects/epigenomics2/users/mbilenky/resources/hg19/ACGT 
    n_all=$($sambamba view $bam -c -F "not (unmapped or duplicate) and mapping_quality >= 5")
    n_peak=$($sambamba view $bam -c -F "not (unmapped or duplicate) and mapping_quality >= 5" -L <(less $dirOut/$sample.FindER2.bed | sed 's/chr//g'))
    echo -e $sample"\t"$(less $dirOut/$sample.FindER2.bed | wc -l)"\t"$(less $dirOut/$sample.FindER2.bed | awk '{s=s+$3-$2}END{print s}')"\t"$n_all"\t"$n_peak | awk '{print $0"\t"int($3/$2)"\t"$5/$4}' >> $dirOut/ER_summary.txt
done
for file in $dirOut/*.FindER2.bed; do
	less $file | awk '{print $1"\t"$2"\t"$3}' > a
	mv a $file
done
export PATH=/home/lli/anaconda2/bin/:$PATH
export PYTHONPATH=/home/lli/anaconda2/lib/python2.7/site-packages
intervene upset -i $dirOut/*.FindER2.bed --project upSet.FindER2 -o $dirOut 

# fractional methylation calls
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
LIB=/home/mbilenky/bin/Solexa_Java/
csizes=/projects/epigenomics2/resources/UCSC_chr/hg19.chrom.sizes
dirw=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/wig/
for region in "CG_25_around_chr" "CG_empty_500_chr"; do
    dirr=/projects/epigenomics/users/mbilenky/CpG/hg19/$region/
    for file in $dirw/*.wig.gz; do
        name=$(basename $file | sed -e 's/.q5.F1028.PET.wig.gz//g')
        out=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/fractional/$region/$name/
        mkdir -p $out
        for chr in {1..22} "X" "Y"; do
            chr="chr"$chr
            mkdir -p $out/$chr/
            echo $region $name $chr
            $JAVA -jar -Xmx15G $LIB/RegionsCoverageFromWigCalculator.jar -w $file -r $dirr/$chr.gz -o $out/$chr/ -c $csizes -n $name
            less $out/$chr/*.coverage | awk '{gsub("chr", "", $1); print $1"_"$2"\t"$4}' > $out/$chr/$chr"."$name.cov
        done
    done
done
dirOut=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/fractional/
mkdir -p $dirOut/CDF_5mC_plots/
mkdir -p $dirOut/CDF_cov_plots/
echo "addpath /home/mbilenky/matlab/dmr -end;
close all; clear all;
set(0,'defaultaxesfontsize',18,'defaultlinelinewidth',2);
names={'MGG_control.24h','MGG_vitc.24h','MGG_vitc.48h','MGG_vitc.72h','MGG_vitc.6d'};
chrs={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'};
for i = 1:5
    name=names{1,i};
    for j = 1:24
        chr=chrs{1,j};
        close all;
        disp([name,' ',chr]);
        [l,cc] = textread(['/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/fractional/CG_25_around_chr/',name,'/',chr,'/',chr,'.',name,'.cov'],'%s %f');
        [c,n,cn]=textread(['/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/fractional/CG_empty_500_chr/',name,'/',chr,'/',chr,'.gz.',name,'.covDist'],'%f %f %f');
        x=c;
        y=cn/max(cn);
        z=cc;
        dip=medip_score2(x,y,z);
        figure('visible','off');box;
        cdfplot(dip);
        xlabel('Fractional methylation');
        ylabel('Fraction of CpGs');
        title(strcat(name,'.',chr));
        dirOut='/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/fractional/CDF_5mC_plots/';
        nameOut=strcat(dirOut, 'CDF_5mC_',name,'.',chr);
        print(gcf, '-dpdf', strcat(nameOut, '.pdf'));
        t=size(dip); n=t(2);
        fileOut = fopen(strcat('/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/fractional/CG_25_around_chr/',name,'/',chr,'/',chr,'.',name,'.dip'),'w');
        for i=1:n
            fprintf(fileOut,'%s\t', l{i});
            fprintf(fileOut,'%7.3f\t',dip(i));
            fprintf(fileOut,'\n');
        end
        fclose(fileOut);
    end
end
exit;" > $dirOut/fractional.m
/gsc/software/linux-x86_64-centos5/matlab-2013a/bin/matlab -nodisplay -nosplash -nodesktop -r "run('/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/fractional/fractional.m')"
for name in "MGG_control.24h" "MGG_vitc.24h" "MGG_vitc.48h" "MGG_vitc.72h" "MGG_vitc.6d"; do
    echo $name
    cat $dirOut/CG_25_around_chr/$name/*/*.dip > $dirOut/$name.dip
done

# check coverage profile and 5mC profile
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/fractional/
echo -e "sample\ttype\tfractional\tN" > $dirIn/qc_5mC_profile.txt 
echo -e "sample\ttype\tmin\tymin\tlower\tmedian\tupper\tymax\tmax" > $dirIn/qc_5mC_quantile.txt #ymin: 10% quantile; ymax: 90% quantile
cd $dirIn
for file in *.dip; do
    lib=$(echo $file | sed 's/.dip//g')
    echo "Processing" $lib
    less $file | awk '{s[int($2*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tgenome\t"i/100"\t"s[i]}}' >> $dirIn/qc_5mC_profile.txt 
    less $file | awk '{gsub("_", "\t"); print $1"\t"$2+23"\t"$2+25"\t"$3}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{n[$8]=n[$8]+1; c[$8]=c[$8]+$4} END{for(i in c){print c[i]/n[i]}}' | awk '{s[int($1*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tCGI\t"i/100"\t"s[i]}}' >> $dirIn/qc_5mC_profile.txt 
    less $file | awk '{print $2}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc_5mC_quantile.txt
    less $file | awk '{gsub("_", "\t"); print $1"\t"$2+23"\t"$2+25"\t"$3}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{n[$8]=n[$8]+1; c[$8]=c[$8]+$4} END{for(i in c){print c[i]/n[i]}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc_5mC_quantile.txt
done

# 5mC matrix
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/fractional/
cd $dirIn
less MGG_control.24h.dip | awk '{print $1}' > x
less /home/lli/hg19/CGI.forProfiles.BED | awk '{print $4}' | sort > a
header="ID"
for file in *.dip; do
    sample=$(echo $file | sed -e 's/.dip//g')
    header=$header" "$sample
    echo $sample
    less $file | awk -F' ' '{print $1" "$2}' | join x - > y
    mv y x
    less $file | awk '{gsub("_", "\t"); print $1"\t"$2+23"\t"$2+25"\t"$3}' | $BEDTOOLS/intersectBed -a /home/lli/hg19/CGI.forProfiles.BED -b stdin -wa -wb | awk '{n[$4]++; c[$4]=c[$4]+$8} END{for(i in c){print i"\t"c[i]/n[i]}}' | sort -k1,1 | join a - > b
    mv b a
done
echo -e $header | cat - x > matrix_genome.5mC
echo -e $header | cat - a > matrix_CGI.5mC
rm x a

# DMR between vitc and control
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/fractional/
dirOut=$dirIn/DMR/
mkdir -p $dirOut
echo -e "name\tm\tdelta\tdm\thyper\thypo" > $dirOut/DM.summary.stats
echo -e "name\tsize\tcut\tlength\tcount\tdmr\thyper\thypo" > $dirOut/DMR.summary.stats
m=0.6; delta=0.4; size=500; cut=4
file1=$dirIn/MGG_control.24h.dip; name1=$(basename $file1 | sed -e 's/.dip//g')
for file2 in $dirIn/MGG_vitc*.dip; do
    name2=$(basename $file2 | sed -e 's/.dip//g'); name=$name1"_"$name2
    echo $name1 $name2
    paste $file1 $file2 | awk -v delta=$delta -v m=$m '{if($1!=$3)print "Bad line", $0; chr="chr"gensub("_[0-9]+", "", "g", $1); start=gensub("[0-9XY]+_", "", "g", $1)+23; end=start+2; if($2-$4 > delta && $2>m)print chr"\t"start"\t"end"\t1\t"$2"\t"$4; else if($2-$4 < -delta && $4>m)print chr"\t"start"\t"end"\t-1\t"$2"\t"$4}' | sort -k1,1 -k2,2n > $dirOut/DM.$name.m$m.d$delta.bed
    less $dirOut/DM.$name.m$m.d$delta.bed | grep 'Bad line'
    Ndm=($(wc -l $dirOut/DM.$name.m$m.d$delta.bed)); Nhyper=($(less $dirOut/DM.$name.m$m.d$delta.bed | awk '{if($4==1){c=c+1}} END{print c}')); Nhypo=($(less $dirOut/DM.$name.m$m.d$delta.bed | awk '{if($4==-1){c=c+1}} END{print c}'))
    echo -e $name"\t"$m"\t"$delta"\t"$Ndm"\t"$Nhyper"\t"$Nhypo >> $dirOut/DM.summary.stats
    /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f DM.$name.m$m.d$delta.bed -n $name -s $size -c $cut
done
export PATH=/home/lli/anaconda2/bin/:$PATH
export PYTHONPATH=/home/lli/anaconda2/lib/python2.7/site-packages
intervene upset -i $dirOut/DM.*.m$m.d$delta.bed --project upSet.DM.m$m.d$delta -o $dirOut --names=MGG_vitc.24h,MGG_vitc.48h,MGG_vitc.6d,MGG_vitc.72h
intervene upset -i $dirOut/DMR.*.hyper.bed --project upSet.DMR.hyper -o $dirOut --names=MGG_vitc.24h,MGG_vitc.48h,MGG_vitc.6d,MGG_vitc.72h
intervene upset -i $dirOut/DMR.*.hypo.bed --project upSet.DMR.hypo -o $dirOut --names=MGG_vitc.24h,MGG_vitc.48h,MGG_vitc.6d,MGG_vitc.72h

# unique ER - FindER2
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
CG=/home/lli/hg19/CG.BED
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/
dirOut=$dirIn/unique/; mkdir -p $dirOut
echo -e "Sample1\Sample2\tN_control\tlen_control\tN_vitc\tlen_vitc\tN_control_unique\tlen_control_unique\tN_vitc_unique\tlen_vitc_unique" > $dirOut/ER.unique.summary
file1=$dirIn/FindER2/MGG_control.24h.FindER2.bed; lib1=$(basename $file1 | sed 's/.FindER2.bed//')
for file2 in $dirIn/FindER2/MGG_vitc.*.bed; do
    lib2=$(basename $file2 | sed 's/.FindER2.bed//')
    echo $lib1 $lib2
    /home/lli/HirstLab/Pipeline/shell/ER.unique.sh -r $file1 -w $dirIn/wig/$lib2.q5.F1028.PET.wig.gz -excl $file2 -o $dirOut -n $lib1-$lib2.$lib1 >> $dirOut/ER.unique.log
    /home/lli/HirstLab/Pipeline/shell/ER.unique.sh -excl $file1 -w $dirIn/wig/$lib1.q5.F1028.PET.wig.gz -r $file2 -o $dirOut -n $lib1-$lib2.$lib2 >> $dirOut/ER.unique.log
    less $dirOut/$lib1-$lib2.$lib1.unique | sed 's/chr//g' | $BEDTOOLS/intersectBed -a stdin -b $CG -c | awk '{if($7>=3)print}' > $dirOut/$lib1-$lib2.$lib1.c3.unique
    less $dirOut/$lib1-$lib2.$lib2.unique | sed 's/chr//g' | $BEDTOOLS/intersectBed -a stdin -b $CG -c | awk '{if($7>=3)print}' > $dirOut/$lib1-$lib2.$lib2.c3.unique
    N_control=$(less $file1 | wc -l)
    len_control=$(less $file1 | awk '{s=s+$3-$2}END{print s}')
    N_vitc=$(less $file2 | wc -l)
    len_vitc=$(less $file2 | awk '{s=s+$3-$2}END{print s}')
    N_control_unique=$(less $dirOut/$lib1-$lib2.$lib1.c3.unique | wc -l)
    len_control_unique=$(less $dirOut/$lib1-$lib2.$lib1.c3.unique | awk '{s=s+$3-$2}END{print s}')
    N_vitc_unique=$(less $dirOut/$lib1-$lib2.$lib2.c3.unique | wc -l)
    len_vitc_unique=$(less $dirOut/$lib1-$lib2.$lib2.c3.unique | awk '{s=s+$3-$2}END{print s}')
    echo -e "$lib1\t$lib2\t$N_control\t$len_control\t$N_vitc\t$len_vitc\t$N_control_unique\t$len_control_unique\t$N_vitc_unique\t$len_vitc_unique" >> $dirOut/ER.unique.summary
done
export PATH=/home/lli/anaconda2/bin/:$PATH
export PYTHONPATH=/home/lli/anaconda2/lib/python2.7/site-packages
intervene upset -i $dirOut/MGG_control.24h-MGG_vitc.*.MGG_control.24h.c3.unique --project upSet.MGG_control.unique -o $dirOut --names=MGG_vitc.24h,MGG_vitc.48h,MGG_vitc.6d,MGG_vitc.72h
intervene upset -i $dirOut/MGG_control.24h-MGG_vitc.*.MGG_vitc.*.c3.unique --project upSet.MGG_vitc.unique -o $dirOut --names=MGG_vitc.24h,MGG_vitc.48h,MGG_vitc.6d,MGG_vitc.72h
cat $dirOut/MGG_control.24h-MGG_vitc.*.MGG_control*c3.unique | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -c 1 -o count | awk '{if($4==4){print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}}' > $dirOut/MGG_control.c3.unique
cat $dirOut/MGG_control.24h-MGG_vitc.*.MGG_vitc*c3.unique | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin -c 1 -o count | awk '{if($4==4){print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}}' > $dirOut/MGG_vitc.c3.unique
intervene upset -i $dirOut/MGG_control.c3.unique $dirOut/MGG_vitc.c3.unique $dirIn/../hMeDIP/FindER2/MGG_control.unique.bed $dirIn/../hMeDIP/FindER2/MGG_vitc.unique.bed -o $dirOut --names=MeDIP.MGG_control.unique,MeDIP.MGG_vitc.unique,hMeDIP.MGG_control.unique,hMeDIP.MGG_vitc.unique
## enhancer/CGI
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
enhancer=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/ChIPseq/FindER2/H3K27ac.MGG.union.bed
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/unique/
echo -e "unique\tDMR\tcategory" > $dirIn/DMR.CGI.category
echo -e "unique\tDMR\tcategory" > $dirIn/DMR.enhancer.category
for file in $dirIn/MGG_vitc.c3.unique $dirIn/MGG_control.c3.unique; do
    DM=$(basename $file | cut -d'.' -f1)
    name=$(echo $sample)
    echo $DM $name
    less /home/lli/hg19/category.CGI.bed | sed 's/chr//g' | $BEDTOOLS/intersectBed -a $file -b stdin -wa -wb | awk '{print "'$DM'""\t"$4"\t"$9}' >> $dirIn/DMR.CGI.category
    less /home/lli/hg19/category.CGI.bed | sed 's/chr//g' | $BEDTOOLS/intersectBed -a $file -b stdin -v | awk '{print "'$DM'""\t"$4"\tnon_CGI"}' >> $dirIn/DMR.CGI.category
    $BEDTOOLS/intersectBed -a $file -b $enhancer -u | $BEDTOOLS/intersectBed -a stdin -b $promoter -u | awk '{print "'$DM'""\t"$4"\tpromoter_enhancer"}' >> $dirIn/DMR.enhancer.category
    $BEDTOOLS/intersectBed -a $file -b $enhancer -u | $BEDTOOLS/intersectBed -a stdin -b $promoter -v | awk '{print "'$DM'""\t"$4"\tdistal_enhancer"}' >> $dirIn/DMR.enhancer.category
    $BEDTOOLS/intersectBed -a $file -b $enhancer -v | awk '{print "'$DM'""\t"$4"\tnon_enhancer"}' >> $dirIn/DMR.enhancer.category
done
## homer
PATH=$PATH:/home/lli/bin/homer/bin/:/home/acarles/weblogo/
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
enhancer=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/ChIPseq/FindER2/H3K27ac.MGG.union.bed
dirIn=/projects/epigenomics3/epigenomics3_results/users/lli/MGG/MeDIP/unique/
mkdir -p $dirIn/homer/
for file in $dirIn/MGG_vitc.c3.unique $dirIn/MGG_control.c3.unique; do
    name=$(basename $file); echo $name
    mkdir -p $dirIn/homer/$name/
    $BEDTOOLS/intersectBed -a $enhancer -b $file | awk '{print "chr"$0"\t"$1":"$2"-"$3}' > $dirIn/enhancer.$name.bed
    /home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/enhancer.$name.bed hg19 $dirIn/homer/$name/ -size 200 -p 12 -len 8
done
