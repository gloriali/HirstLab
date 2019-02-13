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

