#!/bin/sh

# DM CpGs between Cortex vs GE WGBS
dirIn=/projects/epigenomics/users/lli/FetalBrain/WGBS/
dirOut=/projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/
m=0.75 # fractional methylation in at least one sample need to > m  
cd $dirIn
for file in *.sam.bedGraph
do
    /home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i $dirIn -o $dirIn -f $file 
done

/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut -f1 A22475.WGBS.NeurospheresCortex02.sam.bedGraph.combine -f2 A17784-A13819.WGBS.NeurospheresGE02.sam.bedGraph.combine -n Cortex-HuFNSC02_GE-HuFNSC02
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut -f1 A22477.WGBS.NeurospheresCortex04.sam.bedGraph.combine -f2 A22476.WGBS.NeurospheresGE04.sam.bedGraph.combine -n Cortex-HuFNSC04_GE-HuFNSC04

# Tuning parameters
cd $dirOut
> $dirOut/DM.summary.stats  # sample, p, delta, No.of DM CpGs, No.of hypermethylated DM CpGs, No.of hypomethylated DM CpGs    
> $dirOut/DMR.summary.stats # sample, size, cut, Median length of DMRs, Median No.of CpGs per DMR, No.of DMRs, No.of hypermethylated DMRs, No.of hypomethylated DMRs    
for file in *.diff
do
    for pth in 0.05 0.01 0.005 0.001 0.0005
    do
        for delta in 0.5 0.6
        do
            name=$(echo $file | sed -e s/'.diff'//g)
            echo "Processing "$name", p="$pth "delta="$delta 
            less $dirOut/$name.diff | awk 'BEGIN{pth="'$pth'"+0; delta="'$delta'"+0; m="'$m'"+0} {d=$3-$5; chr=gensub(":.*", "", "g", $1); end=gensub(".*-", "", "g", $1); start=end-2; if(1-$6<pth && d>delta && $3>m){print chr"\t"start"\t"end"\t1\t"$3"\t"$5} else if($6<pth && d<-delta && $5>m){print chr"\t"start"\t"end"\t-1\t"$3"\t"$5}}' | sort -k1,1 -k2,2n > $dirOut/DM.$name.m$m.p$pth.d$delta.bed
            dm=($(wc -l $dirOut/DM.$name.m$m.p$pth.d$delta.bed)); hyper=($(less $dirOut/DM.$name.m$m.p$pth.d$delta.bed | awk '{if($4==1){c=c+1}} END{print c}')); hypo=($(less $dirOut/DM.$name.m$m.p$pth.d$delta.bed | awk '{if($4==-1){c=c+1}} END{print c}'))
            echo -e $name"\t"$pth"\t"$delta"\t"$dm"\t"$hyper"\t"$hypo >> $dirOut/DM.summary.stats
        done
    done
done

for file in DM.*.bed
do
    name=$(echo $file | sed -e s/'DM.'//g)
    name=$(echo $name | sed -e s/'.bed'//g)
    for size in 100 200 300 500 
    do
        echo "Processing "$name", size = "$size
        /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f $file -n $name -s $size
    done
done

# set parameters to m=0.75, p=0.005, delta=0.5, size=300
# Bed files for intersect DMR and GREAT  
for f in $dirOut/*.m0.75.p0.005.d0.5.s300*.hyp*
do
    echo "Processing "$f
    less $f | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $f.bed # input for GREAT analysis
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/closestBed -a $f.bed -b /home/lli/hg19/hg19v65_genes.bed -d > $f.closest.gene
done
/home/lli/HirstLab/Pipeline/shell/DMR.intersect.sh -d $dirOut

# validate DMR MeDIP/MRE signal
mkdir -p $dirOut/valid/
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
reg=$dirOut/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed # GE02 UMR
wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2779.bam.q5.F1028.SET_174.wig.gz       # MeDIP Cortex02
name=GE02_UMRs_MeDIP_cortex02
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirOut/valid/ -s $chr -n $name -t Y 
wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2781.bam.q5.F1028.SET_166.wig.gz       # MeDIP GE02
name=GE02_UMRs_MeDIP_GE02
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirOut/valid/ -s $chr -n $name -t Y 
wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2778.bam.q5.F1028.SET_50.wig.gz        # MRE Cortex02
name=GE02_UMRs_MRE_cortex02
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirOut/valid/ -s $chr -n $name -t Y 
wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2780.bam.q5.F1028.SET_50.wig.gz        # MRE GE02
name=GE02_UMRs_MRE_GE02
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirOut/valid/ -s $chr -n $name -t Y 
reg=$dirOut/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed # Cortex02 UMR
wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2779.bam.q5.F1028.SET_174.wig.gz       # MeDIP Cortex02
name=Cortex02_UMRs_MeDIP_cortex02
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirOut/valid/ -s $chr -n $name -t Y 
wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2781.bam.q5.F1028.SET_166.wig.gz       # MeDIP GE02
name=Cortex02_UMRs_MeDIP_GE02
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirOut/valid/ -s $chr -n $name -t Y 
wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2778.bam.q5.F1028.SET_50.wig.gz        # MRE Cortex02
name=Cortex02_UMRs_MRE_cortex02
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirOut/valid/ -s $chr -n $name -t Y 
wig=/home/lli/FetalBrain/MeDIPMRE/wigs/HS2780.bam.q5.F1028.SET_50.wig.gz        # MRE GE02
name=Cortex02_UMRs_MRE_GE02
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $dirOut/valid/ -s $chr -n $name -t Y 

# intersecting with TFBSs
mkdir -p $dirOut/TF/
> $dirOut/TF/DMR.TF.summary
cd $dirOut
for file in DMR.*.bed
do
    name=$(echo $file | sed -e s/'.bed'//g)
    echo "Processing "$name >> $dirOut/TF/DMR.TF.summary
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $file -b /projects/mbilenky/REMC/breast/ENCODE/TFs/wgEncodeRegTfbsClusteredV3.bed.gz -wa -wb | awk '{print $4"\t"$5":"$6"-"$7"\t"$8}' > $dirOut/TF/$name.TF
    less $dirOut/TF/$name.TF | awk '{if(!($1 in h1)){c1=c1+1; h1[$1]=1} if(!($2 in h2)){c2=c2+1; h2[$2]=1} if(!($3 in h3)){c3=c3+1; h3[$3]=1} else{h3[$3]=h3[$3]+1}} END{print "No. of DMRs overlap with TFBSs", c1; print "No. of unique TFBSs", c2; print "No. of unique TFs", c3; for(i in h3){print i"\t"h3[i] >> "'$dirOut'""/TF/""'$name'"".TF.summary"}}' >> $dirOut/TF/DMR.TF.summary
done
awk 'NR==FNR {h[$1]=$2; next} {if($1 in h){print $0"\t"h[$1]"\t"$2/h[$1]}}' $dirOut/TF/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.TF.summary $dirOut/TF/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.TF.summary | sort -k4,4n > $dirOut/TF/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.TF.summary
awk 'NR==FNR {h[$1]=$2; next} {if($1 in h){print $0"\t"h[$1]"\t"$2/h[$1]}}' $dirOut/TF/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.TF.summary $dirOut/TF/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.TF.summary | sort -k4,4n > $dirOut/TF/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.TF.summary
awk 'NR==FNR {cortex[$1]=$2; ge[$1]=$3; ratio[$1]=$4; next} {if($1 in cortex){print $0"\t"cortex[$1]"\t"ge[$1]"\t"ratio[$1]}}' $dirOut/TF/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.TF.summary $dirOut/TF/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.TF.summary | sort -k4,4n > $dirOut/TF/DMR.Cortex_GE.m0.75.p0.005.d0.5.s300.c3.TF.summary

# enrichment at chr ends
cd /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/
for file in DMR.*.bed
do
    echo "Processing "$file
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/hg19.chrom.window.n10.bed -b $file -c > $file.chrom.window
done

##########################################################################
# GW-related UMRs: Cortex02 vs Cortex04; GE02 vs GE04 (17-week vs 13-week)
dirIn=/projects/epigenomics/users/lli/FetalBrain/WGBS/
dirOut=/projects/epigenomics/users/lli/FetalBrain/GW/DMR/
mkdir -p $dirOut
cd $dirOut
> $dirOut/DM.summary.stats  # sample, p, delta, No.of DM CpGs, No.of hypermethylated DM CpGs, No.of hypomethylated DM CpGs    
> $dirOut/DMR.summary.stats # sample, size, cut, Median length of DMRs, Median No.of CpGs per DMR, No.of DMRs, No.of hypermethylated DMRs, No.of hypomethylated DMRs    
m=0.75 
p=0.005
delta=0.5
size=300
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut -f1 A22475.WGBS.NeurospheresCortex02.sam.bedGraph.combine -f2 A22477.WGBS.NeurospheresCortex04.sam.bedGraph.combine -n Cortex-HuFNSC02_Cortex-HuFNSC04 -p $p -d $delta
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut -f1 A17784-A13819.WGBS.NeurospheresGE02.sam.bedGraph.combine -f2 A22476.WGBS.NeurospheresGE04.sam.bedGraph.combine -n GE-HuFNSC02_GE-HuFNSC04 -p $p -d $delta
for file in DM.*.bed
do
    name=$(echo $file | sed -e s/'DM.'//g)
    name=$(echo $name | sed -e s/'.bed'//g)
    echo "Processing "$name", size = "$size
    /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut -o $dirOut -f $file -n $name -s $size
done
for f in $dirOut/*.hyp*
do
    echo "Processing "$f
    less $f | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $f.bed # input for GREAT analysis
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/closestBed -a $f.bed -b /home/lli/hg19/hg19v65_genes.bed -d > $f.closest.gene
done
/home/lli/HirstLab/Pipeline/shell/DMR.intersect.sh -d $dirOut
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed -b DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed > DMR.HuFNSC02_HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed -b DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed > DMR.HuFNSC02_HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed
## intersect with TFBS
mkdir -p $dirOut/TF/
> $dirOut/TF/DMR.TF.summary
cd $dirOut
for file in DMR.*.bed
do
    name=$(echo $file | sed -e s/'.bed'//g)
    echo "Processing "$name >> $dirOut/TF/DMR.TF.summary
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $file -b /projects/mbilenky/REMC/breast/ENCODE/TFs/wgEncodeRegTfbsClusteredV3.bed.gz -wa -wb | awk '{print $4"\t"$5":"$6"-"$7"\t"$8}' > $dirOut/TF/$name.TF
    less $dirOut/TF/$name.TF | awk '{if(!($1 in h1)){c1=c1+1; h1[$1]=1} if(!($2 in h2)){c2=c2+1; h2[$2]=1} if(!($3 in h3)){c3=c3+1; h3[$3]=1} else{h3[$3]=h3[$3]+1}} END{print "No. of DMRs overlap with TFBSs", c1; print "No. of unique TFBSs", c2; print "No. of unique TFs", c3; for(i in h3){print i"\t"h3[i] >> "'$dirOut'""/TF/""'$name'"".TF.summary"}}' >> $dirOut/TF/DMR.TF.summary
done
awk 'NR==FNR {h[$1]=$2; next} {if($1 in h){print $0"\t"h[$1]"\t"$2/h[$1]}}' $dirOut/TF/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.TF.summary $dirOut/TF/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.TF.summary | sort -k4,4n > $dirOut/TF/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.TF.summary
awk 'NR==FNR {h[$1]=$2; next} {if($1 in h){print $0"\t"h[$1]"\t"$2/h[$1]}}' $dirOut/TF/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.TF.summary $dirOut/TF/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.TF.summary | sort -k4,4n > $dirOut/TF/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.TF.summary
awk 'NR==FNR {late[$1]=$2; early[$1]=$3; ratio[$1]=$4; next} {if($1 in ratio){print $0"\t"late[$1]"\t"early[$1]"\t"ratio[$1]}}' $dirOut/TF/DMR.GE-HuFNSC02_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.TF.summary $dirOut/TF/DMR.Cortex-HuFNSC02_Cortex-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.TF.summary | sort -k4,4n > $dirOut/TF/DMR.HuFNSC02_HuFNSC04.m0.75.p0.005.d0.5.s300.c3.TF.summary
## chr End enrichment
cd /projects/epigenomics/users/lli/FetalBrain/GW/DMR/
for file in DMR.*.bed
do
    echo "Processing "$file
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a /home/lli/hg19/hg19.chrom.window.n10.bed -b $file -c > $file.chrom.window
done
## coverage at CGIs
less /projects/epigenomics/resources/UCSC_hg19/CGI/CGI.forProfiles.BED | awk '{chr=gensub("chr", "", "g", $1); print chr"\t"$2"\t"$3"\t"$5}' > /home/lli/hg19/CGI.nochr.BED
region="/home/lli/hg19/CGI.nochr.BED"
dirOut="/projects/epigenomics/users/lli/FetalBrain/WGBS/CGI/"
mkdir -p $dirOut
cd /projects/epigenomics/users/lli/FetalBrain/WGBS/
for file in *.sam.bedGraph
do
    echo "Processing "$file
    /gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $file -b $region -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5+$6"\tchr"$7"\t"$8"\t"$9"\t"$10}' > $dirOut/$file.CGI
    less $dirOut/$file.CGI | awk '{cov[$9]=cov[$9]+$5; count[$9]=count[$9]+1; m[$9]=m[$9]+$4*$5; chr[$9]=$6; start[$9]=$7; end[$9]=$8;} END {for(i in chr){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"m[i]/cov[i]"\t"cov[i]/count[i]}}' | sort -k1,1 -k 2,2n -T /projects/epigenomics/users/lli/tmp/ > $dirOut/$file.CGI.bed
    rm $dirOut/$file.CGI
done

