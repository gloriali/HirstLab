#!/bin/sh

JAVA=/gsc/software/linux-x86_64/jre1.7.0_03/bin/java
LIB=/home/mbilenky/bin/Solexa_Java
csizes=/projects/epigenomics/resources/UCSC_chr/hg19.chrom.sizes
dirr=/projects/epigenomics/mbilenky/CpG/hg19/CG_25_around_chr/
dirw=/projects/mbilenky/REMC/brain/MeDIP/wigs/

for name in "HS2788.MeDIP.Brain01.q5.F1028.SET_174" "HS2790.MeDIP.Brain02.q5.F1028.SET_174" "HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174" "HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174" "HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157" "HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166"
do
    echo "$name"
    out="/projects/epigenomics/lli/FetalBrain/MeDIP/CG_25_around_chr/"$name
    mkdir -p $out
    for chr in {1..22} "X" "Y" 
    do
        chr="chr"$chr
        echo "$chr"
        mkdir -p $out/$chr
        $JAVA -jar -Xmx15G $LIB/RegionsCoverageFromWigCalculator.jar -w $dirw/$name.wig.gz -r $dirr/$chr.gz -o $out/$chr -s $csizes -n $name
        less $out/$chr/*.coverage | awk '{gsub("chr", "", $1); print $1"_"$2"\t"$4}' > $out/$chr/$chr"."$name.cov
    done
done

# check file sizes
for chr in {1..22} "X" "Y" 
do
    chr="chr"$chr
    echo "$chr"
    for name in "HS2788.MeDIP.Brain01.q5.F1028.SET_174" "HS2790.MeDIP.Brain02.q5.F1028.SET_174" "HS2775.MeDIP.NeurospheresCortex01.q5.F1028.SET_174" "HS2779.MeDIP.NeurospheresCortex02.q5.F1028.SET_174" "HS2777.MeDIP.NeurospheresGE01.q5.F1028.SET_157" "HS2781.MeDIP.NeurospheresGE02.q5.F1028.SET_166"
    do
        dir="/projects/epigenomics/lli/FetalBrain/MeDIP/CG_25_around_chr/"$name
        wc -l $dir/$chr/$chr"."$name.cov
    done
done
  
