#!/bin/sh

# RPKM of 5mC modifiers
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
dirRPKM=/projects/epigenomics2/users/lli/glioma/RNAseq/
join <(less $dir5mC/DNAme_regulators.txt | sort -k2,2) $dirRPKM/RPKM/glioma.RPKM -1 2 -2 1 | join - $dirRPKM/NPC_RPKM/NPC.RPKM | sed -e 's/ /\t/g' > $dir5mC/DNAme_regulators.RPKM

# combine strands 5mC for CEMT and NPC
dirIn='/projects/epigenomics2/users/lli/glioma/WGBS/'
cd $dirIn
for file in *.5mC.CpG.gz; do
    lib=$(echo $file | sed -e 's/.5mC.CpG.gz//g')
    echo "Combining strand for" $lib
    /home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i $dirIn -o $dirIn -f $file -n $lib -format novo5mC
done

# Formatting TCGA WGBS
dirIn='/projects/epigenomics2/users/lli/glioma/WGBS/'
cd $dirIn
for file in *TCGA*.bed; do
    lib=$(echo $file | sed -e 's/.bed//g')
    echo "Processing" $lib
    less $file | awk 'NR>1{gsub("chr", ""); if($8>=3){printf "%s\t%d\t%d\t%.0f\t%.0f\t%.4f\n", $1, $2, $3+1, $8-$8*$7/100, $8*$7/100, $7/100}}' | sort -k1,1 -k2,2n > $dirIn/$lib.combine.5mC.CpG
done
rm *TCGA*.bed

# check coverage profile and 5mC profile
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirIn='/projects/epigenomics2/users/lli/glioma/WGBS/'
echo -e "sample\tcoverage\tN" > $dirIn/qc_5mC_coverage.txt
echo -e "sample\ttype\tfractional\tN" > $dirIn/qc_5mC_profile.txt 
echo -e "sample\ttype\tmin\tymin\tlower\tmedian\tupper\tymax\tmax" > $dirIn/qc_5mC_quantile.txt #ymin: 10% quantile; ymax: 90% quantile
cd $dirIn
for file in *.combine.5mC.CpG; do
    lib=$(echo $file | sed -e 's/.5mC.CpG//g' | sed 's/.combine//g')
    echo "Processing" $lib
    less $file | awk '{c = $4 + $5; if(c >= 5000){s[5001]++} else {s[c]++}} END{for(i = 1; i <= 5001; i++){print "'$lib'""\t"i"\t"s[i]}}' >> $dirIn/qc_5mC_coverage.txt
    less $file | awk '{s[int($6*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tgenome\t"i/100"\t"s[i]}}' >> $dirIn/qc_5mC_profile.txt 
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | awk '{s[int($1*100)]++} END{for(i = 0; i<=100; i++){print "'$lib'""\tCGI\t"i/100"\t"s[i]}}' >> $dirIn/qc_5mC_profile.txt 
    less $file | awk '{print $6}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tgenome\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc_5mC_quantile.txt
    less $file | awk '{gsub("chr", ""); print $1"\t"$2"\t"$3"\t"$1":"$2"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{t[$10]=t[$10]+$5; c[$10]=c[$10]+$6} END{for(i in c){print c[i]/(c[i]+t[i])}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\tCGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' >> $dirIn/qc_5mC_quantile.txt
done

# clustering
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirIn='/projects/epigenomics2/users/lli/glioma/WGBS/'
cd $dirIn
less /home/lli/hg19/CG.strand | awk '{if(NR%2){print $2}}' | sort > x
less /home/lli/hg19/CGI.forProfiles.BED | awk '{print $4}' | sort > a
for file in *combine.5mC.CpG; do
    sample=$(echo $file | sed 's/.5mC.CpG//g' | sed 's/.combine//g')
    echo $sample
    less $file | awk '{print $1":"$2"-"$3"\t"$6}' | sort -k1,1 | join x - > y
    mv y x
    less $file | $BEDTOOLS/intersectBed -a /home/lli/hg19/CGI.forProfiles.BED -b stdin -wa -wb | awk '{t[$4]=t[$4]+$8; c[$4]=c[$4]+$9} END{for(i in t){if(t[i]+c[i]>0){print i"\t"c[i]/(c[i]+t[i])}}}' | sort -k1,1 | join a - > b
    mv b a
done
mv x matrix_genome.5mC
mv a matrix_CGI.5mC

# methylation profile around CGI edges
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirIn='/projects/epigenomics2/users/lli/glioma/WGBS/'
dirOut=$dirIn/CGI_edge/
mkdir -p $dirOut
cd $dirIn
echo -e "sample\tedge\tdistance\tfractional" > $dirOut/CGI.edge.profile
for file in *.combine.5mC.CpG; do
    name=$(echo $file | sed -e 's/.5mC.CpG//g' | sed 's/.combine//g')
    echo $name
    less /home/lli/hg19/CGI.edges.bed | awk '$5 ~ /L/ {print $0}' | $BEDTOOLS/closestBed -a $dirIn/$file -b stdin -D a | awk '{if(($12>=-1000)&&($12<=2000)&&($1!="Y")){print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$6"\t"$10"\t"$11"\t"$12"\t""'$name'"}}' > $dirOut/$name.CGI.edge.L
    less /home/lli/hg19/CGI.edges.bed | awk '$5 ~ /R/ {print $0}' | $BEDTOOLS/closestBed -a $dirIn/$file -b stdin -D a | awk '{if(($12>=-2000)&&($12<=1000)&&($1!="Y")){print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$6"\t"$10"\t"$11"\t"$12"\t""'$name'"}}' > $dirOut/$name.CGI.edge.R
    less $dirOut/$name.CGI.edge.L | awk '{s[int($8/50)] = s[int($8/50)] + $5; c[int($8/50)]++} END{for(i in c){print "'$name'""\tL\t"i*50"\t"s[i]/c[i]}}' | sort -k3,3n >> $dirOut/CGI.edge.profile
    less $dirOut/$name.CGI.edge.R | awk '{s[int($8/50)] = s[int($8/50)] + $5; c[int($8/50)]++} END{for(i in c){print "'$name'""\tR\t"i*50"\t"s[i]/c[i]}}' | sort -k3,3n >> $dirOut/CGI.edge.profile
    cat $dirOut/$name.CGI.edge.L $dirOut/$name.CGI.edge.R > $dirOut/$name.CGI.edge
    rm $dirOut/$name.CGI.edge.L $dirOut/$name.CGI.edge.R 
done
cd $dirOut
echo -e "sample1\tsample2\tedge\tdistance\tdelta" > $dirOut/CGI.edge.delta.profile
for f1 in IDH*.CGI.edge; do
    for f2 in NPC*.CGI.edge; do
        s1=$(echo $f1 | sed 's/.CGI.edge//g'); s2=$(echo $f2 | sed 's/.CGI.edge//g');
        echo $s1 $s2;
        join <(less $f1 | awk '{print $4"+"$6"+"$7"\t"$5"\t"$8}' | sort -k1,1) <(less $f2 | awk '{print $4"+"$6"+"$7"\t"$5"\t"$8}' | sort -k1,1) | awk -F' ' '{gsub("+", "\t"); print $1"\t"$2"\t"$3"\t"$5"\t"$4-$6"\t""'$s1'""\t""'$s1'""-""'$s2'"}' > $dirOut/$s1.$s2.CGI.edge.delta
        less $dirOut/$s1.$s2.CGI.edge.delta | awk '{s[int($4/50)] = s[int($4/50)] + $5; c[int($4/50)]++; edge[int($4/50)]=$3} END{for(i in c){print "'$s1'""\t""'$s2'""\t"edge[i]"\t"i*50"\t"s[i]/c[i]}}' | sort -k3,3 -k4,4n >> $dirOut/CGI.edge.delta.profile
    done
done

## enhancer vs 5mC
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
gene=/home/lli/hg19/hg19v69_genes.bed
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
RPKM=/projects/epigenomics2/users/lli/glioma/RNAseq/RPKM.long
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/
dirOut=/projects/epigenomics2/users/lli/glioma/WGBS/H3K27ac/
mkdir -p $dirOut
cd $dirIn
for file in *.FindER.bed.gz; do
	sample=$(echo $file | sed 's/A19308.//g' | sed 's/\..*//g')
	echo $sample
	$BEDTOOLS/intersectBed -a $file -b <(less $promoter | awk '{print "chr"$0}') -u -f 0.5 | awk '{print $0"\t"$1":"$2"-"$3"_promoter"}' > $dirIn/$sample.promoter.enhancer
	$BEDTOOLS/intersectBed -a $file -b ./SE/$sample'_Gateway_SuperEnhancers.bed' -u | $BEDTOOLS/intersectBed -a stdin -b $dirIn/$sample.promoter.enhancer -v | awk '{print $0"\t"$1":"$2"-"$3"_super"}' > $dirIn/$sample.SE.enhancer
	$BEDTOOLS/intersectBed -a $file -b $dirIn/$sample.promoter.enhancer -v | $BEDTOOLS/intersectBed -a stdin -b $dirIn/$sample.SE.enhancer -v | awk '{print $0"\t"$1":"$2"-"$3"_regular"}' > $dirIn/$sample.regular.enhancer
	cat $dirIn/$sample.promoter.enhancer $dirIn/$sample.SE.enhancer $dirIn/$sample.regular.enhancer | sort -k1,1 -k2,2n > $dirIn/$sample.enhancer
	$JAVA -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $(ls ../../wig/H3K27ac/$sample.*wig.gz) -r $dirIn/$sample.enhancer -o $dirOut -c $chr -n $sample > $dirOut/$sample.coverage.log
	less $dirOut/$sample.enhancer.$sample.coverage | sed 's/chr//g' | $BEDTOOLS/intersectBed -a stdin -b $(ls $dir5mC/*$sample*.combine.5mC.CpG) -wa -wb | awk '{t[$5]=t[$5]+$11; c[$5]=c[$5]+$12; chr[$5]=$1; start[$5]=$2; end[$5]=$3; signal[$5]=$6} END{for(i in chr){if(t[i]+c[i]>0){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"signal[i]"\t"c[i]/(c[i]+t[i])"\t""'$sample'"}}}' | sort -k1,1 -k2,2n > $dirOut/$sample.enhancer.5mC
done
cat $dirOut/*.enhancer.5mC > $dirOut/enhancer.5mC
less $gene | awk '$4 ~ /protein_coding/ {gsub("_protein_coding", ""); print $0}' | sort -k1,1 -k2,2n | $BEDTOOLS/closestBed -a <(less $dirOut/enhancer.5mC | sort -k1,1 -k2,2n) -b stdin -wa -wb | awk '{gsub("_", "\t", $4); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"."$11}' > $dirOut/enhancer.5mC.gene
less $dirOut/enhancer.5mC.gene | sort -k8,8 | join - <(sort $RPKM -k1,1) -1 8 -2 1 | awk '{gsub("\\.", "\t", $1); print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$1"\t"$9}' > $dirOut/enhancer.5mC.RPKM
$BEDTOOLS/intersectBed -a $dirOut/enhancer.5mC -b /home/lli/hg19/CG.BED -c | awk '{print $0"\t"$8/($3-$2)*1000}' > $dirOut/enhancer.5mC.CpG
$BEDTOOLS/intersectBed -a $dirOut/enhancer.5mC -b /home/lli/hg19/CGI.forProfiles.BED -u > $dirOut/enhancer.5mC.CGI 
$BEDTOOLS/intersectBed -a $dirOut/enhancer.5mC -b /home/lli/hg19/CGI.forProfiles.BED -v > $dirOut/enhancer.5mC.nonCGI 
### change of state
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
dirDE=/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine/
dirHM=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/
dirOut=/projects/epigenomics2/users/lli/glioma/WGBS/H3K27ac/
cd $dirOut
for file in CEMT*.enhancer.5mC.*.*; do
    sample=$(echo $file | cut -d'.' -f1)
    echo $file $sample
    less $file | sed 's/chr//g' | $BEDTOOLS/intersectBed -a stdin -b $dir5mC/NPC_GE04.combine.5mC.CpG -wa -wb | awk '{t[$4]=t[$4]+$12; c[$4]=c[$4]+$13; chr[$4]=$1; start[$4]=$2; end[$4]=$3; signal[$4]=$5; fractional[$4]=$6; sample[$4]=$7; type[$4]=$8} END{for(i in chr){if(t[i]+c[i]>0){f=c[i]/(c[i]+t[i]); if(f>0.7){category="hyper";}else if(f<0.3){category="hypo";}else{category="median";}; print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"signal[i]"\t"fractional[i]"\t"sample[i]"\t"type[i]"\t"f"\t"category}}}' | sort -k1,1 -k2,2n > $dirOut/$file.states
    $BEDTOOLS/closestBed -a $dirOut/$file.states -b <(less $promoter | sort -k1,1 -k2,2n) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$14}' > $dirOut/$file.states.gene
    cat <(less $dirDE/UP."$sample"_NPC.FDR_0.01.rmin_0.005.Nmin_25 | awk '{print $0"\tUP"}') <(less $dirDE/DN."$sample"_NPC.FDR_0.01.rmin_0.005.Nmin_25 | awk '{print $0"\tDN"}') > $dirDE/DE."$sample"_NPC.FDR_0.01.rmin_0.005.Nmin_25
    awk 'NR==FNR{de[$1]=$9; next} {if(!($11 in de))de[$11]="ST"; print $0"\t"de[$11]}' $dirDE/DE."$sample"_NPC.FDR_0.01.rmin_0.005.Nmin_25 $dirOut/$file.states.gene > $dirOut/$file.states
    for mark in H3K27ac H3K4me1 H3K4me3 H3K27me3; do
        less $(ls $dirHM/$mark/"$sample"*.FindER.bed.gz) | sed 's/chr//g' | $BEDTOOLS/intersectBed -a $dirOut/$file.states -b stdin -u | awk '{print $0"\t""'$mark'""_glioma_T"}' > $dirOut/$file.states.T
        less $(ls $dirHM/$mark/"$sample"*.FindER.bed.gz) | sed 's/chr//g' | $BEDTOOLS/intersectBed -a $dirOut/$file.states -b stdin -v | awk '{print $0"\t""'$mark'""_glioma_F"}' > $dirOut/$file.states.F
        cat $dirOut/$file.states.T $dirOut/$file.states.F | sort -k1,1 -k2,2n > $dirOut/$file.states
        less $(ls $dirHM/$mark/*GE04*.FindER.bed.gz) | sed 's/chr//g' | $BEDTOOLS/intersectBed -a $dirOut/$file.states -b stdin -u | awk '{print $0"\t""'$mark'""_NPC_T"}' > $dirOut/$file.states.T
        less $(ls $dirHM/$mark/*GE04*.FindER.bed.gz) | sed 's/chr//g' | $BEDTOOLS/intersectBed -a $dirOut/$file.states -b stdin -v | awk '{print $0"\t""'$mark'""_NPC_F"}' > $dirOut/$file.states.F
        cat $dirOut/$file.states.T $dirOut/$file.states.F | sort -k1,1 -k2,2n > $dirOut/$file.states
        rm $dirOut/$file.states.T $dirOut/$file.states.F
    done
done
cat CEMT_19.enhancer.5mC.promoter.hyper.states CEMT_22.enhancer.5mC.promoter.hyper.states CEMT_47.enhancer.5mC.promoter.hyper.states > IDHmut.enhancer.5mC.promoter.hyper.states
less IDHmut.enhancer.5mC.promoter.hyper.states | awk '$12 !~ /ST/ {if($10!="hyper"){print $0}}' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$9"\t"$11"\t"$13"\t"$14}' > IDHmut.enhancer.5mC.promoter.hyper.states.DE.DM
less IDHmut.enhancer.5mC.promoter.hyper.states.DE.DM | awk '{print $7}' | sort | uniq -c | awk '$1>1{print $2}' | sort | join - /projects/epigenomics2/resources/Ensembl/hg37v69/hg37v69_genes | awk '{print $1"\t"$7"\t"$8}' > IDHmut.enhancer.5mC.promoter.hyper.states.DE.DM.multi.gene
### Homer
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirOut=/projects/epigenomics2/users/lli/glioma/WGBS/H3K27ac/
cd $dirOut
mkdir -p $dirOut/homer/
for file in *.enhancer.5mC; do
    sample=$(echo $file | sed 's/.enhancer.5mC//g')
    echo $sample
    for category in promoter regular super; do
        less $file | awk '{if($6>0.7 && $4~"'$category'"){print "chr"$0"\thyper"}}' > $file.$category.hyper
        less $file | awk '{if($6<0.3 && $4~"'$category'"){print "chr"$0"\thypo"}}' > $file.$category.hypo
        mkdir -p $dirOut/homer/$file.$category.hyper/
        mkdir -p $dirOut/homer/$file.$category.hypo/
        /home/lli/bin/homer/bin/findMotifsGenome.pl $file.$category.hyper hg19 $dirOut/homer/$file.$category.hyper/ -size 200 -len 8
        /home/lli/bin/homer/bin/findMotifsGenome.pl $file.$category.hypo hg19 $dirOut/homer/$file.$category.hypo/ -size 200 -len 8
    done
done
cd $dirOut/homer/
for category in promoter regular super; do
    echo -e "TF\tmotif\tq\tpercent_with_motif\tsample\tgroup\tcategory" > homer.knownResults.summary.$category
    for lib in CEMT_19 CEMT_22 CEMT_23 CEMT_47 NPC_GE04; do
        for dm in hyper hypo; do
            echo $lib $dm $category
            less ./$lib.enhancer.5mC.$category.$dm/knownResults.txt | awk 'NR > 1 {tf=gensub("/.*", "", "g", $1); gsub("%", ""); print tf"\t"$1"\t"$5"\t"$7"\t""'$lib'""\t""'$dm'""\t""'$category'"}' >> homer.knownResults.summary.$category
        done
    done
done
### TCGA validation
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/
dirOut=/projects/epigenomics2/users/lli/glioma/WGBS/H3K27ac/
cd $dirIn
cat $dirIn/CEMT_23.promoter.enhancer > $dirIn/IDHwt.promoter.enhancer
cat $dirIn/CEMT_23.regular.enhancer > $dirIn/IDHwt.regular.enhancer
cat $dirIn/CEMT_23.SE.enhancer > $dirIn/IDHwt.SE.enhancer
cat $dirIn/CEMT_19.promoter.enhancer $dirIn/CEMT_22.promoter.enhancer $dirIn/CEMT_47.promoter.enhancer | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin | awk '{print $0"\t.\t"$1":"$2"-"$3"_promoter"}' > $dirIn/IDHmut.promoter.enhancer
cat $dirIn/CEMT_19.regular.enhancer $dirIn/CEMT_22.regular.enhancer $dirIn/CEMT_47.regular.enhancer | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin | awk '{print $0"\t.\t"$1":"$2"-"$3"_regular"}' > $dirIn/IDHmut.regular.enhancer
cat $dirIn/CEMT_19.SE.enhancer $dirIn/CEMT_22.SE.enhancer $dirIn/CEMT_47.SE.enhancer | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin | awk '{print $0"\t.\t"$1":"$2"-"$3"_super"}' > $dirIn/IDHmut.SE.enhancer
for file in IDHmut*enhancer; do
    echo $file
    for sample in CEMT_19 CEMT_22 CEMT_47; do
        echo $sample
        $JAVA -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $(ls ../../wig/H3K27ac/$sample.*wig.gz) -r $file -o $dirOut -c $chr -n $sample > $dirOut/IDHmut.$sample.coverage.log
    done
    join <(sort -k5,5 $dirOut/$file.CEMT_19.coverage) <(sort -k5,5 $dirOut/$file.CEMT_22.coverage) -1 5 -2 5 | join - <(sort -k5,5 $dirOut/$file.CEMT_47.coverage) -1 1 -2 5 | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"($6+$12+$18)/3}' > $dirOut/$file.coverage
    > $dirOut/$file.coverage.5mC
    for mC in $dir5mC/IDHmut_TCGA*.combine.5mC.CpG; do
        sample2=$(basename $mC | sed 's/.combine.5mC.CpG//g')
        echo $sample2
        less $dirOut/$file.coverage | sed 's/chr//g' | $BEDTOOLS/intersectBed -a stdin -b $mC -wa -wb | awk '{t[$4]=t[$4]+$9; c[$4]=c[$4]+$10; chr[$4]=$1; start[$4]=$2; end[$4]=$3; signal[$4]=$5} END{for(i in chr){if(t[i]+c[i]>0){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"signal[i]"\t"c[i]/(c[i]+t[i])"\t""'$sample2'"}}}' | sort -k1,1 -k2,2n >> $dirOut/$file.coverage.5mC
    done
done
for file in IDHwt*enhancer; do
    echo $file
    sample=CEMT_23
    $JAVA -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $(ls ../../wig/H3K27ac/$sample.*wig.gz) -r $file -o $dirOut -c $chr -n $sample > $dirOut/IDHwt.$sample.coverage.log
    less $dirOut/$file.CEMT_23.coverage | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6}' > $dirOut/$file.coverage
    > $dirOut/$file.coverage.5mC
    for mC in $dir5mC/IDHwt_TCGA*.combine.5mC.CpG; do
        sample2=$(basename $mC | sed 's/.combine.5mC.CpG//g')
        echo $sample2
        less $dirOut/$file.coverage | sed 's/chr//g' | $BEDTOOLS/intersectBed -a stdin -b $mC -wa -wb | awk '{t[$4]=t[$4]+$9; c[$4]=c[$4]+$10; chr[$4]=$1; start[$4]=$2; end[$4]=$3; signal[$4]=$5} END{for(i in chr){if(t[i]+c[i]>0){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"signal[i]"\t"c[i]/(c[i]+t[i])"\t""'$sample2'"}}}' | sort -k1,1 -k2,2n >> $dirOut/$file.coverage.5mC
    done
done
cat $dirOut/IDH*.enhancer.coverage.5mC > $dirOut/TCGA.enhancer.coverage.5mC
$BEDTOOLS/intersectBed -a $dirOut/TCGA.enhancer.coverage.5mC -b /home/lli/hg19/CG.BED -c | awk '{print $0"\t"$8/($3-$2)*1000}' > $dirOut/TCGA.enhancer.coverage.5mC.CpG
### promoter hyper enhancer in other tissues
dirIn=/projects/edcc_new/reference_epigenomes/
dirOut=/projects/epigenomics2/users/lli/glioma/WGBS/H3K27ac/others/
mkdir -p $dirOut
cd $dirIn
for lib in CEMT*; do
    echo $lib
    if [[ -f $(ls ./$lib/bams/Bisulfite-Seq/DNAMethylation/*.5mC.CpG) && -f $(ls ./$lib/bams/ChIP-Seq/H3K27ac/wig/*.wig.gz) && -f $(ls ./$lib/bams/ChIP-Seq/H3K27ac/FindER.1.0.0b/*.FindER.bed.gz) ]]; then
        echo "Files exist"
        ln -s $dirIn/$lib/bams/Bisulfite-Seq/DNAMethylation/*.5mC.CpG $dirOut/$lib.5mC.CpG
        ln -s $dirIn/$lib/bams/ChIP-Seq/H3K27ac/wig/*.wig.gz $dirOut/$lib.H3K27ac.wig.gz
        ln -s $dirIn/$lib/bams/ChIP-Seq/H3K27ac/FindER.1.0.0b/*.FindER.bed.gz $dirOut/$lib.FindER.bed.gz
    fi
done
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
gene=/home/lli/hg19/hg19v69_genes.bed
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
cd $dirOut
for file in *.FindER.bed.gz; do
    sample=$(echo $file | sed 's/.FindER.bed.gz//g')
    echo $sample
    $BEDTOOLS/intersectBed -a $file -b <(less $promoter | awk '{print "chr"$0}') -u -f 0.5 | awk '{print $0"\t"$1":"$2"-"$3"_promoter"}' > $sample.promoter.enhancer
    $JAVA -jar -Xmx15G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $sample.H3K27ac.wig.gz -r $sample.promoter.enhancer -o $dirOut -c $chr -n $sample > $dirOut/$sample.coverage.log
    less $sample.promoter.enhancer.$sample.coverage | sed 's/chr//g' | $BEDTOOLS/intersectBed -a stdin -b <(less $sample.5mC.CpG | awk '{gsub("chr", ""); print $1"\t"$2"\t"$2+1"\t"$4"\t"$5}') -wa -wb | awk '{if($11+$12 >= 3){t[$5]=t[$5]+$11; c[$5]=c[$5]+$12; chr[$5]=$1; start[$5]=$2; end[$5]=$3; signal[$5]=$6}} END{for(i in chr){if(t[i]+c[i]>0){print chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"signal[i]"\t"c[i]/(c[i]+t[i])"\t""'$sample'"}}}' | sort -k1,1 -k2,2n > $dirOut/$sample.promoter.enhancer.5mC
done
cat $dirOut/*.promoter.enhancer.5mC > $dirOut/promoter.enhancer.5mC
$BEDTOOLS/intersectBed -a $dirOut/promoter.enhancer.5mC -b /home/lli/hg19/CG.BED -c | awk '{print $0"\t"$8/($3-$2)*1000}' > $dirOut/promoter.enhancer.5mC.CpG

# DMR between glioma and NPCs: pairwise between each glioma and all 4 NPCs
dirIn='/projects/epigenomics2/users/lli/glioma/WGBS/'
dirOut='/projects/epigenomics2/users/lli/glioma/WGBS/DMR/'
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
mkdir -p $dirOut/intermediate/
echo -e "sample\tp-value\tdelta\tm\ttotal\thyper\thypo" > $dirOut/intermediate/DM.summary.stats
echo -e "sample\tsize\tcut\tmedian_length\tmedian_N_CpG\ttotal\thyper\thypo" > $dirOut/intermediate/DMR.summary.stats
pth=0.0005
delta=0.6
m=0.75
cov=3
size=500  
cut=3
cd $dirIn
for file1 in IDH*.combine.5mC.CpG; do
    lib1=$(echo $file1 | sed -e 's/.5mC.CpG//g' | sed 's/.combine//g')
    echo -e "\n\n"$lib1
    for file2 in NPC*.combine.5mC.CpG; do
        lib2=$(echo $file2 | sed -e 's/.combine.5mC.CpG//g')
        name=$lib1'_'$lib2
        echo -e "\n"$name
        /home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut/intermediate/ -f1 $file1 -f2 $file2 -n $name -p $pth -d $delta -m $m -c $cov
        /home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut/intermediate/ -o $dirOut/intermediate/ -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
    done
done

## For each glioma sample, take the intersect of DMRs compared to all 4 NPCs
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirOut='/projects/epigenomics2/users/lli/glioma/WGBS/DMR/'
size=500  
cut=3
echo -e "sample\tsize\tcut\thyper\thypo\tlength_hyper\tlength_hypo" > $dirOut/DMR.summary.stats
cd /projects/epigenomics2/users/lli/glioma/WGBS/
for file1 in IDH*.combine.5mC.CpG; do
    lib1=$(echo $file1 | sed -e 's/.5mC.CpG//g' | sed 's/.combine//g')
    echo -e $lib1
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC_Cortex02.s$size.c$cut.hyper.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC_Cortex04.s$size.c$cut.hyper.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC_Cortex.hyper.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC_GE02.s$size.c$cut.hyper.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC_GE04.s$size.c$cut.hyper.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC_GE.hyper.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC_Cortex.hyper.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC_GE.hyper.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t1\t"$3-$2}' > $dirOut/DMR.$lib1'_'NPC.hyper
    less $dirOut/DMR.$lib1'_'NPC.hyper | awk '!/GL/ {print $1"\t"$2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n > $dirOut/DMR.$lib1'_'NPC.hyper.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC_Cortex02.s$size.c$cut.hypo.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC_Cortex04.s$size.c$cut.hypo.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC_Cortex.hypo.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC_GE02.s$size.c$cut.hypo.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC_GE04.s$size.c$cut.hypo.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/intermediate/DMR.$lib1'_'NPC_GE.hypo.bed
    $BEDTOOLS/intersectBed -a $dirOut/intermediate/DMR.$lib1'_'NPC_Cortex.hypo.bed -b $dirOut/intermediate/DMR.$lib1'_'NPC_GE.hypo.bed | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t-1\t"$3-$2}' > $dirOut/DMR.$lib1'_'NPC.hypo
    cat $dirOut/DMR.$lib1'_'NPC.hyper $dirOut/DMR.$lib1'_'NPC.hypo > $dirOut/DMR.$lib1'_'NPC
    less $dirOut/DMR.$lib1'_'NPC.hypo | awk '!/GL/ {print $1"\t"$2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n > $dirOut/DMR.$lib1'_'NPC.hypo.bed
    hyper=($(wc -l $dirOut/DMR.$lib1'_'NPC.hyper))
    hypo=($(wc -l $dirOut/DMR.$lib1'_'NPC.hypo))
    length_hyper=($(less $dirOut/DMR.$lib1'_'NPC.hyper | awk '{len=len+$6}END{print len}'))
    length_hypo=($(less $dirOut/DMR.$lib1'_'NPC.hypo | awk '{len=len+$6}END{print len}'))
    echo -e $lib1"_NPCs\t"$size"\t"$cut"\t"$hyper"\t"$hypo"\t"$length_hyper"\t"$length_hypo >> $dirOut/DMR.summary.stats
done

## DMR enrichment in genomic regions
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
> /projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/glioma_enhancer_all.bed
for file in /projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/CEMT*.bed.gz; do
    less $file >> /projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/glioma_enhancer_all.bed
done
less /projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/glioma_enhancer_all.bed | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin > /projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/glioma_enhancer.bed 
dirOut='/projects/epigenomics2/users/lli/glioma/WGBS/DMR/'
/home/lli/HirstLab/Pipeline/shell/DMR.intersect.sh -d $dirOut -r /projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/glioma_enhancer.bed -n enhancer

## DMR distance to closest CGI (midpoint of DMR to midpoint of CGI)
less /home/lli/hg19/CGI.forProfiles.BED | awk '{mid=int(($2+$3)/2); print $1"\t"mid"\t"mid+1"\t"$4}' > /home/lli/hg19/CGI.midpoint.BED
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirOut='/projects/epigenomics2/users/lli/glioma/WGBS/DMR/'
cd $dirOut
mkdir $dirOut/CGI_dis/
for file in DMR.*.bed; do
    name=$(echo $file | sed -e s/'DMR.'//g | sed -e s/'.bed'//g)
    echo "Processing $name"
    less $file | awk '{gsub("chr", ""); mid=int(($2+$3)/2); print $1"\t"mid"\t"mid+1"\t"$4}' | sort -k1,1 -k 2,2n > $dirOut/CGI_dis/$name.tmp.bed
    $BEDTOOLS/closestBed -a $dirOut/CGI_dis/$name.tmp.bed -b /home/lli/hg19/CGI.midpoint.BED -D b -t first > $dirOut/CGI_dis/DMR.$name.CGI.dis.tmp
    awk 'NR==FNR {len[$4]=$3-$2; next} {print $0"\t"$9/(len[$8]/2)}' /home/lli/hg19/CGI.forProfiles.BED $dirOut/CGI_dis/DMR.$name.CGI.dis.tmp > $dirOut/CGI_dis/DMR.$name.CGI.dis
    rm $dirOut/CGI_dis/$name.tmp.bed $dirOut/CGI_dis/DMR.$name.CGI.dis.tmp 
done

## intersect CEMT_21 with IDH mut: noise or priming events?
dirOut='/projects/epigenomics2/users/lli/glioma/WGBS/DMR/'
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
cd $dirOut
$BEDTOOLS/intersectBed -a DMR.CEMT_21_NPC.hyper.bed -b DMR.CEMT_19_NPC.hyper.bed DMR.CEMT_22_NPC.hyper.bed DMR.CEMT_47_NPC.hyper.bed -wa -u | wc -l
wc -l DMR.CEMT_21_NPC.hyper.bed
$BEDTOOLS/intersectBed -a DMR.CEMT_21_NPC.hypo.bed -b DMR.CEMT_19_NPC.hypo.bed DMR.CEMT_22_NPC.hypo.bed DMR.CEMT_47_NPC.hypo.bed -wa -u | wc -l
wc -l DMR.CEMT_21_NPC.hypo.bed
$BEDTOOLS/intersectBed -a DMR.CEMT_23_NPC.hyper.bed -b DMR.CEMT_19_NPC.hyper.bed DMR.CEMT_22_NPC.hyper.bed DMR.CEMT_47_NPC.hyper.bed -wa -u | wc -l
wc -l DMR.CEMT_23_NPC.hyper.bed
$BEDTOOLS/intersectBed -a DMR.CEMT_23_NPC.hypo.bed -b DMR.CEMT_19_NPC.hypo.bed DMR.CEMT_22_NPC.hypo.bed DMR.CEMT_47_NPC.hypo.bed -wa -u | wc -l
wc -l DMR.CEMT_23_NPC.hypo.bed
$BEDTOOLS/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed -b DMR.CEMT_19_NPC.hyper.bed DMR.CEMT_22_NPC.hyper.bed DMR.CEMT_47_NPC.hyper.bed -wa -u | wc -l
wc -l /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hyper.bed
$BEDTOOLS/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed -b DMR.CEMT_19_NPC.hypo.bed DMR.CEMT_22_NPC.hypo.bed DMR.CEMT_47_NPC.hypo.bed -wa -u | wc -l
wc -l /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC02_GE-HuFNSC02.m0.75.p0.005.d0.5.s300.c3.hypo.bed
$BEDTOOLS/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed -b DMR.CEMT_19_NPC.hyper.bed DMR.CEMT_22_NPC.hyper.bed DMR.CEMT_47_NPC.hyper.bed -wa -u | wc -l
wc -l /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hyper.bed
$BEDTOOLS/intersectBed -a /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed -b DMR.CEMT_19_NPC.hypo.bed DMR.CEMT_22_NPC.hypo.bed DMR.CEMT_47_NPC.hypo.bed -wa -u | wc -l
wc -l /projects/epigenomics/users/lli/FetalBrain/WGBS/DMR/DMR.Cortex-HuFNSC04_GE-HuFNSC04.m0.75.p0.005.d0.5.s300.c3.hypo.bed

## % of hyper CpGs in hyper CGIs
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
CGI=/home/lli/hg19/CGI.forProfiles.BED
dirIn=/projects/epigenomics2/users/lli/glioma/WGBS/DMR/intermediate/
dirOut=/projects/epigenomics2/users/lli/glioma/WGBS/DMR/CGI/
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
cov=3
cd $dirIn
echo -e "chr\tstart\tend\tID\tDM\ttotal\tpercent\tglioma\tNPC\tType" > $dirOut/CGI.DMR.hyper.DM.all
for dm in DM.*.bed; do
    samples=$(echo $dm | sed 's/DM.//g' | sed 's/.m0.75.p0.0005.d0.6.bed//g')
    s1=$(echo $samples | sed 's/_NPC.*//g'); s2=$(echo $samples | sed 's/.*_//g')
    dmr=DMR.$samples.s500.c3.hyper.bed
    echo $samples $s1 $s2
    awk 'NR==FNR {id=$1":"$2; if(($4+$5 >= "'$cov'"+0)&&($4+$5 <= 5000)){h[id]=$6}; next} {id=$1":"$2; if((id in h)&&($4+$5 >= "'$cov'"+0)&&($4+$5 <= 5000)){print "chr"$1"\t"$2"\t"$3"\t"id"\t"$6"\t"h[id]}}' $(ls $dir5mC/*$s2.*combine.5mC.CpG) $(ls $dir5mC/*$s1.*combine.5mC.CpG) > $dirIn/$samples.join
    $BEDTOOLS/intersectBed -a <(less $CGI | awk '{print "chr"$0}') -b $dmr -u | $BEDTOOLS/intersectBed -a stdin -b <(less $dm | awk '$4 ~ /1/ {print "chr"$0}') -c | $BEDTOOLS/intersectBed -a stdin -b $dirIn/$samples.join -c | awk '{if($5>0){print $0"\t"$5/$6"\t""'$s1'""\t""'$s2'""\tDM"}}' >> $dirOut/CGI.DMR.hyper.DM.all
    $BEDTOOLS/intersectBed -a <(less $CGI | awk '{print "chr"$0}') -b $dmr -u | $BEDTOOLS/intersectBed -a stdin -b <(less $dirIn/$samples.join | awk '{if($5-$6>=0.5){print $0}}') -c | $BEDTOOLS/intersectBed -a stdin -b $dirIn/$samples.join -c | awk '{if($5>0){print $0"\t"$5/$6"\t""'$s1'""\t""'$s2'""\tdelta>=0.5"}}' >> $dirOut/CGI.DMR.hyper.DM.all
    $BEDTOOLS/intersectBed -a <(less $CGI | awk '{print "chr"$0}') -b $dmr -u | $BEDTOOLS/intersectBed -a stdin -b <(less $dirIn/$samples.join | awk '{if($5-$6>=0.2){print $0}}') -c | $BEDTOOLS/intersectBed -a stdin -b $dirIn/$samples.join -c | awk '{if($5>0){print $0"\t"$5/$6"\t""'$s1'""\t""'$s2'""\tdelta>=0.2"}}' >> $dirOut/CGI.DMR.hyper.DM.all
done
### hyper CGI with K36
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
CGI=/home/lli/hg19/CGI.forProfiles.BED
gene=/home/lli/hg19/hg19v69_genes.bed
dirER=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/
dirIn=/projects/epigenomics2/users/lli/glioma/WGBS/DMR/
dirOut=$dirIn/CGI/
mkdir -p $dirOut
cd $dirIn
echo -e "Sample\thyper\tCGI\tNon-gene_CGI\tH3K36me3\tNon-gene\tp_Fisher" > $dirOut/CGI.DMR.hyper.summary
for file in DMR.*CEMT*.hyper.bed; do
    sample=$(echo $file | sed 's/DMR.//g' | sed 's/_NPC.hyper.bed//g' | sed 's/IDH.*t_//g')
    echo $sample
    $BEDTOOLS/intersectBed -a <(less $CGI | awk '{print "chr"$0}') -b $file -u > $dirOut/CGI.$file
    $BEDTOOLS/intersectBed -a $dirOut/CGI.$file -b $dirER/H3K36me3/$sample.FDR_0.05.FindER.bed.gz -u > $dirOut/CGI.$file.H3K36me3
    $BEDTOOLS/intersectBed -a $dirOut/CGI.$file.H3K36me3 -b <(less $gene | awk '{print "chr"$0}') -v > $dirOut/CGI.$file.H3K36me3.nongene
    Nhyper=$(less $file | wc -l)
    NCGI=$(less $dirOut/CGI.$file | wc -l)
    NCGI_nongene=$($BEDTOOLS/intersectBed -a $dirOut/CGI.$file -b <(less $gene | awk '{print "chr"$0}') -v | wc -l)
    NK36=$(less $dirOut/CGI.$file.H3K36me3 | wc -l)
    Nnongene=$(less $dirOut/CGI.$file.H3K36me3.nongene | wc -l)
    p=$(echo "phyper($Nnongene, $NCGI_nongene, $NCGI - $NCGI_nongene, $NK36, lower.tail = F)" | $R - | sed -e 's/\[1\] //g')
    echo -e "$sample\t$Nhyper\t$NCGI\t$NCGI_nongene\t$NK36\t$Nnongene\t$p" >> $dirOut/CGI.DMR.hyper.summary
done
cat $dirOut/CGI.DMR.IDHmut_CEMT_19_NPC.hyper.bed.H3K36me3.nongene $dirOut/CGI.DMR.IDHmut_CEMT_22_NPC.hyper.bed.H3K36me3.nongene $dirOut/CGI.DMR.IDHmut_CEMT_47_NPC.hyper.bed.H3K36me3.nongene | awk '{print $4}' | uniq -c
### all methylated CGI
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
CGI=/home/lli/hg19/CGI.forProfiles.BED
gene=/home/lli/hg19/hg19v69_genes.bed
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
dirER=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/
dirbam=/projects/edcc_new/reference_epigenomes/
dirOut=/projects/epigenomics2/users/lli/glioma/WGBS/DMR/CGI/
$BEDTOOLS/intersectBed -a $CGI -b $gene -v > /home/lli/hg19/CGI.nongene.bed
echo -e "Sample\tCGI\tCGI_hyper\tCGI_K36\tCGI_hyper_K36\tp" > $dirOut/CGI.hyper.H3K36me3.summary
cd $dir5mC
for file in *CEMT*.combine.5mC.CpG; do
    sample=$(echo $file | sed 's/.5mC.CpG//g' | sed 's/.combine//g');
    echo $sample;
    $BEDTOOLS/intersectBed -a $file -b /home/lli/hg19/CGI.nongene.bed -wa -wb | awk '{t[$10]=t[$10]+$4; c[$10]=c[$10]+$5; chr[$10]=$7; start[$10]=$8; end[$10]=$9} END{for(i in chr){f=c[i]/(c[i]+t[i]); if(f>=0.75){print "chr"chr[i]"\t"start[i]"\t"end[i]"\t"i"\t"f}}}' | sort -k1,1 -k 2,2n > $dirOut/$sample.CGI.nongene.hyper
    $BEDTOOLS/intersectBed -a $dirOut/$sample.CGI.nongene.hyper -b $dirER/H3K36me3/$sample.FDR_0.05.FindER.bed.gz -u > $dirOut/$sample.CGI.nongene.hyper.H3K36me3
    NCGI=$(less /home/lli/hg19/CGI.nongene.bed | wc -l)
    NCGIK36=$($BEDTOOLS/intersectBed -a <(less /home/lli/hg19/CGI.nongene.bed | awk '{print "chr"$0}') -b $dirER/H3K36me3/$sample.FDR_0.05.FindER.bed.gz -u | wc -l)
    NCGIhyper=$(less $dirOut/$sample.CGI.nongene.hyper | wc -l)
    NCGIK36hyper=$(less $dirOut/$sample.CGI.nongene.hyper.H3K36me3 | wc -l)
    p=$(echo "phyper($NCGIK36hyper, $NCGIK36, $NCGI - $NCGIK36, $NCGIhyper, lower.tail = F)" | $R - | sed -e 's/\[1\] //g')
    echo -e "$sample\t$NCGI\t$NCGIhyper\t$NCGIK36\t$NCGIK36hyper\t$p" >> $dirOut/CGI.hyper.H3K36me3.summary
    N=$($samtools view -q 5 -F 1028 $(ls $dirbam/$sample/bams/RNA-Seq/*.bam) | wc -l)
    $BEDTOOLS/coverageBed -a $dirOut/$sample.CGI.nongene.hyper -b <($samtools view -q 5 -F 1028 -b $(ls $dirbam/$sample/bams/RNA-Seq/*.bam)) -counts
done

## DMR enrichment in chromatin states
dirOut='/projects/epigenomics2/users/lli/glioma/WGBS/DMR/intersect/'
dirState='/projects/epigenomics2/users/lli/glioma/ChIPseq/ChromHMM/'
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
CG='/home/lli/hg19/CG.BED'
cd $dirOut
echo -e "Sample\tDM\tState\tTotal_CpG\tDMR_CpG\tState_CpG\tDMR_State_CpG\tEnrichment" > $dirOut/DMR.chromHMM.enrich
for file in *CEMT*.CpG.bed; do
    lib=$(echo $file | sed -e 's/DMR.//g' | sed -e 's/_NPC.*//g' | sed 's/IDH.*t_//g')
    dm=$(echo $file | sed -e 's/.*NPC.//g' | sed -e 's/.CpG.bed//g')
    echo 'Processing' $lib $dm
    less $dirState/$lib'_18_segments.bed' | awk '{gsub("chr", ""); print $0}' | $BEDTOOLS/intersectBed -a $CG -b stdin -wa -wb | awk '{c[$8]=c[$8]+1}END{for(i in c){print i"\t"c[i]}}' > $lib.state
    total=$(less $lib.state | awk '{c=c+$2}END{print c}')
    less $dirState/$lib'_18_segments.bed' | awk '{gsub("chr", ""); print $0}' | $BEDTOOLS/intersectBed -a $file -b stdin -wa -wb | awk '{c[$8]=c[$8]+1}END{for(i in c){print i"\t"c[i]}}' > $lib.state.dmr
    dmr=$(less $lib.state.dmr | awk '{c=c+$2}END{print c}')
    awk 'NR==FNR {dmr[$1]=$2; next} {print $0"\t"dmr[$1]}' $lib.state.dmr $lib.state | sort -k1,1 | awk '{enrich=($3/$2)/("'$dmr'"/"'$total'"); print "'$lib'""\t""'$dm'""\t"$1"\t"'"$total"'"\t"'"$dmr"'"\t"$2"\t"$3"\t"enrich}' >> $dirOut/DMR.chromHMM.enrich
    rm $lib.state $lib.state.dmr
done
awk 'NR==FNR {name[$1]=$2; next} {print $0"\t"name[$3]}' $dirState/chromatin.states.summary $dirOut/DMR.chromHMM.enrich > $dirOut/DMR.chromHMM.enrich.summary
rm $dirOut/DMR.chromHMM.enrich

## DMR enrichemnt in unique histone modification ERs
dirDMR=/projects/epigenomics2/users/lli/glioma/WGBS/DMR/
dirHM=/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
genome=/home/lli/hg19/hg19.chrom.len
cd $dirDMR
echo -e "Sample\tDMR\tMark\tHM\tN_DMR\tN_HM\tN_intersect\tN_total\tFold\tp-value" > $dirDMR/DMR.uniqueHM.enrichment.summary
for file in DMR.*CEMT*.bed; do
    lib=$(echo $file | sed -e 's/DMR.//g' | sed -e 's/_NPC.*//g' | sed 's/IDH.*t_//g')
    dm=$(echo $file | sed -e 's/.bed//g' | sed -e 's/.*NPC.//g')
    echo $lib $dm
    for mark in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K27ac H3K36me3; do
        for file2 in $dirHM/$mark/$lib*.unique; do
            dhm=$(echo $file2 | sed -e 's/.*.vs.NPC_GE04.//g')
            echo $mark $dhm 
            fc=$(less $file2 | sort -k1,1 -k2,2n | $BEDTOOLS/bedtools fisher -a $dirDMR/$file -b stdin -g $genome | awk 'NR==1{gsub(".*: ", ""); a=$1} NR==2{gsub(".*: ", ""); b=$1} NR==3{gsub(".*: ", ""); i=$1} NR==4{gsub(".*: ", ""); t=$1} END{print a"\t"b"\t"i"\t"t"\t"(i/a)/(b/t)}')
            p=$(less $file2 | sort -k1,1 -k2,2n | $BEDTOOLS/bedtools fisher -a $dirDMR/$file -b stdin -g $genome | awk 'NR==5 {gsub("# ", ""); print}' | $R - | sed -e 's/\[1\] //g')
            echo -e "$lib\t$dm\t$mark\t$dhm\t$fc\t$p" >> $dirDMR/DMR.uniqueHM.enrichment.summary
        done
    done
done

## DMR intersect with DE genes
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
dirDE=/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine/
dirIn=/projects/epigenomics2/users/lli/glioma/WGBS/DMR/
dirOut=$dirIn/DE/
mkdir -p $dirOut
echo -e "Sample\tDM\tDE\tN_DM_promoter\tN_DE\tN_intersect\tp_Fisher\tPercent_intersect" > $dirOut/DMR.DE.summary
n_total=19865 
cd $dirIn
for dmr in DMR.*CEMT*.bed; do
    lib=$(echo $dmr | sed -e 's/DMR.//g' | sed -e 's/_NPC.*.bed//g' | sed 's/IDH.*t_//g')
    dm=$(echo $dmr | sed -e 's/.bed//g' | sed -e 's/.*NPC.//g')
    echo $lib $dm
    less $dmr | awk '{gsub("chr", ""); print}' | $BEDTOOLS/intersectBed -a stdin -b $promoter -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | sort -k5,5 > $dirOut/DMR.$lib.NPC.$dm.promoter.bed
    for file in $dirDE/*$lib*.FDR_0.01.rmin_0.005.Nmin_25; do
        de=$(basename $file | sed -e 's/.CEMT.*//g')
        echo $de
        join $dirOut/DMR.$lib.NPC.$dm.promoter.bed $file -1 5 -2 1 | sed 's/ /\t/g' > $dirOut/DMR.$lib.NPC.$dm.$de
        n_dm=$(less $dirOut/DMR.$lib.NPC.$dm.promoter.bed | wc -l)
        n_de=$(less $file | wc -l)
        n_intersect=$(less $dirOut/DMR.$lib.NPC.$dm.$de | wc -l)
        p=$(echo "phyper($n_intersect, $n_dm, $n_total - $n_dm, $n_de, lower.tail = F)" | $R - | sed -e 's/\[1\] //g')
        echo -e "$lib\t$dm\t$de\t$n_dm\t$n_de\t$n_intersect\t$p" | awk '{print $0"\t"$6/$5}' >> $dirOut/DMR.DE.summary
    done
done
cd $dirOut
echo -e "Sample\tDE_DMR\tCGI_DE_DMR" > DMR.DE.CGI.summary
for file in *CEMT*.UP *CEMT*.DN; do
    lib=$(echo $file | sed -e 's/DMR.//g')
    less $file | awk '{print $2"\t"$3"\t"$4"\t"$1}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb > $file.CGI
    echo -e "$lib\t$(less $file | wc -l)\t$(less $file.CGI | wc -l)" >> DMR.DE.CGI.summary
done
dirHM=/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/
echo -e "Sample\tMark\tDE_DMR\tHMgain_DE_DMR\tHMloss_DE_DMR" > DMR.DE.HM.summary
for mark in H3K27ac H3K27me3; do
    for file in *CEMT*.UP *CEMT*.DN; do
        lib=$(echo $file | sed -e 's/DMR.//g' | sed 's/.NPC.*//g')
        echo $mark $lib
        less $file | awk '{print "chr"$2"\t"$3"\t"$4"\t"$1}' | $BEDTOOLS/intersectBed -a stdin -b $dirHM/$mark/$lib.vs.NPC_GE04.$lib.unique -wa -wb > $file.$mark.gain
        less $file | awk '{print "chr"$2"\t"$3"\t"$4"\t"$1}' | $BEDTOOLS/intersectBed -a stdin -b $dirHM/$mark/$lib.vs.NPC_GE04.NPC_GE04.unique -wa -wb > $file.$mark.loss
        echo -e "$file\t$mark\t$(less $file | wc -l)\t$(less $file.$mark.gain | wc -l)\t$(less $file.$mark.loss | wc -l)" >> DMR.DE.HM.summary
    done
done
cat DMR.CEMT_19.NPC.hyper.UP DMR.CEMT_22.NPC.hyper.UP DMR.CEMT_47.NPC.hyper.UP | awk '{print $1"."$11"."$12}' | sort | uniq | sed -e 's/\./\t/g' > DMR.IDHmut.NPC.hyper.UP
cat DMR.CEMT_19.NPC.hyper.UP.H3K27ac.gain DMR.CEMT_22.NPC.hyper.UP.H3K27ac.gain DMR.CEMT_47.NPC.hyper.UP.H3K27ac.gain | awk '{print $4}' | sort | uniq > DMR.IDHmut.NPC.hyper.UP.H3K27ac.gain 
cat DMR.CEMT_19.NPC.hyper.UP.H3K27me3.loss DMR.CEMT_22.NPC.hyper.UP.H3K27me3.loss DMR.CEMT_47.NPC.hyper.UP.H3K27me3.loss | awk '{print $4}' | sort | uniq > DMR.IDHmut.NPC.hyper.UP.H3K27me3.loss
cat DMR.CEMT_19.NPC.hyper.DN DMR.CEMT_22.NPC.hyper.DN DMR.CEMT_47.NPC.hyper.DN | awk '{print $1"."$11"."$12}' | sort | uniq | sed -e 's/\./\t/g' > DMR.IDHmut.NPC.hyper.DN
cat DMR.CEMT_19.NPC.hypo.UP DMR.CEMT_22.NPC.hypo.UP DMR.CEMT_47.NPC.hypo.UP | awk '{print $1"."$11"."$12}' | sort | uniq | sed -e 's/\./\t/g' > DMR.IDHmut.NPC.hypo.UP
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
dirRPKM=/projects/epigenomics2/users/lli/glioma/RNAseq/
dirHM=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/H3K27ac/
echo -e "Sample\thyper_promoter\tUP_2FC\tH3K27ac" > hyper.UP_2FC.H2K27ac.summary
for file in *hyper.promoter.bed; do
    sample=$(echo $file | cut -d'.' -f2)
    join <(less $file | sort -k5,5) <(less $dirRPKM/RPKM/glioma.RPKM | awk 'NR>1{print $1"\t"($2+$4+$6)/3}') -1 5 -2 1 | join - <(less $dirRPKM/NPC_RPKM/NPC.RPKM | awk 'NR>1{print $1"\t"($4+$5+$10+$11)/4}') | awk '{if(($6+1e-6)/($7+1e-6)>=2){print "chr"$2"\t"$3"\t"$4"\t"$5"\t"$1}}' > DMR.$sample.NPC.hyper.UP.2FC
    join <(less DMR.$sample.NPC.hyper.UP.2FC | sort -k5,5) <(less $promoter | sort -k4,4) -1 5 -2 4 | awk '{print "chr"$6"\t"$7"\t"$8"\t"$1"\t"$5}' | $BEDTOOLS/intersectBed -a stdin -b $dirHM/$sample.FDR_0.05.FindER.bed.gz -wa -wb > DMR.$sample.NPC.hyper.UP.2FC.H3K27ac
    echo -e "$sample\t$(less $file | wc -l)\t$(less DMR.$sample.NPC.hyper.UP.2FC | wc -l)\t$(less DMR.$sample.NPC.hyper.UP.2FC.H3K27ac | wc -l)" >> hyper.UP_2FC.H2K27ac.summary
done
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
for file in DMR.CEMT_19.NPC.hyper.UP.2FC.H3K27ac DMR.CEMT_22.NPC.hyper.UP.2FC.H3K27ac DMR.CEMT_47.NPC.hyper.UP.2FC.H3K27ac; do
    name=$(echo $file | sed 's/DMR.//g')
    mkdir -p $dirOut/Homer/$name/
    /home/lli/bin/homer/bin/findMotifsGenome.pl <(less $file | awk '{print $6"\t"$7"\t"$8"\t"$9}') hg19 $dirOut/Homer/$name/ -size 200 -len 8 
done

