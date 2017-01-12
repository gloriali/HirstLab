#!/bin/sh

# Upgrade GE04 to FindER 1.0.0b PET
dirOut=/projects/epigenomics2/users/lli/glioma/ChIPseq/bam/
/gsc/software/linux-x86_64-centos6/spec-1.3.2/spec2bam --threads 2 --in /home/aldente/private/Projects/REMC/IX0816/AnalyzedData/IX0816-1./Solexa/Data/current/BaseCalls/BWA_2012-09-09/compressed_bams/s_4_GATCAG_paired.dup.sorted.bam.spec --ref /projects/sbs_archive2/spec_ref/9606/hg19/1000genomes/GRCh37-lite.fa.spec.ref --out $dirOut/H3K27ac/A19308.NPC_GE04.bam
/gsc/software/linux-x86_64-centos6/spec-1.3.2/spec2bam --threads 2 --in /home/aldente/private/Projects/REMC/IX0817/AnalyzedData/IX0817-1./Solexa/Data/current/BaseCalls/BWA_2012-09-09/compressed_bams/s_5_CTTGTA_paired.dup.sorted.bam.spec --ref /projects/sbs_archive2/spec_ref/9606/hg19/1000genomes/GRCh37-lite.fa.spec.ref --out $dirOut/H3K27me3/A19306.NPC_GE04.bam
/gsc/software/linux-x86_64-centos6/spec-1.3.2/spec2bam --threads 2 --in /home/aldente/private/Projects/REMC/IX0816/AnalyzedData/IX0816-1./Solexa/Data/current/BaseCalls/BWA_2012-09-09/compressed_bams/s_4_ACTTGA_paired.dup.sorted.bam.spec --ref /projects/sbs_archive2/spec_ref/9606/hg19/1000genomes/GRCh37-lite.fa.spec.ref --out $dirOut/H3K36me3/A19307.NPC_GE04.bam
/gsc/software/linux-x86_64-centos6/spec-1.3.2/spec2bam --threads 2 --in /home/aldente/private/Projects/REMC/IX0817/AnalyzedData/IX0817-1./Solexa/Data/current/BaseCalls/BWA_2012-09-09/compressed_bams/s_5_GATCAG_paired.dup.sorted.bam.spec --ref /projects/sbs_archive2/spec_ref/9606/hg19/1000genomes/GRCh37-lite.fa.spec.ref --out $dirOut/H3K4me1/A19303.NPC_GE04.bam
/gsc/software/linux-x86_64-centos6/spec-1.3.2/spec2bam --threads 2 --in /home/aldente/private/Projects/REMC/IX0817/AnalyzedData/IX0817-1./Solexa/Data/current/BaseCalls/BWA_2012-09-09/compressed_bams/s_5_TAGCTT_paired.dup.sorted.bam.spec --ref /projects/sbs_archive2/spec_ref/9606/hg19/1000genomes/GRCh37-lite.fa.spec.ref --out $dirOut/H3K4me3/A19304.NPC_GE04.bam
/gsc/software/linux-x86_64-centos6/spec-1.3.2/spec2bam --threads 2 --in /home/aldente/private/Projects/REMC/IX0817/AnalyzedData/IX0817-1./Solexa/Data/current/BaseCalls/BWA_2012-09-09/compressed_bams/s_5_GGCTAC_paired.dup.sorted.bam.spec --ref /projects/sbs_archive2/spec_ref/9606/hg19/1000genomes/GRCh37-lite.fa.spec.ref --out $dirOut/H3K9me3/A19305.NPC_GE04.bam
/gsc/software/linux-x86_64-centos6/spec-1.3.2/spec2bam --threads 2 --in /home/aldente/private/Projects/REMC/IX0816/AnalyzedData/IX0816-1./Solexa/Data/current/BaseCalls/BWA_2012-09-09/compressed_bams/s_4_TAGCTT_paired.dup.sorted.bam.spec --ref /projects/sbs_archive2/spec_ref/9606/hg19/1000genomes/GRCh37-lite.fa.spec.ref --out $dirOut/input/A19309.NPC_GE04.bam
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
chr=/projects/epigenomics/resources/UCSC_chr/hg19.bwa2ucsc.names
dirWig=/projects/epigenomics2/users/lli/glioma/ChIPseq/wig/
cd $dirOut
for bam in */*GE04*.bam; do
	name=$(basename $bam | sed -e 's/.bam//g');
	lib=$(echo $name | sed 's/A.*\.//g');
	mark=$(echo $bam | sed 's/\/A.*GE04.bam//g');
	echo "Processing $mark $name $lib";
	$samtools index $bam;
    $samtools flagstat $bam > $dirOut/$mark/$name.flagstat
    $bamstats -g 2864785220 -q 10 -b $bam  > $dirOut/$mark/$name.bamstats
	/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirWig/$mark/ -F:1028,-q:5,-n:$lib,-chr:$chr,-cp
done
JAVA=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/bam/
dirOut=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/
cd $dirIn
for file in H*/*GE04.bam; do
    mark=$(echo $file | sed 's/\/A.*GE04.bam//g')
	file=$(basename $file)
    echo $mark $file
    mkdir -p $dirOut/$mark
    $JAVA -jar -Xmx12G /home/mbilenky/bin/Solexa_Java/FindER.1.0.0b.jar -signalBam $dirIn/$mark/$file -inputBam $dirIn/input/A19309.NPC_GE04.bam -out $dirOut/$mark > $dirOut/$mark/$file.log
done

# Convert pair-end to single-end
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/bam/
dirOut=/projects/epigenomics2/users/lli/glioma/ChIPseq/bam_se/
samtools=/gsc/software/linux-x86_64/samtools-0.1.13/samtools
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
for mark in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3 H3K27ac input; do
	mkdir -p $dirOut/$mark/;
	cd $dirIn/$mark/;
	 
done

# FindER 1.0.0b results
dirIn=/projects/edcc_new/reference_epigenomes/
dirOut=/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/
for lib in 19 21 22 23 47; do
    for mark in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3 H3K27ac; do
        echo $lib $mark
        mkdir -p $dirOut/$mark/
        ln -s $dirIn/CEMT_$lib/bams/ChIP-Seq/$mark/FindER.1.0.0b/*.FDR_0.05.FindER.bed.gz $dirOut/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz
    done
done

echo -e "Mark\tSample\tN_region\tTotal_length" > $dirOut/ER.summary.stats
for mark in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3 H3K27ac; do
    cd $dirOut/$mark/
    for file in *.bed.gz; do
        lib=$(echo $file | sed -e 's/.FDR_0.05.FindER.bed.gz//g' | sed -e 's/.*A19309.//g')
        echo -e $mark"\t"$lib"\t"$(less $file | wc -l)"\t"$(less $file | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER.summary.stats
    done
done

########################################################

# Chromatin states summarize
dirOut=/projects/epigenomics2/users/lli/glioma/ChIPseq/ChromHMM/
cd $dirOut
echo -e "E1\t1_EnhA1
E2\t2_TssA
E3\t3_TssFlnk1
E4\t4_EnhA2
E5\t5_EnhG1
E6\t6_EnhG2
E7\t7_EnhWk
E8\t8_TssFlnkD
E9\t9_TssFlnk2
E10\t10_TssBiv
E11\t11_ReprPC
E12\t12_EnhBiv
E13\t13_Repr
E14\t14_Het
E15\t15_Znf_Rpts
E16\t16_Tx
E17\t17_EnhG3
E18\t18_Quies" | sort -k1,1 | sed '1s/^/State\tName\n/' > $dirOut/chromatin.states.summary

for file in *segments.bed; do
    lib=$(echo $file | sed -e 's/_18_segments.bed//g')
    echo $lib
    less $dirOut/chromatin.states.summary > x
    less $file | awk '{len[$4]=len[$4]+$3-$2}END{for(i in len){print i"\t"len[i]}}' | sort -k1,1 | sed '1s/^/State\t'$lib'\n/' > y
    join -t $'\t' x y > $dirOut/chromatin.states.summary
    rm x y
done

########################################################

# compensate for differences in sequencing depth
## subsampling deeper sequenced library to average depth of all other libraries: iterate 50 times
JAVA=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
dirIn='/projects/epigenomics2/users/lli/glioma/ChIPseq/bam/'
dirOut='/projects/epigenomics2/users/lli/glioma/ChIPseq/FindER/subsample/'
mkdir -p $dirOut
cd $dirIn
for file in ./*/*.bam; do
    if [ ! -e "$file.bai" ]; then
        echo "Indexing" $file
        /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools index $file
    fi
done

cd $dirOut
echo -e "H3K27ac\tCEMT_47\tA37110\t0.3
H3K27me3\tCEMT_47\tA37109\t0.4
H3K4me1\tCEMT_21\tA33467\t0.4
H3K4me1\tCEMT_47\tA37106\t0.4
H3K9me3\tCEMT_47\tA37108\t0.5" | sort -k3,3 > $dirOut/library.adjust
echo -e "library\tmark\tsample\tfrac\tn_adjust\tlen_adjust\tn_original\tlen_original\tn_intersect\tlen_intersect\ti" > $dirOut/ER.summary.adjust

for i in {1..50}; do
    echo $i
    $JAVA -jar -Xmx12G /home/mbilenky/bin/Solexa_Java/FindER.1.0.0b.jar -signalBam $dirIn/H3K27ac/A37110.CEMT_47.bam -inputBam $dirIn/input/A37112.CEMT_47.bam -out $dirOut -fractionS 0.3 > $dirOut/H3K27ac_CEMT_47_f.3.out
    $JAVA -jar -Xmx12G /home/mbilenky/bin/Solexa_Java/FindER.1.0.0b.jar -signalBam $dirIn/H3K4me1/A37106.CEMT_47.bam -inputBam $dirIn/input/A37112.CEMT_47.bam -out $dirOut -fractionS 0.4 > $dirOut/H3K4me1_CEMT_47_f.4.out
    $JAVA -jar -Xmx12G /home/mbilenky/bin/Solexa_Java/FindER.1.0.0b.jar -signalBam $dirIn/H3K27me3/A37109.CEMT_47.bam -inputBam $dirIn/input/A37112.CEMT_47.bam -out $dirOut -fractionS 0.4 > $dirOut/H3K27me3_CEMT_47_f.4.out
    $JAVA -jar -Xmx12G /home/mbilenky/bin/Solexa_Java/FindER.1.0.0b.jar -signalBam $dirIn/H3K9me3/A37108.CEMT_47.bam -inputBam $dirIn/input/A37112.CEMT_47.bam -out $dirOut -fractionS 0.5 > $dirOut/H3K9me3_CEMT_47_f.5.out
    $JAVA -jar -Xmx12G /home/mbilenky/bin/Solexa_Java/FindER.1.0.0b.jar -signalBam $dirIn/H3K4me1/A33467.CEMT_21.bam -inputBam $dirIn/input/A33473.CEMT_21.bam -out $dirOut -fractionS 0.4 > $dirOut/H3K4me1_CEMT_21_f.4.out

    for file in *.bed.gz; do
        name=$(echo $file | sed -e 's/\..*//g')
        mark=$(less $dirOut/library.adjust | awk '$3 ~ /'"$name"'/ {print $1}')
        sample=$(less $dirOut/library.adjust | awk '$3 ~ /'"$name"'/ {print $2}')
        frac=$(less $dirOut/library.adjust | awk '$3 ~ /'"$name"'/ {print $4}')
        echo $name $mark $sample $frac
        n=$(less $file | wc -l)
        len=$(less $file | awk '{s=s+$3-$2}END{print s}')
        n_origin=$(less $dirOut/../$mark/$sample.*.bed.gz | wc -l)
        len_origin=$(less $dirOut/../$mark/$sample.*.bed.gz | awk '{s=s+$3-$2}END{print s}')
        n_intersect=$($BEDTOOLS/intersectBed -a $file -b $dirOut/../$mark/$sample.*.bed.gz | wc -l)
        len_intersect=$($BEDTOOLS/intersectBed -a $file -b $dirOut/../$mark/$sample.*.bed.gz | awk '{s=s+$3-$2}END{print s}')
        echo -e $name"\t"$mark"\t"$sample"\t"$frac"\t"$n"\t"$len"\t"$n_origin"\t"$len_origin"\t"$n_intersect"\t"$len_intersect"\t"$i >> $dirOut/ER.summary.adjust
    done
done

########################################################

# Differentially marked regions
## signal files
dirIn=/projects/edcc_new/reference_epigenomes/
dirOut=/projects/epigenomics2/users/lli/glioma/ChIPseq/wig/
for lib in 19 21 22 23 47; do
    for mark in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3 H3K27ac Input; do
        echo $lib $mark
        mkdir -p $dirOut/$mark/
        ln -s $dirIn/CEMT_$lib/bams/ChIP-Seq/$mark/wig/*.wig.gz $dirOut/$mark/CEMT_$lib.wig.gz
    done
done

## differentially marked regions from FindER peaks and wig singal files
### test: how to set the background coverage cutoff?
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/
dirOut=$dirIn/unique/test/
mkdir -p $dirOut
lib=19
mark=H3K27ac
region=$dirIn/FindER/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz
wig=$dirIn/wig/$mark/NPC_GE04.q5.F1028.wig.gz
excl=$dirIn/FindER/$mark/A19308.H3K27ac.GE04.vs.A19309.Input.GE04.FDR_0.05.FindER.bed.gz
name=CEMT_$lib.vs.NPC_GE04.CEMT_$lib
JAVA=/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java
RegCov=/home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
#### calculate region signal from wig
$JAVA -jar -Xmx15G $RegCov -w $wig -r <(less $region | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}') -o $dirOut -c $chr -n $name.signal > $dirOut/$name.log
for i in {1..10}; do
    #### generate background regions
    $BEDTOOLS/shuffleBed -i $region -g $chr -excl $excl | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"'"$i"'}' > $dirOut/$name.background.$i.bed
    #### calculate coverage for background regions
    $JAVA -jar -Xmx15G $RegCov -w $wig -r $dirOut/$name.background.$i.bed -o $dirOut -c $chr -n $name.background.$i >> $dirOut/$name.log
done

### glioma vs NPC
#### method1: ER and signal files
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/
echo -e "Sample\tMark\tN_glioma\tlen_glioma\tN_NPC\tlen_NPC\tN_glioma_unique\tlen_glioma_unique\tN_NPC_unique\tlen_NPC_unique" > $dirIn/unique/ER.unique.summary
for lib in 19 21 22 23 47; do
    for mark in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3 H3K27ac; do
        echo "CEMT_"$lib $mark
        echo -e "\n\nCEMT_"$lib $mark >> $dirIn/unique/ER.unique.log
        dirOut=$dirIn/unique/$mark/
        mkdir -p $dirOut/
        /home/lli/HirstLab/Pipeline/shell/ER.unique.sh -r $dirIn/FindER/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz -w $dirIn/wig/$mark/NPC_GE04.q5.F1028.wig.gz -excl $(ls $dirIn/FindER/$mark/*GE04.*.FDR_0.05.FindER.bed.gz) -o $dirOut/ -n CEMT_$lib.vs.NPC_GE04.CEMT_$lib >> $dirIn/unique/ER.unique.log
        /home/lli/HirstLab/Pipeline/shell/ER.unique.sh -excl $dirIn/FindER/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz -w $dirIn/wig/$mark/CEMT_$lib.wig.gz -r $(ls $dirIn/FindER/$mark/*GE04.*.FDR_0.05.FindER.bed.gz) -o $dirOut/ -n CEMT_$lib.vs.NPC_GE04.NPC_GE04 >> $dirIn/unique/ER.unique.log
        N_glioma=$(less $dirIn/FindER/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz | wc -l)
        len_glioma=$(less $dirIn/FindER/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
        N_NPC=$(less $dirIn/FindER/$mark/*GE04.*.FDR_0.05.FindER.bed.gz | wc -l)
        len_NPC=$(less $dirIn/FindER/$mark/*GE04.*.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
        N_glioma_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.CEMT_$lib.unique | wc -l)
        len_glioma_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.CEMT_$lib.unique | awk '{s=s+$3-$2}END{print s}')
        N_NPC_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.NPC_GE04.unique | wc -l)
        len_NPC_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.NPC_GE04.unique | awk '{s=s+$3-$2}END{print s}')
        echo -e "$lib\t$mark\t$N_glioma\t$len_glioma\t$N_NPC\t$len_NPC\t$N_glioma_unique\t$len_glioma_unique\t$N_NPC_unique\t$len_NPC_unique" >> $dirIn/unique/ER.unique.summary
    done
done
less ER.unique.log | awk '$1 ~ /^C/ {print $0}' ORS=' ' | sed -e 's/CEMT/\nCEMT/g' | sed -e 's/Coverage cutoff: //g' | sed -e 's/ /\t/g' > $dirIn/unique/ER.unique.cutoff
##### outgroup for H3K36me3
mark=H3K36me3;
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/
dirOut=$dirIn/unique/$mark/
for lib in {4..9} {25..38} {40..45} {50..68} ; do
    echo "CEMT_"$lib $mark
    echo -e "\n\nCEMT_"$lib $mark >> $dirIn/unique/ER.unique.log
    /home/lli/HirstLab/Pipeline/shell/ER.unique.sh -r $(ls /projects/edcc_new/reference_epigenomes/CEMT_$lib/bams/ChIP-Seq/H3K36me3/FindER.1.0.0b/*.FDR_0.05.FindER.bed.gz) -w $dirIn/wig/$mark/NPC_GE04.q5.F1028.wig.gz -excl $(ls $dirIn/FindER/$mark/*GE04.*.FDR_0.05.FindER.bed.gz) -o $dirOut/ -n CEMT_$lib.vs.NPC_GE04.CEMT_$lib >> $dirIn/unique/ER.unique.log
    /home/lli/HirstLab/Pipeline/shell/ER.unique.sh -excl $(ls /projects/edcc_new/reference_epigenomes/CEMT_$lib/bams/ChIP-Seq/H3K36me3/FindER.1.0.0b/*.FDR_0.05.FindER.bed.gz) -w $(ls /projects/edcc_new/reference_epigenomes/CEMT_$lib/bams/ChIP-Seq/H3K36me3/wig/*.wig.gz) -r $(ls $dirIn/FindER/$mark/*GE04.*.FDR_0.05.FindER.bed.gz) -o $dirOut/ -n CEMT_$lib.vs.NPC_GE04.NPC_GE04 >> $dirIn/unique/ER.unique.log
    N_glioma=$(less /projects/edcc_new/reference_epigenomes/CEMT_$lib/bams/ChIP-Seq/H3K36me3/FindER.1.0.0b/*.FDR_0.05.FindER.bed.gz | wc -l)
    len_glioma=$(less /projects/edcc_new/reference_epigenomes/CEMT_$lib/bams/ChIP-Seq/H3K36me3/FindER.1.0.0b/*.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
    N_NPC=$(less $dirIn/FindER/$mark/*GE04.*.FDR_0.05.FindER.bed.gz | wc -l)
    len_NPC=$(less $dirIn/FindER/$mark/*GE04.*.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
    N_glioma_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.CEMT_$lib.unique | wc -l)
    len_glioma_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.CEMT_$lib.unique | awk '{s=s+$3-$2}END{print s}')
    N_NPC_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.NPC_GE04.unique | wc -l)
    len_NPC_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.NPC_GE04.unique | awk '{s=s+$3-$2}END{print s}')
    echo -e "$lib\t$mark\t$N_glioma\t$len_glioma\t$N_NPC\t$len_NPC\t$N_glioma_unique\t$len_glioma_unique\t$N_NPC_unique\t$len_NPC_unique" >> $dirIn/unique/ER.unique.summary
done

#### method2: non-overlapping ER
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
echo -e "Sample\tMark\tN_glioma\tlen_glioma\tN_NPC\tlen_NPC\tN_glioma_unique\tlen_glioma_unique\tN_NPC_unique\tlen_NPC_unique" > $dirIn/unique2/ER.unique.summary
echo -e "Sample\tMark\tN_glioma_unique_method1\tN_glioma_unique_method2\tN_glioma_unique_intersect\tN_NPC_unique_method1\tN_NPC_unique_method2\tN_NPC_unique_intersect" > $dirIn/unique2/ER.unique.compare.methods
for lib in 19 21 22 23 47; do
    for mark in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3 H3K27ac; do
        echo "CEMT_"$lib $mark
        dirOut=$dirIn/unique2/$mark/
        mkdir -p $dirOut/
        $BEDTOOLS/intersectBed -a $dirIn/FindER/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz -b $(ls $dirIn/FindER/$mark/*$mark.GE04.*.FDR_0.05.FindER.bed.gz) -v > $dirOut/CEMT_$lib.vs.NPC_GE04.CEMT_$lib.unique
        $BEDTOOLS/intersectBed -a $(ls $dirIn/FindER/$mark/*$mark.GE04.*.FDR_0.05.FindER.bed.gz) -b $dirIn/FindER/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz -v > $dirOut/CEMT_$lib.vs.NPC_GE04.NPC_GE04.unique
        N_glioma=$(less $dirIn/FindER/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz | wc -l)
        len_glioma=$(less $dirIn/FindER/$mark/CEMT_$lib.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
        N_NPC=$(less $dirIn/FindER/$mark/*$mark.GE04.*.FDR_0.05.FindER.bed.gz | wc -l)
        len_NPC=$(less $dirIn/FindER/$mark/*$mark.GE04.*.FDR_0.05.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
        N_glioma_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.CEMT_$lib.unique | wc -l)
        len_glioma_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.CEMT_$lib.unique | awk '{s=s+$3-$2}END{print s}')
        N_NPC_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.NPC_GE04.unique | wc -l)
        len_NPC_unique=$(less $dirOut/CEMT_$lib.vs.NPC_GE04.NPC_GE04.unique | awk '{s=s+$3-$2}END{print s}')
        echo -e "$lib\t$mark\t$N_glioma\t$len_glioma\t$N_NPC\t$len_NPC\t$N_glioma_unique\t$len_glioma_unique\t$N_NPC_unique\t$len_NPC_unique" >> $dirIn/unique2/ER.unique.summary
        echo -e "$lib\t$mark\t$(less $dirIn/unique/$mark/CEMT_$lib.vs.NPC_GE04.CEMT_$lib.unique | wc -l)\t$N_glioma_unique\t$($BEDTOOLS/intersectBed -a $dirIn/unique/$mark/CEMT_$lib.vs.NPC_GE04.CEMT_$lib.unique -b $dirOut/CEMT_$lib.vs.NPC_GE04.CEMT_$lib.unique | wc -l)\t$(less $dirIn/unique/$mark/CEMT_$lib.vs.NPC_GE04.NPC_GE04.unique | wc -l)\t$N_NPC_unique\t$($BEDTOOLS/intersectBed -a $dirIn/unique/$mark/CEMT_$lib.vs.NPC_GE04.NPC_GE04.unique -b $dirOut/CEMT_$lib.vs.NPC_GE04.NPC_GE04.unique | wc -l)" >> $dirIn/unique2/ER.unique.compare.methods
    done
done

## intersect with DE genes
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
dirDE=/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine/
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/
echo -e "Sample\tMark\tMarked\tDE\tN_DM_promoter\tN_DE\tN_intersect\tp_Fisher\tPercent_intersect" > $dirIn/DHM.DE.summary
n_total=19865
for mark in H3K4me1 H3K4me3 H3K27me3 H3K27ac; do
    cd $dirIn/$mark/
    dirOut=$dirIn/$mark/DE/
    mkdir -p $dirOut
    for dmr in CEMT*.unique; do
        lib=$(echo $dmr | sed -e 's/.vs.*//g')
        dm=$(echo $dmr | cut -d'.' -f 4)
        echo $mark $lib $dm
        less $dmr | awk '{gsub("chr", ""); print}' | $BEDTOOLS/intersectBed -a stdin -b $promoter -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$10}' | sort -k5,5 > $dirOut/$lib.vs.NPC.$dm.unique.promoter.bed
        for file in $dirDE/*$lib*.FDR_0.01.rmin_0.005.Nmin_25; do
            de=$(basename $file | sed -e 's/.CEMT.*//g')
            echo $de
            join $dirOut/$lib.vs.NPC.$dm.unique.promoter.bed $file -1 5 -2 1 | sed 's/ /\t/g' > $dirOut/$lib.vs.NPC.$dm.unique.$de
            n_dm=$(less $dirOut/$lib.vs.NPC.$dm.unique.promoter.bed | wc -l)
            n_de=$(less $file | wc -l)
            n_intersect=$(less $dirOut/$lib.vs.NPC.$dm.unique.$de | wc -l)
            p=$(echo "phyper($n_intersect, $n_dm, $n_total - $n_dm, $n_de, lower.tail = F)" | $R - | sed -e 's/\[1\] //g')
            echo -e "$lib\t$mark\t$dm\t$de\t$n_dm\t$n_de\t$n_intersect\t$p" | awk '{print $0"\t"$7/$6}' >> $dirIn/DHM.DE.summary
        done
    done
done

########################################################

## unique enhancers
### Homer for TFBS motifs
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn='/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/'
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
for mark in H3K27ac H3K4me1; do
    mkdir -p $dirIn/$mark/homer/
    cd $dirIn/$mark/
    $BEDTOOLS/intersectBed -a CEMT_19.vs.NPC_GE04.CEMT_19.unique -b CEMT_22.vs.NPC_GE04.CEMT_22.unique | $BEDTOOLS/intersectBed -a stdin -b CEMT_47.vs.NPC_GE04.CEMT_47.unique | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > IDHmut_NPC.IDHmut.unique
    $BEDTOOLS/intersectBed -a CEMT_19.vs.NPC_GE04.NPC_GE04.unique -b CEMT_22.vs.NPC_GE04.NPC_GE04.unique | $BEDTOOLS/intersectBed -a stdin -b CEMT_47.vs.NPC_GE04.NPC_GE04.unique | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > IDHmut_NPC.NPC.unique
    less CEMT_23.vs.NPC_GE04.CEMT_23.unique | awk '{print $1"\t"$2"\t"$3"\t"$4}' > IDHwt_NPC.IDHwt.unique
    less CEMT_23.vs.NPC_GE04.NPC_GE04.unique | awk '{print $1"\t"$2"\t"$3"\t"$4}' > IDHwt_NPC.NPC.unique
    for file in *.unique; do
        name=$(echo $file | sed -e 's/.unique//g')
        echo "Processing "$mark $name
        dirOut=$dirIn/$mark/homer/$name/
        mkdir -p $dirOut
        /home/lli/bin/homer/bin/findMotifsGenome.pl $file hg19 $dirOut -size 200 -len 8
    done
done
dirRPKM=/projects/epigenomics2/users/lli/glioma/RNAseq/
for mark in H3K27ac H3K4me1; do
    cd $dirIn/$mark/homer/
    less IDHmut_NPC.IDHmut/knownResults.txt | awk '{if(NR==1){print}else{percent=gensub("%", "", "g", $7);if($5<=0.01&&percent>=20){print}}}' > IDHmut_NPC.IDHmut.tf
    less IDHmut_NPC.NPC/knownResults.txt | awk '{if(NR==1){print}else{percent=gensub("%", "", "g", $7);if($5<=0.01&&percent>=20){print}}}' > IDHmut_NPC.NPC.tf
    less IDHmut_NPC.IDHmut.tf | awk 'NR > 1 {print $1}' | sed 's/(.*//g' > IDHmut_NPC.IDHmut.tf.ID
    less IDHmut_NPC.NPC.tf | awk 'NR > 1 {print $1}' | sed 's/(.*//g' > IDHmut_NPC.NPC.tf.ID
done
for mark in H3K27ac H3K4me1; do
    cd $dirIn/$mark/homer/
    less IDHmut_NPC.IDHmut.tf.ID | sort -k2,2 | join - $dirRPKM/NPC_RPKM/NPC.RPKM -1 2 -2 1 | join - $dirRPKM/RPKM/glioma.RPKM | sed 's/ /\t/g' > IDHmut_NPC.IDHmut.tf.RPKM
    less IDHmut_NPC.NPC.tf.ID | sort -k2,2 | join - $dirRPKM/NPC_RPKM/NPC.RPKM -1 2 -2 1 | join - $dirRPKM/RPKM/glioma.RPKM | sed 's/ /\t/g' > IDHmut_NPC.NPC.tf.RPKM
done
#### Ascl1, Olig2 , and HEB in H3K4me1
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn='/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/H3K4me1/'
dirDE='/projects/epigenomics2/users/lli/glioma/RNAseq/DEfine/'
/home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/IDHmut_NPC.IDHmut.unique hg19 $dirIn/homer/ -find $dirIn/homer/IDHmut_NPC.IDHmut/knownResults/known7.motif | awk 'NR>1 {print $1}' | uniq | sed 's/:/\t/g' | sed 's/-/\t/g'  | awk '{print $0"\t"$1":"$2"-"$3}' > $dirIn/homer/IDHmut_NPC.IDHmut.Ascl1.annotate
/home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/IDHmut_NPC.IDHmut.unique hg19 $dirIn/homer/ -find $dirIn/homer/IDHmut_NPC.IDHmut/knownResults/known11.motif | awk 'NR>1 {print $1}' | uniq | sed 's/:/\t/g' | sed 's/-/\t/g'  | awk '{print $0"\t"$1":"$2"-"$3}' > $dirIn/homer/IDHmut_NPC.IDHmut.Olig2.annotate
/home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/IDHmut_NPC.IDHmut.unique hg19 $dirIn/homer/ -find $dirIn/homer/IDHmut_NPC.IDHmut/knownResults/known16.motif | awk 'NR>1 {print $1}' | uniq | sed 's/:/\t/g' | sed 's/-/\t/g'  | awk '{print $0"\t"$1":"$2"-"$3}' > $dirIn/homer/IDHmut_NPC.IDHmut.HEB.annotate
cd $dirIn/homer/
for file in *.annotate; do
    tf=$(echo $file | sed 's/.annotate//g' | sed 's/IDHmut_NPC.IDHmut.//g')
    echo $tf
    less /home/lli/hg19/hg19v69_genes.bed | awk '$4 ~ /protein_coding/ {gsub("_protein_coding", ""); print "chr"$0}' | sort -k1,1 -k2,2n | $BEDTOOLS/closestBed -a <(less $file | sort -k1,1 -k2,2n) -b stdin -D b > $file.closest.gene
    less $file.closest.gene | sort -k8,8 | join - $dirDE/UP.IDHmut.NPC -1 8 | sed 's/ /\t/g' > $file.closest.gene.UP
    less $file.closest.gene | sort -k8,8 | join - $dirDE/DN.IDHmut.NPC -1 8 | sed 's/ /\t/g' > $file.closest.gene.DN
done

### intersect with hypomethylation
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn='/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/'
dirDMR='/projects/epigenomics2/users/lli/glioma/WGBS/DMR/'
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
for mark in H3K27ac H3K4me1; do
    dirOut=$dirIn/$mark/UMR/
    mkdir -p $dirOut
    cd $dirIn/$mark/
    for file in CEMT*.CEMT*.unique; do
        name=$(echo $file | sed -e 's/.unique//g')
        sample=$(echo $file | sed -e 's/.vs.*//g')
        echo "Processing "$mark $sample $name
        $BEDTOOLS/intersectBed -a $file -b $dirDMR/DMR.$sample'_NPC.hypo.bed' | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/$name.unique.UMR.bed
    done
    for file in CEMT*.vs.NPC_GE04.NPC_GE04.unique; do
        name=$(echo $file | sed -e 's/.unique//g')
        sample=$(echo $file | sed -e 's/.vs.*//g')
        echo "Processing "$mark $sample $name
        $BEDTOOLS/intersectBed -a $file -b $dirDMR/DMR.$sample'_NPC.hyper.bed' | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $dirOut/$name.unique.UMR.bed
    done
    cd $dirOut
    $BEDTOOLS/intersectBed -a CEMT_19.vs.NPC_GE04.CEMT_19.unique.UMR.bed -b CEMT_22.vs.NPC_GE04.CEMT_22.unique.UMR.bed > x
    $BEDTOOLS/intersectBed -a CEMT_19.vs.NPC_GE04.CEMT_19.unique.UMR.bed -b CEMT_47.vs.NPC_GE04.CEMT_47.unique.UMR.bed > y
    $BEDTOOLS/intersectBed -a CEMT_22.vs.NPC_GE04.CEMT_22.unique.UMR.bed -b CEMT_47.vs.NPC_GE04.CEMT_47.unique.UMR.bed > z
    cat x y z | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > IDHmut_NPC.IDHmut.unique.UMR.bed
    $BEDTOOLS/intersectBed -a CEMT_19.vs.NPC_GE04.NPC_GE04.unique.UMR.bed -b CEMT_22.vs.NPC_GE04.NPC_GE04.unique.UMR.bed > x
    $BEDTOOLS/intersectBed -a CEMT_19.vs.NPC_GE04.NPC_GE04.unique.UMR.bed -b CEMT_47.vs.NPC_GE04.NPC_GE04.unique.UMR.bed > y
    $BEDTOOLS/intersectBed -a CEMT_22.vs.NPC_GE04.NPC_GE04.unique.UMR.bed -b CEMT_47.vs.NPC_GE04.NPC_GE04.unique.UMR.bed > z
    cat x y z | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > IDHmut_NPC.NPC.unique.UMR.bed
    rm x y z
    cp CEMT_23.vs.NPC_GE04.CEMT_23.unique.UMR.bed IDHwt_NPC.IDHwt.unique.UMR.bed
    cp CEMT_23.vs.NPC_GE04.NPC_GE04.unique.UMR.bed IDHwt_NPC.NPC.unique.UMR.bed
    cd $dirOut
    for file in *.unique.UMR.bed; do
        name=$(echo $file | sed -e 's/.unique.UMR.bed//g')
        echo "Processing homer for "$mark $name
        mkdir -p $dirOut/homer/$name/
        /home/lli/bin/homer/bin/findMotifsGenome.pl $file hg19 $dirOut/homer/$name/ -size 200 -len 8
    done
done
dirRPKM=/projects/epigenomics2/users/lli/glioma/RNAseq/
for mark in H3K27ac H3K4me1; do
    cd $dirIn/$mark/UMR/homer/
    less IDHmut_NPC.IDHmut/knownResults.txt | awk '{if(NR==1){print}else{percent=gensub("%", "", "g", $7);if($5<=0.01&&percent>=20){print}}}' > IDHmut_NPC.IDHmut.tf
    less IDHmut_NPC.NPC/knownResults.txt | awk '{if(NR==1){print}else{percent=gensub("%", "", "g", $7);if($5<=0.01&&percent>=20){print}}}' > IDHmut_NPC.NPC.tf
    less IDHmut_NPC.IDHmut.tf | awk 'NR > 1 {print $1}' | sed 's/(.*//g' > IDHmut_NPC.IDHmut.tf.ID
    less IDHmut_NPC.NPC.tf | awk 'NR > 1 {print $1}' | sed 's/(.*//g' > IDHmut_NPC.NPC.tf.ID
done
for mark in H3K27ac H3K4me1; do
    cd $dirIn/$mark/UMR/homer/
    less IDHmut_NPC.IDHmut.tf.ID | sort -k2,2 | join - $dirRPKM/NPC_RPKM/NPC.RPKM -1 2 -2 1 | join - $dirRPKM/RPKM/glioma.RPKM | sed 's/ /\t/g' > IDHmut_NPC.IDHmut.tf.RPKM
    less IDHmut_NPC.NPC.tf.ID | sort -k2,2 | join - $dirRPKM/NPC_RPKM/NPC.RPKM -1 2 -2 1 | join - $dirRPKM/RPKM/glioma.RPKM | sed 's/ /\t/g' > IDHmut_NPC.NPC.tf.RPKM
done
#### Sox3 in H3K27ac
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn='/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/H3K27ac/'
dirRPKM='/projects/epigenomics2/users/lli/glioma/RNAseq/'
dir5mC='/projects/epigenomics2/users/lli/glioma/WGBS/'
/home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/IDHmut_NPC.NPC.unique hg19 $dirIn/homer/ -find $dirIn/homer/IDHmut_NPC.NPC/knownResults/known13.motif | awk 'NR>1 {print $1}' | uniq | sed 's/:/\t/g' | sed 's/-/\t/g'  | awk '{print $0"\t"$1":"$2"-"$3}' > $dirIn/homer/IDHmut_NPC.NPC.Sox3.annotate
cd $dirIn/homer/
for file in *.annotate; do
    tf=$(echo $file | sed 's/.annotate//g' | sed 's/IDHmut_NPC.NPC.//g')
    echo $tf
    echo -e "ID\tchr\tstart\tend" > x
    less $file | sort -k4,4 | awk '{print $4"\t"$1"\t"$2"\t"$3}' >> x
    for CpG in $dir5mC/*.5mC.CpG.combine.5mC.CpG; do
        sample=$(basename $CpG | sed 's/.5mC.CpG.combine.5mC.CpG//g');
        echo -e "ID\t$sample" > y;
        less $CpG | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a $file -b stdin -wa -wb | awk '{t[$4]=t[$4]+$8; c[$4]=c[$4]+$9} END{for(i in t){if(c[i]+t[i]>0){print i"\t"c[i]/(c[i]+t[i])}}}' | sort -k1,1 >> y;
        join x y | sed 's/ /\t/g' > z;
        mv z x;
    done
    mv x $file.5mC
    echo -e "chr\tstart\tend\tID\tgeneChr\tgeneStart\tgeneEnd\tENSG\tdis" > $file.closest.gene
    less /home/lli/hg19/hg19v69_genes.bed | awk '$4 ~ /protein_coding/ {gsub("_protein_coding", ""); print "chr"$0}' | sort -k1,1 -k2,2n | $BEDTOOLS/closestBed -a <(less $file | sort -k1,1 -k2,2n) -b stdin -D b >> $file.closest.gene
    less $file.closest.gene | sort -k8,8 | join - $dirRPKM/NPC_RPKM/NPC.RPKM -1 8 -2 1 | join - $dirRPKM/RPKM/glioma.RPKM | sed 's/ /\t/g' > $file.closest.gene.RPKM
    join $file.closest.gene.RPKM $dirRPKM/DEfine/UP.IDHmut.NPC | sed 's/ /\t/g' > $file.closest.gene.RPKM.UP
    join $file.closest.gene.RPKM $dirRPKM/DEfine/DN.IDHmut.NPC | sed 's/ /\t/g' > $file.closest.gene.RPKM.DN
done
#### Ascl1, Olig2 and NeuroD1 in H3K4me1
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
PATH=$PATH:/home/lli/bin/homer/.//bin/
PATH=$PATH:/home/acarles/weblogo/
dirIn='/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/H3K4me1/'
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
dirRPKM=/projects/epigenomics2/users/lli/glioma/RNAseq/
/home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/IDHmut_NPC.IDHmut.unique.UMR.bed hg19 $dirIn/homer/ -find $dirIn/homer/IDHmut_NPC.IDHmut/knownResults/known6.motif | awk 'NR>1 {print $1}' | uniq | sed 's/:/\t/g' | sed 's/-/\t/g'  | awk '{print $0"\t"$1":"$2"-"$3}' > $dirIn/homer/IDHmut_NPC.IDHmut.Ascl1.annotate
/home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/IDHmut_NPC.IDHmut.unique.UMR.bed hg19 $dirIn/homer/ -find $dirIn/homer/IDHmut_NPC.IDHmut/knownResults/known14.motif | awk 'NR>1 {print $1}' | uniq | sed 's/:/\t/g' | sed 's/-/\t/g'  | awk '{print $0"\t"$1":"$2"-"$3}' > $dirIn/homer/IDHmut_NPC.IDHmut.Olig2.annotate
/home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/IDHmut_NPC.IDHmut.unique.UMR.bed hg19 $dirIn/homer/ -find $dirIn/homer/IDHmut_NPC.IDHmut/knownResults/known12.motif | awk 'NR>1 {print $1}' | uniq | sed 's/:/\t/g' | sed 's/-/\t/g'  | awk '{print $0"\t"$1":"$2"-"$3}' > $dirIn/homer/IDHmut_NPC.IDHmut.NeuroD1.annotate
/home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/IDHmut_NPC.NPC.unique hg19 $dirIn/homer/ -find $dirIn/homer/IDHmut_NPC.NPC/knownResults/known2.motif | awk 'NR>1 {print $1}' | uniq | sed 's/:/\t/g' | sed 's/-/\t/g'  | awk '{print $0"\t"$1":"$2"-"$3}' > $dirIn/homer/IDHmut_NPC.NPC.Lhx2.annotate
/home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/IDHmut_NPC.NPC.unique hg19 $dirIn/homer/ -find $dirIn/homer/IDHmut_NPC.NPC/knownResults/known14.motif | awk 'NR>1 {print $1}' | uniq | sed 's/:/\t/g' | sed 's/-/\t/g'  | awk '{print $0"\t"$1":"$2"-"$3}' > $dirIn/homer/IDHmut_NPC.NPC.Sox3.annotate
/home/lli/bin/homer/bin/findMotifsGenome.pl $dirIn/IDHmut_NPC.NPC.unique hg19 $dirIn/homer/ -find $dirIn/homer/IDHmut_NPC.NPC/knownResults/known19.motif | awk 'NR>1 {print $1}' | uniq | sed 's/:/\t/g' | sed 's/-/\t/g'  | awk '{print $0"\t"$1":"$2"-"$3}' > $dirIn/homer/IDHmut_NPC.NPC.Sox6.annotate
cd $dirIn/homer/
for file in *.annotate; do
    echo $file
    echo -e "ID\tchr\tstart\tend" > x
    less $file | awk '{print $4"\t"$1"\t"$2"\t"$3}' >> x
    for CpG in $dir5mC/*.5mC.CpG.combine.5mC.CpG; do
        sample=$(basename $CpG | sed 's/.5mC.CpG.combine.5mC.CpG//g');
        echo -e "ID\t$sample" > y;
        less $CpG | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a $file -b stdin -wa -wb | awk '{t[$4]=t[$4]+$8; c[$4]=c[$4]+$9} END{for(i in t){if(c[i]+t[i]>0){print i"\t"c[i]/(c[i]+t[i])}}}' | sort -k1,1 >> y;
        join x y | sed 's/ /\t/g' > z;
        mv z x;
    done
    mv x $file.5mC
    less /home/lli/hg19/hg19v69_genes.bed | awk '$4 ~ /protein_coding/ {gsub("_protein_coding", ""); print "chr"$0}' | sort -k1,1 -k2,2n | $BEDTOOLS/closestBed -a <(less $file | sort -k1,1 -k2,2n) -b stdin -D b > $file.closest.gene
    less $file.closest.gene | sort -k8,8 | join - $dirRPKM/NPC_RPKM/NPC.RPKM -1 8 -2 1 | join - $dirRPKM/RPKM/glioma.RPKM | sed 's/ /\t/g' > $file.closest.gene.RPKM
    join $file.closest.gene.RPKM $dirRPKM/DEfine/UP.IDHmut.NPC | sed 's/ /\t/g' > $file.closest.gene.RPKM.UP
    join $file.closest.gene.RPKM $dirRPKM/DEfine/DN.IDHmut.NPC | sed 's/ /\t/g' > $file.closest.gene.RPKM.DN
done

########################################################

## K27me3 and K36me3
### enrichment in promoters genebody and intergenic regions
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
genome=/home/lli/hg19/hg19.chrom.len
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
gene=/home/lli/hg19/hg19v69_genes.bed
intergenic=/home/lli/hg19/hg19v69_intergenic.bed
#$BEDTOOLS/intersectBed -a $gene -b /home/lli/hg19/hg19.chrlen.autoXY.bed -v | $BEDTOOLS/intersectBed -a $promoter -b stdin -v | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $intergenic
$BEDTOOLS/subtractBed -a /home/lli/hg19/hg19.chrlen.autoXY.bed -b $gene | $BEDTOOLS/subtractBed -a stdin -b $promoter | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > $intergenic
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/
echo -e "Sample\tMark\tMarked\tFC_genebody\tp_genebody\tFC_promoter\tp_promoter\tFC_intergenic\tp_intergenic" > $dirIn/DHM.enrich.summary
for mark in H3K27me3 H3K36me3; do
    cd $dirIn/$mark/
    for file in *.unique; do
        sample=$(echo $file | sed 's/.vs.*//g')
        marked=$(echo $file | sed 's/.unique//g' | sed 's/.*.vs.NPC_GE04.//g')
        echo $sample $mark $marked
        fc_promoter=$(less $file | sort -k1,1 -k2,2n | $BEDTOOLS/bedtools fisher -a stdin -b <(less $promoter | awk '$1 !~ /GL/ {print "chr"$0}' | sort -k1,1 -k2,2n) -g $genome | awk 'NR==1{gsub(".*: ", ""); a=$1} NR==2{gsub(".*: ", ""); b=$1} NR==3{gsub(".*: ", ""); i=$1} NR==4{gsub(".*: ", ""); t=$1} END{print (i/a)/(b/t)}')
        p_promoter=$(less $file | sort -k1,1 -k2,2n | $BEDTOOLS/bedtools fisher -a stdin -b <(less $promoter | awk '$1 !~ /GL/ {print "chr"$0}' | sort -k1,1 -k2,2n) -g $genome | awk 'NR==5 {gsub("# ", ""); print}' | $R - | sed -e 's/\[1\] //g')
        fc_gene=$(less $file | sort -k1,1 -k2,2n | $BEDTOOLS/bedtools fisher -a stdin -b <(less $gene | awk '$1 !~ /GL|name/ {print "chr"$0}' | sort -k1,1 -k2,2n) -g $genome | awk 'NR==1{gsub(".*: ", ""); a=$1} NR==2{gsub(".*: ", ""); b=$1} NR==3{gsub(".*: ", ""); i=$1} NR==4{gsub(".*: ", ""); t=$1} END{print (i/a)/(b/t)}')
        p_gene=$(less $file | sort -k1,1 -k2,2n | $BEDTOOLS/bedtools fisher -a stdin -b <(less $gene | awk '$1 !~ /GL|name/ {print "chr"$0}' | sort -k1,1 -k2,2n) -g $genome | awk 'NR==5 {gsub("# ", ""); print}' | $R - | sed -e 's/\[1\] //g')
        fc_intergenic=$(less $file | sort -k1,1 -k2,2n | $BEDTOOLS/bedtools fisher -a stdin -b <(less $intergenic | awk '$1 !~ /GL/ {print "chr"$0}' | sort -k1,1 -k2,2n) -g $genome | awk 'NR==1{gsub(".*: ", ""); a=$1} NR==2{gsub(".*: ", ""); b=$1} NR==3{gsub(".*: ", ""); i=$1} NR==4{gsub(".*: ", ""); t=$1} END{print (i/a)/(b/t)}')
        p_intergenic=$(less $file | sort -k1,1 -k2,2n | $BEDTOOLS/bedtools fisher -a stdin -b <(less $intergenic | awk '$1 !~ /GL/ {print "chr"$0}' | sort -k1,1 -k2,2n) -g $genome | awk 'NR==5 {gsub("# ", ""); print}' | $R - | sed -e 's/\[1\] //g')
        echo -e "$sample\t$mark\t$marked\t$fc_gene\t$p_gene\t$fc_promoter\t$p_promoter\t$fc_intergenic\t$p_intergenic" >> $dirIn/DHM.enrich.summary
    done
done
### genic/intergenic distribution
samtools=/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools
promoter=/home/lli/hg19/hg19v69_genes_TSS_2000.bed
gene=/home/lli/hg19/hg19v69_genes.bed
intergenic=/home/lli/hg19/hg19v69_intergenic.bed
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/bam/
echo -e "Sample\tMark\tgenebody\tpromoter\tintergenic\tgenic_intergenic" > $dirIn/HM.distribution.summary
for mark in H3K27me3 H3K36me3; do
    cd $dirIn/$mark/
    for bam in *.bam; do
        sample=$(echo $bam | sed 's/.bam//g' | sed 's/A[0-9]*.//g')
        echo $sample $mark
        n_gene=$($samtools view -q 5 -F 1028 -L $gene $bam | wc -l)
        n_promoter=$($samtools view -q 5 -F 1028 -L $promoter $bam | wc -l)
        n_intergenic=$($samtools view -q 5 -F 1028 -L $intergenic $bam | wc -l)
        echo -e "$sample\t$mark\t$n_gene\t$n_promoter\t$n_intergenic" | awk '{print $0"\t"($3+$4)/$5}' >> $dirIn/HM.distribution.summary
    done
done
### Loss of H3K36me3
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
dirIn='/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/H3K36me3/'
dir5mC=/projects/epigenomics2/users/lli/glioma/WGBS/
cd $dirIn
$BEDTOOLS/intersectBed -a CEMT_19.vs.NPC_GE04.NPC_GE04.unique -b CEMT_22.vs.NPC_GE04.NPC_GE04.unique | $BEDTOOLS/intersectBed -a stdin -b CEMT_47.vs.NPC_GE04.NPC_GE04.unique | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > IDHmut_NPC_NPC.unique
cat CEMT_23.vs.NPC_GE04.NPC_GE04.unique | awk '{print $1"\t"$2"\t"$3"\t"$4}' > IDHwt_NPC_NPC.unique
$BEDTOOLS/intersectBed -a IDHmut_NPC_NPC.unique -b IDHwt_NPC_NPC.unique | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > IDH_NPC_NPC.unique
#### significance of intersect
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
genome=/home/lli/hg19/hg19.chrom.len
dirOut=/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/H3K36me3/
lib=( CEMT_19 CEMT_21 CEMT_22 CEMT_23 CEMT_47 )
echo -e "Sample1\tSample2\tN1\tN2\tIntersect\tFC\tPvalue" > $dirOut/H3K36me3.loss.intersect.summary
cd $dirOut; 
for i in {0..3}; do
    for ((j = i+1; j < 5; j++)); do
        sample1=${lib[$i]}; sample2=${lib[$j]};
        file1=$sample1.vs.NPC_GE04.NPC_GE04.unique; file2=$sample2.vs.NPC_GE04.NPC_GE04.unique;
        echo $sample1 $sample2;
        n1=$(less $file1 | wc -l); n2=$(less $file2 | wc -l);
        intersect=$($BEDTOOLS/intersectBed -a $file1 -b $file2 | wc -l);
        fc=$($BEDTOOLS/bedtools fisher -a <(less $file1 | sort -k1,1 -k2,2n) -b <(less $file2| sort -k1,1 -k2,2n) -g $genome | awk 'NR==1{gsub(".*: ", ""); a=$1} NR==2{gsub(".*: ", ""); b=$1} NR==3{gsub(".*: ", ""); i=$1} NR==4{gsub(".*: ", ""); t=$1} END{print (i/a)/(b/t)}')
        p=$($BEDTOOLS/bedtools fisher -a <(less $file1 | sort -k1,1 -k2,2n) -b <(less $file2| sort -k1,1 -k2,2n) -g $genome | awk 'NR==5 {gsub("# ", ""); print}' | $R - | sed -e 's/\[1\] //g')
        echo -e "$sample1\t$sample2\t$n1\t$n2\t$intersect\t$fc\t$p" >> $dirOut/H3K36me3.loss.intersect.summary
    done
done
#### associate with DE
echo -e "File\tN_K36\tN_UP\tN_DN\tN_K36_UP\tN_K36_DN" > $dirIn/H3K36me3_loss_DE.summary
for file in *.vs.NPC_GE04.NPC_GE04.unique; do
    echo $file
    echo -e "ID\tchr\tstart\tend" > x
    less $file | awk '{print $4"\t"$1"\t"$2"\t"$3}' | sort -k1,1 >> x
    for CpG in $dir5mC/*.5mC.CpG.combine.5mC.CpG; do
        sample=$(basename $CpG | sed 's/.5mC.CpG.combine.5mC.CpG//g');
        echo -e "ID\t$sample" > y;
        less $CpG | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a $file -b stdin -wa -wb | awk '{t[$4]=t[$4]+$10; c[$4]=c[$4]+$11} END{for(i in t){if(c[i]+t[i]>0){print i"\t"c[i]/(c[i]+t[i])}}}' | sort -k1,1 >> y;
        join x y | sed 's/ /\t/g' > z;
        mv z x;
    done
    mv x $file.5mC
    less /home/lli/hg19/hg19v69_genes.bed | awk '$4 ~ /protein_coding/ {gsub("_protein_coding", ""); print "chr"$0}' | sort -k1,1 -k2,2n | $BEDTOOLS/intersectBed -a <(less $file | sort -k1,1 -k2,2n) -b stdin -wa -wb > $file.gene
    less $file.gene | sort -k10,10 | join - $dirRPKM/NPC_RPKM/NPC.RPKM -1 10 -2 1 | join - $dirRPKM/RPKM/glioma.RPKM | sed 's/ /\t/g' > $file.gene.RPKM
    join $file.gene.RPKM $dirRPKM/DEfine/UP.IDHmut.NPC | sed 's/ /\t/g' > $file.gene.RPKM.UP
    join $file.gene.RPKM $dirRPKM/DEfine/DN.IDHmut.NPC | sed 's/ /\t/g' > $file.gene.RPKM.DN
    N_K36=$(less $file.gene.RPKM | awk '{print $1}' | uniq | wc -l)
    N_UP=$(less $dirRPKM/DEfine/UP.IDHmut.NPC | wc -l)
    N_DN=$(less $dirRPKM/DEfine/DN.IDHmut.NPC | wc -l)
    N_K36_UP=$(less $file.gene.RPKM.UP | awk '{print $1}' | uniq | wc -l)
    N_K36_DN=$(less $file.gene.RPKM.DN | awk '{print $1}' | uniq | wc -l)
    echo -e "$file\t$N_K36\t$N_UP\t$N_DN\t$N_K36_UP\t$N_K36_DN" >> $dirIn/H3K36me3_loss_DE.summary
done
#### associate with RPKM and 5mC
for file in IDH*.unique; do
    echo $file
    echo -e "ID\tchr\tstart\tend" > x
    less $file | awk '{print $4"\t"$1"\t"$2"\t"$3}' | sort -k1,1 >> x
    for CpG in $dir5mC/*.5mC.CpG.combine.5mC.CpG; do
        sample=$(basename $CpG | sed 's/.5mC.CpG.combine.5mC.CpG//g');
        echo -e "ID\t$sample" > y;
        less $CpG | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5}' | $BEDTOOLS/intersectBed -a $file -b stdin -wa -wb | awk '{t[$4]=t[$4]+$8; c[$4]=c[$4]+$9} END{for(i in t){if(c[i]+t[i]>0){print i"\t"c[i]/(c[i]+t[i])}}}' | sort -k1,1 >> y;
        join x y | sed 's/ /\t/g' > z;
        mv z x;
    done
    mv x $file.5mC
    less /home/lli/hg19/hg19v69_genes.bed | awk '$4 ~ /protein_coding/ {gsub("_protein_coding", ""); print "chr"$0}' | sort -k1,1 -k2,2n | $BEDTOOLS/intersectBed -a <(less $file | sort -k1,1 -k2,2n) -b stdin -wa -wb > $file.gene
    less $file.gene | sort -k8,8 | join - $dirRPKM/NPC_RPKM/NPC.RPKM -1 8 -2 1 | join - $dirRPKM/RPKM/glioma.RPKM | sed 's/ /\t/g' > $file.gene.RPKM
    join $file.gene.RPKM $dirRPKM/DEfine/UP.IDHmut.NPC | sed 's/ /\t/g' > $file.gene.RPKM.UP
    join $file.gene.RPKM $dirRPKM/DEfine/DN.IDHmut.NPC | sed 's/ /\t/g' > $file.gene.RPKM.DN
    N_K36=$(less $file.gene.RPKM | awk '{print $1}' | uniq | wc -l)
    N_UP=$(less $dirRPKM/DEfine/UP.IDHmut.NPC | wc -l)
    N_DN=$(less $dirRPKM/DEfine/DN.IDHmut.NPC | wc -l)
    N_K36_UP=$(less $file.gene.RPKM.UP | awk '{print $1}' | uniq | wc -l)
    N_K36_DN=$(less $file.gene.RPKM.DN | awk '{print $1}' | uniq | wc -l)
    echo -e "$file\t$N_K36\t$N_UP\t$N_DN\t$N_K36_UP\t$N_K36_DN" >> $dirIn/H3K36me3_loss_DE.summary
done
#### assocaite with H3K27me3
dirIn='/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/H3K36me3/'
dirWig=/projects/epigenomics2/users/lli/glioma/ChIPseq/wig/
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
RegCov=/home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar
chr=/home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes
cd $dirIn
for file in *.vs.NPC_GE04.NPC_GE04.unique; do
    sample=$(echo $file | sed 's/.vs.NPC_GE04.NPC_GE04.unique//g')
    echo $sample
    echo -e "ID\tchr\tstart\tend" > a
    less $file | awk '{print $4"\t"$1"\t"$2"\t"$3}' | sort -k1,1 >> a
    cov_glioma=$(less /projects/epigenomics2/users/lli/glioma/ChIPseq/H3K27me3.depth | awk '{if($1=="'$sample'"){print $2}}')
    cov_GE04=35922595
    $JAVA -jar -Xmx15G $RegCov -w $dirWig/H3K27me3/$sample.wig.gz -r <(less $file | awk '{print $1"\t"$2"\t"$3"\t"$4}') -o $dirIn -c $chr -n $file.$sample.H3K27me3 > $dirIn/$file.H3K27me3.log
    $JAVA -jar -Xmx15G $RegCov -w $dirWig/H3K27me3/NPC_GE04.q5.F1028.PET.wig.gz -r <(less $file | awk '{print $1"\t"$2"\t"$3"\t"$4}') -o $dirIn -c $chr -n $file.NPC_GE04.H3K27me3 >> $dirIn/$file.H3K27me3.log
    echo -e "ID\t$sample\tNPC_GE04\tdiff_H3K27me3" > b
    join *$file.$sample.H3K27me3.coverage *$file.NPC_GE04.H3K27me3.coverage -1 4 -2 4 | awk '{print $1"\t"$5/"'$cov_glioma'"*"'$cov_GE04'""\t"$10"\t"$5/"'$cov_glioma'"*"'$cov_GE04'"-$10}' | sort -k1,1 >> b
    join a b | sed 's/ /\t/g' > $file.H3K27me3
    cov_glioma=$(less /projects/epigenomics2/users/lli/glioma/ChIPseq/H3K36me3.depth | awk '{if($1=="'$sample'"){print $2}}')
    cov_GE04=53179871
    $JAVA -jar -Xmx15G $RegCov -w $dirWig/H3K36me3/$sample.wig.gz -r <(less $file | awk '{print $1"\t"$2"\t"$3"\t"$4}') -o $dirIn -c $chr -n $file.$sample.H3K36me3 > $dirIn/$file.H3K36me3.log
    $JAVA -jar -Xmx15G $RegCov -w $dirWig/H3K36me3/NPC_GE04.q5.F1028.PET.wig.gz -r <(less $file | awk '{print $1"\t"$2"\t"$3"\t"$4}') -o $dirIn -c $chr -n $file.NPC_GE04.H3K36me3 >> $dirIn/$file.H3K36me3.log
    echo -e "ID\t$sample\tNPC_GE04\tdiff_H3K36me3" > b
    join *$file.$sample.H3K36me3.coverage *$file.NPC_GE04.H3K36me3.coverage -1 4 -2 4 | awk '{print $1"\t"$5/"'$cov_glioma'"*"'$cov_GE04'""\t"$10"\t"$5/"'$cov_glioma'"*"'$cov_GE04'"-$10}' | sort -k1,1 >> b
    join a b | sed 's/ /\t/g' > $file.H3K36me3
    join $file.H3K27me3 $file.H3K36me3 | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$7"\t"$13}' > $file.H3K27me3.H3K36me3
    rm *$file.$sample.H3K*me3* *$file.NPC_GE04.H3K*me3* a b
done
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
genome=/home/lli/hg19/hg19.chrom.len
gene=/home/lli/hg19/hg19v69_genes.bed
cutoff=2
echo -e "Sample\tH3K27me3\tN_total\tN\tFC_gene\tP_gene" > $dirIn/H3K36me3.loss.H3K27me3.summary
for file in *.H3K27me3.H3K36me3; do
    lib=$(echo $file | sed 's/.vs.NPC_GE04.NPC_GE04.unique.H3K27me3.H3K36me3//g')
    echo $lib
    less $file | awk 'NR>1 {if($5>="'$cutoff'"){print $1"\t"$2"\t"$3"\t"$4}}' > $dirIn/$file.K27gain
    less $file | awk 'NR>1 {if($5<=-"'$cutoff'"){print $1"\t"$2"\t"$3"\t"$4}}' > $dirIn/$file.K27loss
    n_total=$(less $file | wc -l); n_gain=$(less $file.K27gain | wc -l); n_loss=$(less $file.K27loss | wc -l);
    fc_gain=$($BEDTOOLS/bedtools fisher -a <(less $file.K27gain | sort -k1,1 -k2,2n) -b <(less $gene | awk '$1 !~ /GL/ {print "chr"$0}' | sort -k1,1 -k2,2n) -g $genome | awk 'NR==1{gsub(".*: ", ""); a=$1} NR==2{gsub(".*: ", ""); b=$1} NR==3{gsub(".*: ", ""); i=$1} NR==4{gsub(".*: ", ""); t=$1} END{print (i/a)/(b/t)}')
    p_gain=$($BEDTOOLS/bedtools fisher -a <(less $file.K27gain | sort -k1,1 -k2,2n) -b <(less $gene | awk '$1 !~ /GL/ {print "chr"$0}' | sort -k1,1 -k2,2n) -g $genome | awk 'NR==5 {gsub("# ", ""); print}' | $R - | sed -e 's/\[1\] //g')
    fc_loss=$($BEDTOOLS/bedtools fisher -a <(less $file.K27loss | sort -k1,1 -k2,2n) -b <(less $gene | awk '$1 !~ /GL/ {print "chr"$0}' | sort -k1,1 -k2,2n) -g $genome | awk 'NR==1{gsub(".*: ", ""); a=$1} NR==2{gsub(".*: ", ""); b=$1} NR==3{gsub(".*: ", ""); i=$1} NR==4{gsub(".*: ", ""); t=$1} END{print (i/a)/(b/t)}')
    p_loss=$($BEDTOOLS/bedtools fisher -a <(less $file.K27loss | sort -k1,1 -k2,2n) -b <(less $gene | awk '$1 !~ /GL/ {print "chr"$0}' | sort -k1,1 -k2,2n) -g $genome | awk 'NR==5 {gsub("# ", ""); print}' | $R - | sed -e 's/\[1\] //g')
    echo -e "$lib\tgain\t$n_total\t$n_gain\t$fc_gain\t$p_gain" >> $dirIn/H3K36me3.loss.H3K27me3.summary
    echo -e "$lib\tloss\t$n_total\t$n_loss\t$fc_loss\t$p_loss" >> $dirIn/H3K36me3.loss.H3K27me3.summary
done
BEDTOOLS=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/
R=/gsc/software/linux-x86_64-centos5/R-3.1.1/bin/Rscript
genome=/home/lli/hg19/hg19.chrom.len
dirOut=/projects/epigenomics2/users/lli/glioma/ChIPseq/unique/H3K36me3/
lib=( CEMT_19 CEMT_21 CEMT_22 CEMT_23 CEMT_47 )
echo -e "Sample1\tSample2\tK27\tN1\tN2\tIntersect\tFC\tPvalue" > $dirOut/H3K36me3.loss.H3K27me3.intersect.summary
cd $dirOut;
for K27 in gain loss; do
    for i in {0..3}; do
        for ((j = i+1; j < 5; j++)); do
            sample1=${lib[$i]}; sample2=${lib[$j]};
            file1=$sample1'.vs.NPC_GE04.NPC_GE04.unique.H3K27me3.H3K36me3.K27'$K27; file2=$sample2'.vs.NPC_GE04.NPC_GE04.unique.H3K27me3.H3K36me3.K27'$K27;
            echo $sample1 $sample2 $K27;
            n1=$(less $file1 | wc -l); n2=$(less $file2 | wc -l);
            intersect=$($BEDTOOLS/intersectBed -a $file1 -b $file2 | wc -l);
            fc=$($BEDTOOLS/bedtools fisher -a <(less $file1 | sort -k1,1 -k2,2n) -b <(less $file2| sort -k1,1 -k2,2n) -g $genome | awk 'NR==1{gsub(".*: ", ""); a=$1} NR==2{gsub(".*: ", ""); b=$1} NR==3{gsub(".*: ", ""); i=$1} NR==4{gsub(".*: ", ""); t=$1} END{print (i/a)/(b/t)}')
            p=$($BEDTOOLS/bedtools fisher -a <(less $file1 | sort -k1,1 -k2,2n) -b <(less $file2| sort -k1,1 -k2,2n) -g $genome | awk 'NR==5 {gsub("# ", ""); print}' | $R - | sed -e 's/\[1\] //g')
            echo -e "$sample1\t$sample2\t$K27\t$n1\t$n2\t$intersect\t$fc\t$p" >> $dirOut/H3K36me3.loss.H3K27me3.intersect.summary
        done
    done
done
cat *.H3K27me3.H3K36me3.K27gain | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin > IDH_NPC_NPC.unique.H3K27me3.H3K36me3.K27gain
cat *.H3K27me3.H3K36me3.K27loss | sort -k1,1 -k2,2n | $BEDTOOLS/mergeBed -i stdin > IDH_NPC_NPC.unique.H3K27me3.H3K36me3.K27loss

## Histone modifiers expression
cd /projects/epigenomics2/users/lli/glioma/RNAseq/
less Histone.modifiers | sort -k2,2 | join - NPC_RPKM/NPC.RPKM -1 2 -2 1 | join - RPKM/glioma.RPKM | sed 's/ /\t/g' > Histone.modifiers.RPKM
## HOXA/B expression
cd /projects/epigenomics2/users/lli/glioma/RNAseq/
less /home/lli/hg19/HOXA.HOXB.txt | sort -k1,1 | join - NPC_RPKM/NPC.RPKM | join - RPKM/glioma.RPKM | sed 's/ /\t/g' > HOXA.HOXB.RPKM

## Global histone modification changes
JAVA=/home/mbilenky/jdk1.8.0_92/jre/bin/java
RegCov=/home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar
dirOut=/projects/epigenomics2/users/lli/glioma/ChIPseq/global/
mkdir -p $dirOut
### NPC
dirIn=/home/lli/FetalBrain/HisMod/wigs/
libInfo=/projects/epigenomics/users/lli/FetalBrain/FetalBrainLibrariesDetail.tsv
for wig in $dirIn/*.wig.gz; do
	lib=$(basename $wig | cut -d'.' -f1)
	name=$(less $libInfo | awk -F "\t" '$1~"'$lib'" {gsub("HuFNSC", "", $3); gsub("Primary Cell Culture Neurospheres, ", "", $4); gsub(" Derived", "", $4); gsub("Ganglionic Eminence", "GE", $4); print "NPC_"$4$3}')
	mark=$(less $libInfo | awk '$1~"'$lib'" {print $2}')
	echo $lib $mark $name
	mkdir -p $dirOut/$mark/
	$JAVA -jar -Xmx10G $RegCov -w $wig -o $dirOut/$mark/ -s hg19 -n $name -bin 500 -step 500
done
### CEMT
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/wig/
for wig in $dirIn/*/*.wig.gz; do
	name=$(basename $wig | cut -d'.' -f1)
	mark=$(echo $wig | sed "s|$dirIn/||g" | sed "s/\\/.*//g")
	mkdir -p $dirOut/$mark/
	if [ $name != "NPC_GE04" ]; then
		echo $mark $name
		$JAVA -jar -Xmx10G $RegCov -w $wig -o $dirOut/$mark/ -s hg19 -n $name -bin 500 -step 500
	fi
done
### Adult brain
chr=/projects/epigenomics/resources/UCSC_chr/hg19.bwa2ucsc.names
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/NB141/
mkdir -p $dirIn/wigs/
fl=$(echo 175 | awk '{print $1}');
for bam in $dirIn/*.bam; do
	name=$(basename $bam | sed -e 's/.bam//g');
	mark=$(echo $name | sed 's/NB141.//g')
	echo $mark $name;
	mkdir -p $dirOut/$mark/;
	/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $bam $dirIn/wigs/ -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr;
	$JAVA -jar -Xmx10G $RegCov -w $(ls $dirIn/wigs/$name.*.wig.gz) -o $dirOut/$mark/ -s hg19 -n $name -bin 500 -step 500;
done
### Other CEMT
for wig in /projects/edcc_new/reference_epigenomes/CEMT_*/bams/ChIP-Seq/*/wig/*.wig.gz; do
	mark=$(echo $wig | cut -d'/' -f8);
	lib=$(echo $wig | cut -d'/' -f5);
	name="Other_"$lib;
	if [ $lib != "CEMT_19" ] && [ $lib != "CEMT_21" ] && [ $lib != "CEMT_22" ] && [ $lib != "CEMT_23" ] && [ $lib != "CEMT_47" ] && [ $lib != "CEMT_47_PB" ]; then
		echo $mark $name;
		echo $mark $name $(basename $wig) >> $dirOut/signal.log;
		$JAVA -jar -Xmx10G $RegCov -w $wig -o $dirOut/$mark/ -s hg19 -n $name -bin 500 -step 500;
	fi
done
### compare
dirIn=/projects/epigenomics2/users/lli/glioma/ChIPseq/global/;
for mark in H3K27me3 H3K36me3 H3K9me3 H3K27ac H3K4me1 H3K4me3 Input; do
	cd $dirIn/$mark/;
	echo -e "S1\tS2\tMark\tType\tmin\tymin\tbottom\tmiddle\ttop\tymax\tmax" > $dirIn/$mark.glioma_NPC.quantile;
	for file in *.coverage; do
		echo $mark $file
		less $file | awk '{if($4>0){print $1":"$2"-"$3"\t"$4}}' | sort -k2,2nr -T /genesis/scratch/epigenomics/lli/ > $file.sorted;
		less $file.sorted | awk '{if(NR<=400){print $0 >> "'$file'"".sorted.top"}else{print $0 >> "'$file'"".sorted.rest"}}';
	done
	for f1 in *CEMT_19*.coverage *CEMT_21*.coverage *CEMT_22*.coverage *CEMT_23*.coverage *CEMT_47*.coverage *NB141*.coverage; do
		for f2 in *NPC*.coverage *Other*.coverage; do
			s1=$(basename $f1 | cut -d'.' -f4);
			s2=$(basename $f2 | cut -d'.' -f4);
			echo $mark $s1 $s2;
			join <(less $f1.sorted.top | sort -k1,1) <(less $f2.sorted.top | sort -k1,1) -t $'\t' | awk '{print $0"\t"$2/$3}' | sort -k4,4n | awk '{s[NR]=$4} END{print "'$s1'""\t""'$s2'""\t""'$mark'""\ttop\t"s[1]"\t"s[int(NR/10)]"\t"s[int(NR/4)]"\t"s[int(NR/2)]"\t"s[NR-int(NR/4)]"\t"s[NR-int(NR/10)]"\t"s[NR]}' >> $dirIn/$mark.glioma_NPC.quantile;
			join <(less $f1.sorted.rest | sort -k1,1 -T /genesis/scratch/epigenomics/lli/) <(less $f2.sorted.rest | sort -k1,1 -T /genesis/scratch/epigenomics/lli/) -t $'\t' | awk '{print $0"\t"$2/$3}' | sort -k4,4n -T /genesis/scratch/epigenomics/lli/ | awk '{s[NR]=$4} END{print "'$s1'""\t""'$s2'""\t""'$mark'""\tglobal\t"s[1]"\t"s[int(NR/10)]"\t"s[int(NR/4)]"\t"s[int(NR/2)]"\t"s[NR-int(NR/4)]"\t"s[NR-int(NR/10)]"\t"s[NR]}' >> $dirIn/$mark.glioma_NPC.quantile;
		done
	done
        file1=genomeProfile.bin500.step500.CEMT_19.coverage
        file2=genomeProfile.bin500.step500.NPC_GE04.coverage
        join <(less $file1 | awk '{print $1":"$2"-"$3"\t"$4}') <(less $file1 | awk '{print $1":"$2"-"$3"\t"$4}') | awk '{if($2+$3>0){print $1"\t"$2"\t"$3}}' > CEMT_19.NPC_GE04.coverage
done


