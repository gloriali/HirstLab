#!/bin/sh

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
        lib=$(echo $file | sed -e 's/.FDR_0.05.FindER.bed.gz//g' | sed -e 's/.*Input.//g' | sed -e 's/GE04/NPC_GE04/g')
        echo -e $mark"\t"$lib"\t"$(less $file | wc -l)"\t"$(less $file | awk '{s=s+$3-$2}END{print s}') >> $dirOut/ER.summary.stats
    done
done

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



