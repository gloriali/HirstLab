#!/bin/sh

## QC 
bamstats=/gsc/QA-bio/sbs-solexa/opt/linux-x86_64/bwa_stats_0.1.3/bamStats.py
dirIn=/projects/epigenomics2/Brain/NJabado/ChIPseq/bam/
dirOut=/projects/epigenomics2/Brain/NJabado/ChIPseq/QC/
mkdir -p $dirOut
cd $dirIn
for bam in *.bam; do
    name=$(echo $bam | sed -e 's/.bam//g' | sed -e 's/C9E81ANXX_[3-6]_//g')
    echo $name
    $bamstats -g 2864785220 -q 10 -b $dirIn/$bam  > $dirIn/$name.bamstats
    /home/lli/HirstLab/Pipeline/shell/bamstats2report.sh $dirOut $name $dirIn/$name.bamstats
done
/home/lli/HirstLab/Pipeline/shell/bamstats2report.combine.sh $dirOut $dirOut
dirIn=/projects/epigenomics2/Brain/NJabado/ChIPseq/bam/
dirOut=/projects/epigenomics2/Brain/NJabado/ChIPseq/QC/PETlen/
mkdir -p $dirOut
cd $dirIn
for bam in *.bam; do
    name=$(echo $bam | sed -e 's/.bam//g' | sed -e 's/C9E81ANXX_[3-6]_//g')
    echo $name
    /home/mbilenky/bin/PETLengthDist.sh $dirIn/$bam 5 $dirOut 10
done
samtools=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
dirQC=/projects/edcc_new/reference_epigenomes/housekeeping/EDCCProd/resources/ChIPQC1/Brads/
dirIn=/projects/epigenomics2/Brain/NJabado/ChIPseq/bam/
dirOut=/projects/epigenomics2/Brain/NJabado/ChIPseq/QC/
echo -e "Library\tNumber_Target_Regions" > $dirOut/QC.target_regions.txt
cd $dirIn
for bam in *.bam; do
    lib=$(echo $bam | sed -e 's/_C9E81ANXX_[3-6]//g' | sed -e 's/.bam//g');
    sample=$(less ../labels.samples | awk -F "\t" '$5 ~ "'$lib'" {print $1}')
    mark=$(echo $sample | sed 's/PC.//g' | cut -d'_' -f 2 | sed 's/Input[12]/Input/g')
    echo $lib $mark
    if [ $mark == "Input" ]; then
        n=0
    else
        if [ "$mark" == "H3K4me3" ] || [ "$mark" == "H3K27ac" ]; then
            region=$dirQC/ensembl_TSS.uniq.sorted.bed.uniq
        fi
        if [ "$mark" == "H3K4me1" ]; then
            region=$dirQC/encode_ChromHMM_enhancer_state_7.sorted.merged.bed
        fi
        if [ "$mark" == "H3K9me3" ]; then
            region=$dirQC/ensembl_Znf.uniq.sorted.bed
        fi
        if [ "$mark" == "H3K27me3" ]; then
            region=$dirQC/HOX_clusters.sorted.bed
        fi
        if [ "$mark" == "H3K36me3" ]; then
            region=$dirQC/ensembl_genenames.uniq.sorted.bed
        fi
        echo $region
        n=$($samtools view -q 5 -F 1028 -L $region $bam | wc -l)
    fi
    echo -e "$lib\t$n" >> $dirOut/QC.target_regions.txt
done

# FindER
JAVA='/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java'
samtools='/gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools'
FindER='/home/mbilenky/bin/Solexa_Java/FindER.1.0.0b.jar'
dirIn='/projects/epigenomics2/Brain/NJabado/ChIPseq/bam/'
dirOut='/projects/epigenomics2/Brain/NJabado/ChIPseq/FindER/'
mkdir -p $dirOut
echo -e "Sample\tMark\tN_region\tTotal_length" > $dirOut/ER.summary.stats
## indexing
function index {
    bam=$1
    /gsc/software/linux-x86_64-centos5/samtools-0.1.18/bin/samtools index $bam
}
export -f index
cd $dirIn
ls *.bam > $dirIn/BamList.txt
cat BamList.txt | parallel --gnu index
## FindER 
for bam in *.bam; do
    lib=$(echo $bam | sed -e 's/_C9E81ANXX_[3-6]//g' | sed -e 's/.bam//g');
    sample=$(less ../labels.samples | awk -F "\t" '$5 ~ "'$lib'" {print $1}')
    mark=$(echo $sample | sed 's/PC.//g' | cut -d'_' -f 2 | sed 's/Input[12]/Input/g')
    if [ "$mark" != "Input" ]; then
        if [[ "$sample" =~ "Input1" ]]; then
            bam2='PX0378_C9E81ANXX_5_CTATAC.bam'
        elif [[ "$sample" =~ "Input2" ]]; then
            bam2='PX0379_C9E81ANXX_6_ATACGG.bam'
        elif [[ "$sample" =~ "AK1189" ]]; then
            bam2='PX0376_C9E81ANXX_3_CCACGC.bam'
        else
            bam2='PX0377_C9E81ANXX_4_GCAAGG.bam'
        fi
        sample=$(echo $sample | sed 's/PC.//g' | sed 's/_.*//g')
        echo $lib $sample $mark 
        $JAVA -jar -Xmx12G $FindER -signalBam $bam -inputBam $bam2 -out $dirOut > $dirOut/$lib.out
        n=$(less $dirOut/$(echo $bam | sed 's/.bam//g')*.FindER.bed.gz | wc -l)
        len=$(less $dirOut/$(echo $bam | sed 's/.bam//g')*.FindER.bed.gz | awk '{s=s+$3-$2}END{print s}')
        echo -e $sample"\t"$mark"\t"$n"\t"$len >> $dirOut/ER.summary.stats 
    fi
done

