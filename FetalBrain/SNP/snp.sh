#!/bin/sh

# /home/lli/MeDIPMRE/
## brain01 brain02
/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools mpileup -C50 -uf /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa HS2788.bam HS2790.bam | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view -bvcg - > ./SNP/HS2788_HS2790.bcf
/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view HS2788_HS2790.bcf | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/vcfutils.pl varFilter -D100 > HS2788_HS2790.vcf

## cortex01 cortex02
/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools mpileup -C50 -uf /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa HS2775.bam HS2779.bam | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view -bvcg - > ./SNP/HS2775_HS2779.bcf
/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view HS2775_HS2779.bcf | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/vcfutils.pl varFilter -D100 > HS2775_HS2779.vcf

## ge01 ge02
/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools mpileup -C50 -uf /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa HS2777.bam HS2781.bam | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view -bvcg - > ./SNP/HS2777_HS2781.bcf
/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view HS2777_HS2781.bcf | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/vcfutils.pl varFilter -D100 > HS2777_HS2781.vcf

## myo35 myo70
/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools mpileup -C50 -uf /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa HS1394.bam HS2294.bam | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view -bvcg - > HS1394_HS2294.bcf
/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view HS1394_HS2294.bcf | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/vcfutils.pl varFilter -D100 > HS1394_HS2294.vcf

#####################
# use RNA-seq
genome=/home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa
SAMTOOLS=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
BCFTOOLS=/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/
dirIn=/projects/epigenomics/users/lli/FetalBrain/RNAseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/SNP/

cd $dirIn
for cell in "Cortex" "GE"; do
    files=($(ls *$cell*))
    echo ${files[@]}
    for((i=0; i<${#files[@]}; i++ )); do
        for((j=i+1; j<${#files[@]}; j++ )); do
            bam1=${files[i]}
            bam2=${files[j]}
            name=$(echo $bam1 | cut -d'.' -f2)_$(echo $bam2 | cut -d'.' -f2)
            echo $name
            $SAMTOOLS mpileup -C50 -uf $genome $dirIn/$bam1 $dirIn/$bam2 | $BCFTOOLS/bcftools view -bvcg - > $dirOut/$name.bcf
            $BCFTOOLS/bcftools view $dirOut/$name.bcf | $BCFTOOLS/vcfutils.pl varFilter -D100 > $dirOut/$name.vcf
        done
    done
done

cd $dirOut
for file in Cortex*.vcf; do
    file1=$file
    file2=$(echo $file | sed -e 's/Cortex/GE/g')
    name=$(echo $file | sed -e 's/Cortex/HuFNSC/g')
    echo -e $file1"\t"$file2"\t"$name
    awk 'NR==FNR {if(!($1~/^#/)){h[$3]=$0; next}} {if($3 in h){print $0"\t"h[$3]}}' $file2 $file1 > $dirOut/$name
done

