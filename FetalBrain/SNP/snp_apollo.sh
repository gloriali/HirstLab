#!/bin/sh
#$ -S /bin/sh
#$ -m e
#$ -M lli@bcgsc.ca
#$ -N SNP
#$ -V

# SNP calling 
while [ $# -gt 0 ]
do
    case "$1" in
        -i) dirIn="$2"; shift;;
        -o) dirOut="$2"; shift;;
        -b) bamList="$2"; shift;;
        -n) name="$2"; shift;;
    esac
    shift
done

mkdir -p $dirOut
genome=/home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa
SAMTOOLS=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
BCFTOOLS=/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/

$SAMTOOLS mpileup -C50 -uf $genome -b $dirIn/$bamList | $BCFTOOLS/bcftools view -bvcg - > $dirOut/$name.bcf
$BCFTOOLS/bcftools view $dirOut/$name.bcf | $BCFTOOLS/vcfutils.pl varFilter -D100 > $dirOut/$name.vcf
