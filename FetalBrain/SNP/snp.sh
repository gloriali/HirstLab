#!/bin/sh

# use MeDIP/MRE
dirOut=/home/lli/FetalBrain/MeDIPMRE/SNP/
cd $dirOut
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

## discordant SNPs in MeDIP
less HS2775_HS2779.vcf | awk '$1 !~ "#" {gt1=gensub(":.*", "", "g", $10); gt2=gensub(":.*", "", "g", $11); gt3=gensub(":.*", "", "g", $12); gt4=gensub(":.*", "", "g", $13); if(gt1 == gt2 && gt3 == gt4 && gt1 != gt3){print $0}}' > MZ_discordant_SNP_cortex.vcf
less HS2777_HS2781.vcf | awk '$1 !~ "#" {gt1=gensub(":.*", "", "g", $10); gt2=gensub(":.*", "", "g", $11); gt3=gensub(":.*", "", "g", $12); gt4=gensub(":.*", "", "g", $13); if(gt1 == gt2 && gt3 == gt4 && gt1 != gt3){print $0}}' > MZ_discordant_SNP_GE.vcf
awk 'NR==FNR {id=$1":"$2; gt=gensub(":.*", "", "g", $10)"_"gensub(":.*", "", "g", $12); h[id]= gt; next} {id=$1":"$2; gt=gensub(":.*", "", "g", $10)"_"gensub(":.*", "", "g", $12); if(id in h && gt == h[id]){print $0}}' MZ_discordant_SNP_cortex.vcf MZ_discordant_SNP_GE.vcf > MZ_discordant_SNP.vcf

#####################
# use RNA-seq
dirIn=/projects/epigenomics/users/lli/FetalBrain/RNAseq/bam/
dirOut=/projects/epigenomics/users/lli/FetalBrain/SNP/
ls $dirIn/*.bam > $dirIn/bamlist.txt
mkdir -p $dirOut
name=FetalBrain_RNAseq_SNP
genome=/home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa
SAMTOOLS=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
BCFTOOLS=/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/
$SAMTOOLS mpileup -C50 -uf $genome -b $dirIn/bamlist.txt | $BCFTOOLS/bcftools view -bvcg - > $dirOut/$name.bcf
$BCFTOOLS/bcftools view $dirOut/$name.bcf | $BCFTOOLS/vcfutils.pl varFilter -D100 > $dirOut/$name.vcf
less $dirOut/$name.vcf | awk '$1 ~ /#/ {print $0}' > $dirOut/$name.header.vcf
less $dirOut/$name.vcf | awk '$1 !~ /#/ {print $0}' > $dirOut/$name.main.vcf

# filter against dbSNP
awk 'NR==FNR {pos=$2":"$3; id[pos]=$5; next} {pos="chr"$1":"$2; if(pos in id){print id[pos]"\t"$0}}' /projects/epigenomics/resources/UCSC_hg19/dbSNP138/snp138Common.table $dirOut/$name.main.vcf > $dirOut/$name.main.dbSNP.vcf

# 01 vs 02
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $11); gt2=gensub(":.+", "", "g", $13); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}}' > $dirOut/cortex01_cortex02.main.dbSNP.vcf
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $12); gt2=gensub(":.+", "", "g", $14); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$12"\t"$14}}' > $dirOut/GE01_GE02.main.dbSNP.vcf
awk 'NR==FNR {id[$1]=$11"\t"$12; next} {if($1 in id){print $0"\t"id[$1]}}' $dirOut/GE01_GE02.main.dbSNP.vcf $dirOut/cortex01_cortex02.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $11); gt2=gensub(":.+", "", "g", $12); gt3=gensub(":.+", "", "g", $13); gt4=gensub(":.+", "", "g", $14); if(gt1==gt3 && gt2==gt4){print $0}}' > $dirOut/HuFNSC01_HuFNSC02.main.dbSNP.vcf
# 01 vs 03
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $11); gt2=gensub(":.+", "", "g", $16); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$16}}' > $dirOut/cortex01_cortex03.main.dbSNP.vcf
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $12); gt2=gensub(":.+", "", "g", $18); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$12"\t"$18}}' > $dirOut/GE01_GE03.main.dbSNP.vcf
awk 'NR==FNR {id[$1]=$11"\t"$12; next} {if($1 in id){print $0"\t"id[$1]}}' $dirOut/GE01_GE03.main.dbSNP.vcf $dirOut/cortex01_cortex03.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $11); gt2=gensub(":.+", "", "g", $12); gt3=gensub(":.+", "", "g", $13); gt4=gensub(":.+", "", "g", $14); if(gt1==gt3 && gt2==gt4){print $0}}' > $dirOut/HuFNSC01_HuFNSC03.main.dbSNP.vcf
# 01 vs 04
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $11); gt2=gensub(":.+", "", "g", $19); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$19}}' > $dirOut/cortex01_cortex04.main.dbSNP.vcf
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $12); gt2=gensub(":.+", "", "g", $20); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$12"\t"$20}}' > $dirOut/GE01_GE04.main.dbSNP.vcf
awk 'NR==FNR {id[$1]=$11"\t"$12; next} {if($1 in id){print $0"\t"id[$1]}}' $dirOut/GE01_GE04.main.dbSNP.vcf $dirOut/cortex01_cortex04.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $11); gt2=gensub(":.+", "", "g", $12); gt3=gensub(":.+", "", "g", $13); gt4=gensub(":.+", "", "g", $14); if(gt1==gt3 && gt2==gt4){print $0}}' > $dirOut/HuFNSC01_HuFNSC04.main.dbSNP.vcf
# 02 vs 03
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $13); gt2=gensub(":.+", "", "g", $16); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$13"\t"$16}}' > $dirOut/cortex02_cortex03.main.dbSNP.vcf
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $14); gt2=gensub(":.+", "", "g", $18); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$14"\t"$18}}' > $dirOut/GE02_GE03.main.dbSNP.vcf
awk 'NR==FNR {id[$1]=$11"\t"$12; next} {if($1 in id){print $0"\t"id[$1]}}' $dirOut/GE02_GE03.main.dbSNP.vcf $dirOut/cortex02_cortex03.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $11); gt2=gensub(":.+", "", "g", $12); gt3=gensub(":.+", "", "g", $13); gt4=gensub(":.+", "", "g", $14); if(gt1==gt3 && gt2==gt4){print $0}}' > $dirOut/HuFNSC02_HuFNSC03.main.dbSNP.vcf
# 02 vs 04
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $13); gt2=gensub(":.+", "", "g", $19); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$13"\t"$19}}' > $dirOut/cortex02_cortex04.main.dbSNP.vcf
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $14); gt2=gensub(":.+", "", "g", $20); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$14"\t"$20}}' > $dirOut/GE02_GE04.main.dbSNP.vcf
awk 'NR==FNR {id[$1]=$11"\t"$12; next} {if($1 in id){print $0"\t"id[$1]}}' $dirOut/GE02_GE04.main.dbSNP.vcf $dirOut/cortex02_cortex04.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $11); gt2=gensub(":.+", "", "g", $12); gt3=gensub(":.+", "", "g", $13); gt4=gensub(":.+", "", "g", $14); if(gt1==gt3 && gt2==gt4){print $0}}' > $dirOut/HuFNSC02_HuFNSC04.main.dbSNP.vcf
# 03 vs 04
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $16); gt2=gensub(":.+", "", "g", $19); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$16"\t"$19}}' > $dirOut/cortex03_cortex04.main.dbSNP.vcf
less $dirOut/$name.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $18); gt2=gensub(":.+", "", "g", $20); if(gt1!=gt2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$18"\t"$20}}' > $dirOut/GE03_GE04.main.dbSNP.vcf
awk 'NR==FNR {id[$1]=$11"\t"$12; next} {if($1 in id){print $0"\t"id[$1]}}' $dirOut/GE03_GE04.main.dbSNP.vcf $dirOut/cortex03_cortex04.main.dbSNP.vcf | awk '{gt1=gensub(":.+", "", "g", $11); gt2=gensub(":.+", "", "g", $12); gt3=gensub(":.+", "", "g", $13); gt4=gensub(":.+", "", "g", $14); if(gt1==gt3 && gt2==gt4){print $0}}' > $dirOut/HuFNSC03_HuFNSC04.main.dbSNP.vcf

# homozygotic sites
for file in HuFNSC*dbSNP.vcf; do
    name=$(echo $file | sed -e 's/.vcf//g')
    echo $name
    less $file | awk '{if(!(match($0, "0/1") || match($0, "1/0"))){print $0}}' > $name.homo.vcf
done
# summary 
wc -l Hu* | awk '{sample=gensub(".main.+", "", "g", $2); sample=gensub("HuFNSC0", "Subject", "g", sample); if($2 ~ /homo/){homo=$1; value=$1; cat="homozygotic"} else{value=$1-homo; cat="heterozygotic"} print sample"\t"cat"\t"value}' | awk '$1 ~ /Subject/ {print $0}' > SNP.summary

# MZ RNA-seq compare to MeDIP
cd /projects/epigenomics/users/lli/FetalBrain/SNP/
awk 'NR==FNR {id=$1":"$2; gt=gensub(":.*", "", "g", $10)"_"gensub(":.*", "", "g", $12); h[id]= gt; next} {id=$2":"$3; gt=gensub(":.*", "", "g", $11)"_"gensub(":.*", "", "g", $12); if(id in h && gt == h[id]){print $0}}' /home/lli/FetalBrain/MeDIPMRE/SNP/MZ_discordant_SNP.vcf HuFNSC01_HuFNSC02.main.dbSNP.vcf > MZ_discordant_SNP_MeDIP_RNAseq.vcf


