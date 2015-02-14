# /home/lli/MeDIPMRE/
# brain01 brain02
/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools mpileup -C50 -uf /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa HS2788.bam HS2790.bam | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view -bvcg - > ./SNP/HS2788_HS2790.bcf
/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view HS2788_HS2790.bcf | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/vcfutils.pl varFilter -D100 > HS2788_HS2790.vcf

# cortex01 cortex02
/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools mpileup -C50 -uf /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa HS2775.bam HS2779.bam | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view -bvcg - > ./SNP/HS2775_HS2779.bcf
/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view HS2775_HS2779.bcf | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/vcfutils.pl varFilter -D100 > HS2775_HS2779.vcf

# ge01 ge02
/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools mpileup -C50 -uf /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa HS2777.bam HS2781.bam | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view -bvcg - > ./SNP/HS2777_HS2781.bcf
/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view HS2777_HS2781.bcf | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/vcfutils.pl varFilter -D100 > HS2777_HS2781.vcf

# myo35 myo70
/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools mpileup -C50 -uf /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa HS1394.bam HS2294.bam | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view -bvcg - > HS1394_HS2294.bcf
/home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/bcftools view HS1394_HS2294.bcf | /home/pubseq/BioSw/samtools/samtools-0.1.16/bcftools/vcfutils.pl varFilter -D100 > HS1394_HS2294.vcf
