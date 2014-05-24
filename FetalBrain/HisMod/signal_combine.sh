#!/bin/sh

cd /home/lli/FetalBrain/HisMod/signal/
# H3K4me1: promoter
paste hg19v65_genes_TSS_2000.A03485.bam.q5.F1028.SET_190.coverage hg19v65_genes_TSS_2000.A03493.bam.q5.F1028.SET_243.coverage hg19v65_genes_TSS_2000.A03269.bam.q5.F1028.SET_184.coverage hg19v65_genes_TSS_2000.A03281.bam.q5.F1028.SET_240.coverage hg19v65_genes_TSS_2000.A03275.bam.q5.F1028.SET_204.coverage hg19v65_genes_TSS_2000.A03477.bam.q5.F1028.SET_232.coverage | awk '{print $4"\t"$6"\t"$13"\t"$20"\t"$27"\t"$34"\t"$41}' > hg19v65_genes_TSS_2000.H3K4me1
