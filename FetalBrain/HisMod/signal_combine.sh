#!/bin/sh

cd /home/lli/FetalBrain/HisMod/signal/
# H3K4me1: promoter (brain01 brain02 cortex01 cortex02 GE01 GE02)
paste hg19v65_genes_TSS_2000.A03485.bam.q5.F1028.SET_190.coverage hg19v65_genes_TSS_2000.A03493.bam.q5.F1028.SET_243.coverage hg19v65_genes_TSS_2000.A03269.bam.q5.F1028.SET_184.coverage hg19v65_genes_TSS_2000.A03281.bam.q5.F1028.SET_240.coverage hg19v65_genes_TSS_2000.A03275.bam.q5.F1028.SET_204.coverage hg19v65_genes_TSS_2000.A03477.bam.q5.F1028.SET_232.coverage | awk '{print $4"\t"$6"\t"$13"\t"$20"\t"$27"\t"$34"\t"$41}' > hg19v65_genes_TSS_2000.H3K4me1
# H3K4me3: promoter (brain01 brain02 cortex02 GE02)
paste hg19v65_genes_TSS_2000.A03486.bam.q5.F1028.SET_195.coverage hg19v65_genes_TSS_2000.A03494.bam.q5.F1028.SET_257.coverage hg19v65_genes_TSS_2000.A03282.bam.q5.F1028.SET_207.coverage hg19v65_genes_TSS_2000.A03478.bam.q5.F1028.SET_169.coverage | awk '{print $4"\t"$6"\t"$13"\t"$20"\t"$27}' > hg19v65_genes_TSS_2000.H3K4me3
# H3K27me3: promoter (brain01 brain02 cortex01 cortex02 GE01 GE02)
paste hg19v65_genes_TSS_2000.A03488.bam.q5.F1028.SET_172.coverage hg19v65_genes_TSS_2000.A03496.bam.q5.F1028.SET_264.coverage hg19v65_genes_TSS_2000.A03272.bam.q5.F1028.SET_149.coverage hg19v65_genes_TSS_2000.A03284.bam.q5.F1028.SET_194.coverage hg19v65_genes_TSS_2000.A03278.bam.q5.F1028.SET_199.coverage hg19v65_genes_TSS_2000.A03480.bam.q5.F1028.SET_211.coverage | awk '{print $4"\t"$6"\t"$13"\t"$20"\t"$27"\t"$34"\t"$41}' > hg19v65_genes_TSS_2000.H3K27me3
# H3K36me3: pc genebody (brain01 brain02 cortex02 GE02)
paste hg19v65_genes_body.pc.A03489.bam.q5.F1028.SET_175.coverage hg19v65_genes_body.pc.A03497.bam.q5.F1028.SET_233.coverage hg19v65_genes_body.pc.A03273.bam.q5.F1028.SET_219.coverage hg19v65_genes_body.pc.A03285.bam.q5.F1028.SET_238.coverage hg19v65_genes_body.pc.A03279.bam.q5.F1028.SET_229.coverage hg19v65_genes_body.pc.A03481.bam.q5.F1028.SET_215.coverage | awk '{print $4"\t"$5"\t"$11"\t"$17"\t"$23"\t"$29"\t"$35}' > hg19v65_genes_body.H3K36me3.pc
# H3K36me3: nc genebody (brain01 brain02 cortex02 GE02)
paste hg19v65_genes_body.nc.A03489.bam.q5.F1028.SET_175.coverage hg19v65_genes_body.nc.A03497.bam.q5.F1028.SET_233.coverage hg19v65_genes_body.nc.A03273.bam.q5.F1028.SET_219.coverage hg19v65_genes_body.nc.A03285.bam.q5.F1028.SET_238.coverage hg19v65_genes_body.nc.A03279.bam.q5.F1028.SET_229.coverage hg19v65_genes_body.nc.A03481.bam.q5.F1028.SET_215.coverage | awk '{print $4"\t"$5"\t"$11"\t"$17"\t"$23"\t"$29"\t"$35}' > hg19v65_genes_body.H3K36me3.nc
# input: promoter (brain01 brain02 cortex02 GE02)
paste hg19v65_genes_TSS_2000.A03491.bam.q5.F1028.SET_176.coverage hg19v65_genes_TSS_2000.A03499.bam.q5.F1028.SET_130.coverage hg19v65_genes_TSS_2000.A03287.bam.q5.F1028.SET_194.coverage hg19v65_genes_TSS_2000.A03289.bam.q5.F1028.SET_159.coverage hg19v65_genes_TSS_2000.A03288.bam.q5.F1028.SET_174.coverage hg19v65_genes_TSS_2000.A03483.bam.q5.F1028.SET_143.coverage | awk '{print $4"\t"$6"\t"$13"\t"$20"\t"$27"\t"$34"\t"$41}' > hg19v65_genes_TSS_2000.input

# separate pc and nc genes for promoter signal
less hg19v65_genes_TSS_2000.H3K4me1 | awk '/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.H3K4me1.pc
less hg19v65_genes_TSS_2000.H3K4me1 | awk '!/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.H3K4me1.nc
rm hg19v65_genes_TSS_2000.H3K4me1
less hg19v65_genes_TSS_2000.H3K4me3 | awk '/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.H3K4me3.pc
less hg19v65_genes_TSS_2000.H3K4me3 | awk '!/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.H3K4me3.nc
rm hg19v65_genes_TSS_2000.H3K4me3
less hg19v65_genes_TSS_2000.H3K27me3 | awk '/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.H3K27me3.pc
less hg19v65_genes_TSS_2000.H3K27me3 | awk '!/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.H3K27me3.nc
rm hg19v65_genes_TSS_2000.H3K27me3
less hg19v65_genes_TSS_2000.input | awk '/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.input.pc
less hg19v65_genes_TSS_2000.input | awk '!/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.input.nc
rm hg19v65_genes_TSS_2000.input
