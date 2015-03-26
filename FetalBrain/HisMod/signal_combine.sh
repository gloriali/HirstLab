#!/bin/sh

cd /projects/epigenomics/users/lli/FetalBrain/ChIPseq/signal/
# H3K4me3: promoter (brain01 brain02 cortex02 GE02 GE04)
paste hg19v65_genes_TSS_2000.A03486.coverage hg19v65_genes_TSS_2000.A03494.coverage hg19v65_genes_TSS_2000.A03282.coverage hg19v65_genes_TSS_2000.A03478.coverage hg19v65_genes_TSS_2000.A19304.coverage | awk '{print $5"\t"$6"\t"$13"\t"$20"\t"$27"\t"$34}' > hg19v65_genes_TSS_2000.H3K4me3
less hg19v65_genes_TSS_2000.H3K4me3 | awk '/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.H3K4me3.pc
less hg19v65_genes_TSS_2000.H3K4me3 | awk '!/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.H3K4me3.nc
# H3K27me3: promoter (brain01 brain02 cortex01 cortex02 GE01 GE02 GE04)
paste hg19v65_genes_TSS_2000.A03488.coverage hg19v65_genes_TSS_2000.A03496.coverage hg19v65_genes_TSS_2000.A03272.coverage hg19v65_genes_TSS_2000.A03284.coverage hg19v65_genes_TSS_2000.A03278.coverage hg19v65_genes_TSS_2000.A03480.coverage hg19v65_genes_TSS_2000.A19306.coverage | awk '{print $5"\t"$6"\t"$13"\t"$20"\t"$27"\t"$34"\t"$41"\t"$48}' > hg19v65_genes_TSS_2000.H3K27me3
less hg19v65_genes_TSS_2000.H3K27me3 | awk '/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.H3K27me3.pc
less hg19v65_genes_TSS_2000.H3K27me3 | awk '!/protein_coding/ {print $0}' > hg19v65_genes_TSS_2000.H3K27me3.nc
# H3K36me3: genebody (brain01 brain02 cortex02 GE02 GE04)
paste hg19v65_genes_body.pc.A03489.coverage hg19v65_genes_body.pc.A03497.coverage hg19v65_genes_body.pc.A03273.coverage hg19v65_genes_body.pc.A03285.coverage hg19v65_genes_body.pc.A03279.coverage hg19v65_genes_body.pc.A03481.coverage hg19v65_genes_body.pc.A19307.coverage | awk '{print $4"\t"$5"\t"$11"\t"$17"\t"$23"\t"$29"\t"$35"\t"$41}' > hg19v65_genes_body.H3K36me3.pc
paste hg19v65_genes_body.nc.A03489.coverage hg19v65_genes_body.nc.A03497.coverage hg19v65_genes_body.nc.A03273.coverage hg19v65_genes_body.nc.A03285.coverage hg19v65_genes_body.nc.A03279.coverage hg19v65_genes_body.nc.A03481.coverage hg19v65_genes_body.pc.A19307.coverage | awk '{print $4"\t"$5"\t"$11"\t"$17"\t"$23"\t"$29"\t"$35"\t"$41}' > hg19v65_genes_body.H3K36me3.nc

