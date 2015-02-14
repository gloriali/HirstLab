# intersect DMRs with promoters
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_brain.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRbrain_promoter.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_cortex.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRcortex_promoter.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_ge.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRge_promoter.txt

# intersect DMRs with CGIs
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_brain.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRbrain_CGI.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_cortex.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRcortex_CGI.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_ge.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRge_CGI.txt

/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_promoter_pc_brain.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRpromoter_pc_brain_CGI.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_promoter_pc_cortex.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRpromoter_pc_cortex_CGI.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_promoter_pc_ge.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRpromoter_pc_ge_CGI.txt


# intersect DMRs with genebody
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_brain.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRbrain_genebody.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_cortex.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRcortex_genebody.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMRs/DMR_ge.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMRs/intersect/DMRge_genebody.txt


##################################################
# intersect cortex vs GE DMRs with promoter regions
/Applications/bedtools-2.17.0/bin/intersectBed -a ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE.bed -b ~/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_promoter_cortexGE.txt
/Applications/bedtools-2.17.0/bin/intersectBed -a ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE01.bed -b ~/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_promoter_cortexGE01.txt
/Applications/bedtools-2.17.0/bin/intersectBed -a ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE02.bed -b ~/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > ~/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_promoter_cortexGE02.txt

# intersect DMRs with CGIs
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_CGI_cortexGE.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE01.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_CGI_cortexGE01.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE02.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_CGI_cortexGE02.txt

# intersect DMRs with genebody
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_genebody_cortexGE.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE01.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_genebody_cortexGE01.txt
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_cortexGE02.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_cortexVsGE/DMR_genebody_cortexGE02.txt


##################################################
# intersect dynamic growth DMRs brain01 vs 02 with promoter regions
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic.bed -b /home/lli/hg19/hg19v65_genes_TSS_1500.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic_promoter.txt

# genebody
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic.bed -b /home/lli/hg19/hg19v65_genes.bed -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic_genebody.txt

# CGIs
/home/lli/bedtools-2.17.0/bin/intersectBed -a /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic.bed -b /home/lli/hg19/CGI_hg19v65_segregated.txt -wa -wb > /home/lli/FetalBrain/MeDIPMRE/DMR_brain12_dynamic_CGI.txt
