Pipeline shell scripts
=======================
* [WGBS.combine.sh](./WGBS.combine.sh): shell script to combine coverage from both strand for CpGs.
* [methyl_diff.sh](./methyl_diff.sh): shell script to identify DM CpGs for WGBS using methyl_diff.  
* [DMR.dynamic.sh](./DMR.dynamic.sh): shell script to collapse DM CpGs into DMRs.
* [DMR.intersect.sh](./DMR.intersect.sh): intersecting DMRs with genomic regions and report CpG% breakdown in intergenic, intron, exon, gene, promoter regions.     
* [region.intersect.sh](./region.intersect.sh): intersecting region of interest with genomic regions and report bp breakdown in intergenic, intron, exon, gene, promoter regions.     
* [region.mean.sh](./region.mean.sh): shell script to compute regional mean, sum or weighted mean.   
* [ER.unique.sh](./ER.unique.sh): shell script to identify histone modification ChIP-seq unique enriched regions for pairwise comparison.     
* [jaguar.sh](./jaguar.sh): shell script to calculate RPKM from fastq files using Jaguar alignment.      
* [bamstats2report.sh](./bamstats2report.sh): shell script to generate QC summary report from bamstats file.
* [bamstats2report_all.sh](./bamstats2report_all.sh): shell script to generate QC summary report for multiple libraries.
* [bamstats2report.combine.sh](./bamstats2report.combine.sh): shell script to combine multiple QC summary report into a single table.
* [RunB2W.sh](./RunB2W.sh): shell script to generate wig from bam. Run on apollo.
* [RunW2BW.sh](./RunW2BW.sh): shell script to generate BigWig from wig. Run on apollo.
* [snp_apollo.sh](./snp_apollo.sh): shell script to generate VCF files. Run on apollo.
* [SE_ROSE.sh](./SE_ROSE.sh): ROSE rank algorithm for superEnhnacers. 
