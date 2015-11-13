#!/bin/sh

# Total No. of non CpG Cs
## less /projects/production_genomics/bs-seq/a22475/non_CpG_context_A22475_4_lanes_dupsFlagged_readsorted.* | wc -l

# non-zero non CpG methylation
dirOut=/projects/epigenomics/users/lli/FetalBrain/WGBS/nonCG/
mkdir -p $dirOut

less /projects/production_genomics/bs-seq/a22475/non_CpG_context_A22475_4_lanes_dupsFlagged_readsorted.bismark.cov.gz | awk '{if($4!=0){print $0}}' > $dirOut/A22475_Cortex02_nonCG_non0.cov
less /projects/production_genomics/bs-seq/joc163/non_CpG_context_JOC163_4_lanes_dupsFlagged_readsorted.bismark.cov.gz | awk '{if($4!=0){print $0}}' > $dirOut/A17784_A13819_GE02_nonCG_non0.cov
less /projects/production_genomics/bs-seq/a22477/non_CpG_context_A22477_5_lanes_dupsFlagged_readsorted.bismark.cov | awk '{if($4!=0){print $0}}' > $dirOut/A22477_Cortex04_nonCG_non0.cov
less /projects/production_genomics/bs-seq/a22476/non_CpG_context__A22476_4_lanes_dupsFlagged_readsorted.bismark.cov.gz | awk '{if($4!=0){print $0}}' > $dirOut/A22476_GE04_nonCG_non0.cov

