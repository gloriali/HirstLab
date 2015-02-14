REMC Breast Project
=====================
* [BreastRNAseq.tsv](./BreastRNAseq.tsv): breast samples metadata.   
* [isoform](./isoform/): isoform analysis on breast samples.
* [junction](./junction/): isoform validation with junction coverage.
* [epiProfile](./epiProfile/): epigenetic signatures associated with exon/intron transcription.
* [RNA_yield.R](./RNA_yield.R): RNA yield asymmetry between lum and myo.
* [IR_DAVID.R](./IR_DAVID.R): DAVID enrichment for lum and myo intron retention events.
* [GREAT.R](./GREAT.R): GREAT enrichment for lum and myo UMRs.
* [lum_myo.DMR.sh](./lum_myo.DMR.sh): intersecting lum and myo UMRs with genomic regions.
* [TFBS_RPKM.R](./TFBS_RPKM.R): correlation between No. of TFBS and TF RPKM.
* [UMR.network.R](./UMR.network.R): TFBS overlapping UMRs that are associated with DE genes. Input for `cytoscape`.
* [UMR.network.sh](./UMR.network.sh): UMRs overlapping with GATA3, FOXA1, ZNF217, and EGR1, and UMR-associated genes.       
* [IS.DMR.R](./IS.DMR.R): visualization and GREAT enrichment on individual-specific UMRs.
* [IS.DMR.sh](./IS.DMR.sh): intersecting individual-specific UMRs with genomic regions and calculating chrX enrichment.
