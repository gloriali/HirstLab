Analysis pipeline
====================
+ [Path.md](./Path.md): file paths on xhost.
+ [LIMS_API_Sita.md](./LIMS_API_Sita.md): API access to LIMS.
+ [instructions](./instructions): Pipeline instructions.  
    * [DEfine.md](./instructions/DEfine.md): instructions and sample code for using DEfine on genes or exons.  
    * [isoform.md](./instructions/isoform.md): instructions for identifying isoform genes and cassette exons.  
    * [junction.md](./instructions/junction.md): instructions for validating isoforms with junction reads.  
    * [epiProfile.md](./instructions/epiProfile.md): instructions for calculating and plotting epigenetic profiles for isoform and intron retention.       
    * [MeDIP.DMR.md](./instructions/MeDIP.DMR.md): instructions for DMR analysis from MeDIP fractional methylation calls.        
    * [ChromHMM.md](./instructions/ChromHMM.md): instructions on running ChromHMM. _unfinished_.
+ [shell](./shell): shell scripts
    * [WGBS.combine.sh](./shell/WGBS.combine.sh): shell script to combine coverage from both strand for CpGs.
    * [methyl_diff.sh](./shell/methyl_diff.sh): shell script to identify DM CpGs for WGBS using methyl_diff.  
    * [DMR.dynamic.sh](./shell/DMR.dynamic.sh): shell script to collapse DM CpGs into DMRs.
    * [DMR.intersect.sh](./shell/DMR.intersect.sh): intersecting DMRs with genomic regions and report CpG% breakdown in intergenic, intron, exon, gene, promoter regions.     
+ [R](./R): R functions
    * [DMR.figures.R](./R/DMR.figures.R): `DMR_figures` R function for DMR analysis and visualization.
    * [DMR_DE.R](./R/DMR_DE.R): `DMR_DE` R function for identify DE genes for gene-associated DMRs.  
    * [enrich.R](./R/enrich.R): `enrich` R function to plot enriched terms with `ggplot2`.   
    * [enrich_GREAT.R](./R/enrich_GREAT.R): `enrich_GREAT` R function to plot GREAT enriched terms with `ggplot2`.   
    * [isoform.R](./R/isoform.R): `isoform` R function for identifying isoform genes and cassette exons.
    * [junction.R](./R/junction.R): `junction` R function for validating isoform genes and cassette exons with junction reads.   
    * [epiProfile.R](./R/epiProfile.R): `epiProfile` R function for calculating and plotting epigenetic profiles for isoform and intron retention.       
+ [python](./python): python scripts
    * [MeDIP.DMR.py](./python/MeDIP.DMR.py): collapse single DM CpGs to DMRs.
+ [matlab](./matlab): copy of Misha's matlab functions
    * [DEfine.m](./matlab/DEfine.m): copy of Misha's DEfine function.   
    * [myfnorm.m](./matlab/myfnorm.m): copy of Misha's myfnorm function, required for DEfine.
+ [UCSC](./UCSC): script from UCSC genome browser and sample code for building track hub.

