# Fetal Brain - Histone modifications
Gloria Li  
January 27, 2014  

Update Thu Mar 26 12:11:54 2015



## Sanity check   

* No. of peaks, No. of enriched bases, and average peak length seem reasonable except for the unusual high No. of peaks in GE HuFNSC04 input library.              

![](HisMod_files/figure-html/summary-1.png) 

## Correlation with protein-coding gene RPKM 

* Overlapping H3K4me3 and H3K27me3 FindER peaks with protein-coding gene promoters (TSS +/- 1500bp), and overlapping H3K36me3 with genebody.    
* Overall correlations are as expected, with H3K4m3 and H3K36me3 marked genes showing higher RPKM, and H3K27me3 and bivalent promoter marked showing lower RPKM.   

![](HisMod_files/figure-html/RPKM-1.png) 

## Differentially marked genes

* Calculate H3K4me3 and H3K27me3 signal from wig at protein-coding gene promoters (TSS +/- 2000bp) and normalize using quantile normalization.    
* As expected, much less flucturation of No. of differentially marked genes between samples.        
* A significant fraction of DE genes are differentially marked by either H3K4me3 or H3K27me3 in their promoters, more than DMRs. _Results are similar to normalize against total No. of reads although slightly smaller No. of DE genes_      

![](HisMod_files/figure-html/DM-1.png) ![](HisMod_files/figure-html/DM-2.png) 

## Overlapping enhancers with WGBS UMRs

* UMRs between neurospheres (cortex vs GE) are enriched in enhancers (H3K4me1 enriched regions).      
* Between gestational weeks, GW13 UMRs are enriched in enhancers, but GW17 UMRs are not.      
* UMRs overlaped with enhancers are highly enriched for brain development terms. For comparing between neurospheres, GE enhancer UMRs in HuFNSC04 have no significant enrichment. And for comparing between gestational weeks, GW13 enhancer UMRs in Cortex also have no enriched terms.       

![](HisMod_files/figure-html/UMR-1.png) 
![](HisMod_files/figure-html/UMR_GREAT-1.png) ![](HisMod_files/figure-html/UMR_GREAT-2.png) ![](HisMod_files/figure-html/UMR_GREAT-3.png) ![](HisMod_files/figure-html/UMR_GREAT-4.png) ![](HisMod_files/figure-html/UMR_GREAT-5.png) ![](HisMod_files/figure-html/UMR_GREAT-6.png) 


