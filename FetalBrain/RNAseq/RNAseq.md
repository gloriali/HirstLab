Fetal Brain RNA-seq Analysis Summary
========================================================

Gloria Li         
Mon Sep  8 22:39:00 2014 



## Differential gene expression 
### DEfine

  * FDR = 0.01    
  * Minimum sum of RPKM (rmin) = 0.005    
  * Minimum sum of coverage (Nmin) = 25    
  
### Between cortex and GE neurospheres

  * On average, there are __860__ genes differentially expressed between cortex and GE, among them, __454__ are upregulated in cortex, and __406__ are downregulated.    
  * __886__ DE genes are shared by at least two individuals.    
  * DAVID enrichment analysis show significant enrichment in __neuronal development__ and __cell migration__ terms, as well as __EGF-related__ protein domains. 

<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Mon Sep  8 22:39:37 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> UP </TH> <TH> DN </TH> <TH> DE </TH>  </TR>
  <TR> <TD> HuFNSC01 </TD> <TD align="center"> 403 </TD> <TD align="center"> 508 </TD> <TD align="center"> 911 </TD> </TR>
  <TR> <TD> HuFNSC02 </TD> <TD align="center"> 588 </TD> <TD align="center"> 640 </TD> <TD align="center"> 1228 </TD> </TR>
  <TR> <TD> HuFNSC03 </TD> <TD align="center"> 447 </TD> <TD align="center"> 227 </TD> <TD align="center"> 674 </TD> </TR>
  <TR> <TD> HuFNSC04 </TD> <TD align="center"> 378 </TD> <TD align="center"> 249 </TD> <TD align="center"> 627 </TD> </TR>
   </TABLE>
![plot of chunk cortex_GE](./RNAseq_files/figure-html/cortex_GE.png) 

![plot of chunk cortex_GE_enrich](./RNAseq_files/figure-html/cortex_GE_enrich.png) 

### Between MZ twins - HuFNSC01 vs HuFNSC02

  * On average, there are __470__ DE genes across three cells types.   
  * Majority of DE genes are cell type specific, only __98__ are shared between any two cell types.   
  * DAVID enrichment analysis between MZ twins in brain and cortex show similar GO term in __brain development__, but there is no significantly enriched terms in GE.    

<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Mon Sep  8 22:39:38 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> UP </TH> <TH> DN </TH> <TH> DE </TH>  </TR>
  <TR> <TD> brain01_brain02 </TD> <TD align="center"> 461 </TD> <TD align="center"> 181 </TD> <TD align="center"> 642 </TD> </TR>
  <TR> <TD> cortex01_cortex02 </TD> <TD align="center"> 248 </TD> <TD align="center"> 348 </TD> <TD align="center"> 596 </TD> </TR>
  <TR> <TD> GE01_GE02 </TD> <TD align="center">  99 </TD> <TD align="center">  74 </TD> <TD align="center"> 173 </TD> </TR>
  <TR> <TD> cortex03_GE03 </TD> <TD align="center"> 459 </TD> <TD align="center"> 183 </TD> <TD align="center"> 642 </TD> </TR>
  <TR> <TD> cortex04_GE04 </TD> <TD align="center"> 313 </TD> <TD align="center"> 232 </TD> <TD align="center"> 545 </TD> </TR>
   </TABLE>
![plot of chunk individual](./RNAseq_files/figure-html/individual1.png) ![plot of chunk individual](./RNAseq_files/figure-html/individual2.png) 

![plot of chunk individual_enrich](./RNAseq_files/figure-html/individual_enrich1.png) ![plot of chunk individual_enrich](./RNAseq_files/figure-html/individual_enrich2.png) 


## Isoform analysis
### Isoform identification and junction validation  
  * DEfine on exons: FDR = 0.01     
  * Exon expressed in one sample ($\ge$ 10% gene RPKM) and not expressed in the other ($\le$ 1% gene RPKM)   
  * Gene is not DE: DEfine FDR = 0.01
  * Gene is expressed in both samples: gene RPKM > 0.01         
  * Validation: For each isoform exon in the previous pairwise comparison
    + Find junctions associated with this exon with enough coverage, i.e. sum of junction coverage of two samples $\ge$ 1
    + Identify junctions that RPKM change in the same direction as the exon
    + Junction RPKM > 0.1 in one sample and < 0.1 in the other      

#### Between cortex and GE neurospheres

  * On average, __2054__ genes are identified as isoforms between cortex and GE in each individual. __2352__ genes are shared in at least two individuals.       
  * There are more individual-specific isoforms than found in breast cells, although the overlap across individuals are still significant.     
  * Individual specific isoforms between cortex and GE have __no__ significantly enriched terms, suggesting they are more likely random events without biological functions.          
  * Isoforms shared by at least two individuals are enriched in terms related to __cellular signaling__. InterPro protein domain enrichment show enriched terms similar to those observed in breast isoforms.         

<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Mon Sep  8 22:39:40 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> DE_genes </TH> <TH> DE_exons </TH> <TH> with_expressed_genes </TH> <TH> isoform_exons </TH> <TH> exclude_DE_genes </TH> <TH> isoform_genes </TH>  </TR>
  <TR> <TD> cortex01_GE01_summary </TD> <TD align="center"> 911 </TD> <TD align="center"> 32372 </TD> <TD align="center"> 18968 </TD> <TD align="center"> 8440 </TD> <TD align="center"> 7962 </TD> <TD align="center"> 2447 </TD> </TR>
  <TR> <TD> cortex02_GE02_summary </TD> <TD align="center"> 1228 </TD> <TD align="center"> 35196 </TD> <TD align="center"> 21880 </TD> <TD align="center"> 8163 </TD> <TD align="center"> 7374 </TD> <TD align="center"> 2298 </TD> </TR>
  <TR> <TD> cortex03_GE03_summary </TD> <TD align="center"> 674 </TD> <TD align="center"> 29617 </TD> <TD align="center"> 13746 </TD> <TD align="center"> 6401 </TD> <TD align="center"> 6022 </TD> <TD align="center"> 2086 </TD> </TR>
  <TR> <TD> cortex04_GE04_summary </TD> <TD align="center"> 627 </TD> <TD align="center"> 22386 </TD> <TD align="center"> 11253 </TD> <TD align="center"> 4323 </TD> <TD align="center"> 4259 </TD> <TD align="center"> 1386 </TD> </TR>
   </TABLE>
![plot of chunk isoform_cortex_ge](./RNAseq_files/figure-html/isoform_cortex_ge.png) 

![plot of chunk isoform_enrich_cortex_ge](./RNAseq_files/figure-html/isoform_enrich_cortex_ge.png) 

#### Between MZ twins - HuFNSC01 vs HuFNSC02

  * On average, __2617__ genes are identified as isoforms between HuFNSC01 and HuFNSC02 in each cell type. __796__ genes are shared by all three cell types.              
  * On average, __1724__ genes are identified as isoforms between HuFNSC03 and HuFNSC04 in each cell type. __927__ genes are shared between two cell types.            
  * Different regions on the Venn diagram have __no__ significantly enriched terms.     
  * Isoforms between HuFNSC01 and HuFNSC02 in neurospheres show similar terms, related to __cell signaling__, and __blood cell development__ in brain.     

<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Mon Sep  8 22:39:40 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> DE_genes </TH> <TH> DE_exons </TH> <TH> with_expressed_genes </TH> <TH> isoform_exons </TH> <TH> exclude_DE_genes </TH> <TH> isoform_genes </TH>  </TR>
  <TR> <TD> brain01_brain02_summary </TD> <TD align="center"> 642 </TD> <TD align="center"> 32138 </TD> <TD align="center"> 16302 </TD> <TD align="center"> 8980 </TD> <TD align="center"> 8542 </TD> <TD align="center"> 2902 </TD> </TR>
  <TR> <TD> cortex01_cortex02_summary </TD> <TD align="center"> 596 </TD> <TD align="center"> 26983 </TD> <TD align="center"> 15554 </TD> <TD align="center"> 7618 </TD> <TD align="center"> 7445 </TD> <TD align="center"> 2454 </TD> </TR>
  <TR> <TD> GE01_GE02_summary </TD> <TD align="center"> 173 </TD> <TD align="center"> 23810 </TD> <TD align="center"> 12862 </TD> <TD align="center"> 7402 </TD> <TD align="center"> 7351 </TD> <TD align="center"> 2495 </TD> </TR>
  <TR> <TD> cortex03_cortex04_summary </TD> <TD align="center"> 642 </TD> <TD align="center"> 26826 </TD> <TD align="center"> 12185 </TD> <TD align="center"> 5818 </TD> <TD align="center"> 5479 </TD> <TD align="center"> 1994 </TD> </TR>
  <TR> <TD> GE03_GE04_summary </TD> <TD align="center"> 545 </TD> <TD align="center"> 24752 </TD> <TD align="center"> 12223 </TD> <TD align="center"> 4582 </TD> <TD align="center"> 4422 </TD> <TD align="center"> 1454 </TD> </TR>
   </TABLE>
![plot of chunk isoform_HuFNSC01_02](./RNAseq_files/figure-html/isoform_HuFNSC01_021.png) ![plot of chunk isoform_HuFNSC01_02](./RNAseq_files/figure-html/isoform_HuFNSC01_022.png) 

![plot of chunk isoform_enrich_HuFNSC01_02](./RNAseq_files/figure-html/isoform_enrich_HuFNSC01_021.png) ![plot of chunk isoform_enrich_HuFNSC01_02](./RNAseq_files/figure-html/isoform_enrich_HuFNSC01_022.png) ![plot of chunk isoform_enrich_HuFNSC01_02](./RNAseq_files/figure-html/isoform_enrich_HuFNSC01_023.png) 

#### Junction validation
  + For cortex vs GE, on average, __36.4%__ isoform genes have enough junction coverage. Among them, __89.8%__ have support from junction reads.    
  + For comparisons between individuals, on average, __34.4%__ isoform genes have enough junction coverage. Among them, __87.7%__ have support from junction reads.     
  + Between strand specific and non-strand specific libraries, the percentage of isoforms with enough junction coverage are __similar__, however, strand specific libraries have __higher__ percentage of having junction support (> 98% compared to ~ 80%).   

<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Mon Sep  8 22:39:42 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> isoform exons </TH> <TH> isoform genes </TH> <TH> exons with junction coverage </TH> <TH> genes with junction coverage </TH> <TH> exons with junction support </TH> <TH> genes with junction support </TH>  </TR>
  <TR> <TD> cortex01_GE01_summary </TD> <TD align="center"> 7962 </TD> <TD align="center"> 2447 </TD> <TD align="center"> 1900 </TD> <TD align="center"> 832 </TD> <TD align="center"> 1543 </TD> <TD align="center"> 697 </TD> </TR>
  <TR> <TD> cortex02_GE02_summary </TD> <TD align="center"> 7374 </TD> <TD align="center"> 2298 </TD> <TD align="center"> 1962 </TD> <TD align="center"> 801 </TD> <TD align="center"> 1554 </TD> <TD align="center"> 656 </TD> </TR>
  <TR> <TD> cortex03_GE03_summary </TD> <TD align="center"> 6022 </TD> <TD align="center"> 2086 </TD> <TD align="center"> 1972 </TD> <TD align="center"> 872 </TD> <TD align="center"> 1836 </TD> <TD align="center"> 824 </TD> </TR>
  <TR> <TD> cortex04_GE04_summary </TD> <TD align="center"> 4259 </TD> <TD align="center"> 1386 </TD> <TD align="center"> 941 </TD> <TD align="center"> 483 </TD> <TD align="center"> 935 </TD> <TD align="center"> 478 </TD> </TR>
   </TABLE>
<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Mon Sep  8 22:39:42 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> isoform exons </TH> <TH> isoform genes </TH> <TH> exons with junction coverage </TH> <TH> genes with junction coverage </TH> <TH> exons with junction support </TH> <TH> genes with junction support </TH>  </TR>
  <TR> <TD> brain01_brain02_summary </TD> <TD align="center"> 8542 </TD> <TD align="center"> 2902 </TD> <TD align="center"> 1791 </TD> <TD align="center"> 921 </TD> <TD align="center"> 1520 </TD> <TD align="center"> 808 </TD> </TR>
  <TR> <TD> cortex01_cortex02_summary </TD> <TD align="center"> 7445 </TD> <TD align="center"> 2454 </TD> <TD align="center"> 1634 </TD> <TD align="center"> 765 </TD> <TD align="center"> 1221 </TD> <TD align="center"> 623 </TD> </TR>
  <TR> <TD> GE01_GE02_summary </TD> <TD align="center"> 7351 </TD> <TD align="center"> 2495 </TD> <TD align="center"> 1469 </TD> <TD align="center"> 737 </TD> <TD align="center"> 1030 </TD> <TD align="center"> 561 </TD> </TR>
  <TR> <TD> cortex03_cortex04_summary </TD> <TD align="center"> 5479 </TD> <TD align="center"> 1994 </TD> <TD align="center"> 1949 </TD> <TD align="center"> 883 </TD> <TD align="center"> 1825 </TD> <TD align="center"> 830 </TD> </TR>
  <TR> <TD> GE03_GE04_summary </TD> <TD align="center"> 4422 </TD> <TD align="center"> 1454 </TD> <TD align="center"> 1045 </TD> <TD align="center"> 515 </TD> <TD align="center"> 1038 </TD> <TD align="center"> 511 </TD> </TR>
   </TABLE>

### No. of exons for DE genes / isoform genes    
  * DE genes have roughly the same No. of exons as all expressed genes.             
  * Identified isoforms have slightly more No. of exons than DE genes and all expressed genes.  
  * Compared to DE genes, the distribution in No. of exons for isoforms are __much similar__ between different individuals, _not observed in breast libraries_.    
  
![plot of chunk isoform_Nexon](./RNAseq_files/figure-html/isoform_Nexon.png) 

### Position of isoform exons on the gene   
  * In general, there are more alternative spliced exons at the __two ends__ of genes, _similar to observed in breast libraries_.         
  
![plot of chunk isoform_exon_pos](./RNAseq_files/figure-html/isoform_exon_pos.png) 

### Venn Diagram with average expression level, average No. of exons and average exon length   
  * Isoforms have __much lower__ expression level than all expressed genes.          
  * On average, common isoforms between different comparisons have __lower expression level__ than comparison-specific isoforms.                 
  * In general, compared to all isoforms identified, validated isoforms have __lower__ expression levels, _not observed in breast libraries_.     
  * Average No. of exons are very __similar__ in different sections of the Venn diagram, between all, validated isoforms and all expressed genes.        
  * Average length of isoform exons are __shorter__ than all expressed genes. Validated isoform exons are also __shorter__ than all isoforms in general, _not observed in breast libraries_.         
  
![plot of chunk isoform_venn](./RNAseq_files/figure-html/isoform_venn1.png) ![plot of chunk isoform_venn](./RNAseq_files/figure-html/isoform_venn2.png) ![plot of chunk isoform_venn](./RNAseq_files/figure-html/isoform_venn3.png) ![plot of chunk isoform_venn](./RNAseq_files/figure-html/isoform_venn4.png) 


## Epigenetic signature marking 
### Between cortex and GE neurospheres
#### DNA methylation at exon boundaries 
* 5mC at cassette exon boundaries has similar pattern to expressed in both exons.   
* 5mC for cassette exons between cortex and GE neurospheres shows no significant differences.   
* Results from both WGBS (HuFNSC02 & HuFNSC04) and MeDIP (HuFNSC01 & HuFNSC02) support the assumption that 5mC is a stable mark for exon transcription during development.   
* 5mC exon marking is established between neurospheres. _Needs further validation against H1_.   

![plot of chunk epiProfile_5mC_cortexge](./RNAseq_files/figure-html/epiProfile_5mC_cortexge1.png) ![plot of chunk epiProfile_5mC_cortexge](./RNAseq_files/figure-html/epiProfile_5mC_cortexge2.png) ![plot of chunk epiProfile_5mC_cortexge](./RNAseq_files/figure-html/epiProfile_5mC_cortexge3.png) ![plot of chunk epiProfile_5mC_cortexge](./RNAseq_files/figure-html/epiProfile_5mC_cortexge4.png) 

#### H3K36me3 in exon bodies
* H3K36me3 in expressed in both / not expressed exons shows no significant differences between HuFNSC01 and HuFNSC02.    
* H3K36me3 in cassette exons in HuFNSC01 are enriched in GE compared to cortex. However, it is not reproduced in HuFNSC02, where there is no significant differences between cortex and GE. _Not sure what to make of this. Are there any potential bias?_    

![plot of chunk epiProfile_H3K36me3_cortexge](./RNAseq_files/figure-html/epiProfile_H3K36me3_cortexge1.png) ![plot of chunk epiProfile_H3K36me3_cortexge](./RNAseq_files/figure-html/epiProfile_H3K36me3_cortexge2.png) 

### Between MZ twins - HuFNSC01 vs HuFNSC02
#### DNA methylation at exon boundaries 
* 5mC for cassette exons between HuFNSC01 and HuFNSC02 shows no significant differences.   
* There are significant differences in 5mC between HuFNSC01 and HuFNSC02 specific exons in all three cell types.   
* In brain, 5mC in HuFNSC02 specific exons are closer to expressed in both, and HuFNSC01 specific exons are closer to not expressed exons. However, we observe the opposite trend in cortex and GE. _Is this a reflection of developmental stages differences between MZ twins? Are the opposite trends between brain and cortex/GE because of cell culture? Or could all these differences be technical bias / noise?_    

![plot of chunk epiProfile_5mC_MZ](./RNAseq_files/figure-html/epiProfile_5mC_MZ1.png) ![plot of chunk epiProfile_5mC_MZ](./RNAseq_files/figure-html/epiProfile_5mC_MZ2.png) ![plot of chunk epiProfile_5mC_MZ](./RNAseq_files/figure-html/epiProfile_5mC_MZ3.png) 

#### H3K36me3 in exon bodies
* There are some differences between HuFNSC01 and HuFNSC02 in expressed in both and not expressed exons, less so in GE.   
* For cassette exons, HuFNSC02 have enriched H3K36me3 compared to HuFNSC01 in both brain and cortex, but __not__ in GE. _Could this be related to the asymmetry we observed in DNA methylation in brain and cortex?_   

![plot of chunk epiProfile_H3K36me3_MZ](./RNAseq_files/figure-html/epiProfile_H3K36me3_MZ1.png) ![plot of chunk epiProfile_H3K36me3_MZ](./RNAseq_files/figure-html/epiProfile_H3K36me3_MZ2.png) ![plot of chunk epiProfile_H3K36me3_MZ](./RNAseq_files/figure-html/epiProfile_H3K36me3_MZ3.png) 










