REMC isoform analysis - junctions
========================================================

Gloria Li         
Mon May  5 11:01:16 2014 

<!-- re-knit after modify junction.R or junction_valid.R script -->






## Junction coverage and RPKM distribution
<img src="readCount_distribution.png" height="400px" width="400px" />
<img src="RPKM_distribution.png" height="400px" width="400px" />           

## Validate previously identified isoforms with junction RPKM
  * Previous isoform identification with DE exons
    * DE exons by DEfine FDR = 0.015
    * Exon RPKM $\ge$ 10% gene RPKM in one sample & $\le$ 1% in the other
    * Gene RPKM of both samples > 0.005 
    * Exclude DE genes by DEfine FDR = 0.015

  * Validation: For each isoform exon in the previous pairwise comparison
    * Find junctions associated with this exon with enough coverage, i.e. sum of junction coverage of two samples $\ge$ 2
    * Identify junctions that RPKM change in the same direction as the exon
    * Junction RPKM > 0.1 (i.e. coverage $\ge$ 2) in one sample and < 0.1 in the other (i.e. coverage = 0 or 1)
  
## Results:  

  * Only about __40%__ of previously identified exons (__50%__ genes) have junctions with enough coverage for validation. However, __more than 98%__ of them have junction support. 
  
<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Mon May  5 11:01:21 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> No.isoform.exons </TH> <TH> No.isoform.genes </TH> <TH> No.exons.with.junction.cov </TH> <TH> No.genes.with.junction.cov </TH> <TH> No.exons.with.junction.support </TH> <TH> No.genes.with.junction.support </TH>  </TR>
  <TR> <TD> lum084_myo084 </TD> <TD align="center"> 8630 </TD> <TD align="center"> 2381 </TD> <TD align="center"> 3210 </TD> <TD align="center"> 1134 </TD> <TD align="center"> 3166 </TD> <TD align="center"> 1114 </TD> </TR>
  <TR> <TD> lum084_stem084 </TD> <TD align="center"> 8948 </TD> <TD align="center"> 2429 </TD> <TD align="center"> 3793 </TD> <TD align="center"> 1287 </TD> <TD align="center"> 3744 </TD> <TD align="center"> 1270 </TD> </TR>
  <TR> <TD> myo084_stem084 </TD> <TD align="center"> 8619 </TD> <TD align="center"> 2427 </TD> <TD align="center"> 2938 </TD> <TD align="center"> 1138 </TD> <TD align="center"> 2881 </TD> <TD align="center"> 1111 </TD> </TR>
   </TABLE>

  
### No. of exons for DE genes / isoform genes    
  * DE genes have roughly the same No. of exons as all expressed genes.             
  * Identified isoforms have slightly more No. of exons than DE genes and all expressed genes. 
  
![plot of chunk Nexon](figure/Nexon.png) 

  
### Isoform genes in DE genes with small No. of exons

  If a gene has only a few exons ($\le$ 5 exons), and one exon is absent in one sample, i.e. isoform, the absence can bias the overall gene RPKM and this gene may be identified as a DE gene. Check if there are such cases.               
  For DE gene with small No. of exons:           
  * If it is in fact an isoform, not DE gene, only one / few exons in this gene should be differentially expressed.         
  * If it is a truly DE gene, all exons should be differentially expressed.       
                  
  Proportion of DE exons for DE genes with $\le$ 5 exons: (DE exons: fold change $\ge$ 2)      
  * Proportion of DE exons < 1 for the few genes are due to exons not expressed in either samples.         
  * __There is no evidence for isoforms genes identified as DE genes.__                   

![plot of chunk DE_isoform](figure/DE_isoform.png) 

  
### Position of isoform exons on the gene   
  * In general, there are more alternative spliced exons at the __two ends__ of genes.    
  
![plot of chunk exon_pos](figure/exon_pos.png) 

  
### Venn Diagram with average expression level, average No. of exons and average exon length  
  * Common isoforms shared among different comparisons have __much lower ratio of being validated__.        
  * Isoforms have __much lower__ expression level than all expressed genes.          
  * On average, common isoforms between different comparisons have __lower expression level__ than comparison-specific isoforms.                 
  * In general, compared to all isoforms identified, validated isoforms do __not__ have lower expression levels. Our validation approach is not biased towards highly expressed genes.        
  * Average No. of exons are very __similar__ in different sections of the Venn diagram, between all, validated isoforms and all expressed genes.        
  * Average length of isoform exons are __shorter__ than all expressed genes, and validated isoform exons are __sligtly shorter__ than all isoforms in general but the difference is not statistically significant (p = 0.39).         
  
![plot of chunk venn](figure/venn1.png) ![plot of chunk venn](figure/venn2.png) ![plot of chunk venn](figure/venn3.png) ![plot of chunk venn](figure/venn4.png) ![plot of chunk venn](figure/venn5.png) ![plot of chunk venn](figure/venn6.png) ![plot of chunk venn](figure/venn7.png) ![plot of chunk venn](figure/venn8.png) 


### Enrichment of all isoform genes and validated isoform genes for each section on the Venn diagram   
  * Multifunctional correction by ermineJ eliminates general terms and brings more specific terms. Use ermineJ for all GO analysis instead of DAVID.           
  * The enriched terms for validated isoforms are subsets of terms for all isoforms.              

<!--
  * Isoforms identified in lum084 vs myo084 only are enriched in glucuronidation related processes. Isoforms in lum084 vs stem084 only, and myo084 vs stem084 only isoforms do not have significantly enriched terms.         
  * Common isoforms identified in all three comparisons have more enriched terms than comparison-specific isoforms, while isoforms identified in only two comparisons have no significantly enriched terms at all.           
-->

#### lum084 vs myo084

![plot of chunk enrich_lum084_myo084](figure/enrich_lum084_myo0841.png) ![plot of chunk enrich_lum084_myo084](figure/enrich_lum084_myo0842.png) 


#### lum084 vs stem084

![plot of chunk enrich_lum084_stem084](figure/enrich_lum084_stem0841.png) ![plot of chunk enrich_lum084_stem084](figure/enrich_lum084_stem0842.png) 


#### myo084 vs stem084

![plot of chunk enrich_myo084_stem084](figure/enrich_myo084_stem0841.png) ![plot of chunk enrich_myo084_stem084](figure/enrich_myo084_stem0842.png) 


#### Intersecting isoforms in lum084 vs myo084, lum084 vs stem084, and myo084 vs stem084           

<!--
![plot of chunk DAVID_lm_ls_ms](figure/DAVID_lm_ls_ms1.png) ![plot of chunk DAVID_lm_ls_ms](figure/DAVID_lm_ls_ms2.png) 

-->

![plot of chunk erminej_lm_ls_ms](figure/erminej_lm_ls_ms1.png) ![plot of chunk erminej_lm_ls_ms](figure/erminej_lm_ls_ms2.png) 


<!--
#### lum084 vs myo084 only isoforms           
    
![plot of chunk erminej_lm_only](figure/erminej_lm_only1.png) ![plot of chunk erminej_lm_only](figure/erminej_lm_only2.png) 


![plot of chunk DAVID_lm_only](figure/DAVID_lm_only1.png) ![plot of chunk DAVID_lm_only](figure/DAVID_lm_only2.png) 

-->

<!--
#### lum084 vs stem084 only isoforms           
  * No enrichment

```
## Error: Aesthetics must either be length one, or the same length as the
## dataProblems:round(-log10(FDR), 2), Term, -log10(FDR)
```

```
## Error: Aesthetics must either be length one, or the same length as the
## dataProblems:round(-log10(FDR), 2), Term, -log10(FDR)
```

```
## Error: cannot open the connection
```

```
## Error: cannot open the connection
```

```
## Error: object 'valid_ls_only_erminej_BP' not found
```

```
## Error: object 'valid_ls_only_enrich' not found
```

```
## Error: object 'valid_ls_only_enrich' not found
```

```
## Error: object 'valid_ls_only_enrich' not found
```

```
## Error: object 'valid_ls_only_enrich' not found
```

```
## Error: object 'Enrich_valid_ls_only' not found
```


![plot of chunk DAVID_ls_only](figure/DAVID_ls_only1.png) ![plot of chunk DAVID_ls_only](figure/DAVID_ls_only2.png) 


#### myo084 vs stem084 only isoforms           
  * No enrichment

```
## Error: cannot open the connection
```

```
## Error: cannot open the connection
```

```
## Error: object 'all_ms_only_erminej_BP' not found
```

```
## Error: object 'all_ms_only_enrich' not found
```

```
## Error: object 'all_ms_only_enrich' not found
```

```
## Error: object 'all_ms_only_enrich' not found
```

```
## Error: object 'all_ms_only_enrich' not found
```

```
## Error: object 'Enrich_all_ms_only' not found
```

```
## Error: cannot open the connection
```

```
## Error: cannot open the connection
```

```
## Error: object 'valid_ms_only_erminej_BP' not found
```

```
## Error: object 'valid_ms_only_enrich' not found
```

```
## Error: object 'valid_ms_only_enrich' not found
```

```
## Error: object 'valid_ms_only_enrich' not found
```

```
## Error: object 'valid_ms_only_enrich' not found
```

```
## Error: object 'Enrich_valid_ms_only' not found
```


![plot of chunk DAVID_ms_only](figure/DAVID_ms_only1.png) ![plot of chunk DAVID_ms_only](figure/DAVID_ms_only2.png) 

-->

<!--
### Examples for wet lab validation   

#### RM084 lum vs myo    
  * ENSG00000196208: GREB1, growth regulation by estrogen in breast cancer 1          
<img src="./browser/wetlab/ENSG00000196208.png" height="400px" width="400px" />
  * ENSG00000008853: RHOBTB2, Rho-related BTB domain containing 2          
<img src="./browser/wetlab/ENSG00000008853.png" height="400px" width="400px" />
  * ENSG00000108821: COL1A1, collagen, type I, alpha 1            
<img src="./browser/wetlab/ENSG00000108821.png" height="400px" width="400px" />
  * ENSG00000110195: FOLR1, folate receptor 1 (adult)           
<img src="./browser/wetlab/ENSG00000110195.png" height="400px" width="400px" />
  * ENSG00000138795: LEF1, lymphoid enhancer-binding factor 1             
<img src="./browser/wetlab/ENSG00000138795.png" height="400px" width="400px" />
  * ENSG00000170312: CDK1, cyclin-dependent kinase 1               
<img src="./browser/wetlab/ENSG00000170312.png" height="400px" width="400px" />
		
#### RM084 lum vs stem		
  * ENSG00000064787: BCAS1, breast carcinoma amplified sequence 1               
<img src="./browser/wetlab/ENSG00000064787.png" height="400px" width="400px" />
  * ENSG00000127084: FGD3, FYVE RhoGEF and PH domain containing 3             
<img src="./browser/wetlab/ENSG00000127084.png" height="400px" width="400px" />
  * ENSG00000162894: FAIM3, Fas apoptotic inhibitory molecule 3                
<img src="./browser/wetlab/ENSG00000162894.png" height="400px" width="400px" />
  * ENSG00000126217: MCF2L, MCF.2 cell line derived transforming sequence-like             
<img src="./browser/wetlab/ENSG00000126217.png" height="400px" width="400px" />
-->
