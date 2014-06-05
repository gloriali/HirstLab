REMC isoform analysis - junctions
========================================================

Gloria Li         
Thu Jun  5 11:44:23 2014 

<!-- re-knit after modify junction_valid_new.R or junction_other.R script -->






## Validate previously identified isoforms with junction RPKM
  * Previous isoform identification with DE exons
    * DE exons by DEfine FDR = 0.015
    * Exon RPKM $\ge$ 10% gene RPKM in one sample & $\le$ 1% in the other
    * Gene RPKM of both samples > 0.005 
    * Exclude DE genes by DEfine FDR = 0.015

  * Validation: For each isoform exon in the previous pairwise comparison
    * Find junctions associated with this exon with enough coverage, i.e. sum of junction coverage of two samples $\ge$ 1
    * Identify junctions that RPKM change in the same direction as the exon
    * Junction RPKM > 0.1 in one sample and < 0.1 in the other      
  
## Results:  

### Compare strand specific and non-strand specific libraries

  * Non-strand specific libraries (RM080 and RM035) have __much lower__ junction coverage than strand specific libraries (RM084).    
  * Junction RPKM between strand specific and non-strand specific libraries are much more comparable.       

![plot of chunk strand-specifc](figure/strand-specifc1.png) ![plot of chunk strand-specifc](figure/strand-specifc2.png) 


### Junction RPKM clustering

  * Sample clustering on junction RPKM with spearman correlation shows clustering on cell types as well.         
  
![plot of chunk cluster](figure/cluster.png) 


### Isoform junction validation

  * For strand specific libraries, about __40%__ of previously identified exons (__50%__ genes) have junctions with enough coverage for validation. However, __more than 99%__ of them have junction support.          
  * For non-strand specific libraries, only about __20-30%__ isoform genes have enough junction coverage for validation, and among them __75-80%__ genes have support from junction reads.  
  * Venn diagram of isoforms in RM080 is similar to RM084, with the common isoforms shared by different cell types having less ratio of being validated.    
  * Venn diagram of isoforms between same cell types in different individuals show similar patterns. The __majority__ of isoforms are __shared__ among different individual. However, isoforms in strand specific libraries (RM084) have __much higher__ ratio of being validated.       
  
<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Thu Jun  5 11:45:00 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> No.isoform.exons </TH> <TH> No.isoform.genes </TH> <TH> No.exons.with.junction.cov </TH> <TH> No.genes.with.junction.cov </TH> <TH> No.exons.with.junction.support </TH> <TH> No.genes.with.junction.support </TH>  </TR>
  <TR> <TD> lum084_myo084 </TD> <TD align="center"> 8630 </TD> <TD align="center"> 2381 </TD> <TD align="center"> 3618 </TD> <TD align="center"> 1228 </TD> <TD align="center"> 3604 </TD> <TD align="center"> 1217 </TD> </TR>
  <TR> <TD> lum084_stem084 </TD> <TD align="center"> 8948 </TD> <TD align="center"> 2429 </TD> <TD align="center"> 4226 </TD> <TD align="center"> 1411 </TD> <TD align="center"> 4217 </TD> <TD align="center"> 1404 </TD> </TR>
  <TR> <TD> myo084_stem084 </TD> <TD align="center"> 8619 </TD> <TD align="center"> 2427 </TD> <TD align="center"> 3346 </TD> <TD align="center"> 1247 </TD> <TD align="center"> 3333 </TD> <TD align="center"> 1238 </TD> </TR>
  <TR> <TD> lum080_myo080 </TD> <TD align="center"> 12871 </TD> <TD align="center"> 2325 </TD> <TD align="center"> 1687 </TD> <TD align="center"> 509 </TD> <TD align="center"> 1235 </TD> <TD align="center"> 397 </TD> </TR>
  <TR> <TD> lum080_stem080 </TD> <TD align="center"> 8688 </TD> <TD align="center"> 2345 </TD> <TD align="center"> 1559 </TD> <TD align="center"> 686 </TD> <TD align="center"> 1232 </TD> <TD align="center"> 559 </TD> </TR>
  <TR> <TD> myo080_stem080 </TD> <TD align="center"> 8448 </TD> <TD align="center"> 2390 </TD> <TD align="center"> 1105 </TD> <TD align="center"> 548 </TD> <TD align="center"> 791 </TD> <TD align="center"> 427 </TD> </TR>
  <TR> <TD> lum035_myo035 </TD> <TD align="center"> 12911 </TD> <TD align="center"> 2341 </TD> <TD align="center"> 3536 </TD> <TD align="center"> 817 </TD> <TD align="center"> 2324 </TD> <TD align="center"> 619 </TD> </TR>
   </TABLE>
![plot of chunk validate](figure/validate1.png) ![plot of chunk validate](figure/validate2.png) ![plot of chunk validate](figure/validate3.png) ![plot of chunk validate](figure/validate4.png) ![plot of chunk validate](figure/validate5.png) ![plot of chunk validate](figure/validate6.png) ![plot of chunk validate](figure/validate7.png) ![plot of chunk validate](figure/validate8.png) ![plot of chunk validate](figure/validate9.png) ![plot of chunk validate](figure/validate10.png) 

  
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

  * Exon usage along the gene for all multi-transcript genes            
<img src="./figure/exonUsage.png" height="400px" width="400px" />
  
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

<!--
  * The enriched terms for validated isoforms are subsets of terms for all isoforms.              
  * Isoforms identified in lum084 vs myo084 only are enriched in glucuronidation related processes. Isoforms in lum084 vs stem084 only, and myo084 vs stem084 only isoforms do not have significantly enriched terms.         
  * Common isoforms identified in all three comparisons have more enriched terms than comparison-specific isoforms, while isoforms identified in only two comparisons have no significantly enriched terms at all.           
-->

#### lum084 vs myo084

![plot of chunk enrich_lum084_myo084](figure/enrich_lum084_myo0841.png) ![plot of chunk enrich_lum084_myo084](figure/enrich_lum084_myo0842.png) 


#### lum084 vs stem084

![plot of chunk enrich_lum084_stem084](figure/enrich_lum084_stem0841.png) ![plot of chunk enrich_lum084_stem084](figure/enrich_lum084_stem0842.png) 


#### myo084 vs stem084

![plot of chunk enrich_myo084_stem084](figure/enrich_myo084_stem0841.png) ![plot of chunk enrich_myo084_stem084](figure/enrich_myo084_stem0842.png) 


#### Common isoforms shared by lum084 vs myo084, lum084 vs stem084, and myo084 vs stem084           

![plot of chunk erminej_lm_ls_ms](figure/erminej_lm_ls_ms1.png) ![plot of chunk erminej_lm_ls_ms](figure/erminej_lm_ls_ms2.png) 



### Examples for wet lab validation   

#### RM084 lum vs myo    
  * ENSG00000196208: GREB1, growth regulation by estrogen in breast cancer 1          
<img src="./figure/ENSG00000196208.png" height="400px" width="400px" />
  * ENSG00000008853: RHOBTB2, Rho-related BTB domain containing 2          
<img src="./figure/ENSG00000008853.png" height="400px" width="400px" />
  * ENSG00000108821: COL1A1, collagen, type I, alpha 1            
<img src="./figure/ENSG00000108821.png" height="400px" width="400px" />
  * ENSG00000110195: FOLR1, folate receptor 1 (adult)           
<img src="./figure/ENSG00000110195.png" height="400px" width="400px" />
  * ENSG00000138795: LEF1, lymphoid enhancer-binding factor 1             
<img src="./figure/ENSG00000138795.png" height="400px" width="400px" />
  * ENSG00000170312: CDK1, cyclin-dependent kinase 1               
<img src="./figure/ENSG00000170312.png" height="400px" width="400px" />
		
#### RM084 lum vs stem		
  * ENSG00000064787: BCAS1, breast carcinoma amplified sequence 1               
<img src="./figure/ENSG00000064787.png" height="400px" width="400px" />
  * ENSG00000127084: FGD3, FYVE RhoGEF and PH domain containing 3             
<img src="./figure/ENSG00000127084.png" height="400px" width="400px" />
  * ENSG00000162894: FAIM3, Fas apoptotic inhibitory molecule 3                
<img src="./figure/ENSG00000162894.png" height="400px" width="400px" />
  * ENSG00000126217: MCF2L, MCF.2 cell line derived transforming sequence-like             
<img src="./figure/ENSG00000126217.png" height="400px" width="400px" />
