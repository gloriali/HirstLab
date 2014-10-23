Fetal Brain MeDIP Analysis Summary
========================================================

Gloria Li         
Updated: Wed Oct 22 09:49:04 2014 



## DMR analysis from MeDIP fractional calls

  * DM CpG identification: 
    + delta fractional methylation $\ge$ 0.6  
    + fractional methylation of one sample $\ge$ 0.75   
  * Collapse DM CpGs into DMRs:   
    + adjacent CpGs have the same DM status;    
    + distance between adjacent CpGs $\le$ 300bp;   
    + No. of CpGs within each DMR $\ge$ 4.   

## DMRs between neurospheres Cortex and GE derived
### Summary and sanity check  

  * On average, there are __2451__ Cortex UMRs, and __1686__ GE UMRs. There seems to be an asymmetry between Cortex UMRs and GE UMRs.    
  * Median DMR length is __105__, _much smaller than WGBS_. It's similar in all chromosomes.   
  * Median No. of CpGs per DMR is __6__, _larger than WGBS_. Also doesn't fluctuate across chromosomes.  

<!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:05 2014 -->
<TABLE border=1>
<TR> <TH> Sample </TH> <TH> Total.DMR </TH> <TH> Hyper.DMR </TH> <TH> Hypo.DMR </TH>  </TR>
  <TR> <TD align="center"> Cortex-HuFNSC01_GE-HuFNSC01 </TD> <TD align="center"> 3228 </TD> <TD align="center"> 1474 </TD> <TD align="center"> 1754 </TD> </TR>
  <TR> <TD align="center"> Cortex-HuFNSC02_GE-HuFNSC02 </TD> <TD align="center"> 5047 </TD> <TD align="center"> 1899 </TD> <TD align="center"> 3148 </TD> </TR>
   </TABLE>
![](./MeDIP_files/figure-html/MeDIP_neurospheres_sanity-1.png) ![](./MeDIP_files/figure-html/MeDIP_neurospheres_sanity-2.png) 

### Asymmetry between Cortex UMRs and GE UMRs  

  * On average, there are , __1.84__-fold enrichment in total UMR length in Cortex compared to GE, __1.5__ in HuFNSC01, and __2.07__ in HuFNSC02.  
  * The asymmetry appears to be global, in all chromosomes __except for chrX__, and is reproduced in the two individuals.  

![](./MeDIP_files/figure-html/MeDIP_neurospheres_asymmetry-1.png) ![](./MeDIP_files/figure-html/MeDIP_neurospheres_asymmetry-2.png) 

### GREAT analysis on Cortex UMRs and GE UMRs  

  * UMRs in both Cortex and GE in both individuals show enrichment in  __neuron fate commitment__ biological process, __transcriptional regulation__ activities, __Homeobox__ protein domain, and __abnormal brain development__ mouse phenotype.  
  * Cortex UMRs are also enriched in __forebrain regionalization__ and __pattern specification__ processes.   

![](./MeDIP_files/figure-html/MeDIP_neurospheres_GREAT1-1.png) 
![](./MeDIP_files/figure-html/MeDIP_neurospheres_GREAT2-1.png) 
![](./MeDIP_files/figure-html/MeDIP_neurospheres_GREAT3-1.png) 
![](./MeDIP_files/figure-html/MeDIP_neurospheres_GREAT4-1.png) 

### UMR genomic break down  

  + On average, __62.16%__ of CpGs in UMRs overlap with genebody, _similar to WGBS_, and __32.23%__ of CpGs in UMRs overlap with promoters, __2.44-fold__ enriched, _higher than WGBS_. __46.15%__ of CpGs in UMRs overlap with CGIs, __6.23-fold__ than expected by random, _slightly higher than WGBS_.        

<!-- For the entire genome, 3727169 out of 28217448 CpGs overlap with TSS +/- 1500bp promoter regions -->
<!-- For the entire genome, 2089538 out of 28217448 CpGs overlap with CGIs -->

![](./MeDIP_files/figure-html/MeDIP_neurospheres_breakdown-1.png) 

### Proximal UMRs and DE genes

  + On average, there are __490__ UMRs proximally (TSS +/- 1500bp) associated with protein-coding genes, __11.85%__ of all UMRs, _higher than WGBS_.         
  + On average, there are __68__ proximal UMRs associated with DE genes, __13.91%__ of all proximal UMRs, _slightly lower than WGBS_. Among them, there are __57.59%__ unique DE genes change in the same direction as the UMRs, _higher than WGBS_.         
  + The intersect between two individuals are __significant__.  

<!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:11 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> DMRs </TH> <TH> unique.genes </TH> <TH> DE.DMRs </TH> <TH> unique.DE.genes </TH> <TH> same.direction </TH>  </TR>
  <TR> <TD align="center"> GE01.UMRs </TD> <TD align="center"> 419 </TD> <TD align="center"> 466 </TD> <TD align="center">  35 </TD> <TD align="center">  35 </TD> <TD align="center">  20 </TD> </TR>
  <TR> <TD align="center"> Cortex01.UMRs </TD> <TD align="center"> 393 </TD> <TD align="center"> 400 </TD> <TD align="center">  70 </TD> <TD align="center">  67 </TD> <TD align="center">  37 </TD> </TR>
  <TR> <TD align="center"> GE02.UMRs </TD> <TD align="center"> 461 </TD> <TD align="center"> 501 </TD> <TD align="center">  60 </TD> <TD align="center">  58 </TD> <TD align="center">  31 </TD> </TR>
  <TR> <TD align="center"> Cortex02.UMRs </TD> <TD align="center"> 689 </TD> <TD align="center"> 723 </TD> <TD align="center"> 108 </TD> <TD align="center">  97 </TD> <TD align="center">  60 </TD> </TR>
   </TABLE>
![](./MeDIP_files/figure-html/MeDIP_neurospheres_proximal-1.png) ![](./MeDIP_files/figure-html/MeDIP_neurospheres_proximal-2.png) 

#### Cortex UMRs proximal associated DE genes 

<!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:11 2014 -->
<TABLE border=1>
<TR> <TH> name </TH> <TH> description </TH> <TH> DM </TH> <TH> DE </TH>  </TR>
  <TR> <TD align="right"> BMP8B </TD> <TD align="right"> bone_morphogenetic_protein_8b </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> NOS1AP </TD> <TD align="right"> nitric_oxide_synthase_1_(neuronal)_adaptor_protein </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="right"> SFMBT2 </TD> <TD align="right"> Scm-like_with_four_mbt_domains_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> VAX1 </TD> <TD align="right"> ventral_anterior_homeobox_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> RIC3 </TD> <TD align="right"> resistance_to_inhibitors_of_cholinesterase_3_homolog_(C._elegans) </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> TCIRG1 </TD> <TD align="right"> T-cell,_immune_regulator_1,_ATPase,_H+_transporting,_lysosomal_V0_subunit_A3 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="right"> TFCP2 </TD> <TD align="right"> transcription_factor_CP2 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="right"> SIAH3 </TD> <TD align="right"> seven_in_absentia_homolog_3_(Drosophila) </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> NKX2-1 </TD> <TD align="right"> NK2_homeobox_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> ID2 </TD> <TD align="right"> inhibitor_of_DNA_binding_2,_dominant_negative_helix-loop-helix_protein </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="right"> OTX1 </TD> <TD align="right"> orthodenticle_homeobox_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="right"> CYP27A1 </TD> <TD align="right"> cytochrome_P450,_family_27,_subfamily_A,_polypeptide_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="right"> CXCR7 </TD> <TD align="right"> chemokine_(C-X-C_motif)_receptor_7 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> FEZF2 </TD> <TD align="right"> FEZ_family_zinc_finger_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="right"> CYTL1 </TD> <TD align="right"> cytokine-like_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> IRX2 </TD> <TD align="right"> iroquois_homeobox_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> PDGFA </TD> <TD align="right"> platelet-derived_growth_factor_alpha_polypeptide </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> SH3KBP1 </TD> <TD align="right"> SH3-domain_kinase_binding_protein_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
   </TABLE>

#### GE UMRs proximal associated DE genes 

<!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:11 2014 -->
<TABLE border=1>
<TR> <TH> name </TH> <TH> description </TH> <TH> DM </TH> <TH> DE </TH>  </TR>
  <TR> <TD align="right"> C10orf90 </TD> <TD align="right"> chromosome_10_open_reading_frame_90 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> TMEM132B </TD> <TD align="right"> transmembrane_protein_132B </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="right"> PRDM1 </TD> <TD align="right"> PR_domain_containing_1,_with_ZNF_domain </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
   </TABLE>

### Overlap UMRs with TFBSs 

* Overlap UMRs with transcription factor binding sites and count No. of overlapping TFBSs for each TF showed different asymmetry in Cortex and GE.  
* HuFNSC02 show similar trend as WGBS, with TFBSs enriched in Cortex UMRs for most TFs, and top TFs with highest fold cahnge in No. of TFBS between Cortex and GE have overlaps with WGBS, such as GATA2, FOXA1, ZNF217.   
* However, HuFNSC01 showed the opposite trend, with most TFs enriched in GE. _Why?_ In general, the correlation of TFBS Cortex UMR vs GE UMR fold change between the two individual is quite low, 0.21, _similar to WGBS_. Maybe TFBSs trend is individual-specific.   
* With 1.5-fold change cutoff, there are 6 TFs enriched in Cortex UMRs in both individuals, and 4 in GE UMRs.   

![](./MeDIP_files/figure-html/MeDIP_neurospheres_TF-1.png) <!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:11 2014 -->
<TABLE border=1>
<TR> <TH> TF </TH> <TH> Cortex01UMR </TH> <TH> GE01UMR </TH> <TH> Ratio01 </TH> <TH> Cortex02UMR </TH> <TH> GE02UMR </TH> <TH> Ratio02 </TH>  </TR>
  <TR> <TD align="center"> SUZ12 </TD> <TD align="center"> 173 </TD> <TD align="center">  98 </TD> <TD align="center"> 1.77 </TD> <TD align="center"> 256 </TD> <TD align="center"> 160 </TD> <TD align="center"> 1.60 </TD> </TR>
  <TR> <TD align="center"> EZH2 </TD> <TD align="center"> 546 </TD> <TD align="center"> 309 </TD> <TD align="center"> 1.77 </TD> <TD align="center"> 764 </TD> <TD align="center"> 424 </TD> <TD align="center"> 1.80 </TD> </TR>
  <TR> <TD align="center"> ZNF217 </TD> <TD align="center">  13 </TD> <TD align="center">   7 </TD> <TD align="center"> 1.86 </TD> <TD align="center">  22 </TD> <TD align="center">   8 </TD> <TD align="center"> 2.75 </TD> </TR>
  <TR> <TD align="center"> ESR1 </TD> <TD align="center">  47 </TD> <TD align="center">  24 </TD> <TD align="center"> 1.96 </TD> <TD align="center">  91 </TD> <TD align="center">  40 </TD> <TD align="center"> 2.27 </TD> </TR>
  <TR> <TD align="center"> NANOG </TD> <TD align="center">  16 </TD> <TD align="center">   7 </TD> <TD align="center"> 2.29 </TD> <TD align="center">  23 </TD> <TD align="center">  12 </TD> <TD align="center"> 1.92 </TD> </TR>
  <TR> <TD align="center"> ESRRA </TD> <TD align="center">   5 </TD> <TD align="center">   2 </TD> <TD align="center"> 2.50 </TD> <TD align="center">   7 </TD> <TD align="center">   4 </TD> <TD align="center"> 1.75 </TD> </TR>
   </TABLE>
<!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:11 2014 -->
<TABLE border=1>
<TR> <TH> TF </TH> <TH> Cortex01UMR </TH> <TH> GE01UMR </TH> <TH> Ratio01 </TH> <TH> Cortex02UMR </TH> <TH> GE02UMR </TH> <TH> Ratio02 </TH>  </TR>
  <TR> <TD align="center"> SIRT6 </TD> <TD align="center">   2 </TD> <TD align="center">   9 </TD> <TD align="center"> 0.22 </TD> <TD align="center">   3 </TD> <TD align="center">  12 </TD> <TD align="center"> 0.25 </TD> </TR>
  <TR> <TD align="center"> KDM5A </TD> <TD align="center">   4 </TD> <TD align="center">  15 </TD> <TD align="center"> 0.27 </TD> <TD align="center">   6 </TD> <TD align="center">  13 </TD> <TD align="center"> 0.46 </TD> </TR>
  <TR> <TD align="center"> ATF1 </TD> <TD align="center">  16 </TD> <TD align="center">  33 </TD> <TD align="center"> 0.48 </TD> <TD align="center">  22 </TD> <TD align="center">  34 </TD> <TD align="center"> 0.65 </TD> </TR>
  <TR> <TD align="center"> SP2 </TD> <TD align="center">  11 </TD> <TD align="center">  22 </TD> <TD align="center"> 0.50 </TD> <TD align="center">  18 </TD> <TD align="center">  28 </TD> <TD align="center"> 0.64 </TD> </TR>
   </TABLE>


## DMRs between monozygptic twins
### Summary and sanity check

  * On average, there are __3449__ DMR regions identified across three cell types, with __1881__ hypermethylated, and __1568__ hypomethylated.      
  * Median length of all DMRs is __92bp__, _smaller_ than WGBS.    
  * Median No. of CpGs per DMR is __6__, _larger_ than WGBS.        

<!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:11 2014 -->
<TABLE border=1>
<TR> <TH> Sample </TH> <TH> Total.DMR </TH> <TH> Hyper.DMR </TH> <TH> Hypo.DMR </TH>  </TR>
  <TR> <TD> Brain-HuFNSC01_Brain-HuFNSC02 </TD> <TD align="center"> 4472 </TD> <TD align="center"> 2750 </TD> <TD align="center"> 1722 </TD> </TR>
  <TR> <TD> Cortex-HuFNSC01_Cortex-HuFNSC02 </TD> <TD align="center"> 3161 </TD> <TD align="center"> 1758 </TD> <TD align="center"> 1403 </TD> </TR>
  <TR> <TD> GE-HuFNSC01_GE-HuFNSC02 </TD> <TD align="center"> 2716 </TD> <TD align="center"> 1136 </TD> <TD align="center"> 1580 </TD> </TR>
   </TABLE>
![](./MeDIP_files/figure-html/MeDIP_MZ_sanity-1.png) ![](./MeDIP_files/figure-html/MeDIP_MZ_sanity-2.png) 

### Asymmetry between MZ twins 

  + There is an __asymmetry__ between UMRs in HuFNSC01 and HuFNSC02 in the Brain and Cortex neurosphere, but __not__ in GE neurosphere. Fold change in total UMR length HuFNSC02/HuFNSC01 in Brain is __2.33__, in Cortex is __1.76__, and in GE is __0.89__.    
  
![](./MeDIP_files/figure-html/MeDIP_MZ_asymmetry-1.png) ![](./MeDIP_files/figure-html/MeDIP_MZ_asymmetry-2.png) 

### GREAT analysis on MZ UMRs

* __Homeobox__ protein domain is enriched in all UMR lists.   
* __Brain development related__ biological processes are enriched in most lists.   

![](./MeDIP_files/figure-html/MeDIP_MZ_GREAT1-1.png) ![](./MeDIP_files/figure-html/MeDIP_MZ_GREAT1-2.png) ![](./MeDIP_files/figure-html/MeDIP_MZ_GREAT1-3.png) 
![](./MeDIP_files/figure-html/MeDIP_MZ_GREAT2-1.png) ![](./MeDIP_files/figure-html/MeDIP_MZ_GREAT2-2.png) ![](./MeDIP_files/figure-html/MeDIP_MZ_GREAT2-3.png) 

### MZ UMR breakdown

  + Majority of UMR CpGs overlap with genebody, on average __62.59%__. And UMR CpGs are __3-fold__ enriched in promoter regions, with on average __34.29%__ UMR CpGs are in promoter regions, __39.54%__ of CpGs in UMRs overlap with CGIs, __5.34-fold__ than expected by random. 
  + Brain seems to have less CGI/promoter UMRs.   
  
<!-- For the entire genome, 3727169 out of 28217448 CpGs overlap with TSS +/- 1500bp promoter regions -->  
<!-- For the entire genome, 2089538 out of 28217448 CpGs overlap with CGIs -->

![](./MeDIP_files/figure-html/MeDIP_MZ_breakdown-1.png) 

### Proximal UMRs and DE genes  

  + On average, there are __405__ UMRs proximally (TSS +/- 1500bp) associated with protein-coding genes, __11.75%__ of all UMRs, _similar to UMRs between Cortex and GE_.         
  + On average, there are __17__ proximal UMRs associated with DE genes, __4.24%__ of all proximal DMRs, _much less than UMRs between Cortex and GE (less functional?)_. Among them, there are __52.53%__ unique DE genes change in the same direction as the UMRs, _similar to UMRs between Cortex and GE_.         
  + Proximal UMR genes are mostly __cell type specific__, but the overlap is still statistically significant.  
  + No overlap between Brain and Cortex for proximal HuFNSC01 UMR DE genes, only one for HuFNSC02, __neuropeptide_Y__.   

<!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:19 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> DMRs </TH> <TH> unique.genes </TH> <TH> DE.DMRs </TH> <TH> unique.DE.genes </TH> <TH> same.direction </TH>  </TR>
  <TR> <TD> Brain01_Brain02_hyper </TD> <TD align="center"> 432 </TD> <TD align="center"> 473 </TD> <TD align="center">  26 </TD> <TD align="center">  26 </TD> <TD align="center">  10 </TD> </TR>
  <TR> <TD> Brain01_Brain02_hypo </TD> <TD align="center"> 367 </TD> <TD align="center"> 409 </TD> <TD align="center">  20 </TD> <TD align="center">  20 </TD> <TD align="center">  13 </TD> </TR>
  <TR> <TD> Cortex01_Cortex02_hyper </TD> <TD align="center"> 587 </TD> <TD align="center"> 659 </TD> <TD align="center">  30 </TD> <TD align="center">  28 </TD> <TD align="center">  17 </TD> </TR>
  <TR> <TD> Cortex01_Cortex02_hypo </TD> <TD align="center"> 316 </TD> <TD align="center"> 342 </TD> <TD align="center">  24 </TD> <TD align="center">  22 </TD> <TD align="center">  11 </TD> </TR>
  <TR> <TD> GE01_GE02_hyper </TD> <TD align="center"> 369 </TD> <TD align="center"> 402 </TD> <TD align="center">   1 </TD> <TD align="center">   1 </TD> <TD align="center">   1 </TD> </TR>
  <TR> <TD> GE01_GE02_hypo </TD> <TD align="center"> 361 </TD> <TD align="center"> 398 </TD> <TD align="center">   2 </TD> <TD align="center">   2 </TD> <TD align="center">   0 </TD> </TR>
   </TABLE>
![](./MeDIP_files/figure-html/MeDIP_MZ_proximal-1.png) ![](./MeDIP_files/figure-html/MeDIP_MZ_proximal-2.png) 

#### HuFNSC01 UMRs proximal associated DE genes

<!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:19 2014 -->
<TABLE border=1>
<TR> <TH> CellType </TH> <TH> name </TH> <TH> description </TH> <TH> DM </TH> <TH> DE </TH>  </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> SAMD11 </TD> <TD align="right"> sterile_alpha_motif_domain_containing_11 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> ITGA8 </TD> <TD align="right"> integrin,_alpha_8 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> EMX2 </TD> <TD align="right"> empty_spiracles_homeobox_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> ARNTL2 </TD> <TD align="right"> aryl_hydrocarbon_receptor_nuclear_translocator-like_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> CPNE8 </TD> <TD align="right"> copine_VIII </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> PCDH17 </TD> <TD align="right"> protocadherin_17 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right">  </TD> <TD align="right"> Putative_3-phosphoinositide-dependent_protein_kinase_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> CPNE7 </TD> <TD align="right"> copine_VII </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> PKDREJ </TD> <TD align="right"> polycystic_kidney_disease_(polycystin)_and_REJ_homolog_(sperm_receptor_for_egg_jelly_homolog,_sea_urchin) </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> THRB </TD> <TD align="right"> thyroid_hormone_receptor,_beta </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> SFRP2 </TD> <TD align="right"> secreted_frizzled-related_protein_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> IER3 </TD> <TD align="right"> immediate_early_response_3 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> ZBTB12 </TD> <TD align="right"> zinc_finger_and_BTB_domain_containing_12 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> MDGA1 </TD> <TD align="right"> MAM_domain_containing_glycosylphosphatidylinositol_anchor_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> EPHA1 </TD> <TD align="right"> EPH_receptor_A1 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right">  </TD> <TD align="right"> Protein_kinase-like_protein_SgK196 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> TOX </TD> <TD align="right"> thymocyte_selection-associated_high_mobility_group_box </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> NCOA2 </TD> <TD align="right"> nuclear_receptor_coactivator_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> GPR64 </TD> <TD align="right"> G_protein-coupled_receptor_64 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> TRPC5 </TD> <TD align="right"> transient_receptor_potential_cation_channel,_subfamily_C,_member_5 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> ELTD1 </TD> <TD align="right"> EGF,_latrophilin_and_seven_transmembrane_domain_containing_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> CRABP2 </TD> <TD align="right"> cellular_retinoic_acid_binding_protein_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> SPOCK2 </TD> <TD align="right"> sparc/osteonectin,_cwcv_and_kazal-like_domains_proteoglycan_(testican)_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> VAX1 </TD> <TD align="right"> ventral_anterior_homeobox_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> LRRC10B </TD> <TD align="right"> leucine_rich_repeat_containing_10B </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> NDUFA4L2 </TD> <TD align="right"> NADH_dehydrogenase_(ubiquinone)_1_alpha_subcomplex,_4-like_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> FOXO1 </TD> <TD align="right"> forkhead_box_O1 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> C16orf74 </TD> <TD align="right"> chromosome_16_open_reading_frame_74 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> WNT3 </TD> <TD align="right"> wingless-type_MMTV_integration_site_family,_member_3 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> CHST9 </TD> <TD align="right"> carbohydrate_(N-acetylgalactosamine_4-0)_sulfotransferase_9 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> CD97 </TD> <TD align="right"> CD97_molecule </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> RNF144A </TD> <TD align="right"> ring_finger_protein_144A </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> DLX1 </TD> <TD align="right"> distal-less_homeobox_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> OLIG1 </TD> <TD align="right"> oligodendrocyte_transcription_factor_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> SCUBE1 </TD> <TD align="right"> signal_peptide,_CUB_domain,_EGF-like_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> GRIA1 </TD> <TD align="right"> glutamate_receptor,_ionotropic,_AMPA_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> MMD2 </TD> <TD align="right"> monocyte_to_macrophage_differentiation-associated_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> SEMA3E </TD> <TD align="right"> sema_domain,_immunoglobulin_domain_(Ig),_short_basic_domain,_secreted,_(semaphorin)_3E </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> VIPR2 </TD> <TD align="right"> vasoactive_intestinal_peptide_receptor_2 </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> AP1S2 </TD> <TD align="right"> adaptor-related_protein_complex_1,_sigma_2_subunit </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> SH3KBP1 </TD> <TD align="right"> SH3-domain_kinase_binding_protein_1 </TD> <TD align="center"> hypo </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> GE </TD> <TD align="right"> FAM5B </TD> <TD align="right"> family_with_sequence_similarity_5,_member_B </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> GE </TD> <TD align="right"> DLL1 </TD> <TD align="right"> delta-like_1_(Drosophila) </TD> <TD align="center"> hypo </TD> <TD align="center"> DN </TD> </TR>
   </TABLE>

#### HuFNSC02 UMRs proximal associated DE genes

<!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:19 2014 -->
<TABLE border=1>
<TR> <TH> CellType </TH> <TH> name </TH> <TH> description </TH> <TH> DM </TH> <TH> DE </TH>  </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> EPHA8 </TD> <TD align="right"> EPH_receptor_A8 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> IFI6 </TD> <TD align="right"> interferon,_alpha-inducible_protein_6 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> BMP8A </TD> <TD align="right"> bone_morphogenetic_protein_8a </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> ALX3 </TD> <TD align="right"> ALX_homeobox_3 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> NRGN </TD> <TD align="right"> neurogranin_(protein_kinase_C_substrate,_RC3) </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> C1QL4 </TD> <TD align="right"> complement_component_1,_q_subcomponent-like_4 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> CSRP2 </TD> <TD align="right"> cysteine_and_glycine-rich_protein_2 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> TTYH2 </TD> <TD align="right"> tweety_homolog_2_(Drosophila) </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> DTNB </TD> <TD align="right"> dystrobrevin,_beta </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> CXCR7 </TD> <TD align="right"> chemokine_(C-X-C_motif)_receptor_7 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> VHL </TD> <TD align="right"> von_Hippel-Lindau_tumor_suppressor </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> EPHA6 </TD> <TD align="right"> EPH_receptor_A6 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> CAMK2N2 </TD> <TD align="right"> calcium/calmodulin-dependent_protein_kinase_II_inhibitor_2 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> TBC1D1 </TD> <TD align="right"> TBC1_(tre-2/USP6,_BUB2,_cdc16)_domain_family,_member_1 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> NMU </TD> <TD align="right"> neuromedin_U </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> FSTL5 </TD> <TD align="right"> follistatin-like_5 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> CARTPT </TD> <TD align="right"> CART_prepropeptide </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> EDIL3 </TD> <TD align="right"> EGF-like_repeats_and_discoidin_I-like_domains_3 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> CCDC90A </TD> <TD align="right"> coiled-coil_domain_containing_90A </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right">  </TD> <TD align="right"> LOC401296_proteinUncharacterized_protein </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> NPY </TD> <TD align="right"> neuropeptide_Y </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> SFRP1 </TD> <TD align="right"> secreted_frizzled-related_protein_1 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> TOX </TD> <TD align="right"> thymocyte_selection-associated_high_mobility_group_box </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> ENTPD2 </TD> <TD align="right"> ectonucleoside_triphosphate_diphosphohydrolase_2 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> KLHL4 </TD> <TD align="right"> kelch-like_4_(Drosophila) </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Brain </TD> <TD align="right"> SOX3 </TD> <TD align="right"> SRY_(sex_determining_region_Y)-box_3 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> DMRTA2 </TD> <TD align="right"> DMRT-like_family_A2 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> NPR1 </TD> <TD align="right"> natriuretic_peptide_receptor_A/guanylate_cyclase_A_(atrionatriuretic_peptide_receptor_A) </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> PTPRE </TD> <TD align="right"> protein_tyrosine_phosphatase,_receptor_type,_E </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> PDE3B </TD> <TD align="right"> phosphodiesterase_3B,_cGMP-inhibited </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> LDHA </TD> <TD align="right"> lactate_dehydrogenase_A </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> GNG3 </TD> <TD align="right"> guanine_nucleotide_binding_protein_(G_protein),_gamma_3 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> METTL7B </TD> <TD align="right"> methyltransferase_like_7B </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> NME3 </TD> <TD align="right"> non-metastatic_cells_3,_protein_expressed_in </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> GRIN2A </TD> <TD align="right"> glutamate_receptor,_ionotropic,_N-methyl_D-aspartate_2A </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> NOL3 </TD> <TD align="right"> nucleolar_protein_3_(apoptosis_repressor_with_CARD_domain) </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> NR1D1 </TD> <TD align="right"> nuclear_receptor_subfamily_1,_group_D,_member_1 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> IGFBP4 </TD> <TD align="right"> insulin-like_growth_factor_binding_protein_4 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> TBX2 </TD> <TD align="right"> T-box_2 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> LBH </TD> <TD align="right"> limb_bud_and_heart_development_homolog_(mouse) </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> INSIG2 </TD> <TD align="right"> insulin_induced_gene_2 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> WNT7A </TD> <TD align="right"> wingless-type_MMTV_integration_site_family,_member_7A </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> SNCA </TD> <TD align="right"> synuclein,_alpha_(non_A4_component_of_amyloid_precursor) </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> NPNT </TD> <TD align="right"> nephronectin </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> SFRP2 </TD> <TD align="right"> secreted_frizzled-related_protein_2 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> PRR16 </TD> <TD align="right"> proline_rich_16 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> STC2 </TD> <TD align="right"> stanniocalcin_2 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> SNCB </TD> <TD align="right"> synuclein,_beta </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> ITPR3 </TD> <TD align="right"> inositol_1,4,5-trisphosphate_receptor,_type_3 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> NPY </TD> <TD align="right"> neuropeptide_Y </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> SEMA3D </TD> <TD align="right"> sema_domain,_immunoglobulin_domain_(Ig),_short_basic_domain,_secreted,_(semaphorin)_3D </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> VGF </TD> <TD align="right"> VGF_nerve_growth_factor_inducible </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> CAV1 </TD> <TD align="right"> caveolin_1,_caveolae_protein,_22kDa </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> Cortex </TD> <TD align="right"> SHROOM2 </TD> <TD align="right"> shroom_family_member_2 </TD> <TD align="center"> hyper </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="center"> GE </TD> <TD align="right"> DBC1 </TD> <TD align="right"> deleted_in_bladder_cancer_1 </TD> <TD align="center"> hyper </TD> <TD align="center"> DN </TD> </TR>
   </TABLE>

### Overlap UMRs with TFBSs 

* Overlap UMRs with transcription factor binding sites and count No. of overlapping TFBSs for each TF showed similar asymmetry in Brain and Cortex, but is symmetric in GE. The correlation between Brain and Cortex is also very low, __0.15__.    
* With 2-fold change cutoff, there are 18 TFs enriched in HuFNSC02 in both Brain and Cortex.   

![](./MeDIP_files/figure-html/MeDIP_MZ_TF-1.png) <!-- html table generated in R 3.1.1 by xtable 1.7-3 package -->
<!-- Wed Oct 22 09:49:19 2014 -->
<TABLE border=1>
<TR> <TH> TF </TH> <TH> Brain.hypo </TH> <TH> Brain.hyper </TH> <TH> Ratio.Brain </TH> <TH> Cortex.hypo </TH> <TH> Cortex.hyper </TH> <TH> Ratio.Cortex </TH>  </TR>
  <TR> <TD align="center"> ESRRA </TD> <TD align="center">   1 </TD> <TD align="center">   7 </TD> <TD align="center"> 0.14 </TD> <TD align="center">   1 </TD> <TD align="center">   5 </TD> <TD align="center"> 0.20 </TD> </TR>
  <TR> <TD align="center"> TAL1 </TD> <TD align="center">  10 </TD> <TD align="center">  42 </TD> <TD align="center"> 0.24 </TD> <TD align="center">   9 </TD> <TD align="center">  34 </TD> <TD align="center"> 0.26 </TD> </TR>
  <TR> <TD align="center"> HNF4A </TD> <TD align="center">  15 </TD> <TD align="center">  49 </TD> <TD align="center"> 0.31 </TD> <TD align="center">  11 </TD> <TD align="center">  27 </TD> <TD align="center"> 0.41 </TD> </TR>
  <TR> <TD align="center"> HNF4G </TD> <TD align="center">  13 </TD> <TD align="center">  41 </TD> <TD align="center"> 0.32 </TD> <TD align="center">  10 </TD> <TD align="center">  29 </TD> <TD align="center"> 0.34 </TD> </TR>
  <TR> <TD align="center"> ZNF217 </TD> <TD align="center">   6 </TD> <TD align="center">  17 </TD> <TD align="center"> 0.35 </TD> <TD align="center">   6 </TD> <TD align="center">  14 </TD> <TD align="center"> 0.43 </TD> </TR>
  <TR> <TD align="center"> EBF1 </TD> <TD align="center">  49 </TD> <TD align="center"> 130 </TD> <TD align="center"> 0.38 </TD> <TD align="center">  48 </TD> <TD align="center">  98 </TD> <TD align="center"> 0.49 </TD> </TR>
  <TR> <TD align="center"> FOS </TD> <TD align="center">  47 </TD> <TD align="center"> 123 </TD> <TD align="center"> 0.38 </TD> <TD align="center">  43 </TD> <TD align="center"> 118 </TD> <TD align="center"> 0.36 </TD> </TR>
  <TR> <TD align="center"> GATA1 </TD> <TD align="center">  31 </TD> <TD align="center">  73 </TD> <TD align="center"> 0.42 </TD> <TD align="center">  32 </TD> <TD align="center">  80 </TD> <TD align="center"> 0.40 </TD> </TR>
  <TR> <TD align="center"> HSF1 </TD> <TD align="center">   3 </TD> <TD align="center">   7 </TD> <TD align="center"> 0.43 </TD> <TD align="center">   2 </TD> <TD align="center">   6 </TD> <TD align="center"> 0.33 </TD> </TR>
  <TR> <TD align="center"> FOXA2 </TD> <TD align="center">  26 </TD> <TD align="center">  60 </TD> <TD align="center"> 0.43 </TD> <TD align="center">  14 </TD> <TD align="center">  36 </TD> <TD align="center"> 0.39 </TD> </TR>
  <TR> <TD align="center"> GATA3 </TD> <TD align="center">  35 </TD> <TD align="center">  80 </TD> <TD align="center"> 0.44 </TD> <TD align="center">  22 </TD> <TD align="center">  52 </TD> <TD align="center"> 0.42 </TD> </TR>
  <TR> <TD align="center"> GATA2 </TD> <TD align="center">  58 </TD> <TD align="center"> 130 </TD> <TD align="center"> 0.45 </TD> <TD align="center">  39 </TD> <TD align="center">  96 </TD> <TD align="center"> 0.41 </TD> </TR>
  <TR> <TD align="center"> STAT3 </TD> <TD align="center">  41 </TD> <TD align="center">  91 </TD> <TD align="center"> 0.45 </TD> <TD align="center">  26 </TD> <TD align="center">  87 </TD> <TD align="center"> 0.30 </TD> </TR>
  <TR> <TD align="center"> CEBPD </TD> <TD align="center">  22 </TD> <TD align="center">  48 </TD> <TD align="center"> 0.46 </TD> <TD align="center">  23 </TD> <TD align="center">  64 </TD> <TD align="center"> 0.36 </TD> </TR>
  <TR> <TD align="center"> SMC3 </TD> <TD align="center">  71 </TD> <TD align="center"> 150 </TD> <TD align="center"> 0.47 </TD> <TD align="center">  62 </TD> <TD align="center"> 144 </TD> <TD align="center"> 0.43 </TD> </TR>
  <TR> <TD align="center"> NR3C1 </TD> <TD align="center">  46 </TD> <TD align="center">  95 </TD> <TD align="center"> 0.48 </TD> <TD align="center">  39 </TD> <TD align="center">  90 </TD> <TD align="center"> 0.43 </TD> </TR>
  <TR> <TD align="center"> RXRA </TD> <TD align="center">  23 </TD> <TD align="center">  47 </TD> <TD align="center"> 0.49 </TD> <TD align="center">  16 </TD> <TD align="center">  35 </TD> <TD align="center"> 0.46 </TD> </TR>
  <TR> <TD align="center"> HDAC6 </TD> <TD align="center">   1 </TD> <TD align="center">   2 </TD> <TD align="center"> 0.50 </TD> <TD align="center">   2 </TD> <TD align="center">  11 </TD> <TD align="center"> 0.18 </TD> </TR>
   </TABLE>



