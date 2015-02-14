Fetal Brain WGBS Analysis Summary - DMRs between Cortex and GE
========================================================

Gloria Li         
Updated: Mon Dec  1 11:35:29 2014 



## DMR identification with methyl_diff

  * Identify DM CpGs     
    + methyl_diff one-sided p-value $\le$ 0.005  
    + delta fractional methylation $\ge$ 0.5  
    + fractional methylation of one sample $\ge$ 0.75   
  * Collapse DM CpGs into DMRs     
    + adjacent DM CpGs have the same DM status;    
    + distance between adjacent CpGs (size) $\le$ 300bp;   
    + No. of CpGs within each DMR $\ge$ 3.   

## Summary and sanity check  

  * On average, there are __1156__ Cortex UMRs, 179 intersect between two individuals, and __255__ GE UMRs, 10 intersect. The intersect is significant. And there seems to be an asymmetry between Cortex UMRs and GE UMRs.    
  * Median DMR length is __267__, _comparable to breast_. It's similar in all chromosomes in Cortex UMRs, but fluctuate more in GE UMRs, probably due to  small No. of UMRs identified.   
  * Median No. of CpGs per DMR is __5__, _similar to breast_. chr11 and chr13 in GE UMRs have higher No. of CpGs per DMR.  

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec  1 11:35:42 2014 -->
<table border=1>
<tr> <th> Sample </th> <th> Total.DMR </th> <th> Hyper.DMR </th> <th> Hypo.DMR </th>  </tr>
  <tr> <td align="center"> Cortex-HuFNSC02_GE-HuFNSC02 </td> <td align="center"> 2178 </td> <td align="center"> 420 </td> <td align="center"> 1758 </td> </tr>
  <tr> <td align="center"> Cortex-HuFNSC04_GE-HuFNSC04 </td> <td align="center"> 646 </td> <td align="center">  91 </td> <td align="center"> 555 </td> </tr>
   </table>
![](./WGBS_files/figure-html/WGBS_sanity-1.png) ![](./WGBS_files/figure-html/WGBS_sanity-2.png) ![](./WGBS_files/figure-html/WGBS_sanity-3.png) ![](./WGBS_files/figure-html/WGBS_sanity-4.png) ![](./WGBS_files/figure-html/WGBS_sanity-5.png) 

## Asymmetry between Cortex UMRs and GE UMRs  

  * On average, there are , __NaN__-fold enrichment in total UMR length in Cortex compared to GE, __NaN__ in HuFNSC02, and __NaN__ in HuFNSC04.    
  * The asymmetry appears to be global, in all chromosomes. It is reproduced within the two individuals.  
  * __Single CpG level__ differential methylation is __symmetric__, but the asymmetry on UMR level can be reproduced with __different cutoffs__. However, there are __no apparent differences in UMR length__ between Cortex and GE, suggesting that there are more __orphan GE UM CpGs__ that was not able to form UMRs than in Cortex. 

![](./WGBS_files/figure-html/WGBS_asymmetry-1.png) ![](./WGBS_files/figure-html/WGBS_asymmetry-2.png) ![](./WGBS_files/figure-html/WGBS_asymmetry-3.png) ![](./WGBS_files/figure-html/WGBS_asymmetry-4.png) 

## GREAT analysis on Cortex UMRs and GE UMRs  

  * UMRs in both Cortex and GE in both individuals show enrichment in  __transcriptional regulation__ activities.  
  * In HuFNSC02, both Cortex UMRs are enriched in __brain regions development__, and GE UMRs are enriched in __neuron development__.   
  * In HuFNSC04, Cortex UMRs show __abnormal brain development__ in Mouse Phenotype, but are also enriched in __kidney-related processes__.   

![](./WGBS_files/figure-html/WGBS_GREAT1-1.png) 
![](./WGBS_files/figure-html/WGBS_GREAT2-1.png) 
![](./WGBS_files/figure-html/WGBS_GREAT3-1.png) 
![](./WGBS_files/figure-html/WGBS_GREAT4-1.png) 

## UMR genomic break down  

  + On average, __65.94%__ of CpGs in UMRs overlap with genebody, and __15.83%__ of CpGs in UMRs overlap with promoters, not a significant enrichment __(1.2-fold)__. __40.2%__ of CpGs in UMRs overlap with CGIs, __5.43-fold__ than expected by random.        

<!-- For the entire genome, 3727169 out of 28217448 CpGs overlap with TSS +/- 1500bp promoter regions -->
<!-- For the entire genome, 2089538 out of 28217448 CpGs overlap with CGIs -->

![](./WGBS_files/figure-html/WGBS_breakdown-1.png) 

## UMRs intersecting with protein-coding genes and DE genes

  + On average, there are __40__ UMRs proximally associated with protein-coding genes, __2.85%__ of all UMRs.         
  + There are average __8__ proximal UMRs associated with DE genes, __19.88%__ of all proximal UMRs, __much lower than observed in breast__, _maybe more UMRs in enhancers? (ChromHMM)_. Among them, there are __46.67%__ unique DE genes change in the same direction as the UMRs, __also lower than observed in breast__.         
  + The intersect between two individuals are __significant__ in Cortex UMRs. There are 3 DE genes with proximal Cortex UMR that are shared by both individual: __GFAP, NFIX, FEZF1__.     
  + There are no intersect in pc genes with proximal GE UMRs between the two individual.  
  + There are __no significant__ DAVID enrichment terms due to the small number of genes.  

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec  1 11:35:54 2014 -->
<table border=1>
<tr> <th>  </th> <th> pc.Genes </th> <th> unique.Genes </th> <th> pc.Promoters </th> <th> unique.Promoters </th> <th> proximal.DE.Genes </th> <th> same.direction </th> <th> unique.DE.Genes </th>  </tr>
  <tr> <td align="center"> GE_UMRs-HuFNSC02 </td> <td align="center"> 222 </td> <td align="center"> 210 </td> <td align="center">  17 </td> <td align="center">  17 </td> <td align="center">   4 </td> <td align="center">   1 </td> <td align="center">   4 </td> </tr>
  <tr> <td align="center"> Cortex_UMRs-HuFNSC02 </td> <td align="center"> 903 </td> <td align="center"> 712 </td> <td align="center">  53 </td> <td align="center">  52 </td> <td align="center">  15 </td> <td align="center">   8 </td> <td align="center">  14 </td> </tr>
  <tr> <td align="center"> GE_UMRs-HuFNSC04 </td> <td align="center">  47 </td> <td align="center">  46 </td> <td align="center">   9 </td> <td align="center">  11 </td> <td align="center">   3 </td> <td align="center">   2 </td> <td align="center">   3 </td> </tr>
  <tr> <td align="center"> Cortex_UMRs-HuFNSC04 </td> <td align="center"> 286 </td> <td align="center"> 287 </td> <td align="center">  82 </td> <td align="center">  84 </td> <td align="center">  10 </td> <td align="center">   3 </td> <td align="center">   9 </td> </tr>
   </table>
![](./WGBS_files/figure-html/WGBS_proximal-1.png) ![](./WGBS_files/figure-html/WGBS_proximal-2.png) ![](./WGBS_files/figure-html/WGBS_proximal-3.png) ![](./WGBS_files/figure-html/WGBS_proximal-4.png) 

### DE genes with promoter Cortex UMRs  
#### HuFNSC02  
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec  1 11:35:55 2014 -->
<table border=1>
<tr> <th> name </th> <th> description </th> <th> DE </th>  </tr>
  <tr> <td align="right"> CD58 </td> <td align="right"> CD58_molecule </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> RGS10 </td> <td align="right"> regulator_of_G-protein_signaling_10 </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> ADM </td> <td align="right"> adrenomedullin </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> NKX2-1 </td> <td align="right"> NK2_homeobox_1 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> OTX2 </td> <td align="right"> orthodenticle_homeobox_2 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> USP43 </td> <td align="right"> ubiquitin_specific_peptidase_43 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> GFAP </td> <td align="right"> glial_fibrillary_acidic_protein </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> NFIX </td> <td align="right"> nuclear_factor_I/X_(CCAAT-binding_transcription_factor) </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> GAD1 </td> <td align="right"> glutamate_decarboxylase_1_(brain,_67kDa) </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> FZD7 </td> <td align="right"> frizzled_family_receptor_7 </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> FZD5 </td> <td align="right"> frizzled_family_receptor_5 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> CXCR7 </td> <td align="right"> chemokine_(C-X-C_motif)_receptor_7 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> ZAR1 </td> <td align="right"> zygote_arrest_1 </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> FEZF1 </td> <td align="right"> FEZ_family_zinc_finger_1 </td> <td align="center"> DN </td> </tr>
   </table>

#### HuFNSC04  
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec  1 11:35:55 2014 -->
<table border=1>
<tr> <th> name </th> <th> description </th> <th> DE </th>  </tr>
  <tr> <td align="right"> FAM5C </td> <td align="right"> family_with_sequence_similarity_5,_member_C </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> STXBP6 </td> <td align="right"> syntaxin_binding_protein_6_(amisyn) </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> GFAP </td> <td align="right"> glial_fibrillary_acidic_protein </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> NFIX </td> <td align="right"> nuclear_factor_I/X_(CCAAT-binding_transcription_factor) </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> MEIS1 </td> <td align="right"> Meis_homeobox_1 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> FSIP2 </td> <td align="right"> fibrous_sheath_interacting_protein_2 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> INSM1 </td> <td align="right"> insulinoma-associated_1 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> FEZF1 </td> <td align="right"> FEZ_family_zinc_finger_1 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> C9orf172 </td> <td align="right"> chromosome_9_open_reading_frame_172 </td> <td align="center"> DN </td> </tr>
   </table>

### DE genes with promoter GE UMRs  
#### HuFNSC02  
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec  1 11:35:55 2014 -->
<table border=1>
<tr> <th> name </th> <th> description </th> <th> DE </th>  </tr>
  <tr> <td align="right"> PAX6 </td> <td align="right"> paired_box_6 </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> TMEM132B </td> <td align="right"> transmembrane_protein_132B </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> PID1 </td> <td align="right"> phosphotyrosine_interaction_domain_containing_1 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> ZIC3 </td> <td align="right"> Zic_family_member_3 </td> <td align="center"> UP </td> </tr>
   </table>

#### HuFNSC04  
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec  1 11:35:55 2014 -->
<table border=1>
<tr> <th> name </th> <th> description </th> <th> DE </th>  </tr>
  <tr> <td align="right"> MN1 </td> <td align="right"> meningioma_(disrupted_in_balanced_translocation)_1 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> PNCK </td> <td align="right"> pregnancy_up-regulated_non-ubiquitously_expressed_CaM_kinase </td> <td align="center"> UP </td> </tr>
   </table>

## Overlap UMRs with TFBSs 

* Overlap UMRs with transcription factor binding sites and count No. of overlapping TFBSs for each TF showed similar asymmetry between Cortex and GE in both individuals, with TFBSs enriched in Cortex UMRs for most TFs.   
* However, in general, the correlation of TFBS Cortex UMR vs GE UMR fold change between the two individual is quite low, 0.2.  
* There are 15 TFs that are at least 3-fold enriched in TFBSs overlapping Cortex UMRs compared to GE UMRs as shown below.   

![](./WGBS_files/figure-html/WGBS_TF-1.png) <!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec  1 11:35:55 2014 -->
<table border=1>
<tr> <th> TF </th> <th> Cortex02UMR </th> <th> GE02UMR </th> <th> Ratio02 </th> <th> Cortex04UMR </th> <th> GE04UMR </th> <th> Ratio04 </th>  </tr>
  <tr> <td align="center"> TAF7 </td> <td align="center">  12 </td> <td align="center">   4 </td> <td align="center"> 3.00 </td> <td align="center">  10 </td> <td align="center">   1 </td> <td align="center"> 10.00 </td> </tr>
  <tr> <td align="center"> TFAP2A </td> <td align="center">  30 </td> <td align="center">  10 </td> <td align="center"> 3.00 </td> <td align="center">  20 </td> <td align="center">   3 </td> <td align="center"> 6.67 </td> </tr>
  <tr> <td align="center"> USF1 </td> <td align="center">  67 </td> <td align="center">  21 </td> <td align="center"> 3.19 </td> <td align="center">  31 </td> <td align="center">   6 </td> <td align="center"> 5.17 </td> </tr>
  <tr> <td align="center"> ELF1 </td> <td align="center">  45 </td> <td align="center">  14 </td> <td align="center"> 3.21 </td> <td align="center">  39 </td> <td align="center">   5 </td> <td align="center"> 7.80 </td> </tr>
  <tr> <td align="center"> GATA2 </td> <td align="center"> 116 </td> <td align="center">  36 </td> <td align="center"> 3.22 </td> <td align="center">  43 </td> <td align="center">   4 </td> <td align="center"> 10.75 </td> </tr>
  <tr> <td align="center"> ESR1 </td> <td align="center">  55 </td> <td align="center">  17 </td> <td align="center"> 3.24 </td> <td align="center">   9 </td> <td align="center">   3 </td> <td align="center"> 3.00 </td> </tr>
  <tr> <td align="center"> CTCF </td> <td align="center"> 197 </td> <td align="center">  56 </td> <td align="center"> 3.52 </td> <td align="center"> 105 </td> <td align="center">  16 </td> <td align="center"> 6.56 </td> </tr>
  <tr> <td align="center"> TFAP2C </td> <td align="center">  32 </td> <td align="center">   9 </td> <td align="center"> 3.56 </td> <td align="center">  23 </td> <td align="center">   5 </td> <td align="center"> 4.60 </td> </tr>
  <tr> <td align="center"> MYC </td> <td align="center"> 133 </td> <td align="center">  36 </td> <td align="center"> 3.69 </td> <td align="center">  84 </td> <td align="center">   5 </td> <td align="center"> 16.80 </td> </tr>
  <tr> <td align="center"> FOXA1 </td> <td align="center">  89 </td> <td align="center">  24 </td> <td align="center"> 3.71 </td> <td align="center">  26 </td> <td align="center">   5 </td> <td align="center"> 5.20 </td> </tr>
  <tr> <td align="center"> TCF7L2 </td> <td align="center">  82 </td> <td align="center">  22 </td> <td align="center"> 3.73 </td> <td align="center">  47 </td> <td align="center">   4 </td> <td align="center"> 11.75 </td> </tr>
  <tr> <td align="center"> ZNF263 </td> <td align="center">  31 </td> <td align="center">   8 </td> <td align="center"> 3.88 </td> <td align="center">  23 </td> <td align="center">   5 </td> <td align="center"> 4.60 </td> </tr>
  <tr> <td align="center"> E2F1 </td> <td align="center">  35 </td> <td align="center">   9 </td> <td align="center"> 3.89 </td> <td align="center">  38 </td> <td align="center">   5 </td> <td align="center"> 7.60 </td> </tr>
  <tr> <td align="center"> TAL1 </td> <td align="center">  16 </td> <td align="center">   4 </td> <td align="center"> 4.00 </td> <td align="center">   5 </td> <td align="center">   1 </td> <td align="center"> 5.00 </td> </tr>
  <tr> <td align="center"> GATA3 </td> <td align="center">  63 </td> <td align="center">   8 </td> <td align="center"> 7.88 </td> <td align="center">  18 </td> <td align="center">   1 </td> <td align="center"> 18.00 </td> </tr>
   </table>



