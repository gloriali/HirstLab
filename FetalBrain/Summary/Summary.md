# FetalBrain - Summary
Gloria Li  
October 21, 2014  

Updated: Tue Nov  4 17:59:00 2014



## Introduction
### Neurospheres Cortex vs GE

  * Genetic mutations to epigenetic regulators and epigenomic alterations are amongst the earliest events in brain cell transformation. [PMID: 18772396](http://www.ncbi.nlm.nih.gov/pubmed/18772396) 
  * Defining the epigenomic landscape of normal human brain cells is an important first step towards defining the degree of epigenomic deregulation associated with transformation.
  * Neurospheres provide a powerful model to study neural stem cells, which are thought to be population from which malignant clones arise. 
  * Little is known about epigenetic differences that define neurospheres emerged from distinct brain regions. 

### MZ twins 

  * The genomes of monozygotic twins are genetically identical but epigenetically distinct providing evidence for the influence of environment on the phenotype. 
  * When these epigenetic differences arise during development and their consequence is still unknown. 

## Results
### Neurospheres Cortex derived vs GE derived 
#### Asymmetry between Cortex and GE UMRs  

  * On average, there are , __4.31__-fold enrichment in total UMR frequency in Cortex compared to GE, __4.07__ in HuFNSC02, and __5.31__ in HuFNSC04.    
  * The asymmetry appears to be global, in all chromosomes. It is reproduced in the two individuals, and __supported in MeDIP UMRs__.  
  * __Single CpG level__ differential methylation is __symmetric__, but the asymmetry on __UMR level__ can be reproduced with different cutoffs. However, there are __no apparent differences in UMR length__ between Cortex and GE, suggesting that there are more __orphan GE UM CpGs__ that was not able to form UMRs than in Cortex. 

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:18 2014 -->
<table border=1>
<tr> <th> Sample </th> <th> Total.DMR </th> <th> Hyper.DMR </th> <th> Hypo.DMR </th>  </tr>
  <tr> <td align="center"> Cortex-HuFNSC02_GE-HuFNSC02 </td> <td align="center"> 2178 </td> <td align="center"> 420 </td> <td align="center"> 1758 </td> </tr>
  <tr> <td align="center"> Cortex-HuFNSC04_GE-HuFNSC04 </td> <td align="center"> 646 </td> <td align="center">  91 </td> <td align="center"> 555 </td> </tr>
   </table>
![](./Summary_files/figure-html/WGBS_asymmetry-1.png) ![](./Summary_files/figure-html/WGBS_asymmetry-2.png) 

#### GREAT analysis on Cortex and GE UMRs show brain development terms  

  * UMRs in both Cortex and GE show enrichment in __neuron fate commitment__ biological process, __transcriptional regulation__ activities, __Homeobox__ protein domain, and __abnormal brain development__ mouse phenotype.    
  * Cortex UMRs are also enriched in __forebrain regionalization__ and __pattern specification__ processes.   
  * GREAT enriched terms are also __supported by MeDIP__.    
  * In HuFNSC04, Cortex UMRs show __abnormal brain development__ in Mouse Phenotype, but are also enriched in __kidney-related processes__.   

![](./Summary_files/figure-html/WGBS_GREAT1-1.png) 
![](./Summary_files/figure-html/WGBS_GREAT2-1.png) 
![](./Summary_files/figure-html/WGBS_GREAT3-1.png) 
![](./Summary_files/figure-html/WGBS_GREAT4-1.png) 

#### UMR breakdown in chromatin states - _TBC_

  * __Pending: ChromHMM__

#### DE genes with proximal UMRs show key factors in brain development 

  + On average, there are __40__ UMRs proximally associated with protein-coding genes (TSS +/- 1500bp), __2.85%__ of all UMRs.         
  + There are average __8__ proximal UMRs associated with DE genes, __19.88%__ of all proximal UMRs, __much lower than observed in breast and supported by MeDIP__, _maybe more UMRs in enhancers? (ChromHMM)_. Among them, there are __46.67%__ unique DE genes change in the same direction as the UMRs, __also lower than observed in breast and supported by MeDIP__.         
  + The intersect between two individuals are __significant__ in Cortex UMRs. There are 3 DE genes with proximal Cortex UMR that are shared by both individual:       
    * __GFAP__: Glial Fibrillary Acidic Protein. It is used as a marker to distinguish astrocytes from other glial cells during development. Reported associated with many brain disease, including alexanders disease, [PMID: 11567214](http://www.ncbi.nlm.nih.gov/pubmed?cmd=search&term=11567214&dopt=b); gliomas, UP reduce tumor and induce differentiation, [PMID: 15498217](http://www.ncbi.nlm.nih.gov/pubmed?cmd=search&term=15498217&dopt=b), and astrocytoma, tumor suppressor [PMID: 8339269](http://www.ncbi.nlm.nih.gov/pubmed?cmd=search&term=8339269&dopt=b) etc.    
    * __NFIX__: Nuclear Factor I/X (CCAAT-Binding Transcription Factor). _Also supported by MeDIP._ A transcription factor that binds the palindromic sequence in viral and cellular promoters, capable of activating transcription and replication. It is essential for the development of a number of organ systems including brain and bone, e.g. severe neuroanatomical defects (may function in the repression of neural stem cell proliferation or in cell migration) [PMID: 18477394](http://www.ncbi.nlm.nih.gov/pubmed/18477394), [PMID: 19058033](http://www.ncbi.nlm.nih.gov/pubmed/19058033). 
    * __FEZF1__: FEZ Family Zinc Finger 1 (ZNF312B). Transcription repressor. Involved in the axonal projection and proper termination of olfactory sensory neurons. Regulates non-cell-autonomously the layer formation of the olfactory bulb development and the interneurons. May be required for correct rostral migration of the interneuron progenitors. DNA demethylation and histone acetylation in its promoter activates its expression and oncogene effect in gastric cancer [PMID: 19318583](http://www.ncbi.nlm.nih.gov/pubmed?cmd=search&term=19318583&dopt=b). _We observe DN-regulation in Cortex and Cortex UMR though._     
  + There are no intersect in pc genes with proximal GE UMRs between the two individual.  
  + There are __no significant__ DAVID enrichment terms due to the small number of genes.  
  
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:41 2014 -->
<table border=1>
<tr> <th>  </th> <th> proximal.UMRs </th> <th> unique.genes </th> <th> DE.DMRs </th> <th> unique.DE.genes </th> <th> same.direction </th>  </tr>
  <tr> <td align="center"> GE_UMRs-HuFNSC02 </td> <td align="center">  17 </td> <td align="center">  17 </td> <td align="center">   4 </td> <td align="center">   4 </td> <td align="center">   1 </td> </tr>
  <tr> <td align="center"> Cortex_UMRs-HuFNSC02 </td> <td align="center">  53 </td> <td align="center">  52 </td> <td align="center">  15 </td> <td align="center">  14 </td> <td align="center">   8 </td> </tr>
  <tr> <td align="center"> GE_UMRs-HuFNSC04 </td> <td align="center">   9 </td> <td align="center">  11 </td> <td align="center">   3 </td> <td align="center">   3 </td> <td align="center">   2 </td> </tr>
  <tr> <td align="center"> Cortex_UMRs-HuFNSC04 </td> <td align="center">  82 </td> <td align="center">  84 </td> <td align="center">  10 </td> <td align="center">   9 </td> <td align="center">   3 </td> </tr>
   </table>
![](./Summary_files/figure-html/WGBS_proximal-1.png) ![](./Summary_files/figure-html/WGBS_proximal-2.png) 

##### DE genes with proximal Cortex UMRs  

  * HuFNSC02    

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:47 2014 -->
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

  * HuFNSC04    

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:47 2014 -->
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

##### DE genes with proximal GE UMRs  

  * HuFNSC02    

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:47 2014 -->
<table border=1>
<tr> <th> name </th> <th> description </th> <th> DE </th>  </tr>
  <tr> <td align="right"> PAX6 </td> <td align="right"> paired_box_6 </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> TMEM132B </td> <td align="right"> transmembrane_protein_132B </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> PID1 </td> <td align="right"> phosphotyrosine_interaction_domain_containing_1 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> ZIC3 </td> <td align="right"> Zic_family_member_3 </td> <td align="center"> UP </td> </tr>
   </table>

  * HuFNSC04    

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:47 2014 -->
<table border=1>
<tr> <th> name </th> <th> description </th> <th> DE </th>  </tr>
  <tr> <td align="right"> MN1 </td> <td align="right"> meningioma_(disrupted_in_balanced_translocation)_1 </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> PNCK </td> <td align="right"> pregnancy_up-regulated_non-ubiquitously_expressed_CaM_kinase </td> <td align="center"> UP </td> </tr>
   </table>

#### UMR distal associated genes - _TBC_  

  * __PENDING__   

#### Asymmetry in overlapping UMRs with TFBSs 

  * Overlap UMRs with transcription factor binding sites and count No. of overlapping TFBSs for each TF showed similar asymmetry between Cortex and GE in both individuals, with TFBSs enriched in Cortex UMRs for most TFs.  __Opposite trend in HuFNSC01 MeDIP__.     
  * However, in general, the correlation of TFBS Cortex UMR vs GE UMR fold change between the two individual is quite low, only 0.2. __Same in MeDIP__.     
  * There are 15 TFs that are at least 3-fold enriched in TFBSs overlapping Cortex UMRs compared to GE UMRs as shown below.   

![](./Summary_files/figure-html/WGBS_TF-1.png) <!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:48 2014 -->
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

#### DE genes between Cortex and GE are enriched in neuron development and cell migration   

  * On average, there are __860__ genes differentially expressed between cortex and GE, among them, __454__ are upregulated in cortex, and __406__ are downregulated.    
  * __382__ Cortex up-regulated genes, and __456__ GE up-regulated genes are shared by at least two individuals.    
  * DAVID enrichment analysis show significant enrichment in __neuronal development__ and __cell migration__ terms, GE up-regulated genes are enriched in __EGF-related__ protein domains as well. 
  * There are __22__ Cortex up-regulated genes, and __26__ GE up-regulated genes shared by all four individuals.   
  * Brain / neuron development related Cortex up-regulated genes (shared by all):    
    + SLC1A6: Solute Carrier Family 1 Member 6. Glutamate transporter in cerebellar Purkinje neurons, related to traumatic brain injury.    
    + PCDH20: Protocadherin 20. Although its specific function is undetermined, the cadherin-related neuronal receptor is thought to play a role in the establishment and function of specific cell-cell connections in the brain.   
    + __GFAP__: _see above._   
    + __ZIC5__: Zic Family Member 5. Essential for neural crest development, converting cells from an epidermal fate to a neural crest cell fate, associated with meningiomas [PMID: 20199689](http://www.ncbi.nlm.nih.gov/pubmed?cmd=search&term=20199689&dopt=b). 
    + __ZIC2__: Zic Family Member 2. Related to structural anomaly of the forebrain, and neural tube defects.  
    + NEFM: Neurofilament, Medium Polypeptide. Neurofilaments comprise the axoskeleton and functionally maintain neuronal caliber. They may also play a role in intracellular transport to axons and dendrites. This protein
  is commonly used as a biomarker of neuronal damage.         
    + UNC5C: Unc-5 Homolog C. Receptor for netrin required for axon guidance. Mediates axon repulsion of neuronal growth cones in the developing nervous system upon ligand binding.       
    + __ADAM19__: ADAM Metallopeptidase Domain 19. Involved in neurogenesis, cell migration, cell adhesion, cell-cell and cell-matrix interactions, serves as a marker for dendritic cell differentiation. UP in primary brain tumors [PMID: 16772875](http://www.ncbi.nlm.nih.gov/pubmed?cmd=search&term=16772875&dopt=b).       
  * Brain / neuron development related GE up-regulated genes (shared by all):     
    + VAX1: Ventral Anterior Homeobox 1. Transcription factor that may function in dorsoventral specification of the forebrain. Required for axon
  guidance and major tract formation in the developing forebrain.     
    + NRXN3: Neurexin 3. Neuronal cell surface protein that may be involved in cell recognition and cell adhesion.      
    + __FEZF1__: _see above._     
    + SIX3: SIX Homeobox 3. Plays a role in eye development.      
    + ODZ1: Odz, Odd Oz/Ten-M Homolog 1. Involved in neural development, regulating the establishment of proper connectivity within the nervous system.     
    + __EPHA3__: EPH Receptor A3. Involved in short-range contact-mediated axon guidance. UP in glioblastoma [PMID: 23410976](http://www.ncbi.nlm.nih.gov/pubmed?cmd=search&term=23410976&dopt=b).      
    + __OTX2__: Orthodenticle Homeobox 2. Plays a role in the development of the brain and the sense organs, also influences the proliferation and differentiation of dopaminergic neuronal progenitor cells during mitosis. UP in medulloblastomas [PMID: 20028867](http://www.ncbi.nlm.nih.gov/pubmed?cmd=search&term=20028867&dopt=b).     
    + ASTN1: Astrotactin 1. A neuronal adhesion molecule required for glial-guided migration of young postmitotic neuroblasts in cortical regions of developing brain, including cerebrum, hippocampus, cerebellum, and olfactory bulb (Fink et al., 1995).      

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:50 2014 -->
<table border=1>
<tr> <th>  </th> <th> UP </th> <th> DN </th> <th> DE </th>  </tr>
  <tr> <td> HuFNSC01 </td> <td align="center"> 403 </td> <td align="center"> 508 </td> <td align="center"> 911 </td> </tr>
  <tr> <td> HuFNSC02 </td> <td align="center"> 588 </td> <td align="center"> 640 </td> <td align="center"> 1228 </td> </tr>
  <tr> <td> HuFNSC03 </td> <td align="center"> 447 </td> <td align="center"> 227 </td> <td align="center"> 674 </td> </tr>
  <tr> <td> HuFNSC04 </td> <td align="center"> 378 </td> <td align="center"> 249 </td> <td align="center"> 627 </td> </tr>
   </table>
![](./Summary_files/figure-html/DE_Cortex_GE-1.png) ![](./Summary_files/figure-html/DE_Cortex_GE-2.png) <!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:50 2014 -->
<table border=1>
<tr> <th> DE </th> <th> name </th> <th> description </th>  </tr>
  <tr> <td align="center"> UP </td> <td align="right"> SLC1A6 </td> <td align="right"> solute_carrier_family_1_(high_affinity_aspartate/glutamate_transporter),_member_6 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> PCDH20 </td> <td align="right"> protocadherin_20 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> C9orf64 </td> <td align="right"> chromosome_9_open_reading_frame_64 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> GFAP </td> <td align="right"> glial_fibrillary_acidic_protein </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> ZIC5 </td> <td align="right"> Zic_family_member_5 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> NEFM </td> <td align="right"> neurofilament,_medium_polypeptide </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> QPCT </td> <td align="right"> glutaminyl-peptide_cyclotransferase </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> UNC5C </td> <td align="right"> unc-5_homolog_C_(C._elegans) </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> ZIC2 </td> <td align="right"> Zic_family_member_2 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> XKR8 </td> <td align="right"> XK,_Kell_blood_group_complex_subunit-related_family,_member_8 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> PYGL </td> <td align="right"> phosphorylase,_glycogen,_liver </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> SERINC2 </td> <td align="right"> serine_incorporator_2 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> KIAA1239 </td> <td align="right"> KIAA1239 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> B3GALT1 </td> <td align="right"> UDP-Gal:betaGlcNAc_beta_1,3-galactosyltransferase,_polypeptide_1 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> CCDC48 </td> <td align="right"> coiled-coil_domain_containing_48 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> NOS1AP </td> <td align="right"> nitric_oxide_synthase_1_(neuronal)_adaptor_protein </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> ADAM19 </td> <td align="right"> ADAM_metallopeptidase_domain_19 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> H1F0 </td> <td align="right"> H1_histone_family,_member_0 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> ZIC3 </td> <td align="right"> Zic_family_member_3 </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> NT5E </td> <td align="right"> 5'-nucleotidase,_ecto_(CD73) </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> SYNM </td> <td align="right"> synemin,_intermediate_filament_protein </td> </tr>
  <tr> <td align="center"> UP </td> <td align="right"> C1orf226 </td> <td align="right"> chromosome_1_open_reading_frame_226 </td> </tr>
   </table>
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:50 2014 -->
<table border=1>
<tr> <th> DE </th> <th> name </th> <th> description </th>  </tr>
  <tr> <td align="center"> DN </td> <td align="right"> DPPA4 </td> <td align="right"> developmental_pluripotency_associated_4 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> VAX1 </td> <td align="right"> ventral_anterior_homeobox_1 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> SHISA6 </td> <td align="right"> shisa_homolog_6_(Xenopus_laevis) </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> NRXN3 </td> <td align="right"> neurexin_3 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> FGFR2 </td> <td align="right"> fibroblast_growth_factor_receptor_2 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> ST8SIA2 </td> <td align="right"> ST8_alpha-N-acetyl-neuraminide_alpha-2,8-sialyltransferase_2 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> SLC32A1 </td> <td align="right"> solute_carrier_family_32_(GABA_vesicular_transporter),_member_1 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> FEZF1 </td> <td align="right"> FEZ_family_zinc_finger_1 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> AJAP1 </td> <td align="right"> adherens_junctions_associated_protein_1 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> ZNF703 </td> <td align="right"> zinc_finger_protein_703 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> VSIG10L </td> <td align="right"> V-set_and_immunoglobulin_domain_containing_10_like </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> SIX3 </td> <td align="right"> SIX_homeobox_3 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> SORL1 </td> <td align="right"> sortilin-related_receptor,_L(DLR_class)_A_repeats_containing </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> LMO1 </td> <td align="right"> LIM_domain_only_1_(rhombotin_1) </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> ODZ1 </td> <td align="right"> odz,_odd_Oz/ten-m_homolog_1_(Drosophila) </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> EPHA3 </td> <td align="right"> EPH_receptor_A3 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> TIMP3 </td> <td align="right"> TIMP_metallopeptidase_inhibitor_3 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> CAMK2N1 </td> <td align="right"> calcium/calmodulin-dependent_protein_kinase_II_inhibitor_1 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> EEF1A2 </td> <td align="right"> eukaryotic_translation_elongation_factor_1_alpha_2 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> OTX2 </td> <td align="right"> orthodenticle_homeobox_2 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> CHL1 </td> <td align="right"> cell_adhesion_molecule_with_homology_to_L1CAM_(close_homolog_of_L1) </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> PCSK1N </td> <td align="right"> proprotein_convertase_subtilisin/kexin_type_1_inhibitor </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> LMO2 </td> <td align="right"> LIM_domain_only_2_(rhombotin-like_1) </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> ASTN1 </td> <td align="right"> astrotactin_1 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> BASP1 </td> <td align="right"> brain_abundant,_membrane_attached_signal_protein_1 </td> </tr>
  <tr> <td align="center"> DN </td> <td align="right"> ADCYAP1R1 </td> <td align="right"> adenylate_cyclase_activating_polypeptide_1_(pituitary)_receptor_type_I </td> </tr>
   </table>

#### Isoforms between Cortex and GE are enriched in cell signaling proteins

  * On average, __2054__ genes are identified as isoforms between cortex and GE in each individual. __2352__ genes are shared in at least two individuals.       
  * Individual specific isoforms between cortex and GE have __no__ significantly enriched terms, suggesting they are more likely random events without biological functions.          
  * Isoforms shared by at least two individuals are enriched in terms related to __cellular signaling__. InterPro shows enrichment in __EGF protein domain__.         

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:52 2014 -->
<table border=1>
<tr> <th>  </th> <th> DE_genes </th> <th> DE_exons </th> <th> with_expressed_genes </th> <th> isoform_exons </th> <th> exclude_DE_genes </th> <th> isoform_genes </th>  </tr>
  <tr> <td> cortex01_GE01_summary </td> <td align="center"> 911 </td> <td align="center"> 32372 </td> <td align="center"> 18968 </td> <td align="center"> 8440 </td> <td align="center"> 7962 </td> <td align="center"> 2447 </td> </tr>
  <tr> <td> cortex02_GE02_summary </td> <td align="center"> 1228 </td> <td align="center"> 35196 </td> <td align="center"> 21880 </td> <td align="center"> 8163 </td> <td align="center"> 7374 </td> <td align="center"> 2298 </td> </tr>
  <tr> <td> cortex03_GE03_summary </td> <td align="center"> 674 </td> <td align="center"> 29617 </td> <td align="center"> 13746 </td> <td align="center"> 6401 </td> <td align="center"> 6022 </td> <td align="center"> 2086 </td> </tr>
  <tr> <td> cortex04_GE04_summary </td> <td align="center"> 627 </td> <td align="center"> 22386 </td> <td align="center"> 11253 </td> <td align="center"> 4323 </td> <td align="center"> 4259 </td> <td align="center"> 1386 </td> </tr>
   </table>
![](./Summary_files/figure-html/Isoform_Cortex_GE-1.png) 

#### Intron retention between Cortex and GE - _TBC_

  * _PENDING_

#### mCpG is a stable mark for exon transcription during development     

  * mCpG at cassette exon boundaries has similar pattern to expressed in both exons.   
  * mCpG for cassette exons between cortex and GE neurospheres shows no significant differences.   
  * Results from both WGBS (HuFNSC02 & HuFNSC04) and MeDIP (HuFNSC01 & HuFNSC02) support the assumption that mCpG is a stable mark for exon transcription during development.   
  * mCpG exon marking is established between neurospheres. _Needs further validation against H1_.   

![](./Summary_files/figure-html/mCpG_Cortex_GE-1.png) ![](./Summary_files/figure-html/mCpG_Cortex_GE-2.png) ![](./Summary_files/figure-html/mCpG_Cortex_GE-3.png) ![](./Summary_files/figure-html/mCpG_Cortex_GE-4.png) 

#### H3K36me3 in exon bodies - normalization issue? 

  * H3K36me3 in expressed in both / not expressed exons shows no significant differences between HuFNSC01 and HuFNSC02.    
  * H3K36me3 in cassette exons in HuFNSC01 are enriched in GE compared to cortex. However, it is not reproduced in HuFNSC02, where there is no significant differences between cortex and GE. _Not sure what to make of this. Are there any potential bias?_    

![](./Summary_files/figure-html/H3K36me3_Cortex_GE-1.png) ![](./Summary_files/figure-html/H3K36me3_Cortex_GE-2.png) 

### MZ twins  
#### UMR asymmetry between MZ twins in Brain and Cortex 

  + There is an __asymmetry__ between UMRs in HuFNSC01 and HuFNSC02 in the Brain and Cortex neurosphere, but __not__ in GE neurosphere. Fold change in UMR frequency HuFNSC02/HuFNSC01 in Brain is __2.41__, in Cortex is __1.83__, and in GE is __0.899__.    
  
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:00:58 2014 -->
<table border=1>
<tr> <th> Sample </th> <th> Total.DMR </th> <th> Hyper.DMR </th> <th> Hypo.DMR </th>  </tr>
  <tr> <td> Brain-HuFNSC01_Brain-HuFNSC02 </td> <td align="center"> 4472 </td> <td align="center"> 2750 </td> <td align="center"> 1722 </td> </tr>
  <tr> <td> Cortex-HuFNSC01_Cortex-HuFNSC02 </td> <td align="center"> 3161 </td> <td align="center"> 1758 </td> <td align="center"> 1403 </td> </tr>
  <tr> <td> GE-HuFNSC01_GE-HuFNSC02 </td> <td align="center"> 2716 </td> <td align="center"> 1136 </td> <td align="center"> 1580 </td> </tr>
   </table>
![](./Summary_files/figure-html/MeDIP_MZ_asymmetry-1.png) ![](./Summary_files/figure-html/MeDIP_MZ_asymmetry-2.png) 

#### GREAT analysis on MZ UMRs are enriched in Homeobox and brain development  

  * __Homeobox__ protein domain is enriched in all UMR lists.   
  * __Brain development related__ biological processes are enriched in most lists.   

![](./Summary_files/figure-html/MeDIP_MZ_GREAT1-1.png) ![](./Summary_files/figure-html/MeDIP_MZ_GREAT1-2.png) ![](./Summary_files/figure-html/MeDIP_MZ_GREAT1-3.png) 
![](./Summary_files/figure-html/MeDIP_MZ_GREAT2-1.png) ![](./Summary_files/figure-html/MeDIP_MZ_GREAT2-2.png) ![](./Summary_files/figure-html/MeDIP_MZ_GREAT2-3.png) 

#### UMR breakdown in chromatin states - _TBC_

  * __Pending: ChromHMM__

#### MZ UMRs in Brain have less CpGs in CGIs  

  + Majority of UMR CpGs overlap with genebody, on average __62.59%__. And UMR CpGs are __3-fold__ enriched in promoter regions, with on average __34.29%__ UMR CpGs are in promoter regions, __39.54%__ of CpGs in UMRs overlap with CGIs, __5.34-fold__ than expected by random. 
  + Brain seems to have less CGI/promoter UMRs.   
  
<!-- For the entire genome, 3727169 out of 28217448 CpGs overlap with TSS +/- 1500bp promoter regions -->  
<!-- For the entire genome, 2089538 out of 28217448 CpGs overlap with CGIs -->

![](./Summary_files/figure-html/MeDIP_MZ_breakdown-1.png) 

#### DE genes with proximal UMRs between MZ are cell type specific  

  + On average, there are __405__ UMRs proximally (TSS +/- 1500bp) associated with protein-coding genes, __11.75%__ of all UMRs, _similar to UMRs between Cortex and GE_.         
  + On average, there are __17__ proximal UMRs associated with DE genes, __4.24%__ of all proximal DMRs, _much less than UMRs between Cortex and GE (less functional?)_. Among them, there are __52.53%__ unique DE genes change in the same direction as the UMRs, _similar to UMRs between Cortex and GE_.         
  + Proximal UMR genes are mostly __cell type specific__, but the overlap is still statistically significant.  
  + No overlap between Brain and Cortex for proximal HuFNSC01 UMR DE genes, only one for HuFNSC02, __neuropeptide_Y__.   
  + HuFNSC01 UMR DE genes associated with brain development:   
    * __WNT pathway__: SFRP2, WNT3.      
    * EMX2: Empty Spiracles Homeobox 2. Transcription factor. Acts to generate the boundary between the roof and archipallium in the developing brain. May function in combinations with OTX1/2 to specify cell fates in the developing central nervous system.     
    * MDGA1: MAM Domain-Containing Glycosylphosphatidylinositol Anchor Protein 1. Required for radial migration of cortical neurons in the superficial layer of the neocortex.   
    * VAX1: _see above._   
    * DLX1: Distal-Less Homeobox 1. Play a role in the control of craniofacial patterning and the differentiation and survival of inhibitory neurons in the forebrain.    
    * __OLIG1__: Oligodendrocyte Transcription Factor 1. Promotes formation and maturation of oligodendrocytes, especially within the brain. Associated with human glial brain tumors [PMID: 11526205](http://www.ncbi.nlm.nih.gov/pubmed?cmd=search&term=11526205&dopt=b).   
  + HuFNSC02 UMR DE genes associated with brain development:   
    * __WNT pathway__: SFRP1, SFRP2, WNT7A.        
    * EPHA8: EPH Receptor A8. Plays a role in short-range contact-mediated axonal guidance during development of the mammalian nervous system.    
    * NRGN: Neurogranin (Protein Kinase C Substrate, RC3). Acts as a "third messenger" substrate of protein kinase C-mediated molecular cascades during synaptic development and remodeling.    
    * C1QL4: Complement Component 1, Q Subcomponent-Like 4. May regulate the number of excitatory synapses that are formed on hippocampus neurons.       
    * SNCB: Synuclein, Beta. May play a role in neuronal plasticity.       
    * SEMA3D: Sema Domain, Immunoglobulin Domain (Ig), Short Basic Domain, Secreted, (Semaphorin) 3D. Induces the collapse and paralysis of neuronal growth cones. Could potentially act as repulsive cues toward specific neuronal populations.    
    * VGF: May be involved in the regulation of cell-cell interactions or in synatogenesis during the maturation of the nervous system.   

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:01:13 2014 -->
<table border=1>
<tr> <th>  </th> <th> DMRs </th> <th> unique.genes </th> <th> DE.DMRs </th> <th> unique.DE.genes </th> <th> same.direction </th>  </tr>
  <tr> <td> Brain01_Brain02_hyper </td> <td align="center"> 432 </td> <td align="center"> 473 </td> <td align="center">  26 </td> <td align="center">  26 </td> <td align="center">  10 </td> </tr>
  <tr> <td> Brain01_Brain02_hypo </td> <td align="center"> 367 </td> <td align="center"> 409 </td> <td align="center">  20 </td> <td align="center">  20 </td> <td align="center">  13 </td> </tr>
  <tr> <td> Cortex01_Cortex02_hyper </td> <td align="center"> 587 </td> <td align="center"> 659 </td> <td align="center">  30 </td> <td align="center">  28 </td> <td align="center">  17 </td> </tr>
  <tr> <td> Cortex01_Cortex02_hypo </td> <td align="center"> 316 </td> <td align="center"> 342 </td> <td align="center">  24 </td> <td align="center">  22 </td> <td align="center">  11 </td> </tr>
  <tr> <td> GE01_GE02_hyper </td> <td align="center"> 369 </td> <td align="center"> 402 </td> <td align="center">   1 </td> <td align="center">   1 </td> <td align="center">   1 </td> </tr>
  <tr> <td> GE01_GE02_hypo </td> <td align="center"> 361 </td> <td align="center"> 398 </td> <td align="center">   2 </td> <td align="center">   2 </td> <td align="center">   0 </td> </tr>
   </table>
![](./Summary_files/figure-html/MeDIP_MZ_proximal-1.png) ![](./Summary_files/figure-html/MeDIP_MZ_proximal-2.png) 

##### HuFNSC01 UMRs proximal associated DE genes

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:01:13 2014 -->
<table border=1>
<tr> <th> CellType </th> <th> name </th> <th> description </th> <th> DM </th> <th> DE </th>  </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> SAMD11 </td> <td align="right"> sterile_alpha_motif_domain_containing_11 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> ITGA8 </td> <td align="right"> integrin,_alpha_8 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> EMX2 </td> <td align="right"> empty_spiracles_homeobox_2 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> ARNTL2 </td> <td align="right"> aryl_hydrocarbon_receptor_nuclear_translocator-like_2 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> CPNE8 </td> <td align="right"> copine_VIII </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> PCDH17 </td> <td align="right"> protocadherin_17 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right">  </td> <td align="right"> Putative_3-phosphoinositide-dependent_protein_kinase_2 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> CPNE7 </td> <td align="right"> copine_VII </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> PKDREJ </td> <td align="right"> polycystic_kidney_disease_(polycystin)_and_REJ_homolog_(sperm_receptor_for_egg_jelly_homolog,_sea_urchin) </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> THRB </td> <td align="right"> thyroid_hormone_receptor,_beta </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> SFRP2 </td> <td align="right"> secreted_frizzled-related_protein_2 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> IER3 </td> <td align="right"> immediate_early_response_3 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> ZBTB12 </td> <td align="right"> zinc_finger_and_BTB_domain_containing_12 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> MDGA1 </td> <td align="right"> MAM_domain_containing_glycosylphosphatidylinositol_anchor_1 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> EPHA1 </td> <td align="right"> EPH_receptor_A1 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right">  </td> <td align="right"> Protein_kinase-like_protein_SgK196 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> TOX </td> <td align="right"> thymocyte_selection-associated_high_mobility_group_box </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> NCOA2 </td> <td align="right"> nuclear_receptor_coactivator_2 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> GPR64 </td> <td align="right"> G_protein-coupled_receptor_64 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> TRPC5 </td> <td align="right"> transient_receptor_potential_cation_channel,_subfamily_C,_member_5 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> ELTD1 </td> <td align="right"> EGF,_latrophilin_and_seven_transmembrane_domain_containing_1 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> CRABP2 </td> <td align="right"> cellular_retinoic_acid_binding_protein_2 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> SPOCK2 </td> <td align="right"> sparc/osteonectin,_cwcv_and_kazal-like_domains_proteoglycan_(testican)_2 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> VAX1 </td> <td align="right"> ventral_anterior_homeobox_1 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> LRRC10B </td> <td align="right"> leucine_rich_repeat_containing_10B </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> NDUFA4L2 </td> <td align="right"> NADH_dehydrogenase_(ubiquinone)_1_alpha_subcomplex,_4-like_2 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> FOXO1 </td> <td align="right"> forkhead_box_O1 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> C16orf74 </td> <td align="right"> chromosome_16_open_reading_frame_74 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> WNT3 </td> <td align="right"> wingless-type_MMTV_integration_site_family,_member_3 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> CHST9 </td> <td align="right"> carbohydrate_(N-acetylgalactosamine_4-0)_sulfotransferase_9 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> CD97 </td> <td align="right"> CD97_molecule </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> RNF144A </td> <td align="right"> ring_finger_protein_144A </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> DLX1 </td> <td align="right"> distal-less_homeobox_1 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> OLIG1 </td> <td align="right"> oligodendrocyte_transcription_factor_1 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> SCUBE1 </td> <td align="right"> signal_peptide,_CUB_domain,_EGF-like_1 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> GRIA1 </td> <td align="right"> glutamate_receptor,_ionotropic,_AMPA_1 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> MMD2 </td> <td align="right"> monocyte_to_macrophage_differentiation-associated_2 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> SEMA3E </td> <td align="right"> sema_domain,_immunoglobulin_domain_(Ig),_short_basic_domain,_secreted,_(semaphorin)_3E </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> VIPR2 </td> <td align="right"> vasoactive_intestinal_peptide_receptor_2 </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> AP1S2 </td> <td align="right"> adaptor-related_protein_complex_1,_sigma_2_subunit </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> SH3KBP1 </td> <td align="right"> SH3-domain_kinase_binding_protein_1 </td> <td align="center"> hypo </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> GE </td> <td align="right"> FAM5B </td> <td align="right"> family_with_sequence_similarity_5,_member_B </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> GE </td> <td align="right"> DLL1 </td> <td align="right"> delta-like_1_(Drosophila) </td> <td align="center"> hypo </td> <td align="center"> DN </td> </tr>
   </table>

##### HuFNSC02 UMRs proximal associated DE genes

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:01:13 2014 -->
<table border=1>
<tr> <th> CellType </th> <th> name </th> <th> description </th> <th> DM </th> <th> DE </th>  </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> EPHA8 </td> <td align="right"> EPH_receptor_A8 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> IFI6 </td> <td align="right"> interferon,_alpha-inducible_protein_6 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> BMP8A </td> <td align="right"> bone_morphogenetic_protein_8a </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> ALX3 </td> <td align="right"> ALX_homeobox_3 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> NRGN </td> <td align="right"> neurogranin_(protein_kinase_C_substrate,_RC3) </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> C1QL4 </td> <td align="right"> complement_component_1,_q_subcomponent-like_4 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> CSRP2 </td> <td align="right"> cysteine_and_glycine-rich_protein_2 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> TTYH2 </td> <td align="right"> tweety_homolog_2_(Drosophila) </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> DTNB </td> <td align="right"> dystrobrevin,_beta </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> CXCR7 </td> <td align="right"> chemokine_(C-X-C_motif)_receptor_7 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> VHL </td> <td align="right"> von_Hippel-Lindau_tumor_suppressor </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> EPHA6 </td> <td align="right"> EPH_receptor_A6 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> CAMK2N2 </td> <td align="right"> calcium/calmodulin-dependent_protein_kinase_II_inhibitor_2 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> TBC1D1 </td> <td align="right"> TBC1_(tre-2/USP6,_BUB2,_cdc16)_domain_family,_member_1 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> NMU </td> <td align="right"> neuromedin_U </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> FSTL5 </td> <td align="right"> follistatin-like_5 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> CARTPT </td> <td align="right"> CART_prepropeptide </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> EDIL3 </td> <td align="right"> EGF-like_repeats_and_discoidin_I-like_domains_3 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> CCDC90A </td> <td align="right"> coiled-coil_domain_containing_90A </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right">  </td> <td align="right"> LOC401296_proteinUncharacterized_protein </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> NPY </td> <td align="right"> neuropeptide_Y </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> SFRP1 </td> <td align="right"> secreted_frizzled-related_protein_1 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> TOX </td> <td align="right"> thymocyte_selection-associated_high_mobility_group_box </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> ENTPD2 </td> <td align="right"> ectonucleoside_triphosphate_diphosphohydrolase_2 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> KLHL4 </td> <td align="right"> kelch-like_4_(Drosophila) </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Brain </td> <td align="right"> SOX3 </td> <td align="right"> SRY_(sex_determining_region_Y)-box_3 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> DMRTA2 </td> <td align="right"> DMRT-like_family_A2 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> NPR1 </td> <td align="right"> natriuretic_peptide_receptor_A/guanylate_cyclase_A_(atrionatriuretic_peptide_receptor_A) </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> PTPRE </td> <td align="right"> protein_tyrosine_phosphatase,_receptor_type,_E </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> PDE3B </td> <td align="right"> phosphodiesterase_3B,_cGMP-inhibited </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> LDHA </td> <td align="right"> lactate_dehydrogenase_A </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> GNG3 </td> <td align="right"> guanine_nucleotide_binding_protein_(G_protein),_gamma_3 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> METTL7B </td> <td align="right"> methyltransferase_like_7B </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> NME3 </td> <td align="right"> non-metastatic_cells_3,_protein_expressed_in </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> GRIN2A </td> <td align="right"> glutamate_receptor,_ionotropic,_N-methyl_D-aspartate_2A </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> NOL3 </td> <td align="right"> nucleolar_protein_3_(apoptosis_repressor_with_CARD_domain) </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> NR1D1 </td> <td align="right"> nuclear_receptor_subfamily_1,_group_D,_member_1 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> IGFBP4 </td> <td align="right"> insulin-like_growth_factor_binding_protein_4 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> TBX2 </td> <td align="right"> T-box_2 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> LBH </td> <td align="right"> limb_bud_and_heart_development_homolog_(mouse) </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> INSIG2 </td> <td align="right"> insulin_induced_gene_2 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> WNT7A </td> <td align="right"> wingless-type_MMTV_integration_site_family,_member_7A </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> SNCA </td> <td align="right"> synuclein,_alpha_(non_A4_component_of_amyloid_precursor) </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> NPNT </td> <td align="right"> nephronectin </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> SFRP2 </td> <td align="right"> secreted_frizzled-related_protein_2 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> PRR16 </td> <td align="right"> proline_rich_16 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> STC2 </td> <td align="right"> stanniocalcin_2 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> SNCB </td> <td align="right"> synuclein,_beta </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> ITPR3 </td> <td align="right"> inositol_1,4,5-trisphosphate_receptor,_type_3 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> NPY </td> <td align="right"> neuropeptide_Y </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> SEMA3D </td> <td align="right"> sema_domain,_immunoglobulin_domain_(Ig),_short_basic_domain,_secreted,_(semaphorin)_3D </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> VGF </td> <td align="right"> VGF_nerve_growth_factor_inducible </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> CAV1 </td> <td align="right"> caveolin_1,_caveolae_protein,_22kDa </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> Cortex </td> <td align="right"> SHROOM2 </td> <td align="right"> shroom_family_member_2 </td> <td align="center"> hyper </td> <td align="center"> UP </td> </tr>
  <tr> <td align="center"> GE </td> <td align="right"> DBC1 </td> <td align="right"> deleted_in_bladder_cancer_1 </td> <td align="center"> hyper </td> <td align="center"> DN </td> </tr>
   </table>

#### UMR distal associated genes - _TBC_  

  * __PENDING__   

#### Overlapping UMRs with TFBSs show asymmetric in Brain and Cortex  

  * Overlap UMRs with transcription factor binding sites and count No. of overlapping TFBSs for each TF showed similar asymmetry in Brain and Cortex, but is symmetric in GE. The correlation between Brain and Cortex is also very low, __0.15__.    
  * With 2-fold change cutoff, there are 18 TFs enriched in HuFNSC02 in both Brain and Cortex.   

![](./Summary_files/figure-html/MeDIP_MZ_TF-1.png) <!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:01:14 2014 -->
<table border=1>
<tr> <th> TF </th> <th> Brain.hypo </th> <th> Brain.hyper </th> <th> Ratio.Brain </th> <th> Cortex.hypo </th> <th> Cortex.hyper </th> <th> Ratio.Cortex </th>  </tr>
  <tr> <td align="center"> ESRRA </td> <td align="center">   1 </td> <td align="center">   7 </td> <td align="center"> 0.14 </td> <td align="center">   1 </td> <td align="center">   5 </td> <td align="center"> 0.20 </td> </tr>
  <tr> <td align="center"> TAL1 </td> <td align="center">  10 </td> <td align="center">  42 </td> <td align="center"> 0.24 </td> <td align="center">   9 </td> <td align="center">  34 </td> <td align="center"> 0.26 </td> </tr>
  <tr> <td align="center"> HNF4A </td> <td align="center">  15 </td> <td align="center">  49 </td> <td align="center"> 0.31 </td> <td align="center">  11 </td> <td align="center">  27 </td> <td align="center"> 0.41 </td> </tr>
  <tr> <td align="center"> HNF4G </td> <td align="center">  13 </td> <td align="center">  41 </td> <td align="center"> 0.32 </td> <td align="center">  10 </td> <td align="center">  29 </td> <td align="center"> 0.34 </td> </tr>
  <tr> <td align="center"> ZNF217 </td> <td align="center">   6 </td> <td align="center">  17 </td> <td align="center"> 0.35 </td> <td align="center">   6 </td> <td align="center">  14 </td> <td align="center"> 0.43 </td> </tr>
  <tr> <td align="center"> EBF1 </td> <td align="center">  49 </td> <td align="center"> 130 </td> <td align="center"> 0.38 </td> <td align="center">  48 </td> <td align="center">  98 </td> <td align="center"> 0.49 </td> </tr>
  <tr> <td align="center"> FOS </td> <td align="center">  47 </td> <td align="center"> 123 </td> <td align="center"> 0.38 </td> <td align="center">  43 </td> <td align="center"> 118 </td> <td align="center"> 0.36 </td> </tr>
  <tr> <td align="center"> GATA1 </td> <td align="center">  31 </td> <td align="center">  73 </td> <td align="center"> 0.42 </td> <td align="center">  32 </td> <td align="center">  80 </td> <td align="center"> 0.40 </td> </tr>
  <tr> <td align="center"> HSF1 </td> <td align="center">   3 </td> <td align="center">   7 </td> <td align="center"> 0.43 </td> <td align="center">   2 </td> <td align="center">   6 </td> <td align="center"> 0.33 </td> </tr>
  <tr> <td align="center"> FOXA2 </td> <td align="center">  26 </td> <td align="center">  60 </td> <td align="center"> 0.43 </td> <td align="center">  14 </td> <td align="center">  36 </td> <td align="center"> 0.39 </td> </tr>
  <tr> <td align="center"> GATA3 </td> <td align="center">  35 </td> <td align="center">  80 </td> <td align="center"> 0.44 </td> <td align="center">  22 </td> <td align="center">  52 </td> <td align="center"> 0.42 </td> </tr>
  <tr> <td align="center"> GATA2 </td> <td align="center">  58 </td> <td align="center"> 130 </td> <td align="center"> 0.45 </td> <td align="center">  39 </td> <td align="center">  96 </td> <td align="center"> 0.41 </td> </tr>
  <tr> <td align="center"> STAT3 </td> <td align="center">  41 </td> <td align="center">  91 </td> <td align="center"> 0.45 </td> <td align="center">  26 </td> <td align="center">  87 </td> <td align="center"> 0.30 </td> </tr>
  <tr> <td align="center"> CEBPD </td> <td align="center">  22 </td> <td align="center">  48 </td> <td align="center"> 0.46 </td> <td align="center">  23 </td> <td align="center">  64 </td> <td align="center"> 0.36 </td> </tr>
  <tr> <td align="center"> SMC3 </td> <td align="center">  71 </td> <td align="center"> 150 </td> <td align="center"> 0.47 </td> <td align="center">  62 </td> <td align="center"> 144 </td> <td align="center"> 0.43 </td> </tr>
  <tr> <td align="center"> NR3C1 </td> <td align="center">  46 </td> <td align="center">  95 </td> <td align="center"> 0.48 </td> <td align="center">  39 </td> <td align="center">  90 </td> <td align="center"> 0.43 </td> </tr>
  <tr> <td align="center"> RXRA </td> <td align="center">  23 </td> <td align="center">  47 </td> <td align="center"> 0.49 </td> <td align="center">  16 </td> <td align="center">  35 </td> <td align="center"> 0.46 </td> </tr>
  <tr> <td align="center"> HDAC6 </td> <td align="center">   1 </td> <td align="center">   2 </td> <td align="center"> 0.50 </td> <td align="center">   2 </td> <td align="center">  11 </td> <td align="center"> 0.18 </td> </tr>
   </table>

#### DE genes are cell type specific

  * On average, there are __470__ DE genes across three cells types.   
  * Majority of DE genes are cell type specific, only __98__ are shared between any two cell types.   
  * DE genes in Brain is asymmetric, _maybe due to cell heterogenity?_   
  * There are much fewer DE genes in GE.    
  * DAVID enrichment analysis between MZ twins in brain and cortex show similar GO term in __brain development__, but there is no significantly enriched terms in GE.    

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:01:14 2014 -->
<table border=1>
<tr> <th>  </th> <th> UP </th> <th> DN </th> <th> DE </th>  </tr>
  <tr> <td> brain01_brain02 </td> <td align="center"> 461 </td> <td align="center"> 181 </td> <td align="center"> 642 </td> </tr>
  <tr> <td> cortex01_cortex02 </td> <td align="center"> 248 </td> <td align="center"> 348 </td> <td align="center"> 596 </td> </tr>
  <tr> <td> GE01_GE02 </td> <td align="center">  99 </td> <td align="center">  74 </td> <td align="center"> 173 </td> </tr>
   </table>
![](./Summary_files/figure-html/DE_MZ-1.png) 

![](./Summary_files/figure-html/DE_MZ_DAVID-1.png) ![](./Summary_files/figure-html/DE_MZ_DAVID-2.png) 

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:01:15 2014 -->
<table border=1>
<tr> <th> name </th> <th> description </th> <th> Brain </th> <th> Cortex </th> <th> GE </th>  </tr>
  <tr> <td align="right"> LMO1 </td> <td align="right"> LIM_domain_only_1_(rhombotin_1) </td> <td align="center"> UP </td> <td align="center"> DN </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> CXCR7 </td> <td align="right"> chemokine_(C-X-C_motif)_receptor_7 </td> <td align="center"> UP </td> <td align="center"> UP </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> SPP1 </td> <td align="right"> secreted_phosphoprotein_1 </td> <td align="center"> UP </td> <td align="center"> DN </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> NPY </td> <td align="right"> neuropeptide_Y </td> <td align="center"> UP </td> <td align="center"> DN </td> <td align="center"> DN </td> </tr>
  <tr> <td align="right"> COL2A1 </td> <td align="right"> collagen,_type_II,_alpha_1 </td> <td align="center"> DN </td> <td align="center"> UP </td> <td align="center"> UP </td> </tr>
  <tr> <td align="right"> BCL6 </td> <td align="right"> B-cell_CLL/lymphoma_6 </td> <td align="center"> DN </td> <td align="center"> DN </td> <td align="center"> DN </td> </tr>
   </table>

#### Isoforms between MZ are enriched in cell signaling in neurospheres and immune response in Brain

  * On average, __2617__ genes are identified as isoforms between HuFNSC01 and HuFNSC02 in each cell type. __796__ genes are shared by all three cell types.              
  * On average, __1724__ genes are identified as isoforms between HuFNSC03 and HuFNSC04 in each cell type. __927__ genes are shared between two cell types.            
  * Different regions on the Venn diagram have __no__ significantly enriched terms.     
  * Isoforms between HuFNSC01 and HuFNSC02 in neurospheres show similar terms, related to __cell signaling__, and __blood cell development__ in brain.     

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Tue Nov  4 18:01:15 2014 -->
<table border=1>
<tr> <th>  </th> <th> DE_genes </th> <th> DE_exons </th> <th> with_expressed_genes </th> <th> isoform_exons </th> <th> exclude_DE_genes </th> <th> isoform_genes </th>  </tr>
  <tr> <td> brain01_brain02_summary </td> <td align="center"> 642 </td> <td align="center"> 32138 </td> <td align="center"> 16302 </td> <td align="center"> 8980 </td> <td align="center"> 8542 </td> <td align="center"> 2902 </td> </tr>
  <tr> <td> cortex01_cortex02_summary </td> <td align="center"> 596 </td> <td align="center"> 26983 </td> <td align="center"> 15554 </td> <td align="center"> 7618 </td> <td align="center"> 7445 </td> <td align="center"> 2454 </td> </tr>
  <tr> <td> GE01_GE02_summary </td> <td align="center"> 173 </td> <td align="center"> 23810 </td> <td align="center"> 12862 </td> <td align="center"> 7402 </td> <td align="center"> 7351 </td> <td align="center"> 2495 </td> </tr>
  <tr> <td> cortex03_cortex04_summary </td> <td align="center"> 642 </td> <td align="center"> 26826 </td> <td align="center"> 12185 </td> <td align="center"> 5818 </td> <td align="center"> 5479 </td> <td align="center"> 1994 </td> </tr>
  <tr> <td> GE03_GE04_summary </td> <td align="center"> 545 </td> <td align="center"> 24752 </td> <td align="center"> 12223 </td> <td align="center"> 4582 </td> <td align="center"> 4422 </td> <td align="center"> 1454 </td> </tr>
   </table>
![](./Summary_files/figure-html/Isoform_MZ-1.png) ![](./Summary_files/figure-html/Isoform_MZ-2.png) ![](./Summary_files/figure-html/Isoform_MZ-3.png) 

#### Intron retention in MZ twins - _TBC_

  * _PENDING_

## Methods
### DMR identification 
#### WGBS

  * Identify DM CpGs     
    + methyl_diff one-sided p-value $\le$ 0.005  
    + delta fractional methylation $\ge$ 0.5  
    + fractional methylation of one sample $\ge$ 0.75   
  * Collapse DM CpGs into DMRs     
    + adjacent DM CpGs have the same DM status;    
    + distance between adjacent CpGs (size) $\le$ 300bp;   
    + No. of CpGs within each DMR $\ge$ 3.   

#### MeDIP 

  * DM CpG identification: 
    + delta fractional methylation $\ge$ 0.6  
    + fractional methylation of one sample $\ge$ 0.75   
  * Collapse DM CpGs into DMRs:   
    + adjacent CpGs have the same DM status;    
    + distance between adjacent CpGs $\le$ 300bp;   
    + No. of CpGs within each DMR $\ge$ 4.   

### Differential gene expression with DEfine

* FDR = 0.01    
* Minimum sum of RPKM (rmin) = 0.005    
* Minimum sum of coverage (Nmin) = 25    

### Isoform identification and junction validation  

  * DEfine on exons: FDR = 0.01     
  * Exon expressed in one sample ($\ge$ 10% gene RPKM) and not expressed in the other ($\le$ 1% gene RPKM)   
  * Gene is not DE: DEfine FDR = 0.01
  * Gene is expressed in both samples: gene RPKM > 0.01         
  * Validation: For each isoform exon in the previous pairwise comparison
    + Find junctions associated with this exon with enough coverage, i.e. sum of junction coverage of two samples $\ge$ 1
    + Identify junctions that RPKM change in the same direction as the exon
    + Junction RPKM > 0.1 in one sample and < 0.1 in the other      

## Discussions



