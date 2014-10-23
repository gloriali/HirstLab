Fetal Brain RNA-seq - DE genes 
========================================================

Gloria Li         
Updated: Wed Oct 22 21:55:11 2014 





## Differentially expressed genes
### DEfine

  * FDR = 0.01    
  * Minimum sum of RPKM (rmin) = 0.005    
  * Minimum sum of coverage (Nmin) = 25    
  
### Between cortex and GE neurospheres

  * On average, there are __860__ genes differentially expressed between cortex and GE, among them, __454__ are upregulated in cortex, and __406__ are downregulated.    
  * __382__ Cortex up-regulated genes, and __456__ GE up-regulated genes are shared by at least two individuals.    
  * DAVID enrichment analysis show significant enrichment in __neuronal development__ and __cell migration__ terms, GE up-regulated genes are enriched in __EGF-related__ protein domains as well. 

<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Wed Oct 22 21:55:12 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> UP </TH> <TH> DN </TH> <TH> DE </TH>  </TR>
  <TR> <TD> HuFNSC01 </TD> <TD align="center"> 403 </TD> <TD align="center"> 508 </TD> <TD align="center"> 911 </TD> </TR>
  <TR> <TD> HuFNSC02 </TD> <TD align="center"> 588 </TD> <TD align="center"> 640 </TD> <TD align="center"> 1228 </TD> </TR>
  <TR> <TD> HuFNSC03 </TD> <TD align="center"> 447 </TD> <TD align="center"> 227 </TD> <TD align="center"> 674 </TD> </TR>
  <TR> <TD> HuFNSC04 </TD> <TD align="center"> 378 </TD> <TD align="center"> 249 </TD> <TD align="center"> 627 </TD> </TR>
   </TABLE>
![plot of chunk cortex_GE](./DEgene_files/figure-html/cortex_GE1.png) ![plot of chunk cortex_GE](./DEgene_files/figure-html/cortex_GE2.png) 

![plot of chunk cortex_GE_enrich](./DEgene_files/figure-html/cortex_GE_enrich1.png) ![plot of chunk cortex_GE_enrich](./DEgene_files/figure-html/cortex_GE_enrich2.png) 

<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Wed Oct 22 21:55:13 2014 -->
<TABLE border=1>
<TR> <TH> DE </TH> <TH> name </TH> <TH> description </TH>  </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> SLC1A6 </TD> <TD align="right"> solute_carrier_family_1_(high_affinity_aspartate/glutamate_transporter),_member_6 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> PCDH20 </TD> <TD align="right"> protocadherin_20 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> C9orf64 </TD> <TD align="right"> chromosome_9_open_reading_frame_64 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> GFAP </TD> <TD align="right"> glial_fibrillary_acidic_protein </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> ZIC5 </TD> <TD align="right"> Zic_family_member_5 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> NEFM </TD> <TD align="right"> neurofilament,_medium_polypeptide </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> QPCT </TD> <TD align="right"> glutaminyl-peptide_cyclotransferase </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> UNC5C </TD> <TD align="right"> unc-5_homolog_C_(C._elegans) </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> ZIC2 </TD> <TD align="right"> Zic_family_member_2 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> XKR8 </TD> <TD align="right"> XK,_Kell_blood_group_complex_subunit-related_family,_member_8 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> PYGL </TD> <TD align="right"> phosphorylase,_glycogen,_liver </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> SERINC2 </TD> <TD align="right"> serine_incorporator_2 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> KIAA1239 </TD> <TD align="right"> KIAA1239 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> B3GALT1 </TD> <TD align="right"> UDP-Gal:betaGlcNAc_beta_1,3-galactosyltransferase,_polypeptide_1 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> CCDC48 </TD> <TD align="right"> coiled-coil_domain_containing_48 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> NOS1AP </TD> <TD align="right"> nitric_oxide_synthase_1_(neuronal)_adaptor_protein </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> ADAM19 </TD> <TD align="right"> ADAM_metallopeptidase_domain_19 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> H1F0 </TD> <TD align="right"> H1_histone_family,_member_0 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> ZIC3 </TD> <TD align="right"> Zic_family_member_3 </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> NT5E </TD> <TD align="right"> 5'-nucleotidase,_ecto_(CD73) </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> SYNM </TD> <TD align="right"> synemin,_intermediate_filament_protein </TD> </TR>
  <TR> <TD align="center"> UP </TD> <TD align="right"> C1orf226 </TD> <TD align="right"> chromosome_1_open_reading_frame_226 </TD> </TR>
   </TABLE>
<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Wed Oct 22 21:55:13 2014 -->
<TABLE border=1>
<TR> <TH> DE </TH> <TH> name </TH> <TH> description </TH>  </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> DPPA4 </TD> <TD align="right"> developmental_pluripotency_associated_4 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> VAX1 </TD> <TD align="right"> ventral_anterior_homeobox_1 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> SHISA6 </TD> <TD align="right"> shisa_homolog_6_(Xenopus_laevis) </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> NRXN3 </TD> <TD align="right"> neurexin_3 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> FGFR2 </TD> <TD align="right"> fibroblast_growth_factor_receptor_2 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> ST8SIA2 </TD> <TD align="right"> ST8_alpha-N-acetyl-neuraminide_alpha-2,8-sialyltransferase_2 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> SLC32A1 </TD> <TD align="right"> solute_carrier_family_32_(GABA_vesicular_transporter),_member_1 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> FEZF1 </TD> <TD align="right"> FEZ_family_zinc_finger_1 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> AJAP1 </TD> <TD align="right"> adherens_junctions_associated_protein_1 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> ZNF703 </TD> <TD align="right"> zinc_finger_protein_703 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> VSIG10L </TD> <TD align="right"> V-set_and_immunoglobulin_domain_containing_10_like </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> SIX3 </TD> <TD align="right"> SIX_homeobox_3 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> SORL1 </TD> <TD align="right"> sortilin-related_receptor,_L(DLR_class)_A_repeats_containing </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> LMO1 </TD> <TD align="right"> LIM_domain_only_1_(rhombotin_1) </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> ODZ1 </TD> <TD align="right"> odz,_odd_Oz/ten-m_homolog_1_(Drosophila) </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> EPHA3 </TD> <TD align="right"> EPH_receptor_A3 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> TIMP3 </TD> <TD align="right"> TIMP_metallopeptidase_inhibitor_3 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> CAMK2N1 </TD> <TD align="right"> calcium/calmodulin-dependent_protein_kinase_II_inhibitor_1 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> EEF1A2 </TD> <TD align="right"> eukaryotic_translation_elongation_factor_1_alpha_2 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> OTX2 </TD> <TD align="right"> orthodenticle_homeobox_2 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> CHL1 </TD> <TD align="right"> cell_adhesion_molecule_with_homology_to_L1CAM_(close_homolog_of_L1) </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> PCSK1N </TD> <TD align="right"> proprotein_convertase_subtilisin/kexin_type_1_inhibitor </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> LMO2 </TD> <TD align="right"> LIM_domain_only_2_(rhombotin-like_1) </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> ASTN1 </TD> <TD align="right"> astrotactin_1 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> BASP1 </TD> <TD align="right"> brain_abundant,_membrane_attached_signal_protein_1 </TD> </TR>
  <TR> <TD align="center"> DN </TD> <TD align="right"> ADCYAP1R1 </TD> <TD align="right"> adenylate_cyclase_activating_polypeptide_1_(pituitary)_receptor_type_I </TD> </TR>
   </TABLE>

### Between MZ twins - HuFNSC01 vs HuFNSC02

  * On average, there are __470__ DE genes across three cells types.   
  * Majority of DE genes are cell type specific, only __98__ are shared between any two cell types.   
  * DAVID enrichment analysis between MZ twins in brain and cortex show similar GO term in __brain development__, but there is no significantly enriched terms in GE.    

<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Wed Oct 22 21:55:13 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> UP </TH> <TH> DN </TH> <TH> DE </TH>  </TR>
  <TR> <TD> brain01_brain02 </TD> <TD align="center"> 461 </TD> <TD align="center"> 181 </TD> <TD align="center"> 642 </TD> </TR>
  <TR> <TD> cortex01_cortex02 </TD> <TD align="center"> 248 </TD> <TD align="center"> 348 </TD> <TD align="center"> 596 </TD> </TR>
  <TR> <TD> GE01_GE02 </TD> <TD align="center">  99 </TD> <TD align="center">  74 </TD> <TD align="center"> 173 </TD> </TR>
  <TR> <TD> cortex03_GE03 </TD> <TD align="center"> 459 </TD> <TD align="center"> 183 </TD> <TD align="center"> 642 </TD> </TR>
  <TR> <TD> cortex04_GE04 </TD> <TD align="center"> 313 </TD> <TD align="center"> 232 </TD> <TD align="center"> 545 </TD> </TR>
   </TABLE>
![plot of chunk individual](./DEgene_files/figure-html/individual.png) 

![plot of chunk individual_enrich](./DEgene_files/figure-html/individual_enrich1.png) ![plot of chunk individual_enrich](./DEgene_files/figure-html/individual_enrich2.png) 

<!-- html table generated in R 3.0.2 by xtable 1.7-3 package -->
<!-- Wed Oct 22 21:55:14 2014 -->
<TABLE border=1>
<TR> <TH> name </TH> <TH> description </TH> <TH> Brain </TH> <TH> Cortex </TH> <TH> GE </TH>  </TR>
  <TR> <TD align="right"> LMO1 </TD> <TD align="right"> LIM_domain_only_1_(rhombotin_1) </TD> <TD align="center"> UP </TD> <TD align="center"> DN </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> CXCR7 </TD> <TD align="right"> chemokine_(C-X-C_motif)_receptor_7 </TD> <TD align="center"> UP </TD> <TD align="center"> UP </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> SPP1 </TD> <TD align="right"> secreted_phosphoprotein_1 </TD> <TD align="center"> UP </TD> <TD align="center"> DN </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> NPY </TD> <TD align="right"> neuropeptide_Y </TD> <TD align="center"> UP </TD> <TD align="center"> DN </TD> <TD align="center"> DN </TD> </TR>
  <TR> <TD align="right"> COL2A1 </TD> <TD align="right"> collagen,_type_II,_alpha_1 </TD> <TD align="center"> DN </TD> <TD align="center"> UP </TD> <TD align="center"> UP </TD> </TR>
  <TR> <TD align="right"> BCL6 </TD> <TD align="right"> B-cell_CLL/lymphoma_6 </TD> <TD align="center"> DN </TD> <TD align="center"> DN </TD> <TD align="center"> DN </TD> </TR>
   </TABLE>




