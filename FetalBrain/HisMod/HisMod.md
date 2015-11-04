# Fetal Brain - Histone modifications
Gloria Li  
January 27, 2014  

Update Wed Nov  4 13:24:03 2015



## Sanity check   

* No. of peaks, No. of enriched bases, and average peak length seem reasonable except for the unusual high No. of peaks in GE HuFNSC04 input library.              

![](HisMod_files/figure-html/summary-1.png) 

## Correlation with protein-coding gene RPKM 

* Overlapping H3K4me3 and H3K27me3 FindER peaks with protein-coding gene promoters (TSS +/- 1500bp), and overlapping H3K36me3 with genebody.    
* Overall correlations are as expected, with H3K4m3 and H3K36me3 marked genes showing higher RPKM, and H3K27me3 and bivalent promoter marked showing lower RPKM.   

![](HisMod_files/figure-html/RPKM-1.png) 

## Differentially marked genes

* Calculate H3K4me3 and H3K27me3 signal from wig at promoters (TSS +/- 1500bp) and rank all genes (pc + nc).    
* Differential marked genes are defined as rank difference >= 5000.           
* For H3K27me3 differential marked genes, there are opposite trends in pc and nc, _why?_           

### MZ twins

![](HisMod_files/figure-html/DM_MZ1-1.png) ![](HisMod_files/figure-html/DM_MZ1-2.png) ![](HisMod_files/figure-html/DM_MZ1-3.png) ![](HisMod_files/figure-html/DM_MZ1-4.png) 
![](HisMod_files/figure-html/DM_MZ2-1.png) 

```
## [1] "No enrichment for DM_H3K27me3_Brain_Subject2.pc"
```

![](HisMod_files/figure-html/DM_MZ2-2.png) 

```
## [1] "No enrichment for DM_H3K27me3_Cortex_Subject2.pc"
```
![](HisMod_files/figure-html/DM_MZ3-1.png) 
![](HisMod_files/figure-html/DM_MZ4-1.png) ![](HisMod_files/figure-html/DM_MZ4-2.png) 

* Brain01 vs Brain02: DE genes that are DM in both H3K4me3 and H3K27me3         
<table>
 <thead>
  <tr>
   <th style="text-align:right;"> name </th>
   <th style="text-align:right;"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> ARSD </td>
   <td style="text-align:right;"> arylsulfatase_D_[Source:HGNC_Symbol;Acc:717] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ADIPOR2 </td>
   <td style="text-align:right;"> adiponectin_receptor_2_[Source:HGNC_Symbol;Acc:24041] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> KPNA6 </td>
   <td style="text-align:right;"> karyopherin_alpha_6_(importin_alpha_7)_[Source:HGNC_Symbol;Acc:6399] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> CD38 </td>
   <td style="text-align:right;"> CD38_molecule_[Source:HGNC_Symbol;Acc:1667] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> DERA </td>
   <td style="text-align:right;"> deoxyribose-phosphate_aldolase_(putative)_[Source:HGNC_Symbol;Acc:24269] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> HCCS </td>
   <td style="text-align:right;"> holocytochrome_c_synthase_[Source:HGNC_Symbol;Acc:4837] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> RHBDF1 </td>
   <td style="text-align:right;"> rhomboid_5_homolog_1_(Drosophila)_[Source:HGNC_Symbol;Acc:20561] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> C8B </td>
   <td style="text-align:right;"> complement_component_8,_beta_polypeptide_[Source:HGNC_Symbol;Acc:1353] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> DEF6 </td>
   <td style="text-align:right;"> differentially_expressed_in_FDCP_6_homolog_(mouse)_[Source:HGNC_Symbol;Acc:2760] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> VIM </td>
   <td style="text-align:right;"> vimentin_[Source:HGNC_Symbol;Acc:12692] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TYROBP </td>
   <td style="text-align:right;"> TYRO_protein_tyrosine_kinase_binding_protein_[Source:HGNC_Symbol;Acc:12449] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> DYRK4 </td>
   <td style="text-align:right;"> dual-specificity_tyrosine-(Y)-phosphorylation_regulated_kinase_4_[Source:HGNC_Symbol;Acc:3095] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> BIRC3 </td>
   <td style="text-align:right;"> baculoviral_IAP_repeat_containing_3_[Source:HGNC_Symbol;Acc:591] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> IL32 </td>
   <td style="text-align:right;"> interleukin_32_[Source:HGNC_Symbol;Acc:16830] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> PHTF2 </td>
   <td style="text-align:right;"> putative_homeodomain_transcription_factor_2_[Source:HGNC_Symbol;Acc:13411] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SEC62 </td>
   <td style="text-align:right;"> SEC62_homolog_(S._cerevisiae)_[Source:HGNC_Symbol;Acc:11846] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> GPR124 </td>
   <td style="text-align:right;"> G_protein-coupled_receptor_124_[Source:HGNC_Symbol;Acc:17849] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> HHATL </td>
   <td style="text-align:right;"> hedgehog_acyltransferase-like_[Source:HGNC_Symbol;Acc:13242] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> COX15 </td>
   <td style="text-align:right;"> COX15_homolog,_cytochrome_c_oxidase_assembly_protein_(yeast)_[Source:HGNC_Symbol;Acc:2263] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SCYL3 </td>
   <td style="text-align:right;"> SCY1-like_3_(S._cerevisiae)_[Source:HGNC_Symbol;Acc:19285] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> MNAT1 </td>
   <td style="text-align:right;"> menage_a_trois_homolog_1,_cyclin_H_assembly_factor_(Xenopus_laevis)_[Source:HGNC_Symbol;Acc:7181] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SYN1 </td>
   <td style="text-align:right;"> synapsin_I_[Source:HGNC_Symbol;Acc:11494] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> C1orf201 </td>
   <td style="text-align:right;"> chromosome_1_open_reading_frame_201_[Source:HGNC_Symbol;Acc:28070] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> JHDM1D </td>
   <td style="text-align:right;"> jumonji_C_domain_containing_histone_demethylase_1_homolog_D_(S._cerevisiae)_[Source:HGNC_Symbol;Acc:22224] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> LGALS14 </td>
   <td style="text-align:right;"> lectin,_galactoside-binding,_soluble,_14_[Source:HGNC_Symbol;Acc:30054] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SERPINB1 </td>
   <td style="text-align:right;"> serpin_peptidase_inhibitor,_clade_B_(ovalbumin),_member_1_[Source:HGNC_Symbol;Acc:3311] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> XYLT2 </td>
   <td style="text-align:right;"> xylosyltransferase_II_[Source:HGNC_Symbol;Acc:15517] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SLC38A5 </td>
   <td style="text-align:right;"> solute_carrier_family_38,_member_5_[Source:HGNC_Symbol;Acc:18070] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TEAD3 </td>
   <td style="text-align:right;"> TEA_domain_family_member_3_[Source:HGNC_Symbol;Acc:11716] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SNAPC1 </td>
   <td style="text-align:right;"> small_nuclear_RNA_activating_complex,_polypeptide_1,_43kDa_[Source:HGNC_Symbol;Acc:11134] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> BRCA1 </td>
   <td style="text-align:right;"> breast_cancer_1,_early_onset_[Source:HGNC_Symbol;Acc:1100] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> MTMR11 </td>
   <td style="text-align:right;"> myotubularin_related_protein_11_[Source:HGNC_Symbol;Acc:24307] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> IDS </td>
   <td style="text-align:right;"> iduronate_2-sulfatase_[Source:HGNC_Symbol;Acc:5389] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> UBR7 </td>
   <td style="text-align:right;"> ubiquitin_protein_ligase_E3_component_n-recognin_7_(putative)_[Source:HGNC_Symbol;Acc:20344] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> PHF20 </td>
   <td style="text-align:right;"> PHD_finger_protein_20_[Source:HGNC_Symbol;Acc:16098] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> C12orf4 </td>
   <td style="text-align:right;"> chromosome_12_open_reading_frame_4_[Source:HGNC_Symbol;Acc:1184] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> VCL </td>
   <td style="text-align:right;"> vinculin_[Source:HGNC_Symbol;Acc:12665] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ARNTL2 </td>
   <td style="text-align:right;"> aryl_hydrocarbon_receptor_nuclear_translocator-like_2_[Source:HGNC_Symbol;Acc:18984] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> RABEP1 </td>
   <td style="text-align:right;"> rabaptin,_RAB_GTPase_binding_effector_protein_1_[Source:HGNC_Symbol;Acc:17677] </td>
  </tr>
</tbody>
</table>

* Brain01 vs Brain02: DE genes that are DM in both DMR and H3K27me3         

<table>
 <thead>
  <tr>
   <th style="text-align:right;"> name </th>
   <th style="text-align:right;"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> AGK </td>
   <td style="text-align:right;"> acylglycerol_kinase_[Source:HGNC_Symbol;Acc:21869] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> RAD52 </td>
   <td style="text-align:right;"> RAD52_homolog_(S._cerevisiae)_[Source:HGNC_Symbol;Acc:9824] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> PGM3 </td>
   <td style="text-align:right;"> phosphoglucomutase_3_[Source:HGNC_Symbol;Acc:8907] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> CCL18 </td>
   <td style="text-align:right;"> chemokine_(C-C_motif)_ligand_18_(pulmonary_and_activation-regulated)_[Source:HGNC_Symbol;Acc:10616] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TBXA2R </td>
   <td style="text-align:right;"> thromboxane_A2_receptor_[Source:HGNC_Symbol;Acc:11608] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> AGPS </td>
   <td style="text-align:right;"> alkylglycerone_phosphate_synthase_[Source:HGNC_Symbol;Acc:327] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TG </td>
   <td style="text-align:right;"> thyroglobulin_[Source:HGNC_Symbol;Acc:11764] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> XK </td>
   <td style="text-align:right;"> X-linked_Kx_blood_group_(McLeod_syndrome)_[Source:HGNC_Symbol;Acc:12811] </td>
  </tr>
</tbody>
</table>

* Cortex01 vs Cortex02: DE genes that are DM in both DMR and H3K27me3         

<table>
 <thead>
  <tr>
   <th style="text-align:right;"> name </th>
   <th style="text-align:right;"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> JHDM1D </td>
   <td style="text-align:right;"> jumonji_C_domain_containing_histone_demethylase_1_homolog_D_(S._cerevisiae)_[Source:HGNC_Symbol;Acc:22224] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SPATA20 </td>
   <td style="text-align:right;"> spermatogenesis_associated_20_[Source:HGNC_Symbol;Acc:26125] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TNFRSF12A </td>
   <td style="text-align:right;"> tumor_necrosis_factor_receptor_superfamily,_member_12A_[Source:HGNC_Symbol;Acc:18152] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ETV1 </td>
   <td style="text-align:right;"> ets_variant_1_[Source:HGNC_Symbol;Acc:3490] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> CCL3 </td>
   <td style="text-align:right;"> chemokine_(C-C_motif)_ligand_3_[Source:HGNC_Symbol;Acc:10627] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ITGA2B </td>
   <td style="text-align:right;"> integrin,_alpha_2b_(platelet_glycoprotein_IIb_of_IIb/IIIa_complex,_antigen_CD41)_[Source:HGNC_Symbol;Acc:6138] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> CELSR3 </td>
   <td style="text-align:right;"> cadherin,_EGF_LAG_seven-pass_G-type_receptor_3_(flamingo_homolog,_Drosophila)_[Source:HGNC_Symbol;Acc:3230] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> IKZF2 </td>
   <td style="text-align:right;"> IKAROS_family_zinc_finger_2_(Helios)_[Source:HGNC_Symbol;Acc:13177] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> DERA </td>
   <td style="text-align:right;"> deoxyribose-phosphate_aldolase_(putative)_[Source:HGNC_Symbol;Acc:24269] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ZNF806 </td>
   <td style="text-align:right;"> zinc_finger_protein_806_[Source:HGNC_Symbol;Acc:33228] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> BID </td>
   <td style="text-align:right;"> BH3_interacting_domain_death_agonist_[Source:HGNC_Symbol;Acc:1050] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> RAB27B </td>
   <td style="text-align:right;"> RAB27B,_member_RAS_oncogene_family_[Source:HGNC_Symbol;Acc:9767] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> PTBP1 </td>
   <td style="text-align:right;"> polypyrimidine_tract_binding_protein_1_[Source:HGNC_Symbol;Acc:9583] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> PLEKHO1 </td>
   <td style="text-align:right;"> pleckstrin_homology_domain_containing,_family_O_member_1_[Source:HGNC_Symbol;Acc:24310] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ANK1 </td>
   <td style="text-align:right;"> ankyrin_1,_erythrocytic_[Source:HGNC_Symbol;Acc:492] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> RTN4R </td>
   <td style="text-align:right;"> reticulon_4_receptor_[Source:HGNC_Symbol;Acc:18601] </td>
  </tr>
</tbody>
</table>

* GE01 vs GE02: DE genes that are DM in both DMR and H3K27me3         

<table>
 <thead>
  <tr>
   <th style="text-align:right;"> name </th>
   <th style="text-align:right;"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> ADIPOR2 </td>
   <td style="text-align:right;"> adiponectin_receptor_2_[Source:HGNC_Symbol;Acc:24041] </td>
  </tr>
</tbody>
</table>

### NPCs

![](HisMod_files/figure-html/DM_NPC1-1.png) ![](HisMod_files/figure-html/DM_NPC1-2.png) ![](HisMod_files/figure-html/DM_NPC1-3.png) ![](HisMod_files/figure-html/DM_NPC1-4.png) 
![](HisMod_files/figure-html/DM_NPC2-1.png) 
![](HisMod_files/figure-html/DM_NPC3-1.png) ![](HisMod_files/figure-html/DM_NPC3-2.png) ![](HisMod_files/figure-html/DM_NPC3-3.png) ![](HisMod_files/figure-html/DM_NPC3-4.png) ![](HisMod_files/figure-html/DM_NPC3-5.png) ![](HisMod_files/figure-html/DM_NPC3-6.png) ![](HisMod_files/figure-html/DM_NPC3-7.png) 

* Cortex01 vs GE01: DE genes that are DM in both DMR and H3K7me3        

<table>
 <thead>
  <tr>
   <th style="text-align:right;"> name </th>
   <th style="text-align:right;"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> TCIRG1 </td>
   <td style="text-align:right;"> T-cell,_immune_regulator_1,_ATPase,_H+_transporting,_lysosomal_V0_subunit_A3_[Source:HGNC_Symbol;Acc:11647] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> AGT </td>
   <td style="text-align:right;"> angiotensinogen_(serpin_peptidase_inhibitor,_clade_A,_member_8)_[Source:HGNC_Symbol;Acc:333] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> FES </td>
   <td style="text-align:right;"> feline_sarcoma_oncogene_[Source:HGNC_Symbol;Acc:3657] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ECHDC3 </td>
   <td style="text-align:right;"> enoyl_CoA_hydratase_domain_containing_3_[Source:HGNC_Symbol;Acc:23489] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> DOCK9 </td>
   <td style="text-align:right;"> dedicator_of_cytokinesis_9_[Source:HGNC_Symbol;Acc:14132] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ELF4 </td>
   <td style="text-align:right;"> E74-like_factor_4_(ets_domain_transcription_factor)_[Source:HGNC_Symbol;Acc:3319] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> CTSZ </td>
   <td style="text-align:right;"> cathepsin_Z_[Source:HGNC_Symbol;Acc:2547] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> NEK9 </td>
   <td style="text-align:right;"> NIMA_(never_in_mitosis_gene_a)-_related_kinase_9_[Source:HGNC_Symbol;Acc:18591] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> Uncharacterized_protein_[Source:UniProtKB/TrEMBL;Acc:A8MV45] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> KIAA1244 </td>
   <td style="text-align:right;"> KIAA1244_[Source:HGNC_Symbol;Acc:21213] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TMEM132B </td>
   <td style="text-align:right;"> transmembrane_protein_132B_[Source:HGNC_Symbol;Acc:29397] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ITPR3 </td>
   <td style="text-align:right;"> inositol_1,4,5-trisphosphate_receptor,_type_3_[Source:HGNC_Symbol;Acc:6182] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> THSD7A </td>
   <td style="text-align:right;"> thrombospondin,_type_I,_domain_containing_7A_[Source:HGNC_Symbol;Acc:22207] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TSPYL5 </td>
   <td style="text-align:right;"> TSPY-like_5_[Source:HGNC_Symbol;Acc:29367] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> CYP27A1 </td>
   <td style="text-align:right;"> cytochrome_P450,_family_27,_subfamily_A,_polypeptide_1_[Source:HGNC_Symbol;Acc:2605] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> PYGL </td>
   <td style="text-align:right;"> phosphorylase,_glycogen,_liver_[Source:HGNC_Symbol;Acc:9725] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> LIPG </td>
   <td style="text-align:right;"> lipase,_endothelial_[Source:HGNC_Symbol;Acc:6623] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> FUT9 </td>
   <td style="text-align:right;"> fucosyltransferase_9_(alpha_(1,3)_fucosyltransferase)_[Source:HGNC_Symbol;Acc:4020] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> OTX1 </td>
   <td style="text-align:right;"> orthodenticle_homeobox_1_[Source:HGNC_Symbol;Acc:8521] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TFCP2 </td>
   <td style="text-align:right;"> transcription_factor_CP2_[Source:HGNC_Symbol;Acc:11748] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> GSX2 </td>
   <td style="text-align:right;"> GS_homeobox_2_[Source:HGNC_Symbol;Acc:24959] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> COL23A1 </td>
   <td style="text-align:right;"> collagen,_type_XXIII,_alpha_1_[Source:HGNC_Symbol;Acc:22990] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> C13orf15 </td>
   <td style="text-align:right;"> chromosome_13_open_reading_frame_15_[Source:HGNC_Symbol;Acc:20369] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> MDK </td>
   <td style="text-align:right;"> midkine_(neurite_growth-promoting_factor_2)_[Source:HGNC_Symbol;Acc:6972] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TSKU </td>
   <td style="text-align:right;"> tsukushi_small_leucine_rich_proteoglycan_homolog_(Xenopus_laevis)_[Source:HGNC_Symbol;Acc:28850] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> LRRC10B </td>
   <td style="text-align:right;"> leucine_rich_repeat_containing_10B_[Source:HGNC_Symbol;Acc:37215] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ID2 </td>
   <td style="text-align:right;"> inhibitor_of_DNA_binding_2,_dominant_negative_helix-loop-helix_protein_[Source:HGNC_Symbol;Acc:5361] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> APCDD1 </td>
   <td style="text-align:right;"> adenomatosis_polyposis_coli_down-regulated_1_[Source:HGNC_Symbol;Acc:15718] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> COL9A3 </td>
   <td style="text-align:right;"> collagen,_type_IX,_alpha_3_[Source:HGNC_Symbol;Acc:2219] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> APOE </td>
   <td style="text-align:right;"> apolipoprotein_E_[Source:HGNC_Symbol;Acc:613] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> C10orf90 </td>
   <td style="text-align:right;"> chromosome_10_open_reading_frame_90_[Source:HGNC_Symbol;Acc:26563] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> PRDM1 </td>
   <td style="text-align:right;"> PR_domain_containing_1,_with_ZNF_domain_[Source:HGNC_Symbol;Acc:9346] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> LMO3 </td>
   <td style="text-align:right;"> LIM_domain_only_3_(rhombotin-like_2)_[Source:HGNC_Symbol;Acc:6643] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> VAX1 </td>
   <td style="text-align:right;"> ventral_anterior_homeobox_1_[Source:HGNC_Symbol;Acc:12660] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> FZD5 </td>
   <td style="text-align:right;"> frizzled_family_receptor_5_[Source:HGNC_Symbol;Acc:4043] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> CXCR7 </td>
   <td style="text-align:right;"> chemokine_(C-X-C_motif)_receptor_7_[Source:HGNC_Symbol;Acc:23692] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> FEZF1 </td>
   <td style="text-align:right;"> FEZ_family_zinc_finger_1_[Source:HGNC_Symbol;Acc:22788] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> CYTL1 </td>
   <td style="text-align:right;"> cytokine-like_1_[Source:HGNC_Symbol;Acc:24435] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> KCTD12 </td>
   <td style="text-align:right;"> potassium_channel_tetramerisation_domain_containing_12_[Source:HGNC_Symbol;Acc:14678] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> CHST6 </td>
   <td style="text-align:right;"> carbohydrate_(N-acetylglucosamine_6-O)_sulfotransferase_6_[Source:HGNC_Symbol;Acc:6938] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> BRSK2 </td>
   <td style="text-align:right;"> BR_serine/threonine_kinase_2_[Source:HGNC_Symbol;Acc:11405] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> OTX2 </td>
   <td style="text-align:right;"> orthodenticle_homeobox_2_[Source:HGNC_Symbol;Acc:8522] </td>
  </tr>
</tbody>
</table>

* Cortex02 vs GE02: DE genes that are DM in DMR, H3K4me3, and H3K7me3        

<table>
 <thead>
  <tr>
   <th style="text-align:right;"> name </th>
   <th style="text-align:right;"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> PAX6 </td>
   <td style="text-align:right;"> paired_box_6_[Source:HGNC_Symbol;Acc:8620] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> EMX2 </td>
   <td style="text-align:right;"> empty_spiracles_homeobox_2_[Source:HGNC_Symbol;Acc:3341] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TMEM132B </td>
   <td style="text-align:right;"> transmembrane_protein_132B_[Source:HGNC_Symbol;Acc:29397] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TFCP2 </td>
   <td style="text-align:right;"> transcription_factor_CP2_[Source:HGNC_Symbol;Acc:11748] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> HOPX </td>
   <td style="text-align:right;"> HOP_homeobox_[Source:HGNC_Symbol;Acc:24961] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> C10orf90 </td>
   <td style="text-align:right;"> chromosome_10_open_reading_frame_90_[Source:HGNC_Symbol;Acc:26563] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SFTA3 </td>
   <td style="text-align:right;"> surfactant_associated_3_[Source:HGNC_Symbol;Acc:18387] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> VAX1 </td>
   <td style="text-align:right;"> ventral_anterior_homeobox_1_[Source:HGNC_Symbol;Acc:12660] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> TBX2 </td>
   <td style="text-align:right;"> T-box_2_[Source:HGNC_Symbol;Acc:11597] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SIX3 </td>
   <td style="text-align:right;"> SIX_homeobox_3_[Source:HGNC_Symbol;Acc:10889] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SLIT2 </td>
   <td style="text-align:right;"> slit_homolog_2_(Drosophila)_[Source:HGNC_Symbol;Acc:11086] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> NKX2-1 </td>
   <td style="text-align:right;"> NK2_homeobox_1_[Source:HGNC_Symbol;Acc:11825] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> DLX1 </td>
   <td style="text-align:right;"> distal-less_homeobox_1_[Source:HGNC_Symbol;Acc:2914] </td>
  </tr>
</tbody>
</table>

### GW

![](HisMod_files/figure-html/DM_GW1-1.png) ![](HisMod_files/figure-html/DM_GW1-2.png) 
![](HisMod_files/figure-html/DM_GW2-1.png) ![](HisMod_files/figure-html/DM_GW2-2.png) ![](HisMod_files/figure-html/DM_GW2-3.png) ![](HisMod_files/figure-html/DM_GW2-4.png) ![](HisMod_files/figure-html/DM_GW2-5.png) 
![](HisMod_files/figure-html/DM_GW3-1.png) 
![](HisMod_files/figure-html/DM_GW4-1.png) 

* GE02 vs GE04: DE genes that are DM in DMR, H3K4me3, and H3K27me3      

<table>
 <thead>
  <tr>
   <th style="text-align:right;"> name </th>
   <th style="text-align:right;"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> HIST4H4 </td>
   <td style="text-align:right;"> histone_cluster_4,_H4_[Source:HGNC_Symbol;Acc:20510] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> OCIAD2 </td>
   <td style="text-align:right;"> OCIA_domain_containing_2_[Source:HGNC_Symbol;Acc:28685] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> C6orf201 </td>
   <td style="text-align:right;"> chromosome_6_open_reading_frame_201_[Source:HGNC_Symbol;Acc:21620] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> RHCG </td>
   <td style="text-align:right;"> Rh_family,_C_glycoprotein_[Source:HGNC_Symbol;Acc:18140] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ISLR2 </td>
   <td style="text-align:right;"> immunoglobulin_superfamily_containing_leucine-rich_repeat_2_[Source:HGNC_Symbol;Acc:29286] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> EMILIN1 </td>
   <td style="text-align:right;"> elastin_microfibril_interfacer_1_[Source:HGNC_Symbol;Acc:19880] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> EPHA1 </td>
   <td style="text-align:right;"> EPH_receptor_A1_[Source:HGNC_Symbol;Acc:3385] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> CCDC8 </td>
   <td style="text-align:right;"> coiled-coil_domain_containing_8_[Source:HGNC_Symbol;Acc:25367] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> VAX1 </td>
   <td style="text-align:right;"> ventral_anterior_homeobox_1_[Source:HGNC_Symbol;Acc:12660] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SIX3 </td>
   <td style="text-align:right;"> SIX_homeobox_3_[Source:HGNC_Symbol;Acc:10889] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ARAP3 </td>
   <td style="text-align:right;"> ArfGAP_with_RhoGAP_domain,_ankyrin_repeat_and_PH_domain_3_[Source:HGNC_Symbol;Acc:24097] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> NKX2-1 </td>
   <td style="text-align:right;"> NK2_homeobox_1_[Source:HGNC_Symbol;Acc:11825] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> OTX2 </td>
   <td style="text-align:right;"> orthodenticle_homeobox_2_[Source:HGNC_Symbol;Acc:8522] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> OLIG2 </td>
   <td style="text-align:right;"> oligodendrocyte_lineage_transcription_factor_2_[Source:HGNC_Symbol;Acc:9398] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> EDNRB </td>
   <td style="text-align:right;"> endothelin_receptor_type_B_[Source:HGNC_Symbol;Acc:3180] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> ZIC5 </td>
   <td style="text-align:right;"> Zic_family_member_5_[Source:HGNC_Symbol;Acc:20322] </td>
  </tr>
  <tr>
   <td style="text-align:right;">  </td>
   <td style="text-align:right;">  </td>
  </tr>
  <tr>
   <td style="text-align:right;"> SP8 </td>
   <td style="text-align:right;"> Sp8_transcription_factor_[Source:HGNC_Symbol;Acc:19196] </td>
  </tr>
  <tr>
   <td style="text-align:right;"> FZD7 </td>
   <td style="text-align:right;"> frizzled_family_receptor_7_[Source:HGNC_Symbol;Acc:4045] </td>
  </tr>
</tbody>
</table>

## Core enhancers 

* Overlapping all NPC enhancers (Cortex01, Cortex02, GE01, GE02, GE04), in total 79033 regions, average length 2904.36bp.     

### GWAS in core enhancers  

* 1460 enhancers overlap with GWAS, involved in 546 traits.     
* Several traits are related to brain development and function: Cortical structure, Glaucoma (exfoliation), Glioblastoma, Neuranatomic and neurocognitive phenotypes, Neuroblastoma (high-risk), Odorant perception (isobutyraldehyde), Schizophrenia (cytomegalovirus infection interaction), Alzheimer's disease biomarkers.        

### Homer TFBSs

* There are 13 TFs significantly (Benjamini q value < 0.01) enriched in core enhancers and present in > 20% of the core enhancers.    
* Ptf1a: plays an important role in cerebellar and pancreatic development.       
* Isl1: central to the development of pancreatic cell lineages and may also be required for motor neuron generation.        
* Lhx3: involved in the development of interneurons and motor neurons.    
* NF1: negative regulator of the ras signal transduction pathway, associated with neurofibromatosis type 1.     
* Sox3: function as a switch in neuronal development. Keeps neural cells undifferentiated by counteracting the activity of proneural proteins and suppresses neuronal differentiation.       
* Olig2: oligodendrocyte specific marker.      

![](HisMod_files/figure-html/core_enhancer_homer-1.png) 

### Overlapping with WGBS UMRs

* UMRs between neurospheres (cortex vs GE) are enriched in enhancers (H3K4me1 enriched regions).      
* Between gestational weeks, GW13 UMRs are enriched in enhancers, but GW17 UMRs are not.      
* UMRs overlaped with enhancers are highly enriched for brain development terms. For comparing between neurospheres, GE enhancer UMRs in HuFNSC04 have no significant enrichment. And for comparing between gestational weeks, GW13 enhancer UMRs in Cortex also have no enriched terms.       

![](HisMod_files/figure-html/core_enhancer_UMR-1.png) 
![](HisMod_files/figure-html/core_enhancer_UMR_GREAT-1.png) ![](HisMod_files/figure-html/core_enhancer_UMR_GREAT-2.png) ![](HisMod_files/figure-html/core_enhancer_UMR_GREAT-3.png) ![](HisMod_files/figure-html/core_enhancer_UMR_GREAT-4.png) ![](HisMod_files/figure-html/core_enhancer_UMR_GREAT-5.png) ![](HisMod_files/figure-html/core_enhancer_UMR_GREAT-6.png) 

## Unique enhancers 

* Overall, there are about 20-40% enhancers that are unique within each pairwise comparisons.      
* Between MZ twins, HuFNSC01 has more unique enhancers than HuFNSC02 in all three cell types. The asymmetry between two progenitor cell types and between GWs is not supported by both samples.   
* The intersect of two samples between both progenitors and GWs are statistically significant.    

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Comparison </th>
   <th style="text-align:left;"> Samples </th>
   <th style="text-align:center;"> Sample1 </th>
   <th style="text-align:center;"> Sample2 </th>
   <th style="text-align:center;"> Sample1_unique </th>
   <th style="text-align:center;"> Sample2_unique </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MZ </td>
   <td style="text-align:left;"> Brain01_Brain02 </td>
   <td style="text-align:center;"> 69300 </td>
   <td style="text-align:center;"> 58831 </td>
   <td style="text-align:center;"> 30038 </td>
   <td style="text-align:center;"> 18240 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MZ </td>
   <td style="text-align:left;"> Cortex01_Cortex02 </td>
   <td style="text-align:center;"> 125995 </td>
   <td style="text-align:center;"> 105033 </td>
   <td style="text-align:center;"> 35489 </td>
   <td style="text-align:center;"> 5466 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MZ </td>
   <td style="text-align:left;"> GE01_GE02 </td>
   <td style="text-align:center;"> 135489 </td>
   <td style="text-align:center;"> 102950 </td>
   <td style="text-align:center;"> 47186 </td>
   <td style="text-align:center;"> 3765 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Neurospheres </td>
   <td style="text-align:left;"> Cortex01_GE01 </td>
   <td style="text-align:center;"> 125995 </td>
   <td style="text-align:center;"> 135489 </td>
   <td style="text-align:center;"> 14961 </td>
   <td style="text-align:center;"> 24022 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Neurospheres </td>
   <td style="text-align:left;"> Cortex02_GE02 </td>
   <td style="text-align:center;"> 105033 </td>
   <td style="text-align:center;"> 102950 </td>
   <td style="text-align:center;"> 23099 </td>
   <td style="text-align:center;"> 18732 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GW </td>
   <td style="text-align:left;"> GE01_GE04 </td>
   <td style="text-align:center;"> 135489 </td>
   <td style="text-align:center;"> 102477 </td>
   <td style="text-align:center;"> 49603 </td>
   <td style="text-align:center;"> 6361 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GW </td>
   <td style="text-align:left;"> GE02_GE04 </td>
   <td style="text-align:center;"> 102950 </td>
   <td style="text-align:center;"> 102477 </td>
   <td style="text-align:center;"> 22826 </td>
   <td style="text-align:center;"> 22889 </td>
  </tr>
</tbody>
</table>

![](HisMod_files/figure-html/unique_enhancer-1.png) 

### Functional enrichment of unique enhancers

* For MZ twins, no GREAT enriched terms in cortex and GE. In mixed brain, HuFNSC02 unique enhancers show brain organ development terms, WNT pathway, and kidney-related terms.      
* For cortex and GE progenitors, both showed brain-related terms, more so in cortex.     
* For GWs, both GW13 and GW17 unique enhancers showed enrichment for brain organ development terms. GW17 also showed WNT pathway terms.      

![](HisMod_files/figure-html/unique_enhancer_GREAT1-1.png) 
![](HisMod_files/figure-html/unique_enhancer_GREAT2-1.png) 
![](HisMod_files/figure-html/unique_enhancer_GREAT3-1.png) 
![](HisMod_files/figure-html/unique_enhancer_GREAT4-1.png) ![](HisMod_files/figure-html/unique_enhancer_GREAT4-2.png) ![](HisMod_files/figure-html/unique_enhancer_GREAT4-3.png) 

### GWAS in unique enhancers 

* All sets of unique enhancers showed brain or brain disease related GWAS sites, such as Gliomas, Alzheimer's disease, Autism, Schizophrenia or bipolar disorder, Cognitive performance, Neuroblastoma (high-risk), Normalized brain volume, Intelligence.               

### Homer TFBSs
#### Unique enhancers between GW13 and GW17 

* Intersect of unique enhancers between GE01 vs GE04 and GE02 vs GE04. There are 5840 GW13-specific enhancers, and 17909 GW17-specific enhancers.          
* There are 52 TFs significantly (Benjamini q-value < 0.01) enriched in GW13-specific enhancers, and 117 in GW17-specific enhancers. Among them, 42 TFs are in common.        
* In the common TFs, 9 are present in > 20% enhancers in either GW13 or GW17. Their percent occupancy in GW13 and GW17 specific enhancers are similar, and downstream target genes are significantly overlapped (for example, Sox3 hypergeometric test p-value = 0).      
* There are only 10 GW13-specific TFs, and they all present in < 5% of the enhancers.    
* There are 75 GW17-specific TFs, and 5 of them are present in > 20% of the enhancers, including Olig2.      
* GW17-specific enhancers overlapped with Olig2 binding sites are associated with 3 GW-specific DE genes in cortex, and 28 in GE. Both gene lists are enriched for neuron development terms.              

![](HisMod_files/figure-html/unique_enhancer_homer_GW-1.png) ![](HisMod_files/figure-html/unique_enhancer_homer_GW-2.png) ![](HisMod_files/figure-html/unique_enhancer_homer_GW-3.png) ![](HisMod_files/figure-html/unique_enhancer_homer_GW-4.png) ![](HisMod_files/figure-html/unique_enhancer_homer_GW-5.png) ![](HisMod_files/figure-html/unique_enhancer_homer_GW-6.png) ![](HisMod_files/figure-html/unique_enhancer_homer_GW-7.png) 

* Olig2 downstream DE genes in cortex (distance to TSS < 10kb)    
  + EMX2: Empty Spiracles Homeobox 2. acts to generate the boundary between the roof and archipallium in the developing brain. May function in combinations with OTX1/2 to specify cell fates in the developing central nervous system.       
  + RGS10: Regulator Of G-Protein Signaling 10. associated with schizophrenia.     
  + LGI1: Leucine-Rich, Glioma Inactivated 1. This gene is predominantly expressed in neural tissues and its expression is reduced in low grade brain tumors and significantly reduced or absent in malignant gliomas. Mutations in this gene result in autosomal dominant lateral temporal epilepsy. May play a role in the control of neuroblastoma cell survival.     
  
![Olig2_cortex_DE](./HisMod_files/figure-html/homer_unique_enhancer_GW_GW17only_Olig2_DE_cortex_network.png)          

* Olig2 downstream DE genes in GE (distance to TSS < 10kb)      
  + PAX6: Paired Box 6. expressed in the developing nervous system, and in developing eyes. Mutations in this gene are known to cause ocular disorders such as aniridia and Peter's anomaly.     
  + HIVEP2: Human Immunodeficiency Virus Type I Enhancer Binding Protein 2. zinc finger-containing transcription factors. The encoded protein regulates transcription by binding to regulatory regions of various cellular and viral genes that maybe involved in growth, development and metastasis.    
  
![Olig2_GE_DE](./HisMod_files/figure-html/homer_unique_enhancer_GW_GW17only_Olig2_DE_GE_network.png)       

#### Unique enhancers between cortex and GE

* Intersect of unique enhancers between cortex01 vs GE01 and cortex02 vs GE02. There are 3279 cortex-specific enhancers, and 5731 GE-specific enhancers.     
* There are 42 TFs significantly (Benjamini q-value < 0.01) enriched in cortex-specific enhancers, and 78 in GE-specific enhancers. Among them, 32 TFs are in common.        
* In the common TFs, 11 are present in > 20% enhancers in either cortex or GE. Their percent occupancy in cortex and GE specific enhancers are similar, and downstream target genes are significantly overlapped (for example, Lhx3 hypergeometric test p-value = 0).      
* There are only 10 cortex-specific TFs, and they all present in < 20% of the enhancers.    
* There are 46 GW17-specific TFs, and 4 of them are present in > 20% of the enhancers, including Olig2.      
* 252 GE-specific enhancers overlapped with Olig2 binding sites are associated with NPC-specific DE genes, and the DE genes are enriched for neuron development terms.              

![](HisMod_files/figure-html/unique_enhancer_homer_neurospheres-1.png) ![](HisMod_files/figure-html/unique_enhancer_homer_neurospheres-2.png) ![](HisMod_files/figure-html/unique_enhancer_homer_neurospheres-3.png) ![](HisMod_files/figure-html/unique_enhancer_homer_neurospheres-4.png) ![](HisMod_files/figure-html/unique_enhancer_homer_neurospheres-5.png) ![](HisMod_files/figure-html/unique_enhancer_homer_neurospheres-6.png) 


