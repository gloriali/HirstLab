# FetalBrain - Figures
Gloria Li  
November 4, 2014  

Updated: Wed Oct  7 13:53:56 2015



## Figure 1: Epigenetic differences between MZ twins are on the same scale as differences between cell types, and are asymmetric between the twins.       
### Figure 1a: Experimental design overview   

![](./Figures_files/figure-html/design.png)     

### Figure 1b: Comparisons setup  

* Summary of sample information and libraries. Histone modifications include H3K4me1, H3K4me3, H3K9me3, H3K27me3, H3K36me3, and input. * marked samples without H3K4me3. Boxes outline the setup for pairwise comparisons, including comparing between MZ twins, between cortex and GE NPCs, and across three gestational weeks.    

![](./Figures_files/figure-html/comparisons.png)     
    
### Figure 1c: Scale of UMRs and DE genes between MZ twins and NPCs 

* MeDIP UMR frequency (top panel, bp/Mb), No. of unique enhancers (second panel), No. of differential expressed genes (third panel), and No. of isoform genes (bottom panel) between MZ twins. HuFNSC01-specific UMRs and genes upregulated in HuFNSC01 are shown in red, and HuFNSC02-specific UMRs and genes upregulated in HuFNSC02 are shown in blue.        
* _Note the opposite direction of UMRs and unique enhancers, why?_      

![](Figures_files/figure-html/MZ_scale-1.png) 

```
## png 
##   2
```

### Figure S1a: Scale of UMRs and DE genes between MZ twins and NPCs 

* MeDIP UMR frequency (top panel, bp/Mb), No. of unique enhancers (second panel), No. of differential expressed genes (third panel), and No. of isoform genes (bottom panel) between MZ twins (left panel), and between two NPC cell types (right panel). HuFNSC01-specific/cortex-specific UMRs and genes upregulated in HuFNSC01/cortex are shown in red, and HuFNSC02-specific/GE-specific UMRs and genes upregulated in HuFNSC02/GE are shown in blue.        
* _Note the opposite direction of UMRs and unique enhancers, why?_      

![](Figures_files/figure-html/MZ_scale_compare-1.png) 

```
## png 
##   2
```

### Figure S1b: UMR asymmetry
* UMRs frequency (bp/MB) for each chromosome for MeDIP UMRs between MZ twins    

![](Figures_files/figure-html/UMR_asymmetry_MZ-1.png) 

### Figure S1c: UMR locations
* Locations of UMRs along each chromosome for MeDIP UMRs between MZ twins                  

![](Figures_files/figure-html/UMR_pos_MZ-1.png) 

### Figure S1d: UMR genomic breakdown
* Fold enrichment on log2 scale for overlapping genomic regions with MeDIP UMRs between MZ twins    

![](Figures_files/figure-html/genomicBreak_MZ-1.png) 

### Figure S1e: DAVID for DE genes
* DAVID GO biological processes enriched for DE genes between MZ twins        

![](Figures_files/figure-html/DE_DAVID_MZ-1.png) 

### Figure S1f: Venn Diagram of DE genes
* Venn diagrame of No. of differential expressed genes between MZ twins in different cell types     
  
![](Figures_files/figure-html/Venn_DE_MZ-1.png) 

### Figure S1g: miRNA
* Heatmap of differential expressed and highly expressed (RPM > 100) miRNAs between MZ twins        

![](Figures_files/figure-html/miR_MZ-1.png) 

## Figure 2: UMRs between cortex and GE NPCs are more dynamic in GW17 than GW13, and lead to differential expression of key factors in brain development       
### Figure 2a: Scale of UMRs and DE genes between NPCs 

* WGBS UMR frequency (top panel, bp/Mb), fold enrichment of UMRs in enhancers (second panel), No. of differential expressed genes (third panel), and No. of isoform genes (bottom panel) between cortex and GE NPCs in GW13 (left panel), and GW17 (right panel). Cortex-specific UMRs and genes upregulated in cortex are shown in red, and GE-specific UMRs and genes upregulated in GE are shown in blue.        

![](Figures_files/figure-html/NPCs_scale-1.png) 

```
## png 
##   2
```

### Figure 2b: GREAT enrichment for NPCs UMRs

* GREAT GO biological processes terms statistically enriched (region-based binomial and hypergeometric FDR < 0.05) in cortex UMRs (red), and GE UMRs (blue).     

![](Figures_files/figure-html/NPCs_GREAT-1.png) 

### Figure 2c: GFAP and NFIX

* UCSC genome browser shot of examples of cortex UMR activated genes, top panel: Glial Fibrillary Acidic Protein (GFAP), bottom panel: Nuclear Factor I/X (CCAAT-Binding Transcription Factor, NFIX).     

![](./Figures_files/figure-html/GFAP.png)     
![](./Figures_files/figure-html/NFIX.png)   

### Figure 2d: DAVID for DE genes
* DAVID GO biological processes enriched for DE genes between MZ twins (A), between cortex and GE NPCs (B), and between gestational weeks (C).    

![](Figures_files/figure-html/DE_DAVID_NPC-1.png) 

### Figure S2a: UMR asymmetry

![](Figures_files/figure-html/UMR_asymmetry_NPC-1.png) 

### Figure S2b: UMR locations
* Locations of UMRs along each chromosome for WGBS UMRs between cortex and GE NPCs    

![](Figures_files/figure-html/UMR_pos_NPC-1.png) 

### Figure S2c: UMR genomic breakdown
* Fold enrichment on log2 scale for overlapping genomic regions with WGBS UMRs between cortex and GE NPCs    

![](Figures_files/figure-html/genomicBreak_NPC-1.png) 

### Figure S2d: Validate WGBS UMRs with MeDIP/MRE

* (A). For each UMR between NPCs identified by WGBS we calculated the normalized MeDIP-seq (methylated, top panel) and MRE-seq (unmethylated, bottom panel) signal. From this we show boxplot of methylation asymmetry between MeDIP-seq and MRE-seq signals in cortex and GE cells defined as (signal(cortex)-signal(GE))/(signal(cortex)+signal(GE)). (B). UMR frequency (bp/MB) across all chromosomes for MeDIP UMRs between cortex and GE NPCs. (C). GREAT GO biological processes terms enriched (region-based binomial and hypergeometric FDR < 0.05) in MeDIP cortex UMRs (red), and GE UMRs (blue). 

![](Figures_files/figure-html/WGBS_valid-1.png) 

```
## png 
##   2
```

### Figure S2e: DAVID for isoform
* DAVID enriched terms (FDR < 0.01) for isoforms between cortex and GE NPCs         

![](Figures_files/figure-html/isoform_DAVID_NPC-1.png) 

### Figure S2f: miRNA
* (A). Unique miRNA detected (>0.1 reads per million mapped) across all samples. Heatmap of differential expressed and highly expressed (RPM > 100) miRNAs between MZ twins (B), cortex and GE NPCs (C), and gestational weeks (D).         

![](Figures_files/figure-html/miR_NPC-1.png) 

## Figure 3: Differences between GW are asymmetric, with more GW17-specific UMRs and up-regulated genes    
### Figure 3a: UMRs and DE genes between GW 

* WGBS UMR frequency (top panel, bp/Mb), fold enrichment of UMRs in enhancers (second panel), No. of differential expressed genes (third panel), and No. of isoform genes (bottom panel) between GW13 and GW17 in NPCs cortex (left panel), and GE (right panel). GW13-specific UMRs and genes upregulated in GW13 are shown in red, and GW17-specific UMRs and genes upregulated in GW17 are shown in blue.        

![](Figures_files/figure-html/GW_scale-1.png) 

```
## png 
##   2
```

### Figure 3b: GW unique enhancer TFBSs

* Transcription factors exclusively statistically enriched (Benjamini corrected p-value < 0.01, and percent of enhancers with motif > 20%) in GW17 unique enhancers.     

![](Figures_files/figure-html/GW_enhancer_TFBS-1.png) 

```
## png 
##   2
```

### Figure 3c: TF-binding enhancers are hypomethylated

![](Figures_files/figure-html/GW_enhancer_TFBS_5mC-1.png) 

### Figure 3d: TF-binding enhancer targets are up-regulated

![](Figures_files/figure-html/GW_enhancer_TFBS_RPKM-1.png) 

### Figure S3a: UMR asymmetry
* UMRs frequency (bp/MB) for each chromosome for WGBS UMRs between GW13 and GW17     
  
![](Figures_files/figure-html/UMR_asymmetry_GW-1.png) 

### Figure S3b: UMR locations
* Locations of UMRs along each chromosome for WGBS UMRs between GW13 and GW17        

![](Figures_files/figure-html/UMR_pos_GW-1.png) 

### Figure S3c: UMR genomic breakdown
* Fold enrichment on log2 scale for overlapping genomic regions with WGBS UMRs between GW13 and GW17     

![](Figures_files/figure-html/genomicBreak_GW-1.png) 

### Figure S3d: GREAT enrichment for GW UMRs  

* GREAT GO biological processes enriched terms for UMRs between gestational weeks in cortex (A), GE (B), and shared by two cell types (C). GW13 UMRs are shown in red, and GW17 UMRs in blue.       

![](Figures_files/figure-html/GREAT_GW-1.png) 

### Figure 4a: OLIG2

* UCSC genome browser shot of Oligodendrocyte Lineage Transcription Factor 2 (OLIG2) hypomethylated in promoter region and upregulated in GW17.     

![](./Figures_files/figure-html/OLIG2.png)     

### Figure 4b: Heatmap of OLIG2 target genes   

![](Figures_files/figure-html/GW_OLIG2_heatmap-1.png) 

### Figure 4c: DAVID for OLIG2 target genes
* Gene Ontology biological processes enriched in OLIG2 downstream target genes in cortex NPCs (red), and GE NPCs (blue).    

![](Figures_files/figure-html/GW_OLIG2_DAVID-1.png) 

### Figure 4d: Differential expressed OLIG2 target genes   

* Cytoscape network of OLIG2 target genes that are differential expressed.       

![](./Figures_files/figure-html/homer_unique_enhancer_GW_GW17only_Olig2_DE_network.png)     

## Figure 5: Transcriptional activation in neurospheres happens at different gestational stages, with major wave at GW13-GW15 in cortex, and at GW15-GW17 in GE     
### Figure 5a: summary  

* No. of differential expressed genes between different gestational weeks (left panel: GW17 vs GW13, middle panel: GW17 vs GW15, right panel: GW15 vs GW13) in NPCs cortex (red), GE (blue), and shared by two cell types (purple). Bars pointing up shows upregulation in later stages, and bars pointing down shows upregulation in earlier stages.    

![](Figures_files/figure-html/GW_DE_summary-1.png) 

### Figure 5b: Shared by cortex and GE NPCs  

* Patterns of expression for genes differentially expressed between gestational weeks shared by cortex and GE NPCs. Genes are divided into eight expression profile groups, represented by eight different colours. The thickness of the line represents relative No. of genes in the category, and dashed line means no genes are in the category.                

![](Figures_files/figure-html/GW_DE_trend-1.png) 

### Figure 5c: Cortex NPCs  

* Patterns of expression for genes differentially expressed between gestational weeks in cortex NPCs.           

![](Figures_files/figure-html/GW_DE_trend_cortex-1.png) 

### Figure 5d: GE NPCs

* Patterns of expression for genes differentially expressed between gestational weeks in GE NPCs.       

![](Figures_files/figure-html/GW_DE_trend_GE-1.png) 

### Figure 5e,f: Heatmap of stage specific expressed genes  

* Heatmap for RPKM of stage-specific differential expressed genes in Figure 3c and 3d.         

![](Figures_files/figure-html/GW_DE_heatmap-1.png) 

```
## png 
##   2
```

### Figure S5a: Transcriptional clustering

* Unsupervised clustering on gene-level (A), exon-level (B), and miRNA (C) reveals fetal brain cell type relationships.         

![](Figures_files/figure-html/cluster-1.png) 

```
## png 
##   2
```

### Figure S5b: Venn diagram of differential expressed genes between GW 

* Venn diagrame of No. of differential expressed genes between different gestational weeks    
  
![](Figures_files/figure-html/Venn_DE_GW-1.png) 

### Figure S5c: DAVID enrichment for DE genes

* DAVID GO biological processes enriched for DE genes  between gestational weeks    

![](Figures_files/figure-html/DE_DAVID_GW-1.png) 

### Figure S5d: DAVID enrichment for isoform genes   

* DAVID enriched terms (FDR < 0.01) for isoforms across gestational weeks    

![](Figures_files/figure-html/isoform_DAVID_GW-1.png) 

### Figure S5e: miRNA 

* Heatmap of differential expressed and highly expressed (RPM > 100) miRNAs between gestational weeks           

![](Figures_files/figure-html/miR_GW-1.png) 


-----------------

### Supplemental 11: Fraction of DE genes with proximal UMRs

* Fraction of differential expressed genes proximally (TSS +/- 1500bp) associated with MeDIP UMRs between MZ twins (A), WGBS UMRs between cortex and GE NPCs (B), and WGBS UMRs between GW13 and GW17 (C).     

![](Figures_files/figure-html/proximal_DE-1.png) 

```
## png 
##   2
```

### Supplemental 12: TFBS asymmetry between neurosphere UMRs

* Relative abundance of transcription factor binding sites overlapping UMRs in NPCs cortex and GE cells reveals regulatory asymmetry between cell types. Red represents transcription factor more enriched in cortex, and blue represents transcription factor enriched in GE.    

![](Figures_files/figure-html/NPCs_TFBS-1.png) 

### Supplemental 14: Expression of DNA methylation reulators

* Heatmap of RPKM for DNA methylation regulators across all samples.         

![](Figures_files/figure-html/heatmap_5mC_regulator-1.png) 

### Supplemental 16: Correlating histone modifications with transcription

* Divide genes into 8 groups according to whether their genebody are marked by H3K36me3 (top panel: marked, bottom panel: not marked), and whether their promoters are marked by H3K4me3 and H3K27me3 (panels from left to right: marked by H3K4me3 only, both, H3K37me3 only, neither), and plot boxplot of gene RPKM on the log10 scale (RPKM < 0.005 are set to 0.005).   

![](Figures_files/figure-html/histone_RPKM-1.png) 

### Supplemental 17: DNA methylation at exon boundaries  

* Exon-intron junction mCpGs provide an inherited signature of exon expression. Average number of CpGs (black, bottom panel) and average mCpG levels (whole genome bisulfite shotgun, 20bp bins) at exon junctions +/- 200bp in cortex (solid line with round dots) and GE (dashed line with triangles). Exons are divided into four groups: 1) exons expressed in both cell types (exon RPKM > 0.1 in cortex and GE, purple); 2) cortex-specific exons (isoform exons expressed in cortex but not in GE, red); 3) GE-specific exons (isoform exons expressed in GE but not in cortex, green), and 4) exons not expressed in either cell types (all other exons, blue).           

![](Figures_files/figure-html/NPCs_epiProfile_5mC-1.png) 

```
## png 
##   2
```

### Supplemental 18: Neurosphere UMR enrichemnt at chromosome ends  

* Fraction of UMRs along normalized length of the genome in cortex and GE specific UMRs (A), and GW specific UMRs (B).          

![](Figures_files/figure-html/chrEnd_WGBS-1.png) 

```
## png 
##   2
```

### Supplemental 19: TFBS asymmetry in GW UMRs  

* Relative abundance of transcription factor binding sites overlapping UMRs between gestational weeks.       

![](Figures_files/figure-html/GW_TFBS-1.png) 

----------------------------------------
<!--
### Supplemental : GREAT analysis on UMRs between MZ twins - GOBP 


-->




