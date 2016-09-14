# Glioma - Histone modifications ChIP-seq
Gloria Li  
Jun, 19, 2016  

Updated: Wed Sep 14 12:27:03 2016



## FindER enriched regions 
* Some libraries were much more deeply sequenced than others, i.e. CEMT_47, might lead to increased enriched regions due to increased sequencing depth.     
* Subsample bam files of deeply sequenced libraries to the average coverage of other libraries of the same mark (ran 50 times), and re-run FindER:   
	+ CEMT_21 H3K4me1: 0.4
	+ CEMT_47 H3K4me1: 0.4
	+ CEMT_47 H3K9me3: 0.5
	+ CEMT_47 H3K27me3: 0.4
	+ CEMT_47 H3K27ac: 0.3
* Subsampled libraries (panel 2 and 4) had less enriched regions, but were able to retain the majority of enriched regions.    
* FindER v 1.0.0b      

![](ChIPseq_files/figure-html/ER_summary-1.png)

## Unique enriched regions 
* Pairwise comparisons between glioma samples (S1) and NPC GE04 (S2).   
* Challenge: S1 and S2 often have quite different sequencing depth, thus direct comparisons between ER are likely biased.     
* Input: ER in S1, wig in S1, ER in S2, wig in S2.   
* Method: 
	+ Calculate S1 ER signal in S2 wig.    
	+ Randomly generate a set of background regions (not overlapping with S2 ER) with same length as S1 ER.   
	+ Calculate signal for these background regions in S2 wig, use 90% quantile as cutoff. 
	+ S1 unique ER: S1 ER with signal < background cutoff in S2.   
![](ChIPseq_files/figure-html/unique_ER_summary-1.png)

### Histone modification associated DE genes
* As expected, active marks, i.e. H3K27ac and H3K4me3 were significantly associated with transcriptional activation, while H3K27me3 was associated with down regulation.    

![](ChIPseq_files/figure-html/unique_ER_DE-1.png)

### Unique enhancers and transcription factor activities
* Homer transcription factor binding enrichment analysis for unique enhancers.   
	+ q-value < 0.01
	+ Fraction of enhancers with motif > 20%

#### H3K27ac
* Four transcription factors enriched in NPC-specific H3K27ac, only Sox3 was expressed in NPCs.     
	+ Sox3: function as a switch in neuronal development. Keeps neural cells undifferentiated by counteracting the activity of proneural proteins and suppresses neuronal differentiation.      

![](ChIPseq_files/figure-html/unique_H3K27ac_NPC-1.png)![](ChIPseq_files/figure-html/unique_H3K27ac_NPC-2.png)
![](ChIPseq_files/figure-html/unique_H3K27ac_NPC_DAVID-1.png)![](ChIPseq_files/figure-html/unique_H3K27ac_NPC_DAVID-2.png)

* No significantly enriched transcription factor in glioma-specific H3K27ac.   

#### H3K4me1
* Three transcription factors were enirched and expressed in NPC-specific H3K4me1: Sox3, Sox6, Lhx2.     
	+ Sox3: function as a switch in neuronal development. Keeps neural cells undifferentiated by counteracting the activity of proneural proteins and suppresses neuronal differentiation.      
	+ Sox6: plays a key role in several developmental processes, including neurogenesis and skeleton formation.      
	+ Lhx2: acts as a transcriptional activator. Transcriptional regulatory protein involved in the control of cell differentiation in developing lymphoid and neural cell types.        
	
![](ChIPseq_files/figure-html/unique_H3K4me1_NPC-1.png)![](ChIPseq_files/figure-html/unique_H3K4me1_NPC-2.png)
![](ChIPseq_files/figure-html/unique_H3K4me1_NPC_DAVID-1.png)![](ChIPseq_files/figure-html/unique_H3K4me1_NPC_DAVID-2.png)![](ChIPseq_files/figure-html/unique_H3K4me1_NPC_DAVID-3.png)

* Three transcription factors were enriched and expressed in glioma-specific H3K4me1: HEB (TCF12), Ascl1, Olig2.      
	+ HEB: involved in the initiation of neuronal differentiation.     
	+ Ascl1: plays a role in the neuronal commitment and differentiation and in the generation of olfactory and autonomic neurons.      
	+ Olig2: is an essential regulator of ventral neuroectodermal progenitor cell fate. Required for oligodendrocyte and motor neuron specification in the spinal cord, as well as for the development of somatic motor neurons in the hindbrain.      
![](ChIPseq_files/figure-html/unique_H3K4me1_IDHmut-1.png)![](ChIPseq_files/figure-html/unique_H3K4me1_IDHmut-2.png)![](ChIPseq_files/figure-html/unique_H3K4me1_IDHmut-3.png)![](ChIPseq_files/figure-html/unique_H3K4me1_IDHmut-4.png)![](ChIPseq_files/figure-html/unique_H3K4me1_IDHmut-5.png)

## Chromatin states
* ChromHMM for 6 core histone marks in glioma and GE04.    
![](/projects/epigenomics2/users/lli/glioma/ChIPseq/ChromHMM/emission.pdf)

![](ChIPseq_files/figure-html/chromHMM-1.png)

