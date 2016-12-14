# Glioma - Histone modifications ChIP-seq
Gloria Li  
Jun, 19, 2016  

Updated: Wed Dec 14 05:36:10 2016



## Global comparison
* Calculate an average signal of 500bp bins for both samples.           
* Look into only bins that that have non zero signal in both.       
* Find bins that are common in top 400 bins in both samples, likely blacklist regions, e.g. tandem repeats.       
* Plot box plot for the ratio log2(signal2/signal1) in top bins.       
* Compare it the same ratio for the rest of the bins.       
* Results and conclusion: inconsistency in input samples suggest that this method is not accurate when comparing samples from different biological origins, maybe due to genetic differences in these blacklist regions, e.g. different number of copies of the repeats. Thus these top ranked regions cannot be used as normalization standard.   

![](ChIPseq_files/figure-html/global-1.png)<!-- -->![](ChIPseq_files/figure-html/global-2.png)<!-- -->![](ChIPseq_files/figure-html/global-3.png)<!-- -->

## FindER enriched regions 
* FindER v 1.0.0b      

![](ChIPseq_files/figure-html/ER_summary-1.png)<!-- -->

## Unique enriched regions 
* Pairwise comparisons between glioma samples (S1) and NPC GE04 (S2).   
* Challenge: S1 and S2 often have quite different sequencing depth, thus direct comparisons between ER are likely biased.     
* Input: ER in S1, wig in S1, ER in S2, wig in S2.   
* Method: 
	+ Calculate S1 ER signal in S2 wig.    
	+ Randomly generate a set of background regions (not overlapping with S2 ER) with same length as S1 ER.   
	+ Calculate signal for these background regions in S2 wig, use 90% quantile as cutoff. 
	+ S1 unique ER: S1 ER with signal < background cutoff in S2.   
	
![](ChIPseq_files/figure-html/unique_ER_summary-1.png)<!-- -->

### Loss of H3K36me3      
* Global loss of H3K36me3 in both IDH mutant and wt gliomas.   
* H3K36me3 loss regions are enriched in genebody but not intergenic regions, while H3K36me3 gain regions are also enriched in intergenic regions.      

![](ChIPseq_files/figure-html/K36_loss-1.png)<!-- -->![](ChIPseq_files/figure-html/K36_loss-2.png)<!-- -->

#### H3K36me3 loss regions are consistent in all glioma samples

![](ChIPseq_files/figure-html/K36_intersect-1.png)<!-- -->![](ChIPseq_files/figure-html/K36_intersect-2.png)<!-- --><table>
 <thead>
  <tr>
   <th style="text-align:center;"> Sample1 </th>
   <th style="text-align:center;"> Sample2 </th>
   <th style="text-align:center;"> N1 </th>
   <th style="text-align:center;"> N2 </th>
   <th style="text-align:center;"> Intersect </th>
   <th style="text-align:center;"> FC </th>
   <th style="text-align:center;"> Pvalue </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> 121025 </td>
   <td style="text-align:center;"> 86378 </td>
   <td style="text-align:center;"> 64807 </td>
   <td style="text-align:center;"> 29.2078 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> 121025 </td>
   <td style="text-align:center;"> 134598 </td>
   <td style="text-align:center;"> 101956 </td>
   <td style="text-align:center;"> 26.8387 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> 121025 </td>
   <td style="text-align:center;"> 126327 </td>
   <td style="text-align:center;"> 94837 </td>
   <td style="text-align:center;"> 27.7623 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> 121025 </td>
   <td style="text-align:center;"> 122223 </td>
   <td style="text-align:center;"> 79591 </td>
   <td style="text-align:center;"> 23.5931 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> 86378 </td>
   <td style="text-align:center;"> 134598 </td>
   <td style="text-align:center;"> 64738 </td>
   <td style="text-align:center;"> 24.9022 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> 86378 </td>
   <td style="text-align:center;"> 126327 </td>
   <td style="text-align:center;"> 62212 </td>
   <td style="text-align:center;"> 26.6623 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> 86378 </td>
   <td style="text-align:center;"> 122223 </td>
   <td style="text-align:center;"> 63933 </td>
   <td style="text-align:center;"> 27.7202 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> 134598 </td>
   <td style="text-align:center;"> 126327 </td>
   <td style="text-align:center;"> 98600 </td>
   <td style="text-align:center;"> 24.6980 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> 134598 </td>
   <td style="text-align:center;"> 122223 </td>
   <td style="text-align:center;"> 84817 </td>
   <td style="text-align:center;"> 21.5346 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> 126327 </td>
   <td style="text-align:center;"> 122223 </td>
   <td style="text-align:center;"> 74938 </td>
   <td style="text-align:center;"> 21.1345 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
</tbody>
</table>

#### H3K36me3 loss regions are associated with cell cycle regulation

![](ChIPseq_files/figure-html/K36_loss_GREAT-1.png)<!-- -->![](ChIPseq_files/figure-html/K36_loss_GREAT-2.png)<!-- -->

#### H3K36me3 loss regions are associated with transcriptional downregulation

![](ChIPseq_files/figure-html/K36_loss_RPKM-1.png)<!-- --><table>
 <thead>
  <tr>
   <th style="text-align:center;"> File </th>
   <th style="text-align:center;"> N_K36 </th>
   <th style="text-align:center;"> N_UP </th>
   <th style="text-align:center;"> N_DN </th>
   <th style="text-align:center;"> N_K36_UP </th>
   <th style="text-align:center;"> N_K36_DN </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> 8205 </td>
   <td style="text-align:center;"> 485 </td>
   <td style="text-align:center;"> 297 </td>
   <td style="text-align:center;"> 78 </td>
   <td style="text-align:center;"> 191 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> 7787 </td>
   <td style="text-align:center;"> 485 </td>
   <td style="text-align:center;"> 297 </td>
   <td style="text-align:center;"> 78 </td>
   <td style="text-align:center;"> 204 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> 8700 </td>
   <td style="text-align:center;"> 485 </td>
   <td style="text-align:center;"> 297 </td>
   <td style="text-align:center;"> 85 </td>
   <td style="text-align:center;"> 196 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> 8868 </td>
   <td style="text-align:center;"> 485 </td>
   <td style="text-align:center;"> 297 </td>
   <td style="text-align:center;"> 111 </td>
   <td style="text-align:center;"> 182 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> 9106 </td>
   <td style="text-align:center;"> 485 </td>
   <td style="text-align:center;"> 297 </td>
   <td style="text-align:center;"> 78 </td>
   <td style="text-align:center;"> 216 </td>
  </tr>
</tbody>
</table>

#### Subtle gain of H3K27me3   
* No significant correlation between H3K36me3 and H3K27me3 signals.    
* The intersect of H3K27me3 gain/loss regions between different samples are still significant.     
* Both H3K27me3 gain and loss regions are enriched in genebody.       
* The functions of genes associated with H3K27me3 gain/loss regions are different.      

![](ChIPseq_files/figure-html/unique_ER_K36_K27-1.png)<!-- --><table>
 <thead>
  <tr>
   <th style="text-align:center;"> Sample1 </th>
   <th style="text-align:center;"> Sample2 </th>
   <th style="text-align:center;"> K27 </th>
   <th style="text-align:center;"> N1 </th>
   <th style="text-align:center;"> N2 </th>
   <th style="text-align:center;"> Intersect </th>
   <th style="text-align:center;"> FC </th>
   <th style="text-align:center;"> Pvalue </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> gain </td>
   <td style="text-align:center;"> 12114 </td>
   <td style="text-align:center;"> 11376 </td>
   <td style="text-align:center;"> 3171 </td>
   <td style="text-align:center;"> 111.3230 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> gain </td>
   <td style="text-align:center;"> 12114 </td>
   <td style="text-align:center;"> 12604 </td>
   <td style="text-align:center;"> 3893 </td>
   <td style="text-align:center;"> 121.8710 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> gain </td>
   <td style="text-align:center;"> 12114 </td>
   <td style="text-align:center;"> 11066 </td>
   <td style="text-align:center;"> 2565 </td>
   <td style="text-align:center;"> 85.7780 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> gain </td>
   <td style="text-align:center;"> 12114 </td>
   <td style="text-align:center;"> 7314 </td>
   <td style="text-align:center;"> 2416 </td>
   <td style="text-align:center;"> 117.4300 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> gain </td>
   <td style="text-align:center;"> 11376 </td>
   <td style="text-align:center;"> 12604 </td>
   <td style="text-align:center;"> 3375 </td>
   <td style="text-align:center;"> 114.2300 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> gain </td>
   <td style="text-align:center;"> 11376 </td>
   <td style="text-align:center;"> 11066 </td>
   <td style="text-align:center;"> 2022 </td>
   <td style="text-align:center;"> 73.0374 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> gain </td>
   <td style="text-align:center;"> 11376 </td>
   <td style="text-align:center;"> 7314 </td>
   <td style="text-align:center;"> 2380 </td>
   <td style="text-align:center;"> 124.8790 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> gain </td>
   <td style="text-align:center;"> 12604 </td>
   <td style="text-align:center;"> 11066 </td>
   <td style="text-align:center;"> 2609 </td>
   <td style="text-align:center;"> 84.0968 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> gain </td>
   <td style="text-align:center;"> 12604 </td>
   <td style="text-align:center;"> 7314 </td>
   <td style="text-align:center;"> 2509 </td>
   <td style="text-align:center;"> 117.5300 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> gain </td>
   <td style="text-align:center;"> 11066 </td>
   <td style="text-align:center;"> 7314 </td>
   <td style="text-align:center;"> 1343 </td>
   <td style="text-align:center;"> 67.6098 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> loss </td>
   <td style="text-align:center;"> 17238 </td>
   <td style="text-align:center;"> 11725 </td>
   <td style="text-align:center;"> 6180 </td>
   <td style="text-align:center;"> 179.2170 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> loss </td>
   <td style="text-align:center;"> 17238 </td>
   <td style="text-align:center;"> 20772 </td>
   <td style="text-align:center;"> 12592 </td>
   <td style="text-align:center;"> 193.0650 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> loss </td>
   <td style="text-align:center;"> 17238 </td>
   <td style="text-align:center;"> 22221 </td>
   <td style="text-align:center;"> 12401 </td>
   <td style="text-align:center;"> 183.8320 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_19 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> loss </td>
   <td style="text-align:center;"> 17238 </td>
   <td style="text-align:center;"> 18407 </td>
   <td style="text-align:center;"> 9116 </td>
   <td style="text-align:center;"> 161.8510 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> loss </td>
   <td style="text-align:center;"> 11725 </td>
   <td style="text-align:center;"> 20772 </td>
   <td style="text-align:center;"> 7637 </td>
   <td style="text-align:center;"> 175.0760 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> loss </td>
   <td style="text-align:center;"> 11725 </td>
   <td style="text-align:center;"> 22221 </td>
   <td style="text-align:center;"> 6610 </td>
   <td style="text-align:center;"> 146.5940 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_21 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> loss </td>
   <td style="text-align:center;"> 11725 </td>
   <td style="text-align:center;"> 18407 </td>
   <td style="text-align:center;"> 8126 </td>
   <td style="text-align:center;"> 215.8130 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> loss </td>
   <td style="text-align:center;"> 20772 </td>
   <td style="text-align:center;"> 22221 </td>
   <td style="text-align:center;"> 13828 </td>
   <td style="text-align:center;"> 162.2870 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_22 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> loss </td>
   <td style="text-align:center;"> 20772 </td>
   <td style="text-align:center;"> 18407 </td>
   <td style="text-align:center;"> 11328 </td>
   <td style="text-align:center;"> 159.2870 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> CEMT_23 </td>
   <td style="text-align:center;"> CEMT_47 </td>
   <td style="text-align:center;"> loss </td>
   <td style="text-align:center;"> 22221 </td>
   <td style="text-align:center;"> 18407 </td>
   <td style="text-align:center;"> 9639 </td>
   <td style="text-align:center;"> 130.9510 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
</tbody>
</table>

![](ChIPseq_files/figure-html/unique_ER_K36_K27-2.png)<!-- -->![](ChIPseq_files/figure-html/unique_ER_K36_K27-3.png)<!-- -->![](ChIPseq_files/figure-html/unique_ER_K36_K27-4.png)<!-- -->

#### No redistribution of H3K27me3      

![](ChIPseq_files/figure-html/K27me3_redistribution-1.png)<!-- -->

* [Nada's results](http://science.sciencemag.org/content/352/6287/844.figures-only)   

### Histone modification associated DE genes
* As expected, active marks, i.e. H3K27ac and H3K4me3 were significantly associated with transcriptional activation, while H3K27me3 was associated with down regulation.    

![](ChIPseq_files/figure-html/unique_ER_DE-1.png)<!-- -->

### Unique enhancers and transcription factor activities
* Homer transcription factor binding enrichment analysis for unique enhancers.   
	+ q-value < 0.01
	+ Fraction of enhancers with motif > 20%

#### H3K27ac
* Four transcription factors enriched in NPC-specific H3K27ac, only Sox3 was expressed in NPCs.     
	+ Sox3: function as a switch in neuronal development. Keeps neural cells undifferentiated by counteracting the activity of proneural proteins and suppresses neuronal differentiation.      

![](ChIPseq_files/figure-html/unique_H3K27ac_NPC-1.png)<!-- -->![](ChIPseq_files/figure-html/unique_H3K27ac_NPC-2.png)<!-- --><table>
 <thead>
  <tr>
   <th style="text-align:left;"> ENSG </th>
   <th style="text-align:left;"> Name </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000049323 </td>
   <td style="text-align:left;"> LTBP1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000075213 </td>
   <td style="text-align:left;"> SEMA3A </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000079931 </td>
   <td style="text-align:left;"> MOXD1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000112319 </td>
   <td style="text-align:left;"> EYA4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000115457 </td>
   <td style="text-align:left;"> IGFBP2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000116194 </td>
   <td style="text-align:left;"> ANGPTL1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000120738 </td>
   <td style="text-align:left;"> EGR1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000135862 </td>
   <td style="text-align:left;"> LAMC1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000137872 </td>
   <td style="text-align:left;"> SEMA6D </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000138193 </td>
   <td style="text-align:left;"> PLCE1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000138771 </td>
   <td style="text-align:left;"> SHROOM3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000147145 </td>
   <td style="text-align:left;"> LPAR4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000156140 </td>
   <td style="text-align:left;"> ADAMTS3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000158966 </td>
   <td style="text-align:left;"> CACHD1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000165891 </td>
   <td style="text-align:left;"> E2F7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000169436 </td>
   <td style="text-align:left;"> COL22A1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000171004 </td>
   <td style="text-align:left;"> HS6ST2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000177283 </td>
   <td style="text-align:left;"> FZD8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000185008 </td>
   <td style="text-align:left;"> ROBO2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000188467 </td>
   <td style="text-align:left;"> SLC24A5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000196353 </td>
   <td style="text-align:left;"> CPNE4 </td>
  </tr>
</tbody>
</table>
![](ChIPseq_files/figure-html/unique_H3K27ac_NPC_DAVID-1.png)<!-- -->![](ChIPseq_files/figure-html/unique_H3K27ac_NPC_DAVID-2.png)<!-- -->

* No significantly enriched transcription factor in glioma-specific H3K27ac.   

#### H3K4me1
* Three transcription factors were enirched and expressed in NPC-specific H3K4me1: Sox3, Sox6, Lhx2.     
	+ Sox3: function as a switch in neuronal development. Keeps neural cells undifferentiated by counteracting the activity of proneural proteins and suppresses neuronal differentiation.      
	+ Sox6: plays a key role in several developmental processes, including neurogenesis and skeleton formation.      
	+ Lhx2: acts as a transcriptional activator. Transcriptional regulatory protein involved in the control of cell differentiation in developing lymphoid and neural cell types.        
	
![](ChIPseq_files/figure-html/unique_H3K4me1_NPC-1.png)<!-- -->![](ChIPseq_files/figure-html/unique_H3K4me1_NPC-2.png)<!-- -->
![](ChIPseq_files/figure-html/unique_H3K4me1_NPC_DAVID-1.png)<!-- -->![](ChIPseq_files/figure-html/unique_H3K4me1_NPC_DAVID-2.png)<!-- -->![](ChIPseq_files/figure-html/unique_H3K4me1_NPC_DAVID-3.png)<!-- -->

* Three transcription factors were enriched and expressed in glioma-specific H3K4me1: HEB (TCF12), Ascl1, Olig2.      
	+ HEB: involved in the initiation of neuronal differentiation.     
	+ Ascl1: plays a role in the neuronal commitment and differentiation and in the generation of olfactory and autonomic neurons.      
	+ Olig2: is an essential regulator of ventral neuroectodermal progenitor cell fate. Required for oligodendrocyte and motor neuron specification in the spinal cord, as well as for the development of somatic motor neurons in the hindbrain.      
![](ChIPseq_files/figure-html/unique_H3K4me1_IDHmut-1.png)<!-- -->![](ChIPseq_files/figure-html/unique_H3K4me1_IDHmut-2.png)<!-- -->![](ChIPseq_files/figure-html/unique_H3K4me1_IDHmut-3.png)<!-- -->![](ChIPseq_files/figure-html/unique_H3K4me1_IDHmut-4.png)<!-- -->![](ChIPseq_files/figure-html/unique_H3K4me1_IDHmut-5.png)<!-- -->

## Chromatin states
* ChromHMM for 6 core histone marks in glioma and GE04.    
![](/projects/epigenomics2/users/lli/glioma/ChIPseq/ChromHMM/emission.pdf)

![](ChIPseq_files/figure-html/chromHMM-1.png)<!-- -->

