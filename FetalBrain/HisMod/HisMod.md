# Fetal Brain - Histone modifications
Gloria Li  
January 27, 2014  

Update Wed Jan 28 20:59:29 2015



## Sanity check   

* No. of peaks, No. of enriched bases, and average peak length seem reasonable except for the unusual high No. of peaks in GE HuFNSC04 input library.              

![](HisMod_files/figure-html/summary-1.png) 

## Correlation with protein-coding gene RPKM 

* Overlapping H3K4me3 and H3K27me3 with protein-coding gene promoters (TSS +/- 1500bp), and overlapping H3K36me3 with genebody.    
* Overall correlations are as expected, with H3K4m3 and H3K36me3 marked genes showing higher RPKM, and H3K27me3 and bivalent promoter marked showing lower RPKM.   
* H3K4me3 unmarked genes in GE HuFNSC04 shows significantly higher RPKM than any other samples, and No. of genes marked by H3K4me3 are much lower than other samples. Biology or bad quality library?    

![](HisMod_files/figure-html/RPKM-1.png) 

## Overlapping enhancers with WGBS UMRs

* UMRs between neurospheres (cortex vs GE) are enriched in enhancers (H3K4me1 enriched regions).      
* Between gestational weeks, GW13 UMRs are enriched in enhancers, but GW17 UMRs are not.      

<table>
 <thead>
  <tr>
   <th style="text-align:right;"> Sample </th>
   <th style="text-align:right;"> UMR </th>
   <th style="text-align:right;"> No.enhancers </th>
   <th style="text-align:right;"> No.UMR </th>
   <th style="text-align:right;"> No.enhancerUMR </th>
   <th style="text-align:right;"> percent </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> neurosphere02 </td>
   <td style="text-align:right;"> hyper </td>
   <td style="text-align:right;"> 79512 </td>
   <td style="text-align:right;"> 420 </td>
   <td style="text-align:right;"> 157 </td>
   <td style="text-align:right;"> 0.37 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> neurosphere02 </td>
   <td style="text-align:right;"> hypo </td>
   <td style="text-align:right;"> 79512 </td>
   <td style="text-align:right;"> 1758 </td>
   <td style="text-align:right;"> 631 </td>
   <td style="text-align:right;"> 0.36 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> neurosphere04 </td>
   <td style="text-align:right;"> hyper </td>
   <td style="text-align:right;"> 79512 </td>
   <td style="text-align:right;"> 91 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 0.33 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> neurosphere04 </td>
   <td style="text-align:right;"> hypo </td>
   <td style="text-align:right;"> 79512 </td>
   <td style="text-align:right;"> 555 </td>
   <td style="text-align:right;"> 231 </td>
   <td style="text-align:right;"> 0.42 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> GW_Cortex </td>
   <td style="text-align:right;"> hyper </td>
   <td style="text-align:right;"> 79512 </td>
   <td style="text-align:right;"> 179 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 0.15 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> GW_Cortex </td>
   <td style="text-align:right;"> hypo </td>
   <td style="text-align:right;"> 79512 </td>
   <td style="text-align:right;"> 1104 </td>
   <td style="text-align:right;"> 570 </td>
   <td style="text-align:right;"> 0.52 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> GW_GE </td>
   <td style="text-align:right;"> hyper </td>
   <td style="text-align:right;"> 79512 </td>
   <td style="text-align:right;"> 722 </td>
   <td style="text-align:right;"> 163 </td>
   <td style="text-align:right;"> 0.23 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> GW_GE </td>
   <td style="text-align:right;"> hypo </td>
   <td style="text-align:right;"> 79512 </td>
   <td style="text-align:right;"> 2280 </td>
   <td style="text-align:right;"> 1166 </td>
   <td style="text-align:right;"> 0.51 </td>
  </tr>
</tbody>
</table>
![](HisMod_files/figure-html/UMR-1.png) 


