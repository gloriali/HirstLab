Fetal Brain - MeDIP vs MethylCRF vs WGBS
========================================================

Gloria Li         
Tue Mar 11 21:48:00 2014 
<!-- re-knit after modify compare_bin200.R WGBS.R and MeDIP.R -->

## WGBS   
* Is methylation on two strands symmetric?   
    Seems the differences between strands are most likely result from technical variations.          
<img src="../WGBS/strand.symmetry.cortex02.png" height="400px" width="400px" />
<img src="../WGBS/strand.symmetry.distribution.png" height="400px" width="400px" />    
* Take the weighted average on coverage for fractional calls on both strand.     

## MeDIP
* MeDIP fractional calls for 28M CpGs generated for 6 brain (brain, cortex, GE) samples from HuFNSC01 and HuFNSC02 with Misha's algorithm.   

## MeDIP vs MethylCRF vs WGBS  
### Distribution and global patterns
* WGBS fractional calls from different samples seem more stable after taking the weighted average.     
* MeDIP calls vary more between different samples than other two methods. (not observed in breast datasets) 
  * Possible reason: sequencing fragment length in the brain samples varies more. But wig files are generated with average fragment length, so regional coverage calculated from the wig files are more unstable. -- Misha is looking into this.   
<img src="../MeDIP/ECDF.MeDIP.MethylCRF.WGBS.png" height="400px" width="400px" />          
<img src="../MeDIP/cortex02.3methods_boxplot.png" height="200px" width="600px" />            
<img src="../MeDIP/ge02.3methods_boxplot.png" height="200px" width="600px" />           



