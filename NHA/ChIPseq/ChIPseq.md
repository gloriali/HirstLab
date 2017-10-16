# VitC - ChIPseq
Gloria Li  
Sept 6, 2017  

Updated: Mon Oct 16 10:49:27 2017



## QC  
* Sequencing depth are even across samples, and sufficient (50M for narrow marks, 100M for broad marks).     
* All QC metric looks good.       
* All results suggest that NHAR VitC H3K9me3 sample is actually a H3K4me3 sample, probably due to using the wrong antibody in IP: remove H3K9me3 from further analysis.          

![](ChIPseq_files/figure-html/QC-1.png)<!-- -->![](ChIPseq_files/figure-html/QC-2.png)<!-- -->![](ChIPseq_files/figure-html/QC-3.png)<!-- -->

## Enrich regions
* From [UCSC tracks](http://www.bcgsc.ca/downloads/mb/VitC_glioma/HistoneHub/hub.txt), FindER results make more sense than MACS2.      

![](ChIPseq_files/figure-html/ER_summary-1.png)<!-- -->

## Unique enrich regions
* Unique enrich regions were identified by non-overlapping regions in pairwise comparisons.    

![](ChIPseq_files/figure-html/ER_unique_summary-1.png)<!-- -->

### Homer for unique H3K27ac 
* Results correlate those in CEMT glioma.    

![](ChIPseq_files/figure-html/homer_K27ac1-1.png)<!-- -->
![](ChIPseq_files/figure-html/homer_K27ac2-1.png)<!-- -->
![](ChIPseq_files/figure-html/homer_K27ac3-1.png)<!-- -->
![](ChIPseq_files/figure-html/homer_K27ac4-1.png)<!-- -->

