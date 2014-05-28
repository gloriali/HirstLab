## Step 1: signal calculation with FindER
* Run on xhost             
* Chromosome sizes: /home/mbilenky/UCSC_chr/hg19_auto_XY.chrom.sizes               
* Mappability: /projects/epigenomics/resources/UCSC_hg19/mappability/mappability_lt_1             
* `/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx30G /home/mbilenky/bin/Solexa_Java/FindER.jar -chr <chrom sizes> -i <input wig file> -o <output directory> -n <name prefix> -v -m <mappability file> -print -bin 200 -step 200`          
* e.g. `../FetalBrain/HisMod/FindER.sh`

### FindER output
* .signal.bedGraph.gz: (for ChromHMM use) < chr >    < start > < end >   < -log10(p-value) >   
* .log10pval12.0.bedGraph.gz: (thresholded peaks) < chr >   < start >   < end > < # >    
* .info: sanity check for Misha

### peak file sanity check
* No. of enriched regions: `for file in *.log10pval12.0.bedGraph.gz; do less $file | wc -l; done`          
* No. of total bases enriched: `for file in *.log10pval12.0.bedGraph.gz; do less $file | awk '{s+=$3-$2} END {print s}'; done`
* Average peak length: `for file in *.log10pval12.0.bedGraph.gz; do less $file | awk '{s+=$3-$2; n++} END {print s/n}'; done`      
```{r}
(Npeak <- ggplot(finder, aes(x = sample, y = Npeak_chromHMM, group = factor(Mark))) + 
   geom_line(aes(color = factor(Mark))) + 
   theme_bw())
ggsave(Npeak, file = "FindER_Npeak.pdf")
(Nbase <- ggplot(finder, aes(x = sample, y = Nbase_chromHMM, group = factor(Mark))) + 
   geom_line(aes(color = factor(Mark))) + 
   theme_bw())
ggsave(Nbase, file = "FindER_Nbase.pdf")
```
