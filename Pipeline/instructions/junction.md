## Step 1. Generate junction coverage files   
+ Check strand specificity: `less <library.report> | grep 'Percent of junctions with <1% leakage'`: ~100% for strand specific libraries.     
+ For non-strand specific libraries: `/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh <RNA-seq bam file> <output directory>`   
+ For strand specific libraries: `/home/mbilenky/bin/Solexa_Shell/RunJunctions.sh <RNA-seq bam file> <output directory> S`   

## Step 2. Isoform validation       
+ [junction.R](./junction.R)
    
### Criteria:    
+ Enough junction coverage: sum of junction coverage in two libraries > `covcut`
+ Junction RPKM changes in the same direction as exon RPKM  
+ Junction RPKM > `jcut` in one sample and < `jcut` in the other
    
### Parameters:       
* lib1, lib2: library IDs   
* cell1, cell2: cell types   
* donor1, donor2: individuals   
* Nread1, Nread2: `total number of exonic reads for RPKM (protein coding; no MT, no ribo proteins, top expressed 0.005 exons excluded)` from RNA-seq library report   
* [read length1, read length2]: read length from samtools bam files for both libraries, default to `75`  
* [dirIn]: input directory for junction coverage files, default to `<current wd>`   
* [dirIsoform]: directory for isoform results, default to `<current wd>../isoform/`   
* [dirOut]: output directory, default to `<current wd>`   
* [jcut]: cutoff for junction RPKM, i.e.one sample > `jcut`, the other < `jcut`, default to `0.1`    
* [covcut]: cutoff for sum junction coverage of two samples to be considered enough junction coverage, default to `1`
    
### Required input       
* annotation files: `/home/lli/hg19/hg19v65_genes.Rdata` and `/home/lli/hg19/hg19v65_junctions_vs_exons_for_genes.relations.unique`   
* junction coverage file for both libraries: `<lib>.q5.F516.jb.bedGraph.gz`   
* RNA-seq library report: `<lib>.report`   
* isoform exon results: `<cell1>-<donor1>_<cell2>-<donor2>_isoform.txt`
    
### Output:        
* validated exons: `<dirOut><cell1>-<donor1>_<cell2>-<donor2>_isoform_valid.txt`  
* validated genes: `<dirOut><cell1>-<donor1>_<cell2>-<donor2>_isoform_valid_gene.txt`  
* return list:  
   + `summary`: `isoform exons, isoform_genes, exons_with_junction_cov, genes_with_junction_cov, exons_with_junction_support, genes with junction support`  
   + `isoform_valid_exon`: validated exons    
   + `isoform_valid_gene`: validated genes   

### Usage:           
* source on xhost: `source("/home/lli/HirstLab/Pipeline/R/junction.R")`          
* example:  

```
lib1='A03473'; cell1='Cortex'; donor1='HuFNSC01'; Nread1=Ncortex01;
lib2='A03474'; cell2='GE'; donor2='HuFNSC01'; Nread2=Nge01;
cortex01_GE01 <- junction(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2, Nread1 = Nread1, Nread2 = Nread2)
cortex01_GE01_summary <- cortex01_GE01$summary
cortex01_GE01_isoform_valid_exon <- cortex01_GE01$isoform_valid_exon
cortex01_GE01_isoform_valid_gene <- cortex01_GE01$isoform_valid_gene
rm(cortex01_GE01, lib1, lib2, cell1, cell2, donor1, donor2)
```


