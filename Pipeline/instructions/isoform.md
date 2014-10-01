## Step1 DEfine on genes and exons       
  * Run DEfine on both coding genes and coding exons with __same__ parameters; see [DEfine.md](./DEfine.md)       
  * Result filename format: `UP/DN.<cell1>-<donor1>_<cell2>-<donor2>.FDR_<fdr cutoff>.rmin_<rmin>.Nmin_<Nmin>`           
  * Result file format: `Gene id, rpkm1, rpkm2, p-value, Multiple testing corrected p-value`     

## Step2 isoform identification 
  + [isoform.R](./isoform.R)      

### Criteria: 
  + DEfine on exons 
  + gene is expressed in both samples (`gene RPKM > RPKMmin`)
  + `exon RPKM < cutoff*gene RPKM` for one sample & `> cutoff2*gene RPKM` for the other 
  + gene is not DE         

### Parameters: 
  + lib1, lib2: library IDs
  + cell1, cell2: cell types
  + donor1, donor2: individuals
  + [cutoff, cutoff2]: cutoff for exons expressed or not. `Exon RPKM < cutoff*geneRPKM` are considered not expressed, default set to 1%; `exon RPKM > cutoff2*geneRPKM` are considered expressed, default set to 10%
  + [RPKMmin]: min RPKM for gene to be considered expressed, default set to 0.01
  + [fdr, rmin, Nmin]: FDR, rmin and Nmin used in DEfine, should be the same for both exons and genes, default set to 0.01, 0.005, 25
  + [dirExon]: directory to DE exon results from DEfine, default to `<current working directory>../DEfine/exon/`
  + [dirGene]: directory to DE gene results from DEfine, default to `<current working directory>../DEfine/gene/`
  + [dirOut]: output directory, default to `current working directory`
  + [dirIn1, dirIn2]: path to all libraries, default set to the same; default to `/projects/epigenomics/ep50/internal/jqc.1.7.6/`
  + [RPKM1, RPKM2]: gene RPKM files, default to `<dirIn><lib>/coverage/<lib1>.G.A.rpkm.pc`
  + [geneUP, geneDN]: gene DEfine results, default to `<dirGene>UP/DN.<cell1>-<donor1>_<cell2>-<donor2>.FDR_<fdr>.rmin_<rmin>.Nmin_<Nmin>`
  + [exonUP, exonDN]: exon DEfine results, default to `<dirExon>UP/DN.<cell1>-<donor1>_<cell2>-<donor2>.FDR_<fdr>.rmin_<rmin>.Nmin_<Nmin>`   

### Required input: 
  + DEfine on genes and exons: run with __same__ parameters
  + gene RPKM file for both library: `<dirIn><lib>/coverage/<lib>.G.A.rpkm.pc`; file format: `ENSG_id Nreads  gene_RPKM average_RPKM  min_RPKM  max_RPKM`      

### Output:
  + cassette exons: `<dirOut><cell1>-<donor1>_<cell2>-<donor2>_isoform.txt`
  + isoform genes: `<dirOut><cell1>-<donor1>_<cell2>-<donor2>_isoform_gene.txt`
  + return list:
    * summary: No. of DE_genes, DE_exons, with_expressed_genes, isoform_exons, exclude_DE_genes, isoform_genes
    * isoform_exon: cassette exons 
    * isoform_gene: genes identified as isoforms

### Usage: 
  + source on xhost: `source("/home/lli/bin/R-3.0.2/isoform.R")`          
  + example:   
```
lib1='A03473'; cell1='Cortex'; donor1='HuFNSC01';
lib2='A03474'; cell2='GE'; donor2='HuFNSC01';
cortex01_GE01 <- isoform(lib1 = lib1, lib2 = lib2, cell1 = cell1, cell2 = cell2, donor1 = donor1, donor2 = donor2)
cortex01_GE01_summary <- cortex01_GE01$summary
cortex01_GE01_isoform_exon <- cortex01_GE01$isoform_exon
cortex01_GE01_isoform_gene <- cortex01_GE01$isoform_gene
rm(cortex01_GE01)
``` 

  
  
