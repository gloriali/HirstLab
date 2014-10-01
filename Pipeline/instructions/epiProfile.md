## Step 1. Calculate signal profiles
* Coordinates and boundaries files
    + exons: `exons=/home/lli/hg19/hg19v65_exons_for_genes`   
    + exon boundaries: `exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique` and `exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique`
    + exon boundary bin 20bp: `exons3p=/home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique_bin_20` and `exons5p=/home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique_bin_20`   
    + introns: `introns=/home/lli/hg19/hg19v65_introns_for_genes`   
    + intron boundaries: `introns3p=/home/lli/hg19/hg19v65_introns_for_genes.3prime_200` and `introns5p=/home/lli/hg19/hg19v65_introns_for_genes.5prime_200`
    + intron boundary bin 20bp: `introns3p=/home/lli/hg19/hg19v65_introns_for_genes.3prime_200_bin_20` and `introns5p=/home/lli/hg19/hg19v65_introns_for_genes.5prime_200_bin_20`   
* WGBS @ exon boundaries         
```
bed=<fractional methylation call bedGraph file>; name=<cell_type><donor>_WGBS;  
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique  
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique  
out=<output directory>/exons3p_200/; reg=$exons3p  
#out=<output directory>/exons5p_200/; reg=$exons5p  
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromBEDCalculator.jar -b $bed -r $reg -o $out -s hg19 -n $name -bin 20 -t Y -trueBED  
```
* H3K36me3 @ exon body          
```
wig=<H3K36me3 wig file>; name=<cell_type><donor>_H3K36me3;  
exons=/home/lli/hg19/hg19v65_exons_for_genes  
out=<output directory>/exons/; reg=$exons  
chr=/projects/epigenomics/resources/UCSC_chr/hg19.chrom.sizes  
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsCoverageFromWigCalculator.jar -w $wig -r $reg -o $out -s $chr -n $name  
```
* Other marks @ exon boundaries: MeDIP, H3K4me1, H3K4me3, H3K9me3, H3K27me3  
```
wig=<mark wig file>; name=<cell_type><donor>_<mark>;  
exons3p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.3prime_200.unique
exons5p=/projects/epigenomics/resources/Ensembl/hg19v65/hg19v65_exons_for_genes.5prime_200.unique
out=<output directory>/exons3p_200/; reg=$exons3p  
#out=<output directory>/exons5p_200/; reg=$exons5p  
/gsc/software/linux-x86_64-centos5/java-1.7.0-u13/bin/java -jar -Xmx80G /home/mbilenky/bin/Solexa_Java/RegionsProfileFromWigCalculator.jar -w $wig -r $reg -o $out -s hg19 -n $name -bin 20 -t Y 
```  
* No. of CpGs @ exon boundaries   
```
exons3p=/home/lli/hg19/hg19v65_exons_for_genes.3prime_200.unique_bin_20
exons5p=/home/lli/hg19/hg19v65_exons_for_genes.5prime_200.unique_bin_20
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $exons3p -b /home/lli/hg19/CG.BED -c > <output directory>/exons3p_200/CpG.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt   
/gsc/software/linux-x86_64-centos5/bedtools-2.17.0/bin/intersectBed -a $exons5p -b /home/lli/hg19/CG.BED -c > <output directory>/exons5p_200/CpG.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt   
less <output directory>/exons3p_200/CpG.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt | awk 'BEGIN {id=""; record=""} {if ($4 == id){record=(record" "$5)} else {print record; id=$4; record=(id" "$5)}}' > <output directory>/exons3p_200/CpG.hg19v65_exons_for_genes.3prime_200.unique
rm <output directory>/exons3p_200/CpG.hg19v65_exons_for_genes.3prime_200.unique_bin_20.txt
less <output directory>/exons5p_200/CpG.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt | awk 'BEGIN {id=""; record=""} {if ($4 == id){record=(record" "$5)} else {print record; id=$4; record=(id" "$5)}}' > <output directory>/exons5p_200/CpG.hg19v65_exons_for_genes.5prime_200.unique
rm <output directory>/exons5p_200/CpG.hg19v65_exons_for_genes.5prime_200.unique_bin_20.txt
```

## Step 2. Calculate and plot average profile for isoform exons / intron retention
* [epiProfile.R](./epiProfile.R)

### Required input:   
* H3K36me3: exon coverage results from RegionsCoverageFromWigCalculator.jar, name <cell><donor>_H3K36me3       
* WGBS: exon boundary profile results from RegionsProfileFromBEDCalculator.jar, name <cell><donor>_WGBS       
* other marks: exon boundary profile results from RegionsProfileFromWigCalculator.jar, name <cell><donor>_mark        
* isoform results from isoform function between two cell types        
* exon ID in exon profiles and isoforms should match format

### Parameters:   
* mark: name of the epigenetic mark  
* cell1, cell2: name of the two cell types used for epigenetic profile    
* donor1, donor2: name of the two individuals used for epigenetic profile    
* dirIn: input directory to exon signal coverage and exon boundary profiles   
* dirOut: output directory, default to current working directory  
* both, neither, cell1_specific, cell2_specific: exon IDs for exons expressed in both, not expressed, expressed only in cell1, and expressed only in cell2  
* geneRPKM: required for H3K36me3, list of exon IDs for different gene RPKM groups

### Output:   
* pdf figure in dirOut  
* return exon info and summary stats for H3K36me3   
* return exon boundaries profile summary for all other marks  

### Usage:
* source on xhost: `source("/home/lli/bin/R-3.0.2/epiProfile.R")`          
* source locally: `source("~/HirstLab/Pipeline/epiProfile.R")`    
* example: `WGBS_boundaries <- epiProfile(mark = "WGBS", cell1 = "lum", cell2 = "myo", donor1 = "RM066", donor2 = "RM045", dirIn = "~/REMC/epiProfile/", both = both, neither = neither, cell1_specific = lum_specific, cell2_specific = myo_specific)`   
