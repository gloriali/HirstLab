load("~/hg19/hg19v65_genes.Rdata")
junction_exon <- read.delim("~/hg19/hg19v65_junctions_vs_exons_for_genes.relations.unique", head = T, as.is = T)
high_fold <- function(exon, up, covcut, junctions, junction_exons){
  e <- 1e-6
  junctions <- na.omit(junctions[junction_exons[junction_exons$exonID == exon, "junctionID"], ])
  junctions <- junctions[(junctions[, "N1"] + junctions[, "N2"]) > covcut,]
  if(nrow(junctions) == 0){
    return(c(NA, NA, NA)) # not enough coverage for junctionss of this gene
  }else{
    junctions$fold <- log2((junctions[, "rpkm1"] + e) / (junctions[, "rpkm2"] + e))
    if(up > 0){
      highest <- junctions[which.max(junctions$fold),]
    }else{
      highest <- junctions[which.min(junctions$fold),]
    }
    return(list(junctionID = as.character(highest[, "junctionID"]), junction1 = highest[, "rpkm1"], junction2 = highest[, "rpkm2"]))
  }
}
Vhigh_fold <- Vectorize(high_fold, vectorize.args = c("exon", "up"), SIMPLIFY = T)

junction <- function(lib1, lib2, cell1, cell2, donor1, donor2, Nread1, Nread2, read_length1 = 75, read_length2 = 75, jcut = 0.1, covcut = 1, 
                     dirIn = "", dirIsoform = "../isoform/", dirOut = ""){
  junction1 <- read.delim(paste0(dirIn, lib1, ".q5.F516.jb.bedGraph.gz"), as.is =T, head = F)
  junction1 <- junction1[!duplicated(junction1$V1), ]
  rownames(junction1) <- gsub("chr", "", junction1$V1)
  junction2 <- read.delim(paste0(dirIn, lib2, ".q5.F516.jb.bedGraph.gz"), as.is =T, head = F)
  junction2 <- junction2[!duplicated(junction2$V1), ]
  rownames(junction2) <- gsub("chr", "", junction2$V1)
  junction <- data.frame(junctionID = unique(c(rownames(junction1), rownames(junction2))), N1 = 0, N2 = 0)
  rownames(junction) <- junction$junctionID
  junction[rownames(junction1), "N1"] <- junction1$V3
  junction[rownames(junction2), "N2"] <- junction2$V3
  junction$rpkm1 <- junction$N1*10^3*10^6 / (Nread1*read_length1) 
  junction$rpkm2 <- junction$N2*10^3*10^6 / (Nread2*read_length2) 
  isoform_valid <- read.delim(paste0(dirIsoform, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform.txt"), as.is =T)
  isoform_exons <- nrow(isoform_valid)
  isoform_genes <- length(unique(isoform_valid$id))
  isoform_valid <- data.frame(exonID = isoform_valid$V1, geneID = isoform_valid$id, DEfine.p = isoform_valid$V5, gene1 = isoform_valid$gene1, gene2 = isoform_valid$gene2, exon1 = isoform_valid$V2, exon2 = isoform_valid$V3, 
                              t(Vhigh_fold(exon = as.character(isoform_valid$V1), up = isoform_valid$V2 - isoform_valid$V3, covcut = covcut, junctions = junction, junction_exons = junction_exon)))
  colnames(isoform_valid)[(ncol(isoform_valid) - 2):ncol(isoform_valid)] <- c("junctionID", "junction1", "junction2")
  isoform_valid$junctionID <- as.character(isoform_valid$junctionID)
  isoform_valid$junction1 <- as.numeric(isoform_valid$junction1)
  isoform_valid$junction2 <- as.numeric(isoform_valid$junction2)
  # No. of exons/genes with enough junction coverage 
  isoform_valid <- na.omit(isoform_valid)
  exons_with_junction_cov <- nrow(isoform_valid)
  genes_with_junction_cov <- length(unique(isoform_valid$geneID))
  # No. of exons/genes with junction support
  isoform_valid <- isoform_valid[(isoform_valid$junction1 >= jcut & isoform_valid$junction2 < jcut) | (isoform_valid$junction1 < jcut & isoform_valid$junction2 >= jcut), ]
  isoform_valid <- isoform_valid[order((isoform_valid$gene1 + isoform_valid$gene2), decreasing = T), ]
  isoform_valid[,c("coord", "chr", "start", "end", "strand", "type", "name", "description")] <- droplevels(ensembl[as.character(isoform_valid$geneID), c("coord", "chr", "start", "end", "strand", "type", "name", "description")])
  exons_with_junction_support <- nrow(isoform_valid)
  genes_with_junction_support <- length(unique(isoform_valid$geneID))
  write.table(isoform_valid, file = paste0(dirOut, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform_valid.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  isoform_valid_gene <- isoform_valid[!duplicated(isoform_valid$geneID), ]
  write.table(isoform_valid_gene, file = paste0(dirOut, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform_valid_gene.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  summary = c("isoform exons" = isoform_exons, "isoform genes" = isoform_genes, "exons with junction coverage" = exons_with_junction_cov, "genes with junction coverage" = genes_with_junction_cov, "exons with junction support" = exons_with_junction_support, "genes with junction support" = genes_with_junction_support)
  return(list(summary = summary, isoform_valid_exon = isoform_valid, isoform_valid_gene = isoform_valid_gene))
}
# isoform validation with junction reads
# Criteria: 
#   enough junction coverage: sum of junction coverage in two libraries > `covcut`
#   junction RPKM changes in the same direction as exon RPKM  
#   junction RPKM cutoff: junction RPKM > `jcut` in one sample and < `jcut` in the other   
# Parameters: 
#   lib1, lib2: library IDs   
#   cell1, cell2: cell types   
#   donor1, donor2: individuals   
#   Nread1, Nread2: `total number of exonic reads for RPKM (protein coding; no MT, no ribo proteins, top expressed 0.005 exons excluded)` from RNA-seq library report   
#   [read_length1, read_length2]: read length from samtools bam files for both libraries, default to 75   
#   [dirIn]: input directory for junction coverage files, default to `<current wd>`   
#   [dirIsoform]: directory for isoform results, default to `<current wd>../isoform/`   
#   [dirOut]: output directory, default to `<current wd>`   
#   [jcut]: cutoff for junction RPKM, i.e.one sample > jcut, the other < jcut, default to 0.1    
#   [covcut]: cutoff for sum junction coverage of two samples to be considered enough junction coverage, default to 1    
# Required input
#   annotation files: `~/hg19/hg19v65_genes.Rdata` and `~/hg19/hg19v65_junctions_vs_exons_for_genes.relations.unique`   
#   junction coverage file for both libraries: `<lib>.q5.F516.jb.bedGraph.gz`   
#   RNA-seq library report: `<lib>.report`   
#   isoform exon results: `<cell1>-<donor1>_<cell2>-<donor2>_isoform.txt`
# Output:
#  validated exons: `<dirOut><cell1>-<donor1>_<cell2>-<donor2>_isoform_valid.txt`  
#  validated genes: `<dirOut><cell1>-<donor1>_<cell2>-<donor2>_isoform_valid_gene.txt`  
#  return list: 
#   summary: isoform exons, isoform_genes, exons_with_junction_cov, genes_with_junction_cov, exons_with_junction_support, genes with junction support
#   isoform_valid_exon: validated exons   
#   isoform_valid_gene: validatedgenes   


