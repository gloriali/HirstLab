DMR_DE <- function(DMR_gene, DE, DM, name){
  load("~/hg19/hg19v65_genes.Rdata")
  DMR_gene$DMR <- gsub("_.*", "", DMR_gene$V4)
  DMR_gene <- DMR_gene[!duplicated(DMR_gene[, c("DMR", "V5")]),]
  DMR_gene <- data.frame(DMR = DMR_gene$DMR, ensembl[DMR_gene$V5, c("id", "name", "description", "coord", "strand")], stringsAsFactors = F)
  DMR_gene$DM <- DM
  DMR_gene_DE <- DMR_gene[DMR_gene$id %in% as.character(DE$V1), ]
  DMR_gene_DE$DE <- DE[DMR_gene_DE$id, "DE"]
  summary <- c(proximal.DMRs = length(unique(DMR_gene$DMR)), unique.genes = length(unique(DMR_gene$id)), 
               DE.DMRs = length(unique(DMR_gene_DE$DMR)), unique.DE.genes = length(unique(DMR_gene_DE$id)), 
               same.direction = sum((DMR_gene_DE$DM == "hyper" & DMR_gene_DE$DE == "DN")|(DMR_gene_DE$DM == "hypo" & DMR_gene_DE$DE == "UP")))
  write.table(DMR_gene, file = paste0("DMR_gene.", name, ".txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  write.table(DMR_gene_DE, file = paste0("DMR_gene_DE.", name, ".txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  return(list(summary = summary, DMR_gene = DMR_gene, DMR_gene_DE = DMR_gene_DE))
}
# DE genes for gene-associated DMRs 
# Required input: 
#   DMR_gene association list: output from DMR.intersect.sh;  
#   DE gene list: output from DEfine;  
# Parameters: 
#   DMR_gene: DMR_gene association list, format: chr, start, end, DMR.ID_CpG.ID, gene.ID
#   DE: DE gene list, format: gene.ID, rpkm1, rpkm2, p-value, adjusted.p-value, DE (UP/DN)
#   DM: hyper / hypo 
#   name: name for output files
# Output: 
#   DMR_gene.<name>.txt: unique DMR gene association list, format: DMR.ID, gene.ID, name, description, coord, strand, DM
#   DMR_gene_DE.<name>.txt: DE genes with gene-associated DMRs, format:  DMR.ID, gene.ID, name, description, coord, strand, DM, DE  
#   return list: 
#     summary: summary statistics, c(DMRs, unique.genes, DE.DMRs, unique.DE.genes, same.direction)
#     DMR_gene, DMR_gene_DE: same as output file
