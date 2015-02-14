isoform <- function(lib1, lib2, cell1, cell2, donor1, donor2, cutoff = 0.01, cutoff2 = 0.1, RPKMmin = 0.01, fdr = 0.01, rmin = 0.005, Nmin = 25, 
                    dirExon = "../DEfine/exon/", dirGene = "../DEfine/gene/", dirOut = "", dirIn1 = "/projects/epigenomics/ep50/internal/jqc.1.7.6/", dirIn2 = dirIn1, 
                    RPKM1 = paste0(dirIn1, lib1, "/coverage/", lib1, ".G.A.rpkm.pc"), RPKM2 = paste0(dirIn2, lib2, "/coverage/", lib2, ".G.A.rpkm.pc"), 
                    geneUP = paste0(dirGene, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_", fdr, ".rmin_", rmin, ".Nmin_", Nmin), geneDN = paste0(dirGene, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_", fdr, ".rmin_", rmin, ".Nmin_", Nmin),
                    exonUP = paste0(dirExon, "UP.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_", fdr, ".rmin_", rmin, ".Nmin_", Nmin), exonDN = paste0(dirExon, "DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_", fdr, ".rmin_", rmin, ".Nmin_", Nmin)){
  # read gene RPKM files
  pc_lib1 <- read.delim(RPKM1, row.names = 1, as.is = T, head = F)
  pc_lib2 <- read.delim(RPKM2, row.names = 1, as.is = T, head = F)
  # read DE exon results
  exon_up <- read.delim(exonUP, as.is = T, head = F)
  exon_dn <- read.delim(exonDN, as.is = T, head = F)
  # read DE exon results
  gene_up <- read.delim(geneUP, as.is = T, head = F)
  gene_dn <- read.delim(geneDN, as.is = T, head = F)
  gene <- c(gene_up$V1, gene_dn$V1)
  (DE_genes <- length(gene))
  # get gene RPKM for DE exons
  exon <- rbind(exon_up, exon_dn)
  exon$id <- unlist(strsplit(as.character(exon$V1), '_'))[2*(1:nrow(exon))]
  exon$gene1 <- pc_lib1[exon$id, "V3"]
  exon$gene2 <- pc_lib2[exon$id, "V3"]
  (DE_exons <- nrow(exon))
  # exclude exons in not expressed genes
  exon <- na.omit(exon)
  exon <- exon[exon$gene1 > RPKMmin & exon$gene2 > RPKMmin, ]
  (with_expressed_genes <- nrow(exon))
  # identify isoform exons
  exon <-exon[(exon$V2 <= cutoff*exon$gene1 & exon$V3 >= cutoff2*exon$gene2)|(exon$V2 >= cutoff2*exon$gene1 & exon$V3 <= cutoff*exon$gene2),]
  (isoform_exons <- nrow(exon))
  # exclude DE genes
  exon <- exon[!(exon$id %in% gene), ]
  exon <- exon[order(exon$V5), ]
  (exclude_DE_genes <- nrow(exon))
  # isoform genes
  isoform_gene <- exon[!duplicated(exon$id), ]
  (isoform_genes <- nrow(isoform_gene))
  write.table(exon, file = paste0(dirOut, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform.txt"), sep = "\t", quote = F, row.names = F)
  write.table(isoform_gene, file = paste0(dirOut, cell1, "-", donor1, "_", cell2, "-", donor2, "_isoform_gene.txt"), sep = "\t", quote = F, row.names = F)
  summary <- c("Sample" = paste0(cell1, "-", donor1, "_", cell2, "-", donor2), "DE_genes" = DE_genes, "DE_exons" = DE_exons, "with_expressed_genes" = with_expressed_genes, "isoform_exons" = isoform_exons, "exclude_DE_genes" = exclude_DE_genes, "isoform_genes" = isoform_genes)
  return(list(summary = summary, isoform_exon = exon, isoform_gene = isoform_gene))
}
# isoform (cassete exon) identification 
# Criteria: 
#  DEfine on exons 
#  gene is expressed in both samples (gene RPKM > RPKMmin)
#  exon RPKM < cutoff * gene RPKM for one sample & > cutoff2 for the other 
#  gene is not DE
# Parameters: 
#  lib1, lib2: library IDs
#  cell1, cell2: cell types
#  donor1, donor2: individuals
#  [cutoff, cutoff2]: cutoff for exons expressed or not. Exon RPKM < cutoff*geneRPKM are considered not expressed, default set to 1%; exon RPKM > cutoff2*geneRPKM are considered expressed, default set to 10%
#  [RPKMmin]: min RPKM for gene to be considered expressed, default set to 0.01
#  [fdr, rmin, Nmin]: FDR, rmin and Nmin used in DEfine, should be the same for both exons and genes, default set to 0.01, 0.005, 25
#  [dirExon]: directory to DE exon results from DEfine, default to <current working directory>../DEfine/exon/
#  [dirGene]: directory to DE gene results from DEfine, default to <current working directory>../DEfine/gene/
#  [dirOut]: output directory, default to current working directory
#  [dirIn1, dirIn2]: path to all libraries, default set to the same; default to /projects/epigenomics/ep50/internal/jqc.1.7.6/
#  [RPKM1, RPKM2]: gene RPKM files, default to <dirIn><lib>/coverage/<lib1>.G.A.rpkm.pc
#  [geneUP, geneDN]: gene DEfine results, default to <dirGene>UP/DN.<cell1>-<donor1>_<cell2>-<donor2>.FDR_<fdr>.rmin_<rmin>.Nmin_<Nmin>
#  [exonUP, exonDN]: exon DEfine results, default to <dirExon>UP/DN.<cell1>-<donor1>_<cell2>-<donor2>.FDR_<fdr>.rmin_<rmin>.Nmin_<Nmin>
# Required input: 
#  DEfine on genes and exons: run with same parameters; result filename format: "UP/DN.", cell1, "-", donor1, "_", cell2, "-", donor2, ".FDR_", fdr, ".rmin_", rmin, ".Nmin_", Nmin
#  gene RPKM file for both library: dirIn, lib, "/coverage/", lib, ".G.A.rpkm.pc"; file format: ENSG_id \t #reads \t gene_RPKM \t average_RPKM \t min_RPKM \t max_RPKM
# Output:
#  cassete exons: <dirOut><cell1>"-"<donor1>"_"<cell2>"-"<donor2>_isoform.txt
#  isoform genes: <dirOut><cell1>"-"<donor1>"_"<cell2>"-"<donor2>_isoform_gene.txt
#  return list: 
#   summary: No. of DE_genes, DE_exons, with_expressed_genes, isoform_exons, exclude_DE_genes, isoform_genes
#   isoform_exon: exons identified as isoforms 
#   isoform_gene: genes identified as isoforms

