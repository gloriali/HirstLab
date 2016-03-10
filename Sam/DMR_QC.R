library(ggplot2)
library(ggbio)
library(GenomicRanges)

# Take input directory and filename pattern, produce dataframes of all that match.
# Assumes that the files will contain a DMR.[...]-[...].methyl pattern
# Be careful when entering patterns that you're only including files that you want

list.DMR <- function(pattern, path = "") {
  
  files = list.files("DMR/")[grep(pattern, list.files("DMR/"))]
  tables = list()
  
  for(i in files) {
    filename = sub("DMR\\.(.*)-(.*)\\.methyl.*", "\\1 - \\2", i)
    assign(filename, read.delim(paste0(getwd(), "/", path, i), header = F))
    tables[[filename]] = get(filename)
    names(tables[[filename]])[1:7] <- c("Chromosome", "Start", "End", "ID", "Hypo/Hyper", "No. of CpGs per DMR", "DMR Length (bp)")
    
    # Assumes chromosomes are misordered from bash (eg. chr1, chr10, chr11, ..., chr19, chr2, chr 20, etc.) and reorders them
    tables[[filename]]$Chromosome = factor(tables[[filename]]$Chromosome, levels = levels(factor(tables[[filename]]$Chromosome))[c(1, 12, 16:22, 2:11, 13:15, 23, 24)])
  }
  tables
}


DMR.len.plot <- function(y) {
  len.plot = function(x) {
    titlename = paste0(names(y)[x], "\nDMR Lengths by Chromosome")
    ggplot(data = y[[x]], aes(x = `Chromosome`, y = `DMR Length (bp)`, fill = `Chromosome`)) +
      geom_boxplot() + 
      coord_cartesian(ylim = c(0,1500)) + 
      ggtitle(titlename)
  }
  lapply(seq_along(y), len.plot)
}

DMR.cpg.plot <- function(y) {
  cpg.plot = function(x) {
    titlename = paste0(names(y)[x], "\nNo. of CpGs per DMR by Chromosome")
    ggplot(data = y[[x]], aes(x = `Chromosome`, y = `No. of CpGs per DMR`, fill = `Chromosome`)) +
      geom_boxplot() +
      coord_cartesian(ylim = c(4,15)) + 
      ggtitle(titlename)
  }
  lapply(seq_along(y), cpg.plot)
}

data(hg19IdeogramCyto, package = "biovizBase")   # Mixed chr order
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))
DMR_pos_WGBS_gr <- keepSeqlevels(as(DMRs[[1]], "GRanges"), paste0("chr", c(1:22, "X", "Y")))
ggplot(hg19) + 
  layout_karyogram(cytoband = TRUE) +
#  layout_karyogram(data = DMR_pos_WGBS_gr, geom = "point", aes(y = y, x = pos, color = factor(`Hypo/Hyper`, levels = c("1", "-1"))), position = position_jitter(height = 2), size = 1.5, alpha = 0.5) +
  xlab("Position of DMRs on the chromosome") + 
  coord_cartesian(ylim = c(-15, 25)) #+
#  facet_grid(seqnames ~ donor)   Need to read up on facets

  
