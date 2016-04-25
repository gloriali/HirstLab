library(ggplot2)
library(ggbio)
library(GenomicRanges)
source("/projects/epigenomics/users/smcconnell/glioma/WGBS/comparisons/DMR/DMR.figures.R")
source("/home/lli/HirstLab/Pipeline/R/enrich_GREAT.R")

# Take input directory and filename pattern, produce dataframes of all that match.
# Assumes that the files will contain a DMR.[...]-[...].methyl pattern
# Be careful when entering patterns that you're only including files that you want

list.DMR <- function(pattern, path = "") {
  # Given a pattern describing 1 or more filenames, produce a list of bedgraphs
  
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

DMR.avg.len <- function(dmrs) {
  avg.len = function(dmrName) {
    mean(dmrs[[dmrName]][,7])
  }
  lapply(seq_along(dmrs), avg.len)
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

DMR.freq.plot <- function(y) {
  freq.plot <- function(x) {
    DMR_figures(DMR = y[[x]],
                sample1 = strsplit(names(y)[x], split = " ")[[1]][1],
                sample2 = strsplit(names(y)[x], split = " ")[[1]][3],
                figures = "frequency")
  }
  lapply(seq_along(y), freq.plot)
}

DMR.GREAT.plot <- function(path, filename) {
  
  great.list = list.files(path, filename)
  great.sub = function(files) {
    subbed = sub(x = files, pattern = "^GREAT.*_", replacement = "")
    subbed = sub(x = subbed, pattern = "\\.tsv$", replacement = "")
  }
  great.list = unique(unlist(lapply(great.list, great.sub)))
  
  enrich.simple = function(greatfile) {
    enrich_GREAT(file = greatfile, 
                 name = greatfile,
                 dirIn = paste0(getwd(), "/great/Hyper-Hypo_Split/"),
                 dirOut = paste0(getwd(), "/great/Hyper-Hypo_Split/"))
  }
  lapply(great.list, enrich.simple)
}

DMR.stats <- function(hypo, hyper) {
  # Given two equal length lists of hypo and hyper DMRs, produce dataframes of their summary statistics
  
  if (length(hypo) != length(hyper)) {
    return("Error: hypo and hyper are not of equal length")
  }
  if (!all(names(hyperDMRs) == names(hypoDMRs))) {
    return("Error: Elements of the same index in both hypo and hyper must have the same sample name")
  }
  findStats <- function(num) {
    hypoDMR  = hypo[[num]]
    hyperDMR = hyper[[num]]
    lengths = c(length(hypoDMR[,7]), length(hyperDMR[,7]))
    sums = c(sum(hypoDMR[,7]), sum(hyperDMR[,7]))
    df = data.frame(`Number of UMRs` = lengths, `Avg UMR Length` = sums)
    names(df) <- gsub(x = names(df),
                      pattern = "\\.",
                      replacement = " ")
    row.names(df) <- c("Hypo", "Hyper")
    df
  }
  result <- lapply(seq_along(hypo), findStats)
  names(result) <- names(hypo)
  result
}


