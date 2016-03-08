library(ggplot2)

# bed file of DMRs -> data frame containing chr's and their average DMR lengths

DMR_averages <- function(dmr) {
  dmrl = list()
  chr = unique(dmr$V1)
  for(i in chr)
  {
    dmrl = rbind(dmrl, c(i, mean(dmr$V7[dmr$V1 == i])))
  }
  dmrl = data.frame(dmrl)
  `colnames<-`(dmrl, c("chr", "avg"))
}

# Take input directory and filename pattern, produce dataframes of all that match.
# Assumes that the files will contain a DMR.[...]-[...].methyl pattern
# Be careful when entering patterns that you're only including files that you want

list.DMR <- function(pattern, path = "") {
  
  files = list.files("DMR/")[grep(pattern, list.files("DMR/"))]
  tables = list()
  
  for(i in files) {
    filename = sub("DMR\\.(.*)-(.*)\\.methyl.*", "\\1_\\2", i)
    assign(filename, read.delim(paste0(getwd(), "/", path, i), header = F))
    tables[[filename]] = get(filename)
    
    # Assumes chromosomes are misordered from bash (eg. chr1, chr10, chr11, ..., chr19, chr2, chr 20, etc.) and reorders them
    tables[[filename]]$V1 = factor(tables[[filename]]$V1, levels = levels(factor(tables[[filename]]$V1))[c(1, 12, 16:22, 2:11, 13:15, 23, 24)])
  }
  tables
}

# Take input directory and filename pattern, produce dataframes of chr's and average DMR lengths
# Assumes that the files will contain a DMR.[...]-[...].methyl pattern
# Be careful when entering patterns that you're only including files that you want

DMR.lengths <- function(pattern, path = "") {
  lapply(list.DMR(pattern, path), DMR_averages)
}

# Given a list of BED formatted DMRs, y: create a boxplot of average DMR length for each chromosome for every element of y
# ~~~~ Make axis labels next ~~~~~

DMR.len.plot <- function(y) {
  lapply(y, function(i) ggplot(data = i, aes(x = V1, y = V7, fill = V1)) + geom_boxplot() + coord_cartesian(ylim = c(0,1500)))
}


