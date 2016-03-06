
# produce uniquely named data.frames for each input

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

list.DMR <- function(pattern, path = "") {
  
  files = list.files("DMR/")[grep(pattern, list.files("DMR/"))]
  tables = list()
  
  for(i in files) {
    filename = sub("DMR\\.(.*)-(.*)\\.methyl.*", "\\1_\\2", i)
    assign(filename, read.delim(paste0(getwd(), "/", path, i), header = F))
    tables[[filename]] = get(filename)
  }
  tables
}

# Take input directory and filename pattern, produce dataframes of chr's and average DMR lengths

DMR.lengths <- function(pattern, path = "") {
  lapply(list.DMR(pattern, path), DMR_averages)
}

