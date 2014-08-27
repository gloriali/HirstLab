# TFBS overlapping UMR with DE genes network

RPKM <- read.csv("~/快盘/Publications/breast/revision/sup/TableS12.csv", head = T, as.is = T, row.names = 1)
UP <- read.csv("~/快盘/Publications/breast/revision/sup/UP.lum_myo.csv", head = T, as.is = T, row.names = NULL)
DN <- read.csv("~/快盘/Publications/breast/revision/sup/DN.lum_myo.csv", head = T, as.is = T, row.names = NULL)
setwd("~/快盘/REMC/UMR")
lum.UMR.genes <- read.delim("lum.UMR.genes", head = F, as.is = T)
myo.UMR.genes <- read.delim("myo.UMR.genes", head = F, as.is = T)

# GATA3 downstream DE genes in lum 
GATA3 <- read.delim("GATA3.lum.UMR", head = F, as.is = T)
GATA3_UP <- intersect(UP$Ensembl.ID, lum.UMR.genes[lum.UMR.genes$V1 %in% GATA3$V1, "V2"])
GATA3_UP <- rbind(RPKM[RPKM$HUGO %in% "GATA3", ], RPKM[GATA3_UP, ])
GATA3_UP$lum_myo <- apply(GATA3_UP[, grep("lum", colnames(GATA3_UP), value = T)], 1, mean) / apply(GATA3_UP[, grep("myo", colnames(GATA3_UP), value = T)], 1, mean)
GATA3_UP <- data.frame(Source = "GATA3", Interaction = "Activate", Ensembl = row.names(GATA3_UP), HUGO = GATA3_UP$HUGO, log2FC = log2(GATA3_UP$lum_myo))
write.table(GATA3_UP[2:nrow(GATA3_UP),], file = "GATA3_UP.txt", sep = "\t", quote = F, row.names = F)
GATA3_UP_node <- c("log2FC", paste(GATA3_UP$HUGO, "=", GATA3_UP$log2FC))
write.table(GATA3_UP[,c("HUGO", "log2FC")], file = "GATA3_UP_node.txt", sep = "\t", quote = F, row.names = F)


