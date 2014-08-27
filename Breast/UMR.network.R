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
# write.table(GATA3_UP[2:nrow(GATA3_UP),], file = "GATA3_UP.txt", sep = "\t", quote = F, row.names = F)
# write.table(GATA3_UP[,c("HUGO", "log2FC")], file = "GATA3_UP_node.txt", sep = "\t", quote = F, row.names = F)

# FOXA1 downstream DE genes in lum 
FOXA1 <- read.delim("FOXA1.lum.UMR", head = F, as.is = T)
FOXA1_UP <- intersect(UP$Ensembl.ID, lum.UMR.genes[lum.UMR.genes$V1 %in% FOXA1$V1, "V2"])
FOXA1_UP <- rbind(RPKM[RPKM$HUGO %in% "FOXA1", ], RPKM[FOXA1_UP, ])
FOXA1_UP$lum_myo <- apply(FOXA1_UP[, grep("lum", colnames(FOXA1_UP), value = T)], 1, mean) / apply(FOXA1_UP[, grep("myo", colnames(FOXA1_UP), value = T)], 1, mean)
FOXA1_UP <- data.frame(Source = "FOXA1", Interaction = "Activate", Ensembl = row.names(FOXA1_UP), HUGO = FOXA1_UP$HUGO, log2FC = log2(FOXA1_UP$lum_myo))

# ZNF217 downstream DE genes in lum 
ZNF217 <- read.delim("ZNF217.lum.UMR", head = F, as.is = T)
ZNF217_UP <- intersect(UP$Ensembl.ID, lum.UMR.genes[lum.UMR.genes$V1 %in% ZNF217$V1, "V2"])
ZNF217_UP <- rbind(RPKM[RPKM$HUGO %in% "ZNF217", ], RPKM[ZNF217_UP, ])
ZNF217_UP$lum_myo <- apply(ZNF217_UP[, grep("lum", colnames(ZNF217_UP), value = T)], 1, mean) / apply(ZNF217_UP[, grep("myo", colnames(ZNF217_UP), value = T)], 1, mean)
ZNF217_UP <- data.frame(Source = "ZNF217", Interaction = "Activate", Ensembl = row.names(ZNF217_UP), HUGO = ZNF217_UP$HUGO, log2FC = log2(ZNF217_UP$lum_myo))

network <- rbind(GATA3_UP[2:nrow(GATA3_UP),], FOXA1_UP[2:nrow(FOXA1_UP),], ZNF217_UP[2:nrow(ZNF217_UP),])
network_node <- rbind(GATA3_UP[,c("HUGO", "log2FC")], FOXA1_UP[,c("HUGO", "log2FC")], ZNF217_UP[,c("HUGO", "log2FC")])
write.table(network, file = "network.txt", sep = "\t", quote = F, row.names = F)
write.table(network_node, file = "network_node.txt", sep = "\t", quote = F, row.names = F)


