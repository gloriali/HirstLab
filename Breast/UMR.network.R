# TFBS overlapping UMR with DE genes network

RPKM <- read.csv("~/快盘/Publications/breast/revision/sup/TableS12.csv", head = T, as.is = T, row.names = 1)
UP <- read.csv("~/快盘/Publications/breast/revision/sup/UP.lum_myo.csv", head = T, as.is = T, row.names = NULL)
DN <- read.csv("~/快盘/Publications/breast/revision/sup/DN.lum_myo.csv", head = T, as.is = T, row.names = NULL)
DE <- rbind(UP, DN)
setwd("~/快盘/REMC/UMR")
lum.UMR.genes <- read.delim("lum.UMR.genes", head = F, as.is = T)
myo.UMR.genes <- read.delim("myo.UMR.genes", head = F, as.is = T)
lum.UMR.TFs <- read.delim("DMRs.-1_TFs", head = F, as.is = T)
myo.UMR.TFs <- read.delim("DMRs.1_TFs", head = F, as.is = T)

# lum UMR TF activators 
# GATA3 downstream DE genes 
GATA3 <- read.delim("GATA3.lum.UMR", head = F, as.is = T)
GATA3_DE <- intersect(DE$Ensembl.ID, lum.UMR.genes[lum.UMR.genes$V1 %in% GATA3$V1, "V2"])
GATA3_DE <- rbind(RPKM[RPKM$HUGO %in% "GATA3", ], RPKM[GATA3_DE, ])
GATA3_DE$lum_myo <- apply(GATA3_DE[, grep("lum", colnames(GATA3_DE), value = T)], 1, mean) / apply(GATA3_DE[, grep("myo", colnames(GATA3_DE), value = T)], 1, mean)
GATA3_DE <- data.frame(Source = "GATA3", Interaction = "Activate", Ensembl = row.names(GATA3_DE), HUGO = GATA3_DE$HUGO, log2FC = log2(GATA3_DE$lum_myo))
# FOXA1 downstream DE genes in lum 
FOXA1 <- read.delim("FOXA1.lum.UMR", head = F, as.is = T)
FOXA1_DE <- intersect(DE$Ensembl.ID, lum.UMR.genes[lum.UMR.genes$V1 %in% FOXA1$V1, "V2"])
FOXA1_DE <- rbind(RPKM[RPKM$HUGO %in% "FOXA1", ], RPKM[FOXA1_DE, ])
FOXA1_DE$lum_myo <- apply(FOXA1_DE[, grep("lum", colnames(FOXA1_DE), value = T)], 1, mean) / apply(FOXA1_DE[, grep("myo", colnames(FOXA1_DE), value = T)], 1, mean)
FOXA1_DE <- data.frame(Source = "FOXA1", Interaction = "Activate", Ensembl = row.names(FOXA1_DE), HUGO = FOXA1_DE$HUGO, log2FC = log2(FOXA1_DE$lum_myo))
# ZNF217 downstream DE genes in lum 
ZNF217 <- read.delim("ZNF217.lum.UMR", head = F, as.is = T)
ZNF217_DE <- intersect(DE$Ensembl.ID, lum.UMR.genes[lum.UMR.genes$V1 %in% ZNF217$V1, "V2"])
ZNF217_DE <- rbind(RPKM[RPKM$HUGO %in% "ZNF217", ], RPKM[ZNF217_DE, ])
ZNF217_DE$lum_myo <- apply(ZNF217_DE[, grep("lum", colnames(ZNF217_DE), value = T)], 1, mean) / apply(ZNF217_DE[, grep("myo", colnames(ZNF217_DE), value = T)], 1, mean)
ZNF217_DE <- data.frame(Source = "ZNF217", Interaction = "Activate", Ensembl = row.names(ZNF217_DE), HUGO = ZNF217_DE$HUGO, log2FC = log2(ZNF217_DE$lum_myo))
# EGR1 downstream DE genes 
EGR1 <- read.delim("EGR1.myo.UMR", head = F, as.is = T)
EGR1_DE <- intersect(DE$Ensembl.ID, myo.UMR.genes[myo.UMR.genes$V1 %in% EGR1$V1, "V2"])
EGR1_DE <- rbind(RPKM[RPKM$HUGO %in% "EGR1", ], RPKM[EGR1_DE, ])
EGR1_DE$lum_myo <- apply(EGR1_DE[, grep("lum", colnames(EGR1_DE), value = T)], 1, mean) / apply(EGR1_DE[, grep("myo", colnames(EGR1_DE), value = T)], 1, mean)
EGR1_DE <- data.frame(Source = "EGR1", Interaction = "Activate", Ensembl = row.names(EGR1_DE), HUGO = EGR1_DE$HUGO, log2FC = log2(EGR1_DE$lum_myo))
# GATA3 + FOXA1 + ZNF217: top 3 lum UMR TF activators 
network <- rbind(GATA3_DE[2:nrow(GATA3_DE),], FOXA1_DE[2:nrow(FOXA1_DE),], ZNF217_DE[2:nrow(ZNF217_DE),], EGR1_DE[2:nrow(EGR1_DE),])
network_node <- rbind(GATA3_DE[,c("HUGO", "log2FC")], FOXA1_DE[,c("HUGO", "log2FC")], ZNF217_DE[,c("HUGO", "log2FC")], EGR1_DE[,c("HUGO", "log2FC")])
network_node$UP_lum <- network_node$log2FC > 0
write.table(network, file = "network.txt", sep = "\t", quote = F, row.names = F)
write.table(network_node, file = "network_node.txt", sep = "\t", quote = F, row.names = F)

# lum UMR DN genes TF
lum.UMR.DN.TF <- unique(unlist(strsplit(lum.UMR.TFs[lum.UMR.TFs$V1 %in% lum.UMR.genes[lum.UMR.genes$V2 %in% DN$Ensembl.ID, "V1"], "V2"], ",")))
# myo UMR UP genes TF
myo.UMR.UP.TF <- unique(unlist(strsplit(myo.UMR.TFs[myo.UMR.TFs$V1 %in% myo.UMR.genes[myo.UMR.genes$V2 %in% UP$Ensembl.ID, "V1"], "V2"], ",")))


