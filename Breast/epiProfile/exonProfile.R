# Epigenetic profile at exon boundaries
setwd("~/快盘/REMC/epiProfile/")
# get lum-specific, myo-specific, expressed in both cell types exons  
isoform <- read.delim("lum084_myo084_isoform_all.txt", head = T, as.is = T)
isoform$coord <- gsub("<[-]*1_ENSG[0-9]+", "", isoform$V1)
isoform$coord  <- gsub("chr", "", isoform$coord)
isoform$coord  <- gsub(":", "_", isoform$coord)
isoform$coord  <- gsub("-", "_", isoform$coord)
isoform$coord  <- paste0(isoform$id, "_", isoform$coord)
lum_specific <- isoform[isoform$V2 < isoform$V3, "coord"]
myo_specific <- isoform[isoform$V2 > isoform$V3, "coord"]
lum084 <- read.delim("A17918.G.exn.A.rpkm", head = F, as.is = T)
lum084$coord <- gsub("<[-]*1", "", lum084$V1)
lum084$coord  <- gsub("chr", "", lum084$coord)
lum084$coord  <- gsub(":", "_", lum084$coord)
lum084$coord  <- gsub("-", "_", lum084$coord)
lum084$coord  <- paste0(lum084$V2, "_", lum084$coord)
myo084 <- read.delim("A17919.G.exn.A.rpkm", head = F, as.is = T)
myo084$coord <- gsub("<[-]*1", "", myo084$V1)
myo084$coord  <- gsub("chr", "", myo084$coord)
myo084$coord  <- gsub(":", "_", myo084$coord)
myo084$coord  <- gsub("-", "_", myo084$coord)
myo084$coord  <- paste0(myo084$V2, "_", myo084$coord)
both <- intersect(lum084[lum084$V4 > 0.1, "coord"], myo084[myo084$V4 > 0.1, "coord"])
neither <- setdiff(lum084$coord, c(lum_specific, myo_specific, both))
rm(lum084, myo084, isoform)

# DNA methylation profile 
lum_bismark_3p <- read.table("./exons3p_200/lumRM066_bismark.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_bismark_3p <- read.table("./exons3p_200/myoRM045_bismark.hg19v65_exons_for_genes.3prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
lum_bismark_5p <- read.table("./exons5p_200/lumRM066_bismark.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
myo_bismark_5p <- read.table("./exons5p_200/myoRM045_bismark.hg19v65_exons_for_genes.5prime_200.unique.profile", sep = " ", head = F, as.is = T, row.names = 1)
methyl_3p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), methyl = -1)
methyl_3p$group <- interaction(methyl_3p$Cell_type, methyl_3p$Expression)
methyl_3p[methyl_3p$group == "lum.lum-specific", "methyl"] <- colMeans(lum_bismark_3p[lum_specific,], na.rm = T)
methyl_3p[methyl_3p$group == "lum.myo-specific", "methyl"] <- colMeans(lum_bismark_3p[myo_specific,], na.rm = T)
methyl_3p[methyl_3p$group == "lum.expressed_in_both", "methyl"] <- colMeans(lum_bismark_3p[both,], na.rm = T)
methyl_3p[methyl_3p$group == "lum.not_expressed", "methyl"] <- colMeans(lum_bismark_3p[neither,], na.rm = T)
methyl_3p[methyl_3p$group == "myo.lum-specific", "methyl"] <- colMeans(myo_bismark_3p[lum_specific,], na.rm = T)
methyl_3p[methyl_3p$group == "myo.myo-specific", "methyl"] <- colMeans(myo_bismark_3p[myo_specific,], na.rm = T)
methyl_3p[methyl_3p$group == "myo.expressed_in_both", "methyl"] <- colMeans(myo_bismark_3p[both,], na.rm = T)
methyl_3p[methyl_3p$group == "myo.not_expressed", "methyl"] <- colMeans(myo_bismark_3p[neither,], na.rm = T)
methyl_5p <- data.frame(Cell_type = rep(c("lum", "myo"), each = 20*4), Expression = rep(c(rep(c("lum-specific", "myo-specific", "expressed_in_both", "not_expressed"), each = 20)), 2), Position = rep(seq(-190, 190, by = 20), times = 8), methyl = -1)
methyl_5p$group <- interaction(methyl_5p$Cell_type, methyl_5p$Expression)
methyl_5p[methyl_5p$group == "lum.lum-specific", "methyl"] <- colMeans(lum_bismark_5p[lum_specific,], na.rm = T)
methyl_5p[methyl_5p$group == "lum.myo-specific", "methyl"] <- colMeans(lum_bismark_5p[myo_specific,], na.rm = T)
methyl_5p[methyl_5p$group == "lum.expressed_in_both", "methyl"] <- colMeans(lum_bismark_5p[both,], na.rm = T)
methyl_5p[methyl_5p$group == "lum.not_expressed", "methyl"] <- colMeans(lum_bismark_5p[neither,], na.rm = T)
methyl_5p[methyl_5p$group == "myo.lum-specific", "methyl"] <- colMeans(myo_bismark_5p[lum_specific,], na.rm = T)
methyl_5p[methyl_5p$group == "myo.myo-specific", "methyl"] <- colMeans(myo_bismark_5p[myo_specific,], na.rm = T)
methyl_5p[methyl_5p$group == "myo.expressed_in_both", "methyl"] <- colMeans(myo_bismark_5p[both,], na.rm = T)
methyl_5p[methyl_5p$group == "myo.not_expressed", "methyl"] <- colMeans(myo_bismark_5p[neither,], na.rm = T)
save(methyl_3p, methyl_5p, file = "exonProfile_methyl.Rdata")
library(ggplot2)
(methyl_3p_profile <- ggplot(methyl_3p, aes(x = Position, y = methyl, group = group)) + 
   geom_line(aes(color = Expression)) + 
   geom_point(aes(color = Expression, shape = Cell_type)) + 
   ggtitle("DNA methylation profile around 3-prime exon boundaries") + 
   ylab("Average DNA methylation level") + 
   theme_bw())
ggsave(methyl_3p_profile, file = "methyl_3p_profile.pdf")
(methyl_5p_profile <- ggplot(methyl_5p, aes(x = Position, y = methyl, group = group)) + 
   geom_line(aes(color = Expression)) + 
   geom_point(aes(color = Expression, shape = Cell_type)) + 
   ggtitle("DNA methylation profile around 5-prime exon boundaries") + 
   ylab("Average DNA methylation level") + 
   theme_bw())
ggsave(methyl_5p_profile, file = "methyl_5p_profile.pdf")
