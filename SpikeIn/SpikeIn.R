# Spike-in and RNA yield normalization for IL10 expression in naive, primary, memory, and secondary T-cells
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
setwd("/projects/edcc_prj2/RNAseq/PX0348/")

## Spike-in control: RPKM ~ concentration
spikein <- read.delim("./SpikeIn/SpikeIn.concentration.RPKM", head = F, as.is = T, col.names = c("Species", "Count", "Mix", "Concentration", "Length", "RPKM", "Index")) %>% 
	mutate(Index = revalue(Index, c(AACCCC = "Primary", ACCCAG = "Secondary", AGCGCT = "Memory", CAAAAG = "Naive")))
for(i in levels(factor(spikein$Index))){
	assign(paste0("lm_", i), lm(log2(Concentration) ~ log2(RPKM), spikein %>% filter(Index == i)))
}
(spikein_figure <- ggplot(spikein, aes(log2(RPKM), log2(Concentration), color = Index)) + 
	geom_point() + 
	geom_smooth(method = "lm", se=FALSE) + 
	scale_color_manual(name = "", values = c("red", "blue", "purple", "orange")) + 
	facet_wrap(~ Index, nrow = 2) + 
	theme_bw())
spikein_figure_gt <- ggplot_gtable(ggplot_build(spikein_figure))
pdf("./SpikeIn/SpikeIn.pdf")
grid.draw(spikein_figure_gt)
grid.text(substitute(italic(y) == a + b %.% italic(x), list(a = format(coef(lm_Memory)[1], digits = 2), b = format(coef(lm_Memory)[2], digits = 2))), x = unit(0.25, "npc"), y = unit(0.9, "npc"))
grid.text(substitute(italic(r)^2~"="~r2, list(r2 = format(summary(lm_Memory)$r.squared, digits = 3))), x = unit(0.25, "npc"), y = unit(0.88, "npc"))
grid.text(substitute(italic(y) == a + b %.% italic(x), list(a = format(coef(lm_Naive)[1], digits = 2), b = format(coef(lm_Naive)[2], digits = 2))), x = unit(0.55, "npc"), y = unit(0.9, "npc"))
grid.text(substitute(italic(r)^2~"="~r2, list(r2 = format(summary(lm_Naive)$r.squared, digits = 3))), x = unit(0.55, "npc"), y = unit(0.88, "npc"))
grid.text(substitute(italic(y) == a + b %.% italic(x), list(a = format(coef(lm_Primary)[1], digits = 2), b = format(coef(lm_Primary)[2], digits = 2))), x = unit(0.25, "npc"), y = unit(0.45, "npc"))
grid.text(substitute(italic(r)^2~"="~r2, list(r2 = format(summary(lm_Primary)$r.squared, digits = 3))), x = unit(0.25, "npc"), y = unit(0.43, "npc"))
grid.text(substitute(italic(y) == a + b %.% italic(x), list(a = format(coef(lm_Secondary)[1], digits = 2), b = format(coef(lm_Secondary)[2], digits = 2))), x = unit(0.55, "npc"), y = unit(0.45, "npc"))
grid.text(substitute(italic(r)^2~"="~r2, list(r2 = format(summary(lm_Secondary)$r.squared, digits = 3))), x = unit(0.55, "npc"), y = unit(0.43, "npc"))
dev.off()

## RNA yield
yield <- data.frame(Index = c("AACCCC", "ACCCAG", "AGCGCT", "CAAAAG"), 
										Sample = c("Primary", "Secondary", "Memory", "Naive"), 
										Yield = c(45, 15, 100, 60))
rownames(yield) <- yield$Index

## Noramlization
e <- 10^-6
for(i in yield$Index){
	assign(paste0("RPKM_", yield[i, "Sample"]), read.delim(paste0("./mm10v71/PX0348_", i, "/coverage/PX0348_", i, ".G.A.rpkm.pc"), head = F, as.is = T, col.names = c("ENSG", "count", "RPKM", "RPKM_mean", "RPKM_min", "RPKM_max")) %>% mutate(RPKM = RPKM + e))
	assign(paste0("Expression_", yield[i, "Sample"]), data.frame(ENSG = get(paste0("RPKM_", yield[i, "Sample"]))$ENSG, RPKM = get(paste0("RPKM_", yield[i, "Sample"]))$RPKM) %>% 
				 	mutate(Concentration = 2^(predict(get(paste0("lm_", yield[i, "Sample"])), get(paste0("RPKM_", yield[i, "Sample"])))), Expression = Concentration * yield[i, "Yield"]/15))
	write.table(get(paste0("Expression_", yield[i, "Sample"])), file = paste0("./mm10v71/PX0348_", i, ".expression.pc"), sep = "\t", row.names = F, quote = F)
}
expression_all <- data.frame(ENSG = Expression_Memory$ENSG, Naive = Expression_Naive$Expression, Primary = Expression_Primary$Expression, Memory = Expression_Memory$Expression, Secondary = Expression_Secondary$Expression)
rownames(expression_all) <- expression_all$ENSG
write.table(expression_all, file = "./mm10v71/Expression.normalized.matrix", sep = "\t", row.names = F, quote = F)

