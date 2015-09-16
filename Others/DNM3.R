library(dplyr)
library(reshape2)
library(ggplot2)
setwd("~/Desktop/DNM3")
dnm3 <- data.frame()
for(file in list.files(pattern = "ENST.*txt")){
	id=gsub("\\.txt", "", file)
	print(id)
	enst <- read.delim(file, head = T) %>% filter(grepl("ENSE", Exon...Intron)) %>% 
		mutate(Start = gsub(",", "", Start), End = gsub(",", "", End), id = id, chr = "chr1") %>% 
		select(chr, Start, End, Exon...Intron, id) 
	dnm3 <- rbind(dnm3, enst)
}
dnm3 <- dnm3 %>% arrange(Start)
write.table(dnm3, file = "DNM3_ENST.bed", sep = "\t", quote = F, col.names = F, row.names = F)
dnm3_rpkm <- read.delim("DNM3.txt", head = T, as.is = T) %>% select(Transcript.ID, RPKM001, RPKM002)
dnm3_rpkm_tall <- melt(dnm3_rpkm, id = "Transcript.ID")
(dnm3_rpkm_figure <- ggplot(dnm3_rpkm_tall, aes(x = Transcript.ID, y = value, fill = variable)) + 
 	geom_bar(stat = "identity", width = 0.5, position = position_dodge()) + 
 	scale_fill_manual(name = "Sample", labels = c("001", "002"), values = c("RPKM001" = "red", "RPKM002" = "blue")) + 
 	coord_flip() + 
 	xlab("") + 
 	ylab("RPKM") + 
 	theme_bw())
ggsave(dnm3_rpkm_figure, file = "DNM3_RPKM_figure.pdf")
