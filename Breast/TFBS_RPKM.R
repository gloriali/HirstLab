setwd("~/快盘/Publications/breast/revision/sup/")
TF1 <- read.csv("TableS11.csv", head = T, as.is = T)
e <- 1e-6
TF <- data.frame(Donor = rep(c("RM035", "RM080", "RM084"), each = nrow(TF1)), 
                 HUGO = rep(TF1$HUGO, 3), 
                 FC_BS = rep(TF1$Luminal.Myoepithelial, 3), 
                 FC_RPKM = c((TF1$lum.RM035 + e)/(TF1$myo.RM035 + e), (TF1$lum.RM080 + e)/(TF1$myo.RM080 + e), (TF1$lum.RM084 + e)/(TF1$myo.RM084 + e)))
library(ggplot2)
(cor_plot <- ggplot(data = TF, aes(x = log2(FC_BS), y = log2(FC_RPKM), color = Donor)) + 
   geom_point() + 
   geom_smooth(method = lm, se = F) +
   xlab("log2 TFBS FC lum/myo") + 
   ylab("log2 RPKM FC lum/myo") + 
   theme_bw())
ggsave(cor_plot, file = "TFBS_RPKM.pdf")

