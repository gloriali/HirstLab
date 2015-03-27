setwd("~/快盘/IDH")
library(ggplot2)
library(dplyr)

# ============= IDH prevalence =========
prevalence <- read.delim("IDH_prevalence.txt", as.is = T) %>% mutate(Subtype = factor(Subtype, levels = rev(Subtype)), Type = factor(Type, levels = c("Glioma", "Brain other", "Leukemia", "Blood other")))
(prevalence_figure <- ggplot(prevalence, aes(x = Subtype, y = Percent, fill = Type)) + 
   geom_bar(stat = "identity", width = .5) + 
   geom_text(aes(label = Total), hjust = -.5, size = 6) + 
   coord_flip(ylim = c(0, 100)) + 
   xlab("") + 
   ylab("Percentage") + 
   theme_bw() + 
   theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 0)))
ggsave(prevalence_figure, file = "prevalence_figure.pdf")
