#script for plotting pairwise W&C Fst and DEGs for pairs of clones
#Only controls

library(ggplot2)
library(viridis)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(tidyverse)


wd <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/updated_figures_02-23-25/pairwise_control_DEGs_Fst_corr/")
setwd(wd)
list.files()

d <- read_csv("./pairwise_fst_control_degs.csv")


#
corr <- cor.test(x=d$degs, y=d$weighted, method = 'pearson')
corr


r <- 0.009020235 * 0.009020235 

#
nice_layout <- theme_cowplot()+
  theme(plot.background = element_rect("white"))
#



p1 <- ggplot(d, aes(x = d$degs, y = d$weighted, label = pair, col = d$weighted)) +
  geom_smooth(method = "lm") +
  geom_point(size = 3) +
  nice_layout +
  xlab("pairwise DEGs") +
  ylab("weighted W&C Fst") +
  scale_color_viridis(option = "plasma", discrete = FALSE, alpha = 1) 



p1

p1a <- p1 + geom_text_repel(mapping = aes(x = d$degs, y = d$weighted, label = pair),inherit.aes = FALSE, max.overlaps = 100, box.padding = 0.5) 

p1a


wdII <-("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/MS_v1_figures_tables_submit/Sup_Figs/")
ggsave(paste(wdII,"Sup_Fig6_pairwise_wc_Fst-control_DEGs_corr.png"), plot = p1, dpi = 1000, height = 6, width = 10)
ggsave(paste(wdII,"Sup_Fig6_pairwise_wc_Fst-control_DEGs_corr.svg"), plot = p1, dpi = 1000, height = 6, width = 10)

###
ggsave(paste(wdII,"daphnia_pairwise_wc_Fst-control_DEGs_corr_ggrepel.png"), plot = p1a, dpi = 1000, height = 6, width = 10)
ggsave(paste(wdII,"daphnia_pairwise_wc_Fst-control_DEGs_corr_ggrepel.svg"), plot = p1a, dpi = 1000, height = 6, width = 10)
