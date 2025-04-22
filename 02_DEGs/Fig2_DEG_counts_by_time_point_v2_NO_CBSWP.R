#plot summary of differentially expressed genes for Daphnia RNAseq

library(ggplot2)
library(viridis)
library(cowplot)
library(tidyverse)
library(forcats)

######################

gctemp <- data.frame(
  count = c(80,167,226,517),
  day = c(96,96,168,168),
  log = c('z_down', 'up', 'z_down', 'up')
)

gctemp$log <- ordered(gctemp$log)

nice_layout <- theme_cowplot()+
  panel_border(color = "grey85", size = 1, linetype = 1,
               remove = FALSE, "black") 


p2 <- ggplot(gctemp, aes(x=as.factor(day), y=count, fill = log)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("steelblue3","midnightblue")) +
  #ggtitle(label = "DEG count (genotype_clone + temp) by sampling time") +
  xlab("sampling hour") +
  ylab("DEG count") +
  nice_layout

p2

wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/Fig2_build/")
#ggsave(paste(wdII, "DEG_count_genoclone_temp_v2.png", sep =""), plot = p2, dpi = 1000, height = 6, width = 6)
#ggsave(paste(wdII, "DEG_count_genoclone_temp_v2.svg", sep =""), plot = p2, dpi = 1000, height = 6, width = 6)
#WORKING PLOT SAVED BELOW
ggsave(paste(wdII, "DEG_count_NO_CBSWP_genoclone_temp_v3.svg", sep =""), plot = p2, dpi = 1000, height = 4, width = 6)


#Add in barplot of common gene counts (Fig2A inset)

shared <- data.frame(
  count = c(41,32,1),
  log2FC = c("up-regulated", "down-regulated", "n=1")
)


p3 <- ggplot(shared, aes(x=forcats::fct_reorder(log2FC, desc(count)), y=count, fill = log2FC)) +
  geom_col() +
  scale_fill_manual(values = c("lightgreen", "darkgreen", "darkseagreen")) +
  #ggtitle(label = "DEG count (genotype_clone + temp) by sampling time") +
  xlab("") +
  #ylab("DEG count") +
  nice_layout +
  theme(legend.position = "none")


p3

wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/Fig2_build/")
ggsave(paste(wdII, "shared_NO_CBSWP_logFC_count_v2.svg", sep =""), plot = p3, dpi = 1000, height = 4, width = 5)

