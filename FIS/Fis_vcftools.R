#script for plotting He/Ho in Daphnia SNPs called from RNAseq data

library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(ggsignif)
library(viridis)

wd <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/updated_figures_02-23-25/FIS/")
setwd(wd)
list.files()

#read in vcftools --het output
het <- read.table("./daphnia_FIS.het ", header = TRUE)

#read in sample data
samples <- read.csv("../mito_PCA/sample.info.all.modified.csv")

#prune samples
samples <- samples[samples$rename_id != "s180", ]

#merge
het2 <- full_join(het, samples, by = c("INDV" = "rename_id"))

#plot
nice_layout <- theme_cowplot()+
  theme(plot.background = element_rect("white")) 

colors <- viridis(6, begin = 0.5)

p <- ggplot(het2, aes(genotype_clone, F, fill = genotype_clone, color = genotype_clone)) +
  geom_violin(scale = "width") +
  geom_point(aes(color = genotype_clone), size = 1, alpha = 0.9, position=position_jitter(width=0.1, height=0.1)) +
  xlab("clone")+
  ylab("FIS") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = c("#000000", "#000000", "#000000", "#000000", "#000000", "#000000")) + 
  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  nice_layout +
  theme(legend.position="none") 
  
p

wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/MS_v1_figures_tables_submit/Sup_Figs/")
ggsave(paste(wdII,"Sup_Fig4_Fis_plot_by_clone_v3.svg", sep = ""), p, dpi = 1000, height = 8, width = 6)

