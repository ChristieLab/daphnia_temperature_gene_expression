library(ggplot2)
library(tidyverse)
library(cowplot)
library(viridis)


wd <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/updated_figures_02-23-25/read_stats/")
setwd(wd)

reads <- read_csv("./read_stats_by_sample.csv")
head(reads)
samples <- read_csv("sample.info.all.modified.csv")
head(samples)

cat <- left_join(samples, reads, by = c("rename_id" = "sample"))

#calculate average numbers of reads
seq_avg <- mean(cat$total_seqs)
#calculated average length across all samples
len_avg <- mean(cat$total_len)

#daphnia pulex genome size
cat$dpulex <- 133196385 
#calculate coverage for all samples
cat <- mutate(cat, cov = total_len/dpulex)

len_avg


nice_layout <- theme_cowplot()+
  theme(plot.background = element_rect("white"))



p <- ggplot(cat, aes(rename_id, total_len))+
  geom_col(aes(fill=genotype_clone)) +
  geom_hline(yintercept = len_avg, linetype='dashed') +
  xlab("Sample")+
  ylab("Total Sequence Length") +
  scale_fill_viridis(discrete = TRUE) +
  nice_layout +
  ggtitle("Sequenced Bases per Sample") +
  theme(axis.text.x.bottom = element_blank())+ 
  guides(fill=guide_legend(title = "Clone"))
p

wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/MS_v1_figures_tables_submit/Sup_Figs/")
ggsave(paste(wdII, "Sup_Fig3-total_sequence_length.svg"), plot = p, dpi = 1000, height = 6, width = 10)
ggsave(paste(wdII, "Sup_Fig3-total_sequence_length.png"), plot = p, dpi = 1000, height = 6, width = 10)


