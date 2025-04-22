#plot summary of shared deferentially expressed genes for Daphnia RNAseq

library(ggplot2)
library(viridis)
library(cowplot)

wd <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/")
setwd(wd)
# DEGs data
# first column 'X' is the gene name
d4 <- read.csv("./96_results/p_05_degs_96hrs_genotype_clone-temp.csv")
d7 <- read.csv("./168_results/p_05_degs_168hrs_genotype_clone-temp.csv")

# find shared genes
shared_genes <- intersect(d4$X, d7$X)
length(shared_genes)

# Make dataframes of shared genes
d4.1 <- d4[d4$X %in% shared_genes, ]
d7.1 <- d7[d7$X %in% shared_genes, ]

combined_data <- rbind(
  cbind(d4.1, day = "Day 4"),
  cbind(d7.1, day = "Day 7")
)

# Convert Day column to factor for correct ordering on x-axis
combined_data$day <- factor(combined_data$day, levels = c("Day 4", "Day 7"))
combined_data$sign <- factor(sign(combined_data$log2FoldChange), labels = c("negative", "positive"))

nice_layout <- theme_cowplot()+
  panel_border(color = "grey85", size = 1, linetype = 1,
               remove = FALSE, "black") 

p1 <- ggplot(combined_data, aes(x = day, y = log2FoldChange, color = day)) +
  geom_point() +
  labs(x = "Day", y = "log2FoldChange", 
       title = "Log2FoldChange for Shared Genes between Day 4 and Day 7") +
  nice_layout
p1


# all shared genes
p2 <- ggplot(combined_data, aes(x = day, y = log2FoldChange, color = sign)) +
  geom_point(size = 2) +
  geom_line(aes(group = X), linewidth = 0.5, color = "lightgrey") +
  scale_color_manual(values = c("#f8766d", "#03bfc4"), labels = c("Negative", "Positive")) +
  labs(x = "Day", y = "log2FoldChange", title = "25 Shared DEGS")+
  theme_classic() +
  theme(legend.position = "none")
  #scale_y_continuous(limits = c(-30, 30), breaks = seq(-30, 30, by = 10))

p2

#filter for log2FC greater than 0
upreg <- combined_data %>%
  group_by(X) %>%
  filter(log2FoldChange > 0)

p3 <- ggplot(upreg, aes(x = day, y = log2FoldChange, color = day)) +
  geom_point(size = 2) +
  geom_line(aes(group = X), linewidth = 0.5, color = "lightgrey") +
  scale_color_manual(values = c("#03bfc4", "#f8766d"), labels = c("Negative", "Positive")) +
  labs(x = "Day", y = "log2FoldChange") +
  nice_layout +
  theme(legend.position = "none") 
  #scale_y_continuous(limits = c(-30, 30), breaks = seq(-30, 30, by = 10))

p3
wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/Fig2_build/")
ggsave(paste(wdII, "shared_DEG_log2FC_compare_genoclone_temp_v1.svg", sep =""), plot = p3, dpi = 1000, height = 4, width = 6)

# Pull out DEGS from day 4, that are NOT on day 7 <-- continue from here


