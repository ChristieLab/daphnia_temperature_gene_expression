library(EnhancedVolcano)
library(DESeq2)
library(cowplot)


wd <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/assembled_transcripts_final/")
setwd(wd)

#Read in data
############## 96 ##################
#Read in Deseq matrix
dds <- readRDS("./DESeq_96_NO_CBSWP.genotype_clone.temp.rds")

#calcluate significant results
res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
summary(res)

length(which(res[, 6] < 0.05))
summary(res)

min(res$log2FoldChange, na.rm = TRUE)
max(res$log2FoldChange, na.rm = TRUE)


############ 168 #################
#Read in Deseq matrix
dds2 <- readRDS("./DESeq_168_NO_CBSWP.genotype_clone.temp.rds")

#calcluate significant results
res2 <- results(dds2, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
summary(res2)

length(which(res2[, 6] < 0.05))
summary(res2)

min(res2$log2FoldChange, na.rm = TRUE)
max(res2$log2FoldChange, na.rm = TRUE)


####################################
#update for review with new color scheme.
#make new data frame from res results
# Convert DESeqResults object to a data frame for easier manipulation with ggplot
# Create the volcano plot for 96
res_df <- as.data.frame(res)

# Add a new column to categorize genes
res_df$diff_expressed <- "Not significant"
res_df$diff_expressed[res_df$diff_expressed == "Not significant" & res_df$log2FoldChange > 0 & res_df$padj < 0.05] <- "Up-regulated"
res_df$diff_expressed[res_df$diff_expressed == "Not significant" & res_df$log2FoldChange < 0 & res_df$padj < 0.05] <- "Down-regulated"


# Create breaks for x-axis
break_start <- 10
break_end <- 19
break_length <- break_end - break_start

# Create a new, adjusted x-axis column
res_df$log2FoldChange_adj <- res_df$log2FoldChange
res_df$log2FoldChange_adj[res_df$log2FoldChange > break_end] <- res_df$log2FoldChange_adj[res_df$log2FoldChange > break_end] - break_length

# First plot: data from the left side of the break
p_left <- ggplot(subset(res_df, log2FoldChange < break_start),
                 aes(x = log2FoldChange_adj, y = -log10(padj), color = diff_expressed)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "steelblue3", "Down-regulated" = "midnightblue", "Not significant" = "gray19")) +
  labs(y = "-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.position = "none",
    axis.title.x = element_blank(),
  ) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10), limits = c(-12, break_start))

# Create a dummy data frame for the right plot
dummy_df <- data.frame(
  log2FoldChange_adj = c(break_end),
  padj = c(0.5), 
  diff_expressed = "Not significant"
)

# Second plot: data from the right side of the break
p_right <- ggplot(dummy_df, aes(x = log2FoldChange_adj, y = -log10(padj), color = diff_expressed)) +
  geom_point(alpha = 0.0, size = 0) + # Use a transparent point
  scale_color_manual(values = c("Up-regulated" = "steelblue3", "Down-regulated" = "midnightblue", "Not significant" = "gray19")) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color = "black"),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_x_continuous(
    breaks = c(20),
    labels = c("20"),
    limits = c(break_start, 20) # Set a fixed limit slightly above the break
  )

# Combine the plots
combined_plot <- plot_grid(p_left, p_right, align = "h", rel_widths = c(1, 0.2))

# Add the final title and break symbol
v_96 <- ggdraw(combined_plot) +
  draw_label("96 hours", fontface = 'bold', x = 0.5, y = 0.95) +
  draw_label("/", x = 0.83, y = 0.065, size = 12) +
  draw_label("Log2 Fold Change", x = 0.5, y = 0.02) # Add a single x-axis label for the combined plot

v_96

wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/Fig3_MS/")
ggsave(paste(wdII, "96_NO_CBSWP_p0.05_log2FC5_volcano_panel_v2.svg", sep = ""), v_96, dpi = 1000, height = 4, width = 8.5)

#######now for the 168-hour plot
res2_df <- as.data.frame(res2)

# add a new column to categorize genes
res2_df$diff_expressed <- "Not significant"
res2_df$diff_expressed[res2_df$log2FoldChange > 0 & res2_df$padj < 0.05] <- "Up-regulated"
res2_df$diff_expressed[res2_df$log2FoldChange < 0 & res2_df$padj < 0.05] <- "Down-regulated"


# Create a new, adjusted x-axis column
res2_df$log2FoldChange_adj <- res2_df$log2FoldChange
res2_df$log2FoldChange_adj[res2_df$log2FoldChange > break_end] <- res2_df$log2FoldChange_adj[res2_df$log2FoldChange > break_end] - break_length

# First plot: data from the left side of the break
p_left2 <- ggplot(subset(res2_df, log2FoldChange < break_start),
                 aes(x = log2FoldChange_adj, y = -log10(padj), color = diff_expressed)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "steelblue3", "Down-regulated" = "midnightblue", "Not significant" = "gray19")) +
  labs(y = "-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.position = "none",
    axis.title.x = element_blank(),
  ) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10), limits = c(-12, break_start))

# Second plot: data from the right side of the break
p_right2 <- ggplot(subset(res2_df, log2FoldChange >= break_end),
                  aes(x = log2FoldChange_adj, y = -log10(padj), color = diff_expressed)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "steelblue3", "Down-regulated" = "midnightblue", "Not significant" = "gray19")) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color = "black"),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  ) +
  scale_x_continuous(
    breaks = c(20),
    labels = c("20"),
    limits = c(break_start, max(res2_df$log2FoldChange_adj, na.rm = TRUE))
  )

# Combine the plots
combined_plot2 <- plot_grid(p_left2, p_right2, align = "h", rel_widths = c(1, 0.2))

# Add the final title and break symbol
v_168 <- ggdraw(combined_plot2) +
  draw_label("/", x = 0.83, y = 0.065, size = 12) 

v_168

wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/Fig3_MS/")
ggsave(paste(wdII, "168_NO_CBSWP_p0.05_log2FC5_volcano_panel_v2.svg", sep = ""), v_168, dpi = 1000, height = 4, width = 8.5)
