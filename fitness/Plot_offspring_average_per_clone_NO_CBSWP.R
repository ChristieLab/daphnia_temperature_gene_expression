#=============================================================================================================#
# Script created by Nathan Backenstose, nbackens@purdue.edu
# Script created in version R 4.0.4 on 09/08/2024
# This script:Plots fitness data for Daphnia 
# Usage notes: 
#============================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
wd <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/updated_figures_02-23-25/")
setwd(wd)

library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)
library(stats)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(viridis)


#Read in data
fit <- read_csv("./fitness_data/Fitness_data_simplified.csv")
samples <- read_csv("../02_degs/assembled_transcripts_final/sample.info.all.modified.csv")

#Clean up the clone names (some spelling errors)
fit$Clone[fit$Clone == 'SWP4'] <- 'CBSWP4'
fit$Clone[fit$Clone == 'NWP7'] <- 'CBNWP7'
fit$Clone[fit$Clone == 'ARLCO10'] <- 'ARLC10'
fit$Clone[fit$Clone == 'LPI'] <- 'LP1'

#Remove CBSWP4
fit <- fit[fit$Clone != 'CBSWP4', ]

#Convert Clone to Genotype Clone
fit$Clone[fit$Clone == 'ARLC10'] <- 'C1'
fit$Clone[fit$Clone == 'BING2'] <- 'C2'
fit$Clone[fit$Clone == 'CBNWP7'] <- 'C3'
fit$Clone[fit$Clone == 'CF1'] <- 'C4'
fit$Clone[fit$Clone == 'LF8'] <- 'C4'
fit$Clone[fit$Clone == 'LP1'] <- 'C5'
fit$Clone[fit$Clone == 'PWA10'] <- 'C5'
fit$Clone[fit$Clone == 'PWA7'] <- 'C6'

fit$Temperature <- as.factor(fit$Temperature)

# Step 1: Calculate Cumulative Offspring
fit <- fit %>%
  group_by(Unique_ID) %>%
  arrange(Exp_Hour) %>%
  mutate(Cumulative_Babies = cumsum(Babies)) %>%
  ungroup()

# Step 2: Calculate Average Cumulative Offspring per Clone and Treatment
avg_cumulative <- fit %>%
  group_by(Exp_Hour, Clone, Temperature) %>%
  summarize(Mean_Cumulative_Babies = mean(Cumulative_Babies, na.rm = TRUE),
            SD_Cumulative_Babies = sd(Cumulative_Babies, na.rm = TRUE)) %>%
  ungroup()

# Step 3: Plot Average Cumulative Offspring with Error Bars
# Plot theme
nice_layout <- theme_cowplot() +
  theme(plot.background = element_rect("white"))

# Define color palette
colors <- viridis(length(unique(fit$Clone)), option = "plasma", direction = -1, begin = 0)

# Plotting function with aesthetics
plot_data <- function(clone_name, color) {
  clone_data <- avg_cumulative %>% filter(Clone == clone_name)
  
  p <- ggplot(clone_data, aes(x = Exp_Hour, y = Mean_Cumulative_Babies, color = Clone, shape = Temperature)) +
    geom_errorbar(aes(ymin = Mean_Cumulative_Babies - SD_Cumulative_Babies, 
                      ymax = Mean_Cumulative_Babies + SD_Cumulative_Babies), 
                  width = 2.5, colour = "black") +
    geom_line(aes(linetype = Temperature), linewidth = 1) +
    geom_point(size = 3.5) +
    scale_x_continuous(breaks = unique(fit$Exp_Hour)) +
    scale_color_manual(values = color) +
    scale_y_continuous(limits = c(-5, max(avg_cumulative$Mean_Cumulative_Babies) * 1.4)) + #Adjusted Y axis limit.
    nice_layout +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  return(p)
}

# Generate plots for each clone
unique_clones <- unique(fit$Clone)
plots <- list()

for (i in 1:length(unique_clones)) {
  plots[[i]] <- plot_data(unique_clones[i], colors[i])
}

# Arrange plots in a grid
p6 <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]],
                nrow = 2,
                ncol = 3,
                common.legend = TRUE)

print(p6)
#Save plot
wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/MS_v1_figures_tables_submit/Fig2_MS/")
ggsave(paste(wdII, "plot_average_cumulative_offspring_perGenotypeClone_perDay_updated_genotype_clone_split_grid_v6.svg", sep = ""), p6, dpi = 1000, height =7 , width = 12)

