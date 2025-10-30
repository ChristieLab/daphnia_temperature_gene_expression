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
library(broom)
library(purrr)


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

#Calculate Cumulative Offspring
fit <- fit %>%
  group_by(Unique_ID) %>%
  arrange(Exp_Hour) %>%
  mutate(Cumulative_Babies = cumsum(Babies)) %>%
  ungroup()

#Calculate Average Cumulative Offspring per Clone and Treatment
avg_cumulative <- fit %>%
  group_by(Exp_Hour, Clone, Temperature) %>%
  summarize(Mean_Cumulative_Babies = mean(Cumulative_Babies, na.rm = TRUE),
            SD_Cumulative_Babies = sd(Cumulative_Babies, na.rm = TRUE)) %>%
  ungroup()

#write table of data for Sup data.
fit_selection <- fit[, c(3,5,6,8,9,11,15)]
fit_selection <- fit_selection %>%
  rename(Offspring = Babies,
         Cumulative_Offspring = Cumulative_Babies)

write.csv(fit_selection, "~/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/Supplemental_data_1.csv", row.names = FALSE)

f #Plot Average Cumulative Offspring with Error Bars
# Plot theme
nice_layout <- theme_cowplot() +
  theme(plot.background = element_rect("white"))

# Correct order for the clones to match the viridis palette
clone_color_order <- c('C1', 'C2', 'C3', 'C4', 'C5', 'C6')

# Define the exact order you want the plots to appear in the grid
# This is the key to fixing the layout issue
plot_grid_order <- c('C1', 'C2', 'C3', 'C4', 'C5', 'C6')

# Generate a palette with the desired number of colors
color_palette <- viridis(length(clone_color_order), option = "plasma", direction = -1, begin = 0)

# Create a named vector to map colors to clone names
clone_colors <- setNames(color_palette, clone_color_order)

# Plotting function
plot_data <- function(clone_name) {
  clone_data <- avg_cumulative %>% filter(Clone == clone_name)
  
  p <- ggplot(clone_data, aes(x = Exp_Hour, y = Mean_Cumulative_Babies, color = Clone, shape = Temperature)) +
    geom_errorbar(aes(ymin = Mean_Cumulative_Babies - SD_Cumulative_Babies,
                      ymax = Mean_Cumulative_Babies + SD_Cumulative_Babies),
                  width = 2.5, colour = "black") +
    geom_line(aes(linetype = Temperature), linewidth = 1) +
    geom_point(size = 3.5) +
    scale_x_continuous(breaks = unique(fit$Exp_Hour)) +
    scale_color_manual(values = clone_colors) + # Use the named color vector
    scale_y_continuous(limits = c(-5, max(avg_cumulative$Mean_Cumulative_Babies) * 1.4)) +
    nice_layout +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  return(p)
}

# Generate plots and store them in a named list for easy access
plots <- list()
for (i in 1:length(plot_grid_order)) {
  plots[[plot_grid_order[i]]] <- plot_data(plot_grid_order[i])
}

# Arrange plots in a grid using the `plotlist` argument
# This ensures the order is correct
p6 <- ggarrange(plotlist = plots[plot_grid_order],
                nrow = 2,
                ncol = 3,
                common.legend = TRUE)

print(p6)

plots[[6]]

#Save plot
wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/Fig2_MS/")
#ggsave(paste(wdII, "plot_average_cumulative_offspring_perGenotypeClone_perDay_updated_genotype_clone_split_grid_v7.svg", sep = ""), p6, dpi = 1000, height =7 , width = 12)


#Review update stats on cumulative fecundity
#############################################
#Find the right test
data_subset <- fit %>%
  filter(Clone == "C1", Exp_Hour == 168, Temperature == 20)

shapiro.test(data_subset$Cumulative_Babies)
qqline(data_subset$Cumulative_Babies, col = "red")
# A mixed bag, going with the non-parametric test

#stats
wilcox_results_df <- fit %>%
  filter(Temperature %in% c("20", "25")) %>%
  group_by(Clone, Exp_Hour) %>%
  do(tidy(wilcox.test(Cumulative_Babies ~ Temperature, data = .))) %>%
  ungroup() %>%
  dplyr::select(Clone, Exp_Hour, p.value)

wilcox_results_dfII <- wilcox_results_df %>%
  group_by(Clone) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  ungroup()


# determine stats significance for each comparison (corrected)
sig_dataII <- wilcox_results_dfII %>%
  mutate(
    group1 = "20",
    group2 = "25",
    p.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE            ~ "ns" #ns = not significant
    )
  ) 

pvalue_table <- sig_dataII[, 1:4]
write.csv(pvalue_table, "~/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/Sup_tables/Supplemental_Table_2.csv", row.names = FALSE)

# determine stats significance for each comparison
sig_data <- wilcox_results_df %>%
  mutate(
    group1 = "20",
    group2 = "25",
    p.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns" #ns = not significant
    )
  ) %>%
  # We need to find the maximum y-value for each comparison to place the annotation.
  left_join(avg_cumulative, by = c("Clone", "Exp_Hour")) %>%
  group_by(Clone, Exp_Hour) %>%
  mutate(y.position = max(Mean_Cumulative_Babies) * 1.05) %>%
  ungroup() %>%
  dplyr::select(Clone, Exp_Hour, p.signif, y.position, group1, group2) %>%
  distinct()

#Apply to the plot


#Review update II (clutch size)
####################################

ggplot(fit, aes(x = Exp_Hour, y = Babies, group = Unique_ID, linetype = Temperature)) +
  geom_line(alpha = 1) +
  scale_linetype_manual(
    values = c("20" = "solid", "25" = "dotted"),
    name = "Temperature (\u00B0C)"
  ) +
  scale_x_continuous(breaks = c(24, 48, 72, 96, 120, 144, 168)) +
  labs(
    title = "Individual Daphnia Clutch Size Over Experimental Hours",
    x = "Experimental Hour",
    y = "Clutch Size (Number of Babies)"
  ) +
  facet_wrap(~ Clone, scales = "fixed") +
  nice_layout +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#time between clutches
time_between_clutches_df <- fit %>%
  filter(Babies > 0) %>%
  group_by(Unique_ID) %>%
  arrange(Exp_Hour) %>%
  mutate(time_between_clutches = Exp_Hour - lag(Exp_Hour)) %>%
  ungroup()

time_between_clutches_df$Temperature <- as.factor(time_between_clutches_df$Temperature)

p7 <- ggplot(time_between_clutches_df, aes(x = Clone, y = time_between_clutches)) +
  geom_boxplot(aes(middle = mean(time_between_clutches), fill = Temperature), outlier.shape = NA) +
  geom_jitter(aes(color = Temperature), width = 0.2, height = 0, alpha = 0.6) + 
  scale_fill_manual(
    values = c("20" = "lightblue", "25" = "lightcoral"), 
    labels = c("20" = "Control (20\u00B0C)", "25" = "Experimental (25\u00B0C)"), # Custom legend labels
    name = "Treatment"
  ) +
  scale_color_manual(
    values = c("20" = "darkblue", "25" = "darkred"), # Custom color for jittered points
    labels = c("20" = "Control (20\u00B0C)", "25" = "Experimental (25\u00B0C)"),
    name = "Treatment"
  ) +
  labs(
    title = "Time Between Clutches",
    x = "Clone",
    y = "Time Between Clutches (Hours)"
  ) +
  facet_wrap(~ Temperature, scales = "free_x") + # Separate plots for each temperature
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.background = element_rect("white")
  )


p7


wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/")
#ggsave(paste(wdII, "plot_time_between_clutches_perGenotypeClone_perDay_updated_genotype_clone.svg", sep = ""), p7, dpi = 1000, height =6 , width = 8)
#ggsave(paste(wdII, "plot_time_between_clutches_perGenotypeClone_perDay_updated_genotype_clone.png", sep = ""), p7, dpi = 1000, height =6 , width = 8)



#mean time between clutches
mean_time_between_clutches_summary <- time_between_clutches_df %>%
  group_by(Clone, Temperature) %>%
  summarise(
    Mean_Time_Between_Clutches = mean(time_between_clutches, na.rm = TRUE),
    SD_Time_Between_Clutches = sd(time_between_clutches, na.rm = TRUE),
    SE_Time_Between_Clutches = SD_Time_Between_Clutches / sqrt(n())
  ) %>%
  ungroup()

# Plot the mean time between clutches with error bars
p8 <- ggplot(mean_time_between_clutches_summary, aes(x = Clone, y = Mean_Time_Between_Clutches, color = Temperature)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) + # Plot points for means, dodged
  geom_errorbar(
    aes(ymin = Mean_Time_Between_Clutches - SE_Time_Between_Clutches,
        ymax = Mean_Time_Between_Clutches + SE_Time_Between_Clutches),
    width = 0.2, # Width of the error bar ends
    position = position_dodge(width = 0.5) # Dodge error bars to match points
  ) +
  scale_color_manual(
    values = c("20" = "darkblue", "25" = "darkred"), # Consistent colors with previous plot
    labels = c("20" = "Control (20\u00B0C)", "25" = "Experimental (25\u00B0C)"),
    name = "Treatment"
  ) +
  labs(
    title = "Mean Time Between Clutches by Clone and Treatment",
    x = "Clone",
    y = "Mean Time Between Clutches (Hours)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.background = element_rect("white")
  )

p8

wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/")
ggsave(paste(wdII, "plot_mean_time_between_clutches_perGenotypeClone_perDay_updated_genotype_clone.svg", sep = ""), p8, dpi = 1000, height =6 , width = 8)
ggsave(paste(wdII, "plot_mean_time_between_clutches_perGenotypeClone_perDay_updated_genotype_clone.png", sep = ""), p8, dpi = 1000, height =6 , width = 8)


#review update cont.
#normalize by number of clutches per treatment
######################################################


clutch_summary <- fit %>%
  group_by(Unique_ID, Temperature) %>%
  summarise(
    Total_Offspring = sum(Babies),
    Total_Clutches = sum(Babies > 0),
    .groups = 'drop'
  )

# Calculate the average clutch size for each individual.
# Handle cases where Total_Clutches is 0 to avoid division by zero.
clutch_summary <- clutch_summary %>%
  mutate(
    Avg_Clutch_Size = ifelse(Total_Clutches > 0, Total_Offspring / Total_Clutches, 0)
  )

# Separate the data into the two treatment groups.
temp_20 <- clutch_summary %>% filter(Temperature == 20)
temp_25 <- clutch_summary %>% filter(Temperature == 25)

# Perform the Mann-Whitney U test (wilcox.test in R)
mann_whitney_result <- wilcox.test(temp_20$Avg_Clutch_Size, temp_25$Avg_Clutch_Size,
                                   alternative = "two.sided")

print(mann_whitney_result)


p9 <- ggplot(clutch_summary, aes(x = Temperature, y = Avg_Clutch_Size, fill = Temperature)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = max(clutch_summary$Avg_Clutch_Size) * 1.1) +
  scale_fill_manual(
    values = c("20" = "lightblue", "25" = "lightcoral")
  ) +
  labs(
    title = "Average Clutch Size by Temperature Treatment",
    x = "Temperature Treatment (\u00B0C)",
    y = "Average Clutch Size"
  ) +
  nice_layout +
  theme(
    legend.position = "none"
  )
print(p9)

wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/")
#ggsave(paste(wdII, "plot_average_clutch_size_per_treatment.svg", sep = ""), p9, dpi = 1000, height =6 , width = 8)
#ggsave(paste(wdII, "plot_average_clutch_size_per_treatment.png", sep = ""), p9, dpi = 1000, height =6 , width = 8)


## The goal is to perform a Mann-Whitney U test to compare the distribution
# of time between clutches for the 20°C and 25°C treatments.

# Separate the 'time_between_clutches' data into the two treatment groups.
clutch_timing_20 <- time_between_clutches_df %>% filter(Temperature == 20)
clutch_timing_25 <- time_between_clutches_df %>% filter(Temperature == 25)

mean_time_between_clutches_summary_20 <- clutch_timing_20 %>%
  group_by(Clone, Temperature) %>%
  summarise(
    Mean_Time_Between_Clutches = mean(time_between_clutches, na.rm = TRUE),
    SD_Time_Between_Clutches = sd(time_between_clutches, na.rm = TRUE),
    SE_Time_Between_Clutches = SD_Time_Between_Clutches / sqrt(n())
  ) %>%
  ungroup()
  


# Perform the Mann-Whitney U test (wilcox.test in R).
mann_whitney_clutch_timing_result <- wilcox.test(clutch_timing_20$time_between_clutches, 
                                                 clutch_timing_25$time_between_clutches,
                                                 alternative = "two.sided")

print(mann_whitney_clutch_timing_result)


# Calculate the max y-position for the labels
label_y_pos <- max(time_between_clutches_df$time_between_clutches, na.rm = TRUE) * 1.05

# Perform Wilcoxon test for each Clone and apply Bonferroni correction
# We use `unnest` to get the p-values from the list format
clone_p_values <- time_between_clutches_df %>%
  group_by(Clone) %>%
  do(p = wilcox.test(time_between_clutches ~ Temperature, data = .)$p.value) %>%
  ungroup() %>%
  unnest(p)

# Apply the Bonferroni correction
clone_p_values$p.bonferroni <- p.adjust(clone_p_values$p, method = "bonferroni")
clone_p_values$p.holm <- p.adjust(clone_p_values$p, method = "holm")

# Add the y-position to the data frame for plotting
clone_p_values$y.position <- label_y_pos

# Modify the existing plot 'p7' to add the p-value.

p7_updated <- ggplot(time_between_clutches_df, aes(x = Clone, y = time_between_clutches, color = Clone)) +
  geom_boxplot(aes(fill = Temperature)) +
  geom_jitter(
    aes(color = Temperature),
    position = position_jitterdodge(
      jitter.width = 0.25,
      jitter.height = 0,
      dodge.width = 0.75
    ),
    alpha = 0.4
  ) +
  stat_summary(
    aes(group = interaction(Clone, Temperature), color = Temperature),
    fun = mean,
    geom = "crossbar",
    size = 0.5,
    width = 0.75, # Controls the length of the crossbar
    position = position_dodge(width = 0.75)
  ) +
  #stat_compare_means(
  #  aes(group = Temperature), 
  #  method = "wilcox.test", 
  #  label = "p.format",
  #  p.adjust.method = "bonferroni",
  #  label.y = max(time_between_clutches_df$time_between_clutches, na.rm = TRUE) * 1.05,
  #  size = 3
  #) +
  geom_text(
    data = clone_p_values,
    aes(x = Clone, y = y.position, label = paste0("p = ", format.pval(p.holm, digits = 4))),
    inherit.aes = FALSE, # This is important to not inherit main plot aesthetics
    size = 3
  ) +
  scale_fill_manual(
    values = c("20" = "#D5FFFF", "25" = "#FF7F7F"),
    labels = c("20" = "Control (20\u00B0C)", "25" = "Experimental (25\u00B0C)"),
    name = "Treatment"
  ) +
  scale_color_manual(
    values = c("20" = "darkblue", "25" = "darkred"),
    labels = c("20" = "Control (20\u00B0C)", "25" = "Experimental (25\u00B0C)"),
    name = "Treatment"
  ) +
  scale_y_continuous(breaks = seq(0, max(time_between_clutches_df$time_between_clutches, na.rm = TRUE) + 24, by = 24)) +
  labs(
    title = "Time Between Clutches",
    x = "Clone",
    y = "Time Between Clutches (Hours)"
  ) +
  facet_wrap(~ Clone, scales = "free_x") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.background = element_rect("white")
  )

print(p7_updated)
# 
wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/")
ggsave(paste(wdII, "plot_time_between_clutches_perGenotypeClone_perDay_updated_genotype_clone_with_test.svg", sep = ""), p7_updated, dpi = 1000, height =6 , width = 8)
ggsave(paste(wdII, "plot_time_between_clutches_perGenotypeClone_perDay_updated_genotype_clone_with_test.png", sep = ""), p7_updated, dpi = 1000, height =6 , width = 8)

