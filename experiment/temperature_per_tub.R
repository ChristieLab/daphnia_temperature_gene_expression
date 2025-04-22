#D. pulex experiment temperature plotting
library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(scales)
library(patchwork)

wd <- ("Purdue/daphnia_RNAseq/results/experiment_data/")
setwd(wd)

#read in temperature data
d <- read.csv("./complete_temp_dp.csv", header = TRUE)

#filter data for start of experiment
d <- d[d$date_time >= as.POSIXct("2021-10-19 10:00:00"), ]


#set theme
nice_layout <- theme_cowplot()+
  panel_border(color = "grey85", size = 1, linetype = 1,
               remove = FALSE, "black") 

#adjust format of column
d$date_time <- as.POSIXct(d$date_time, format = "%Y-%m-%d %H:%M:%S")

# Calculate minimum and maximum values for x and y axes
x_min <- min(as.Date(d$date_time, format = "%Y-%m-%d %H:%M:%S"))
x_max <- max(as.Date(d$date_time, format = "%Y-%m-%d %H:%M:%S"))
y_min <- min(d$temp_c)
y_max <- max(d$temp_c)

#Boxplots
# Create the normal boxplot using ggplot2
 p <- ggplot(d, aes(x = tub, y = temp_c, fill = measured)) +
  geom_boxplot() +
  scale_fill_manual(values = c("hobo" = "lightblue", "manual" = "purple2")) +
  labs(title = "Temperature by Tub and Measurement Type",
       x = "Tub",
       y = "Temperature (°C)",
       fill = "Measurement Type") +
  nice_layout

p

#Create Standard Error Boxplot




############################
#Create a function to generate the plot for a specific tub
##Update with real times
first_sample	<- c("2021-10-23 12:02:00")
first_sample <- as.POSIXct(first_sample, format = "%Y-%m-%d %H:%M:%S")
#max_96	16:15:00
last_sample <- c("2021-10-26 11:00:00")
last_sample <- as.POSIXct(last_sample, format = "%Y-%m-%d %H:%M:%S") 
#max_168	14:06:00




plot_tub <- function(tub_value) {
  ggplot(d[d$tub == tub_value, ], aes(x = date_time, y = temp_c, color = measured, fill = measured)) +
    geom_point(aes(shape = measured), size = 1) +
    geom_vline(xintercept = first_sample, color = "grey", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = last_sample, color = "grey", linetype = "dashed", linewidth = 1) +
    scale_x_datetime(limits = c(as.POSIXct("2021-10-19 00:00:00"), as.POSIXct("2021-10-26 14:00:00")), breaks = "1 day", date_labels = "%d-%m") +
    scale_y_continuous(limits = c(y_min, y_max)) +
    scale_shape_manual(values = c("manual" = 2, "hobo" = 16)) +
    scale_color_manual(values = c("#0096FF", "#8B0000")) +
    labs(title = paste("Tub", tub_value), x = "Day", y = "Temperature (°C)") +
    nice_layout
    
}



# Generate plots for each tub
tub_plots <- lapply(unique(d$tub), plot_tub)

#check plot aesthetics
tub_plots[[1]]

# Arrange plots in a grid
grid <- tub_plots[[1]] + 
  tub_plots[[2]] +
  tub_plots[[3]] + 
  tub_plots[[4]] + 
  tub_plots[[5]] + 
  tub_plots[[6]] + 
  tub_plots[[7]] + 
  tub_plots[[8]] + 
  tub_plots[[9]] + 
  tub_plots[[10]] & theme(legend.position = "bottom")
grid

grid2 <- grid + plot_layout(guides = "collect", ncol = 2, axis_title = "collect")
grid2

#save plot
wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/experiment/")
ggsave(paste(wdII, "Temperature_of_tubs_boxplot_normal_v3.svg", sep = ""), p, height = 6.5, width = 9, dpi = 1000)
#ggsave(paste(wdII, "Temperature_of_tubs_boxplot_SE_v3.svg", sep = ""), p2, height = 6.5, width = 9, dpi = 1000)
ggsave(paste(wdII, "Temperature_of_tubs_v3.svg", sep = ""), grid2, height = 14, width = 12, dpi = 1000)



###Find min and max for each condition



control <- filter(d, tub %in% c("a", "c", "e", "g", "j"))
unique(control$tub)
min_control <- min(control$temp_c)
max_control <- max(control$temp_c)

experimental <- filter(d, tub %in% c("b", "d", "f", "h", "i"))
min_exp <- min(experimental$temp_c)
max_exp <- max(experimental$temp_c)
