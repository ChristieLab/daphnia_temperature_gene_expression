#Plot PCA results from PLINK2
#Author: NJCB
#Date: 02/23/2025

#load packages
library(tidyverse)
library(viridis)
library(cowplot)
library(ggsci)
library(ggrepel)
library(grid)
library(gridExtra)
library(Cairo)
library(cowplot)
library(vcfR)
library(ggConvexHull)
library(svglite)


#set working directory and check
setwd("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/updated_figures_02-23-25/PCA/")
getwd()

#read in data
pca <- read_table("./all_samples_NO_BQSR_NO_CBSWP_SNPS_filterOut_clean_maf0.05_hwe0_geno0.1_mind0.2_recode_4.2.eigenvec", col_names = TRUE)
eigenval <- scan("./all_samples_NO_BQSR_NO_CBSWP_SNPS_filterOut_clean_maf0.05_hwe0_geno0.1_mind0.2_recode_4.2.eigenval")

#sort out the pca data
#remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
#read in sample data
samples <- read.csv("./sample.info.all.modified.csv")
#had to remov a sample ("s180") that was filtered out for missing data 
#[may need to do a join with other plots in the future] 
samples <- head(samples, -1)

#remake data frame
pca <- as_tibble(data.frame(pca, samples$rename_id, samples$genotype_clone, samples$time, samples$temp, samples$site))

#convert to percentage variance explained.
#Note: use the number of PCs you calculated in pipeline
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

#plot percent variance
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
wdII <- ("~/Purdue/daphnia_RNAseq/figures_tables/MS_figures/MS_v1_figures_tables_submit/Sup_Figs/")
#ggsave(paste(wdII, "mito_skree.png", sep = ""), plot = a, height = 5, width = 5, dpi = 1000, units = "in")


#calculate cuulative sum of the percentage variance explained
cumsum(pve$pve)

#create theme for plot
nice_layout <- theme_cowplot() +
  theme(plot.background = element_rect("white"),
  panel_border(color = "grey85", size = 1, linetype = 1,
              remove = FALSE, "black")) 
  


#plot PCA(s)
p1_2 <- ggplot(pca, aes(PC1, PC2, color = samples.genotype_clone, shape = samples.genotype_clone)) + 
  geom_jitter(size = 3, width = .005, height = .005, alpha = 0.8) +
  scale_shape_manual(values = c(16, 16, 16, 16, 16, 16)) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  #coord_equal() +
  nice_layout +
  scale_color_viridis(discrete = TRUE, option = "plasma", direction = -1, begin = 0)

p1_2

p1_2nl <- p1_2+ theme(legend.position = "none")
p1_2nl



p1_3 <- ggplot(pca, aes(PC1, PC3, color = samples.genotype_clone, shape = samples.genotype_clone)) + 
  geom_jitter(size = 3, width = .005, height = .005, alpha = 0.8) +
  scale_shape_manual(values = c(16, 16, 16, 16, 16, 16)) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  #coord_equal() +
  nice_layout +
  scale_color_viridis(discrete = TRUE) +
  theme(legend.position = "none")
p1_3

p1_4 <- ggplot(pca, aes(PC1, PC4, color = samples.genotype_clone, shape = samples.genotype_clone)) + 
  geom_jitter(size = 3, width = .005, height = .005, alpha = 0.8) +
  scale_shape_manual(values = c(16, 16, 16, 16, 16, 16)) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)")) +
  #coord_equal() +
  nice_layout +
  scale_color_viridis(discrete = TRUE) +
  theme(legend.position = "none")

p1_4

p2_3 <- ggplot(pca, aes(PC2, PC3, color = samples.genotype_clone, shape = samples.genotype_clone)) + 
  geom_jitter(size = 3, width = .005, height = .005, alpha = 0.8) +
  scale_shape_manual(values = c(16, 16, 16, 16, 16, 16)) +
  xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  #coord_equal() +
  nice_layout +
  scale_color_viridis(discrete = TRUE) +
  theme(legend.position = "none")
p2_3

p2_4 <- ggplot(pca, aes(PC2, PC4, color = samples.genotype_clone, shape = samples.genotype_clone)) + 
  geom_jitter(size = 3, width = .005, height = .005, alpha = 0.8) +
  scale_shape_manual(values = c(16, 16, 16, 16, 16, 16)) +
  xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)")) +
  #coord_equal() +
  nice_layout +
  scale_color_viridis(discrete = TRUE) +
  theme(legend.position = "none")
p2_4

p3_4 <- ggplot(pca, aes(PC3, PC4, color = samples.genotype_clone, shape = samples.genotype_clone)) + 
  geom_jitter(size = 3, width = .005, height = .005, alpha = 0.8) +
  scale_shape_manual(values = c(16, 16, 16, 16, 16, 16)) +
  xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)")) +
  #coord_equal() +
  nice_layout +
  scale_color_viridis(discrete = TRUE) +
  theme(legend.position = "none")
p3_4


#capture legend
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

leg <- g_legend(p1_2) 






#plot grid
setwd("~/Purdue/daphnia_RNAseq/figures_tables/MS_figures/MS_v1_figures_tables_submit/Fig1_MS/")
wdII <- ("~/Purdue/daphnia_RNAseq/figures_tables/MS_figures/MS_v1_figures_tables_submit/Fig1_MS/")
ggsave(paste(wdII, "Fig1C_PCA_no_legend.svg", sep = ""), p1_2nl, units = "in", dpi=1000, height=6, width=10)
ggsave(paste(wdII, "Fig1C_PCA_no_legend.png", sep = ""), p1_2nl, units = "in", dpi=1000, height=6, width=10)
ggsave(paste(wdII, "Fig1C_PCA.svg", sep = ""), p1_2, units = "in", dpi=1000, height=6, width=10)
ggsave(paste(wdII, "Fig1C_PCA.png", sep = ""), p1_2, units = "in", dpi=1000, height=6, width=10)




#svg("Fig1C_PCA.svg", height=10,width=12)
#grid.newpage()
#pushViewport(viewport(layout = grid.layout(3,3))) # specify blank 3x3 grid with set widths
#grid.draw(leg)
#vplayout <- function(x,y) viewport(layout.pos.row = x, layout.pos.col = y)
#print(p1_2, vp = vplayout(1,1)) # plot legend
#print(p1_2nl, vp = vplayout(1,1))
#print(p1_3, vp = vplayout(2,1))
#print(p1_4, vp = vplayout(3,1))
#print(p2_3, vp = vplayout(2,2))
#print(p2_4, vp = vplayout(3,2))
#print(p3_4, vp = vplayout(3,3))

popViewport()

dev.off()
while (!is.null(dev.list()))  dev.off()

#save

ggsave(paste("./PCA_MITO_LEGEND.svg", sep=""), plot=leg, units = "in", dpi=1000, height=6, width=8)

#################################
#MC main PCA edit

#set working directory and check
setwd("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/updated_figures_02-23-25/PCA/")
getwd()

pca <- read.table("./pca_output.txt", header = TRUE)


# sort out the individual species and pops
#read in sample data
samples <- read.csv("./sample.info.all.modified.csv")
#had to remov a sample ("s180") that was filtered out for missing data 
#[may need to do a join with other plots in the future] 
samples <- head(samples, -1)


pca_plot_data <- cbind(pca, samples)
pca_plot_data <- as_tibble(pca_plot_data)

pve <- data.frame(PC = 1:2, pve = c(58.5, 14.3)) # Replace with your actual values

nice_layout <- theme_cowplot() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "grey85", size = 1, linetype = 1, fill = NA))

p1_2 <- ggplot(pca_plot_data, aes(x = Dim1, y = Dim2, color = genotype_clone)) +
  geom_jitter(size = 3, width = 0.005, height = 0.005, alpha = 0.8) +
  #scale_shape_manual(values = c(16, 16, 16)) + # Adjust the number of values to match your number of groups
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_cowplot() +
  scale_color_viridis_d(option = "plasma", direction = -1, begin = 0) +
  theme(
    panel.grid.major = element_line(color = "lightgrey"),
    panel.border = element_rect(color = "black"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Display the plot with the legend
print(p1_2)

# Display the plot without the legend
p1_2nl <- p1_2 + theme(legend.position = "none")
print(p1_2nl)


wdII <- ("~/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/Fig1_MS/")
ggsave(paste(wdII, "Fig1b_PCA_no_legend.svg", sep = ""), p1_2nl, units = "in", dpi=1000, height=6, width=8)
ggsave(paste(wdII, "Fig1b_PCA_no_legend.png", sep = ""), p1_2nl, units = "in", dpi=1000, height=6, width=8)

