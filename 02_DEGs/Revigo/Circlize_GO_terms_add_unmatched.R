library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(randomForest)
library(purrr)
library(cluster)
library(factoextra)
#install.packages("NbClust")
library(NbClust)
library(writexl)
library(stringr)
library(tidytext)
library(dplyr)
library(circlize)


wdIII <- c("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/TableX-Revigo_cluster/")
dat_96 <- read.csv(paste(wdIII, "96_Revigo_cluster_table_v2.csv", sep =""))
dat_168 <- read.csv(paste(wdIII, "168_Revigo_cluster_table.csv", sep =""))

# Get unique clusters from both data frames

unique_clusters_96 <- unique(dat_96$Cluster)
unique_clusters_168 <- unique(dat_168$Cluster)

# Create the matrix
comparison_matrix <- matrix(0, nrow = length(unique_clusters_168), ncol = length(unique_clusters_96),
                   dimnames = list(unique_clusters_168, unique_clusters_96))

# Populate the matrix
for (cluster_168 in unique_clusters_168) {
  go_terms_168 <- dat_168$GO[dat_168$Cluster == cluster_168]
  
  # Iterate through each cluster in dat_96
  for (cluster_96 in unique_clusters_96) {
    go_terms_96 <- dat_96$GO[dat_96$Cluster == cluster_96]
    
    # Count the number of shared GO terms
    shared_go_terms <- length(intersect(go_terms_168, go_terms_96))
    comparison_matrix[cluster_168, cluster_96] <- shared_go_terms
  }
}

# Add a new row and column for unmatched terms
comparison_matrix <- rbind(comparison_matrix, rep(0, ncol(comparison_matrix)))
comparison_matrix <- cbind(comparison_matrix, rep(0, nrow(comparison_matrix)))

# Create a set of all GO terms in dat_168
all_go_terms_168 <- unique(dat_168$GO)

# Count unmatched terms
unmatched_count <- 0
for (i in seq_along(unique_clusters_96)) {
  unique_cluster_96 <- unique_clusters_96[i]
  go_terms_96 <- dat_96$GO[dat_96$Cluster == unique_cluster_96]
  
  unmatched_count <- unmatched_count + length(setdiff(go_terms_96, all_go_terms_168))
}

# Add unmatched count to the new row and column
comparison_matrix[nrow(comparison_matrix), ncol(comparison_matrix)] <- unmatched_count

# Create a set of all GO terms in dat_96
all_go_terms_96 <- unique(dat_96$GO)

# Count unmatched terms in dat_168
unmatched_count_168 <- 0
for (i in seq_along(unique_clusters_168)) {
  unique_cluster_168 <- unique_clusters_168[i]
  go_terms_168 <- dat_168$GO[dat_168$Cluster == unique_cluster_168]
  
  unmatched_count_168 <- unmatched_count_168 + length(setdiff(go_terms_168, all_go_terms_96))
}

# Add a new row and column for unmatched terms in dat_168
comparison_matrix <- rbind(comparison_matrix, rep(0, ncol(comparison_matrix)))
comparison_matrix <- cbind(comparison_matrix, rep(0, nrow(comparison_matrix)))

# Add unmatched count for dat_168 to the new row and column
comparison_matrix[nrow(comparison_matrix), ncol(comparison_matrix)] <- unmatched_count_168


#rename columns
colnames(comparison_matrix) <- c("96_1", "96_2", "96_3", "96_4", "96_5", "96_6", "96_7", "96_8", "96_9", "96_96_1", "168_168_1")
rownames(comparison_matrix) <- c("168_1", "168_2", "168_3", "168_4", "168_5", "168_6", "168_7", "168_8", "168_9", "96_96_2", "168_168_2")

###Plot chord diagram#### 
#set colors for clusters
grid.col = c("96_1" = "#E41A1C", 
             "96_2" = "#377EB8", 
             "96_3" = "#4DAF4A", 
             "96_4" = "#984EA3", 
             "96_5" = "#FF7F00", 
             "96_6" = "#FFFF33", 
             "96_7" = "#A65628", 
             "96_8" = "#F781BF", 
             "96_9" = "#999999", 
             "168_1" = "#E41A1C", 
             "168_2" = "#377EB8", 
             "168_3" = "#4DAF4A", 
             "168_4" = "#984EA3", 
             "168_5" = "#FF7F00", 
             "168_6" = "#FFFF33", 
             "168_7" = "#A65628", 
             "168_8" = "#F781BF", 
             "168_9" = "#999999", 
             "96_96_1" = "#B2BEB5", 
             "96_96_2" = "#B2BEB5",
             "168_168_1" = "#71797E",
             "168_168_2" = "#71797E")

circos.par(start.degree = 90)

circos <- chordDiagram(comparison_matrix,
                       preAllocateTracks = 1,
                       annotationTrack = c("grid"),
                       annotationTrackHeight = c(0.08),
                       #order = c("168_3", "168_6", "168_8", "168_9", "168_7", "168_1", "168_2", "168_4", "168_5", "168_168_2", "96_96_1", "96_96_2", "168_168_1", "96_1", "96_2", "96_3", "96_4", "96_5", "96_6", "96_7", "96_8", "96_9"),
                       order = c("168_5", "168_2", "168_4", "168_1", "168_7", "168_9", "168_8", "168_6", "168_3", "168_168_2", "96_96_1", "96_96_2", "168_168_1", "96_9", "96_8", "96_7", "96_6", "96_5", "96_4", "96_3", "96_2", "96_1"),
                       big.gap = 30,
                       grid.col = grid.col)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

circos.clear()

###########
#unique GO terms data

unique_to_96 <- dat_96[!dat_96$GO %in% dat_168$GO,]
unique_to_168 <- dat_168[!dat_168$GO %in% dat_96$GO,]





