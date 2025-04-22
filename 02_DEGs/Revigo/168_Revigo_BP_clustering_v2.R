#=============================================================================================================#
# Script created by Mark Christie, all rights reserved, contact at markchristie1500@gmail.com
# Script edited by Nathan Backenstose, nbackens@purdue.edu
# Script created in version R 4.2.2 on xx/xx/2023
# This script:
# Usage notes: 
#============================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/Revigo_NO_CBSWP/")

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

display.brewer.all()
display.brewer.pal(n = 8, name = 'Dark2')
colors <- brewer.pal(n = 8, name = "Dark2")

list.files()

#used Drosophila for Revigo; medium setting
dat <- read.table("Revigo_168_results/Revigo_BP_Scatterplot.tsv", header=TRUE, sep="\t")
head(dat)

go <- read.table("./go_168_NO_CBSWP.txt", header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
head(go)

# only write column names once if appending to file
#write.table(output.reduced, "FST_output_reduced.txt", col.names=!file.exists("FST_output_reduced.txt"), sep="\t", append = TRUE, row.names = FALSE)
#=============================================================================================================#
head(dat)
table(dat$Value)

dat2 <- dat[order(dat$Value, decreasing = TRUE), ]
head(dat2)
nrow(dat2)

# matching to input go terms and counts (are the counts correlated)
m1   <- match(dat2[, 1], go[, 1])
dat3 <- cbind(go[m1,], dat2)
plot(dat3[, 2], dat3$Value) # yes they are correlated

plot(dat3$PC_0, dat3$PC_1) # lots of null values!

dat4 <- dat3[-which(dat3$PC_0 == "null"), ]

size.key <- dat4$Value
keymax   <- max(size.key)
cex.max  <- 4
cex.min  <- 0.65
ranges   <- (cex.max-cex.min)/keymax
values   <- seq(from = 0.5, to = cex.max, by = ranges)
key <- cbind(1:keymax, values[-1])
m1  <- match(dat4$Value, key[, 1])
cexs<- key[m1, ]
cexs <- cexs*1.5
dat4 <- cbind(dat4, cexs[, 2])
dat4$row_id <- 1:nrow(dat4)


#############STOP HERE IF CLUSTERS AND PLOTS WERE ALREADY MADE #####################
#############                 SKIP TO LINE ~250                #####################  

#K means clustering

#determine best k (elbow method)
#in the plot, look for the sharpest decline

dat4$PC_0 <- as.numeric(dat4$PC_0)
dat4$PC_1 <- as.numeric(dat4$PC_1)
scaled_data <- as.data.frame(scale(dat4[, 10:11]))

#calculate within-cluster sum of squares across 10 k's
tot_withinss <- map_dbl(1:10, function(k) {
  model <- kmeans(x = scaled_data, centers = k)
  model$tot.withinss
})

cluster_df <- data.frame(
  k = 1:10,
  tot_withinss = tot_withinss
)

d <- ggplot(cluster_df, aes(x = k, y = tot_withinss)) +
  geom_line() + geom_point()+
  scale_x_continuous(breaks = 1:10) 
#nice_layout

d


#save
wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/Fig3_build/Revigo_NO_CBSWP/")
#ggsave()

#Trying a silouhette analysis as another test
#Values close to 1 suggest that the observation is well matched to the assigned cluster
#Values close to 0 suggest that the observation is borderline matched between two clusters
#Values close to -1 suggest that the observations may be assigned to the wrong cluster

sil_width <- map_dbl(2:10,  function(k){
  model <- pam(x = scaled_data, k = k)
  model$silinfo$avg.width
})

sil_df <- data.frame(
  k = 2:10,
  sil_width = sil_width
)

d2 <- ggplot(sil_df, aes(x = k, y = sil_width)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 2:10)

d2
#7 is the best

##Final test, gap statistic
set.seed(123)
gap_stat <- clusGap(scaled_data, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)
#9 for the gap statistics

###########another k means testing software#######
# Elbow method
fviz_nbclust(scaled_data, kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2)+
  labs(subtitle = "Elbow method")
# Silhouette method
fviz_nbclust(scaled_data, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(scaled_data, kmeans, nstart = 25,  method = "gap_stat", nboot = 100)+
  labs(subtitle = "Gap statistic method")

########################################################
####### APPLY K MEANS CLUSTERS AND CREATE PLOT #########
########################################################
dev.off()

#Create SVG
wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/Fig3_build/Revigo_NO_CBSWP/")
svg(paste(wdII, "168_Revigo_BP_kmeans_clustering_v2.svg", sep = ""), height = 6.5, width = 11)

#Assign clusters to GO terms
km.cluster <- kmeans(dat4[,10:11], centers = 9, iter.max = 10, nstart =2)
dat4$kmeans.cluster <- km.cluster$cluster
table(dat4$Name, km.cluster$cluster)

#scale points based on the number of times a GO term was counted
scaled_values <- 1 + 4 * (dat4$Value - min(dat4$Value)) / (max(dat4$Value) - min(dat4$Value))

color_palette <- brewer.pal(9, "Set1")

plot(dat4$PC_0, dat4$PC_1, 
     cex = scaled_values, 
     xlab = "PC 1", 
     ylab = "PC 2", 
     bty="l", 
     cex.axis = 1.5, 
     col= "black", 
     pch = 21,
     bg = color_palette[dat4$kmeans.cluster])

#Find the highest GO term frequency for the each kmeans cluster and report row_id for next step
max_values_per_cluster <- dat4 %>%
  group_by(kmeans.cluster) %>%
  summarise(Value = max(Value), row_id = row_id[which.max(Value)])

# Display the results
print(max_values_per_cluster)


# Each row includes the index of the point to annotate, and the new x and y coordinates for the label.
# use the row.id column in dat4 to choose the rows for the index. R row values are misleading.
annotations <- data.frame(
  index = c(26, 3, 6, 85, 13, 1, 5, 14, 33),  # Indices of points to annotate 
  label_x = c(0, 0, 0, 0, 0, 0, 0, 0, 0),  # x coordinates for the labels
  label_y = c(5, 4, 3, 2, 1, 0, -1, -2, -3)  # y coordinates for the labels
)

# Annotate the specified points
for (i in 1:nrow(annotations)) {
  point_to_annotate <- dat4[annotations$index[i], ]
  
  # Extract the label text from the "Name" column
  label_text <- dat4$Name[annotations$index[i]]
  
  # Add the text annotation at the specified location
  text(annotations$label_x[i], annotations$label_y[i], 
       labels = label_text, 
       pos = 4,  # Position text to the right of the point
       cex = 0.8)
  
  # Add the line from the point to the annotation
  segments(point_to_annotate$PC_0, point_to_annotate$PC_1, 
           annotations$label_x[i], annotations$label_y[i])  # Line from point to annotation
}


legend_values <- c(5, 150, 315)  # Example values from the "Value" column
legend_sizes <- 1 + 4 * (legend_values - min(dat4$Value)) / (max(dat4$Value) - min(dat4$Value))

# Add the legend
legend("topright", legend = legend_values, pt.cex = legend_sizes, pch = 21, col = "black", pt.bg = "gray", title = "GO term frequency")


#save plot

while (!is.null(dev.list()))  dev.off()

###############################
######EXPORT GO TABLE #########
###############################

#make new data frame to export
dat5 <- data.frame(
  Cluster = dat4$kmeans.cluster,
  GO = dat4$TermID,
  Frequency = dat4$Value,
  Description = dat4$Name
)

#order by kmeans cluster, then by value
dat5 <- arrange(dat5, Cluster, desc(Frequency))

#write to csv and excel spreadsheet
wdIII <- c("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/TableX-Revigo_cluster/")
write.csv(dat5, paste(wdIII, "168_Revigo_cluster_table.csv", sep=""), col.names = TRUE, row.names = FALSE, quote = FALSE)

write_xlsx(dat5, paste(wdIII, "168_Revigo_cluster_table.xlsx", sep=""))

color_palette


#####################################################################
#######Identify most common words in each cluster ###################
#####################################################################

#read in data with clustering already in place, so new clusters and colors are not assigned.
wdIII <- c("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/TableX-Revigo_cluster/")
dat5 <- read.csv(paste(wdIII, "168_Revigo_cluster_table.csv", sep =""))

#join with dataframe containing Revigo PCA, etc
dat6 <- dat4 %>%
  full_join(dat5, by = c("TermID" = "GO"))

#select just the columns I need
dat6 <- dat6 %>%
  select(3, 4, 5, 10, 11, 14, 15) %>%
  arrange(Cluster, desc(Value))

dat6$PC_0 <- as.numeric(dat6$PC_0)
dat6$PC_1 <- as.numeric(dat6$PC_1)


#Unnest the words from the "Name" column for each row
words_df <- dat6 %>%
  unnest_tokens(word, Name)  # Split the "Name" column into individual words

#Count the frequency of each word within each cluster
word_counts <- words_df %>%
  count(Cluster, word, sort = TRUE)

#Identify the top 10 words for each cluster
top_words_per_cluster <- word_counts %>%
  group_by(Cluster) %>%
  slice_max(n, n = 10, with_ties = FALSE)  # Select the top 10 words per cluster

#View
top_words_per_cluster
print(top_words_per_cluster, n = 90)


####################Plot inset/close up #####################

dev.off()

#Create SVG
wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/Fig3_build/Revigo_NO_CBSWP/")
svg(paste(wdII, "168_Revigo_BP_kmeans_clustering_v2_cluster_inset_v2.svg", sep = ""), height = 6.5, width = 11)

dat6_filtered <- dat6 %>%
  filter(Cluster %in% c(5, 7))

# Scale values for point size
scaled_values <- 1 + 4 * (dat6_filtered$Value - min(dat6_filtered$Value)) / 
  (max(dat6_filtered$Value) - min(dat6_filtered$Value))

# Color palette
color_palette <- brewer.pal(9, "Set1")

# Determine the range of the data
x_range <- range(dat6_filtered$PC_0)
y_range <- range(dat6_filtered$PC_1)

# Expand the range slightly (e.g., 10% beyond the range), this is to accomodate larger points associated with GO terms. Some were cut off.
x_margin <- diff(x_range) * 0.25
y_margin <- diff(y_range) * 0.25

# Set the new limits
xlim <- c(x_range[1] - x_margin, x_range[2] + x_margin)
ylim <- c(y_range[1] - y_margin, y_range[2] + y_margin)


plot(dat6_filtered$PC_0, dat6_filtered$PC_1, 
     cex = scaled_values, 
     xlab = "PC 1", 
     ylab = "PC 2", 
     bty="l", 
     cex.axis = 1.5, 
     col= "black", 
     pch = 21,
     bg = color_palette[dat6_filtered$Cluster],
     xlim = xlim,  # Apply the adjusted x limits
     ylim = ylim)  # Apply the adjusted y limits

# Find rows with "reproductive" or "signaling" in the "Value" column
#dat6_filtered$Value <- as.character(dat6_filtered$Value)

# Filter rows with "reproductive" or "signaling" in the "Value" column
annotations_filtered <- dat6_filtered %>%
  filter(str_detect(Name, regex("development|process", ignore_case = TRUE)))

# Split data into two clusters based on "Name" and "Cluster"
cluster_a_data <- annotations_filtered %>%
  filter(Cluster == "5" & str_detect(Name, "development"))

cluster_b_data <- annotations_filtered %>%
  filter(Cluster == "7" & str_detect(Name, "process"))

#get top twelve
cluster_a_data <- cluster_a_data %>%
  slice(1:12)

cluster_b_data <- cluster_b_data %>%
  slice(1:12)

# Combine the filtered data from both clusters
annotations_filtered <- bind_rows(cluster_a_data, cluster_b_data)

#change to numeric values
annotations_filtered$PC_0 <- as.numeric(annotations_filtered$PC_0)
annotations_filtered$PC_1 <- as.numeric(annotations_filtered$PC_1)

# Annotate the points on the plot
for (i in 1:nrow(annotations_filtered)) {
  point_to_annotate <- annotations_filtered[i, ]
  
  # Add the text annotation at the specified location
  text(point_to_annotate$PC_0, point_to_annotate$PC_1, 
       labels = point_to_annotate$Name, 
       pos = 3, # Position text to the right/left [4,2] of the point
       offset = 0 ,
       cex = 0.5)
  
  # Add the line from the point to the annotation
  segments(point_to_annotate$PC_0, point_to_annotate$PC_1, 
           point_to_annotate$PC_0 , point_to_annotate$PC_1 + 1)  # Adjust the coordinates as needed
}


#########testing############
#make non-annotated points grey
non_annotated <- dat6_filtered %>%
  filter(!(row_number() %in% annotations_filtered$row_number))  # Assuming row_number() or similar identifier is used

# Plot
points(non_annotated$PC_0, non_annotated$PC_1, 
       cex = scaled_values[!dat6_filtered$row_number %in% annotations_filtered$row_number], 
       pch = 21, 
       bg = "#808080")  # Grey color



#save plot

while (!is.null(dev.list()))  dev.off()



###############################################################
#######EXTRA STUFF NOT INCLUDED IN ANALYSIS ###################
###############################################################


########Clustering_test with random forest ########################################
#rf.fit <- randomForest(x=dat4[,10:11], y=NULL, ntree = 50000, proximity = TRUE, oob.prox = TRUE)
#hclust.rf <- hclust(as.dist(1-rf.fit$proximity), method = "ward.D2")
#rf.cluster = cutree(hclust.rf, k=6)
#dat4$rf.clusters <- rf.cluster
#table(rf.cluster, dat4$Name)

#plot(dat4$PC_0, dat4$PC_1, cex = 1, xlab = "PC 1", ylab = "PC 2", bty="l", cex.axis = 1.5, col=dat4$rf.clusters)

#nice_layout <- theme_cowplot()+
#  panel_border(color = "grey85", size = 1, linetype = 1,
#               remove = FALSE, "black") 

#p <- ggplot() +
#geom_point(aes(x=as.factor(dat4$PC_0), y=as.factor(dat4$PC_1)), color = dat4$rf.clusters) +
#xlab("PC 1") +
#ylab("PC 2") +
#nice_layout 

#p +
#scale_x_discrete(breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8), labels = c(-8, -6, -4, -2, 0, 2, 4, 6, 8)) +
#scale_y_discrete(breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8), labels = c(-8, -6, -4, -2, 0, 2, 4, 6, 8))
#theme(axis.text.x = element_blank(), axis.text.y = element_blank())

#p
#########################################
##ggplot2 test###

#nice_layout <- theme_cowplot()+
#anel_border(color = "grey85", size = 1, linetype = 1,
#            remove = FALSE, "black") 

#p <- ggplot() +
#geom_point(aes(dat4$PC_0, dat4$PC_1)) +
#xlab("PC 1") +
#ylab("PC 2") +
#nice_layout +
#scale_x_discrete(breaks=-6:6, labels=-6:6)
#theme(axis.text.x = element_blank(), axis.text.y = element_blank())

#p

###





###############################################################
#######EXTRA STUFF NOT INCLUDED IN ANALYSIS ###################
###############################################################


########Clustering_test with random forest ########################################
#rf.fit <- randomForest(x=dat4[,10:11], y=NULL, ntree = 50000, proximity = TRUE, oob.prox = TRUE)
#hclust.rf <- hclust(as.dist(1-rf.fit$proximity), method = "ward.D2")
#rf.cluster = cutree(hclust.rf, k=6)
#dat4$rf.clusters <- rf.cluster
#table(rf.cluster, dat4$Name)

#plot(dat4$PC_0, dat4$PC_1, cex = 1, xlab = "PC 1", ylab = "PC 2", bty="l", cex.axis = 1.5, col=dat4$rf.clusters)

#nice_layout <- theme_cowplot()+
#  panel_border(color = "grey85", size = 1, linetype = 1,
#               remove = FALSE, "black") 

#p <- ggplot() +
#geom_point(aes(x=as.factor(dat4$PC_0), y=as.factor(dat4$PC_1)), color = dat4$rf.clusters) +
#xlab("PC 1") +
#ylab("PC 2") +
#nice_layout 

#p +
#scale_x_discrete(breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8), labels = c(-8, -6, -4, -2, 0, 2, 4, 6, 8)) +
#scale_y_discrete(breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8), labels = c(-8, -6, -4, -2, 0, 2, 4, 6, 8))
#theme(axis.text.x = element_blank(), axis.text.y = element_blank())

#p
#########################################
##ggplot2 test###

#nice_layout <- theme_cowplot()+
#anel_border(color = "grey85", size = 1, linetype = 1,
#            remove = FALSE, "black") 

#p <- ggplot() +
#geom_point(aes(dat4$PC_0, dat4$PC_1)) +
#xlab("PC 1") +
#ylab("PC 2") +
#nice_layout +
#scale_x_discrete(breaks=-6:6, labels=-6:6)
#theme(axis.text.x = element_blank(), axis.text.y = element_blank())

#p

###

