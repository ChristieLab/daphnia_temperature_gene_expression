#Fig 2CD

###################################
###### UPDATE 10/03/24   ###########
### Only plot shared DEGs       ###
### pair up and down regulated  ###
### only report log2FC          ###
###################################


wd <- ("/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/Fig2cd_NO_CBSWP/")
setwd(wd)

library('DESeq2')
library("rtracklayer")
library("ggplot2")
library("gplots")
library("dplyr")
library("tidyverse")
library(cowplot)

list.files()

###Now do the same thing, but with the log2FC

rm(list = ls())
wd <- ("/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/Fig2cd_NO_CBSWP/")
setwd(wd)

library('DESeq2')
library("rtracklayer")
library("ggplot2")
library("gplots")
library("dplyr")
library("tidyverse")
library(cowplot)

list.files()

#read in list of genes that were differentially expressed on day 4 but not on day 7
list <- read.csv("./96_DEGs_not_DEGs_at_168_NO_CBSWP.csv")
overlap_list <- read.csv("./96_DEGs_overlap_DEGs_at_168_NO_CBSWP.csv")

#combine these and add a column to determine if they're from the overlap or not
listb <- bind_rows(list, overlap_list) 

#read in DESeq2 results from 96 and 168 to get log2FC for the genes in list above
dds <- readRDS("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/assembled_transcripts_final/DESeq_96_NO_CBSWP.genotype_clone.temp.rds")
dds2 <- readRDS("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/assembled_transcripts_final/DESeq_168_NO_CBSWP.genotype_clone.temp.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
res2 <- results(dds2, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)

#Double check DEGs
summary(res)
summary(res2)

#grab information from the 168 DESeq2 results and put them into a data frame.
rownames(res2)
p <- rownames(res2)
p <- data.frame(p)
p$log2FoldChange_168 <- res2$log2FoldChange
head(p)

#parse 168 results for just genes found on the overlap list and add in 
out <- left_join(overlap_list, p, by = c("list" = "p"))

#Column name change
colnames(out)[which(names(out) == "log2FoldChange")] <- "log2FoldChange_96"
colnames(out)[which(names(out) == "list")] <- "gene"

# Create a new column "expression" and initialize it with an empty string
out$expression <- ""

# Use if-else conditions to assign values based on if expression
out$expression[out$log2FoldChange_96 > 0 & out$log2FoldChange_168 > 0] <- "up"
out$expression[out$log2FoldChange_96 < 0 & out$log2FoldChange_168 < 0] <- "down"
out$expression[out$log2FoldChange_96 > 0 & out$log2FoldChange_168 < 0] <- "up_to_down"
out$expression[out$log2FoldChange_96 < 0 & out$log2FoldChange_168 > 0] <- "down_to_up"

#check to see if the numbers match up
table(out$expression)
#We're good, lets plot

#reorder columns
out <- out[, c("gene", "expression", "log2FoldChange_96", "log2FoldChange_168")]

head(out)


#combine into data frame for plotting purposes. experimental group only.
d <- cbind(out[1], out[2], stack(out[3:4]))
head(d)


#change column names
colnames(d) <- c("gene", "change", "expression", "day")

#replace some values
d <- d %>%
  mutate(day = case_when(
    day == "log2FoldChange_96" ~ "96",
    day == "log2FoldChange_168" ~ "168",
    TRUE ~ day
  ))

#Convert to factor
d$day <- factor(d$day, levels = c("96", "168"))

head(d)

#plot
nice_layout <- theme_cowplot()+
  panel_border(color = "grey85", size = 1, linetype = 1,
               remove = FALSE, "black") 

p3 <- ggplot(d, aes(x = day, y = expression, color = change)) +
  geom_line(aes(group = gene), linewidth = 0.5, color = "lightgrey") +
  geom_point(size = 2) +
  geom_hline(yintercept=0, linetype="dashed", color = "cadetblue") +
  scale_color_manual(values = c(up = "darkseagreen", down = "lightgreen", down_to_up = "darkgreen"), labels = c("Shared Positive", "Shared Negative", "Down_to_Up")) +
  labs(x = "Day", y = "Log2 Fold Change") +
  nice_layout +
  theme(legend.position = "none") 
p3


wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/MS_figures/Fig2_build/Fig2CD_NO_CBSWP/")
ggsave(paste(wdII,"96_sig_DEGs_SHARED_at_168__NO_CBSWP_adj_expresssion_up_and_down_v6_LOG2FC.svg"), plot = p3, dpi = 1000, height = 4, width = 6)
#ggsave(paste(wdII,"96_sig_DEGs_log2FC_at_168_ylim_v2.svg"), plot = p3_2, dpi = 1000, height = 4, width = 6)

#plot up and down separately
d_up <- filter(d, change %in% c("up"))

d_down <- filter(d, change %in% c("down", "down_to_up"))


#plot upregulated only
p3_up <- ggplot(d_up, aes(x = day, y = expression, color = change)) +
  geom_line(aes(group = gene), linewidth = 0.5, color = "lightgrey") +
  geom_point(size = 2) +
  scale_color_manual(values = c(up = "darkseagreen", down = "lightgreen", down_to_up = "darkgreen"), labels = c("Shared Positive", "Shared Negative", "Down_to_Up")) +
  labs(x = "Day", y = "Log2 Fold Change") +
  nice_layout +
  theme(legend.position = "none") 
p3_up

p3_down <- ggplot(d_down, aes(x = day, y = expression, color = change)) +
  geom_line(aes(group = gene), linewidth = 0.5, color = "lightgrey") +
  geom_point(size = 2) +
  scale_color_manual(values = c(up = "darkseagreen", down = "lightgreen", down_to_up = "darkgreen"), labels = c("Shared Positive", "Shared Negative", "Down_to_Up")) +
  labs(x = "Day", y = "Log2 Fold Change") +
  nice_layout +
  theme(legend.position = "none") 
p3_down

#save
ggsave(paste(wdII,"96_sig_DEGs_SHARED_at_168_NO_CBSWP_upregulated_v5_LOG2FC.svg"), plot = p3_up, dpi = 1000, height = 4, width = 6)
ggsave(paste(wdII,"96_sig_DEGs_SHARED_at_168_NO_CBSWP_downregulated_v5_LOG2FC.svg"), plot = p3_down, dpi = 1000, height = 4, width = 6)



