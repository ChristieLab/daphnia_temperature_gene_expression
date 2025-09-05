#Fig 3CD

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
library(readxl)

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

#update for reviewer comment (8/20/25)
##########################################
#identify how stable functional categories of the shared genes
#subtract log2fc values at 96 from 168

magnitude_of_change <- d %>%
  pivot_wider(
    id_cols = gene,
    names_from = day,
    values_from = expression
  ) %>%
  # Calculate the difference (Day 168 - Day 96)
  mutate(
    expression_difference = `168` - `96`
  ) %>%
  select(
    gene,
    expression_difference
  )

#plot
p4 <- ggplot(magnitude_of_change, aes(x = expression_difference)) +
  geom_histogram(binwidth = 0.25, fill = "lightblue", color = "black") +
  labs(
    title = "Expression Difference Between 96h and 168h",
    x = "Log2FC Difference (168h - 96h)",
    y = "Number of Genes"
  ) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60)) +
  nice_layout +
  theme(plot.background = element_rect("white"))
  

p4

wdII <- "~/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/"
ggsave(paste(wdII,"shared_DEGs_histo_expression_difference.svg"), plot = p4, dpi = 1000, height = 4, width = 6)
ggsave(paste(wdII,"shared_DEGs_histo_expression_difference.png"), plot = p4, dpi = 1000, height = 4, width = 6)



#how many change more than 0.25?
genes_with_significant_change <- magnitude_of_change %>%
  filter(abs(expression_difference) > 0.25)
#9 total

genes_without_significant_change <- magnitude_of_change %>%
  filter(abs(expression_difference) < 0.25)
#65 total

#now lets use the annotation to try and merge these two groups with functional annot
#read in the annots (just 96, doesn't matter because they are shared)

annot <- read_tsv("~/Purdue/daphnia_RNAseq/results/02_degs/merged_deg_annotations_NO_CBSWP/96.NO_CBSWP.genotype_clone-temp.merged_annotations.tsv")

# Create a new column with the stringtie id
magnitude_of_change_with_id <- magnitude_of_change %>%
  mutate(MSTRG_id = str_extract(gene, "MSTRG\\.\\d+"))

#merge
merged_data <- left_join(magnitude_of_change_with_id, annot, by = c("MSTRG_id" = "gene_id"))


gene_functional_data <- left_join(magnitude_of_change_with_id, annot, by = c("MSTRG_id" = "gene_id")) %>%
  select(
    gene,
    expression_difference,
    Description,
    GOs
  )

gos_96 <- read_tsv("~/Purdue/daphnia_RNAseq/results/02_degs/Revigo_NO_CBSWP/Revigo_96_results/Revigo_BP_Table.tsv")

unpacked_genes <- gene_functional_data %>%
  separate_rows(GOs, sep = ",") %>%
  mutate(GOs = str_trim(GOs))

merged_data_with_names <- unpacked_genes %>%
  left_join(gos_96, by = c("GOs" = "TermID")) %>%
  dplyr::select(
    gene,
    expression_difference,
    GO_ID = GOs,
    GO_Name = Name
  ) %>%
  drop_na(GO_Name)


common_words_to_remove <- c("of", "in", "and", "or", "to", "for", "a", "the", "with")

most_prevalent_categories <- merged_data_with_names %>%
  drop_na(GO_Name) %>%
  group_by(gene, expression_difference) %>%
  summarise(
    all_go_names = paste(GO_Name, collapse = " "),
    .groups = 'drop_last' 
  ) %>%
  rowwise() %>%
  mutate(
    words = str_split(all_go_names, "\\s+")
  ) %>%
  ungroup() %>%
  unnest(words) %>%
  mutate(
    words = str_to_lower(words),
    words = str_remove_all(words, "[[:punct:]0-9]")
  ) %>%
  filter(!words %in% common_words_to_remove, words != "") %>%
  group_by(gene) %>%
  dplyr::count(words, sort = TRUE) %>%
  slice_max(order_by = n, n = 1) %>%
  ungroup() %>%
  left_join(
    distinct(merged_data_with_names, gene, expression_difference),
    by = "gene"
  ) %>%
  arrange(gene, desc(n))

print(most_prevalent_categories)

####Create table for SI
#read in merged data from annotations
annot_96 <- read_excel("~/Purdue/daphnia_RNAseq/figures_tables/MS_figures/MS_v1_figures_tables_submit/Sup_Tables/Sup_Table3_96_DEGs_with_annot_full_table.xlsx")

# Get IDs of the 74 shared DEGs and extract MSTRG.id and chromosome
shared_deg_base_info <- data.frame(full_gene_id = unique(out$gene)) %>%
  mutate(
    MSTRG_id = str_extract(full_gene_id, "MSTRG\\.\\d+"),
    chromosome = str_extract(full_gene_id, "^[^|]+") # Extracts "NC_060017.1"
  ) %>%
  dplyr::select(MSTRG_id, chromosome) %>%
  distinct() # Ensure unique MSTRG_id and chromosome combinations


shared_deg_coords_from_annot_tsv <- shared_deg_base_info %>%
  left_join(annot %>% dplyr::select(gene_id, start, end),
            by = c("MSTRG_id" = "gene_id")) %>%
  dplyr::rename(start = start, end = end)


#left join to annot merge script object to get full coordinates
all_coords <- shared_deg_coords_from_annot_tsv %>%
  left_join(degs_df %>% dplyr::select(seqid, start, end, gene_id),
            by = c("MSTRG_id" = "gene_id")) %>%
  dplyr::select(chromosome, start.y, end.y, MSTRG_id)


all_coords <- all_coords %>%
  rename(start = start.y, end = end.y, gene_id = MSTRG_id)


final_deg_table <- all_coords %>%
  left_join(annot_96 %>% dplyr::select(chromosome, start, end, eggNog_gene_name, eggNog_gene_description),
            by = c("chromosome", "start", "end")) %>%
  dplyr::select(
    chromosome,
    start,
    end = end, 
    gene_id = gene_id, 
    eggNog_gene_name,
    eggNog_gene_description
  )

head(final_deg_table)

output_dir_table <- "~/Purdue/daphnia_RNAseq/MS_v2/overall_response_to_temperature/figures_and_tables/Sup_tables/"

if (!dir.exists(output_dir_table)) {
  dir.create(output_dir_table, recursive = TRUE)
}

# Write the table to a CSV file.
write.csv(final_deg_table,
          file = paste0(output_dir_table, "Sup_Table3_shared_DEGs_detailed_annotation_table.csv"),
          row.names = FALSE)
