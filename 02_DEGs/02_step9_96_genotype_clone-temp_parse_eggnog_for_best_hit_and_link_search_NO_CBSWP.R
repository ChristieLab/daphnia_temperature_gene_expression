setwd("/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/assembled_transcripts_final/")
list.files()
library("DESeq2")
library("rtracklayer")
library("ggplot2")
library("gplots")
library("seqinr")
library("dplyr")
library("stringr")
library("tidyverse")



dds <- readRDS("./DESeq_96_NO_CBSWP.genotype_clone.temp.rds")
assembly <- readGFF("stringtie_all_merged.gtf")

#--<>--<>--<>--<>-- Extract and examine results from DESeq analysis --<>--<>--<>--<>--
res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)

################### ID OF DEGs ############################################
#--<>--<>--<>--<>-- How many and which genes are sig DE? --<>--<>--<>--<>--
sum(res$padj < 0.05, na.rm=T)            
sum(res$padj < 0.05 & abs(res$log2FoldChange) >= 0, na.rm=T)
degs         <- res[which(res$padj < 0.05),]   ## put all genes w/ padj < 0.05 in resSig

#PLAN == match back to gtf, find location, extract sequences, blastx
#Extract locus names 
deg.names  <- rownames(degs)
deg.names2 <- unlist(strsplit(deg.names, "|", fixed = TRUE))   
trans.id   <- deg.names2[seq(from = 3, to = length(deg.names2), by = 3)]
gene.pos   <- deg.names2[seq(from = 3, to = length(deg.names2), by = 3)]  # use to double check matching

#Match to assembly
m1 <- match(trans.id, assembly$transcript_id)
m2 <- match(assembly$transcript_id, trans.id) # only 1 match? No; isolate all relevant matches
genes <- assembly[m1, ] # genes object has correct chromosome and start locations for all genes (to extract from genome for blast)

#match to GTF file; used in check below
gtf   <- readGFF("/Users/NJCB/Documents/Purdue/daphnia_RNAseq/reference/new.gtf") ## read in daphnia GTF

setwd("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/degs_annotations_NO_CBSWP/96_genotype-clone.temp_NO_CBSWP_emapper_search_default/")
list.files()

#read in files from eggnog (pre-parsed with regex to split first column into gene chrom+location and then a new column for hit number)
dii <- read_tsv("./96.NO_CBSWP.genotype-clone.temp.annot.emapper.annotations", comment = "#", col_names = TRUE)

#aggregate using deplyr for each query and the lowest evalue

eg <- dii %>%
  group_by(query) %>%
  slice(which.min(evalue)) %>%
  ungroup()


#######OLD METHOD
#pick best hit from e value and retain columns query, evalue, preferred gene name, description, and GO terms
#best <- aggregate(dii$evalue ~ dii$query+dii$seed_ortholog+dii$score+dii$eggNOG_OGs+dii$max_annot_lvl+dii$COG_category+dii$Description+dii$Preferred_name+dii$GOs+dii$EC+dii$KEGG_ko+dii$KEGG_Pathway+dii$KEGG_Module+dii$KEGG_Reaction+dii$KEGG_rclass+dii$BRITE+dii$KEGG_TC+dii$CAZy+dii$BiGG_Reaction+dii$PFAMs, FUN = min)
#best <- aggregate(dii, dii$evalue ~ dii$query+dii$score+dii$Description+dii$Preferred_name+dii$GOs, FUN = min)
#######



eg[c('seqid', 'start')] <- str_split_fixed(eg$query, ':', 2)
eg <- eg[c("seqid", "start", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")]
eg[c('start', 'end')] <- str_split_fixed(eg$start, '-', 2)
eg <- eg[c("seqid", "start", "end", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")]
#eg[c('end', 'number')] <- str_split_fixed(eg$end, '_', 2)
eg_df <- eg[c("seqid", "start", "end", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")]
eg_df$start = as.integer(as.character(eg_df$start))
eg_df$end = as.integer(as.character(eg_df$end))

#combine with padj <0.05 degs
degs_data <- as.data.frame(degs)
com1 <- cbind(genes, degs_data)
degs_df <- com1[c("seqid", "start", "end", "strand", "gene_id", "transcript_id", "gene_name", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

merged_annot <- degs_df %>% inner_join(eg_df, by=c('seqid'='seqid', 'start'='start', 'end'='end'))
merged_annot <- merged_annot[c("seqid", "start", "end", "strand", "gene_id", "transcript_id", "gene_name", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Preferred_name", "Description", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")]
final_merged <- merged_annot[order(-rank(merged_annot$log2FoldChange), decreasing = FALSE), ]

###########SAVE LINKED SPREADSHEETS [CHANGE]###################
setwd("/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/merged_deg_annotations_NO_CBSWP/")
write.table(final_merged, "96.NO_CBSWP.genotype_clone-temp.merged_annotations.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t", append = FALSE)


