#Pairwise DEG comparison between Daphnia clones
#Author: NJCB - 02/24/25

library('DESeq2')
library("rtracklayer")
library("ggplot2")
library("gplots")
#library("beyonce")
library(dplyr)

wd <- ("/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/updated_figures_02-23-25/pairwise_DEGs/")
setwd(wd)
list.files()

#C1_C2
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C1")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C2")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C1 and C2

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C1_C2_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C1_C2_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C1_C3
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C1")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C3")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C1 and C3

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C1_C3_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C1_C3_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C1_C4
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C1")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C4")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C1 and C4

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

# Subsampling for C4
c4_indices <- which(col.data$genotype_clone == "C4")
if (length(c4_indices) > 10) {
  set.seed(123)
  selected_c4_indices <- sample(c4_indices, 10)
  c1_indices <- which(col.data$genotype_clone == "C1")
  final_indices <- c(c1_indices, selected_c4_indices)
  gene.counts <- gene.counts[, final_indices]
  col.data <- col.data[final_indices, ]
}


# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C1_C4_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C1_C4_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C1_C5
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C1")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C5")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C1 and C5

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

# Subsampling for C5
c5_indices <- which(col.data$genotype_clone == "C5")
if (length(c5_indices) > 10) {
  set.seed(123)
  selected_c5_indices <- sample(c5_indices, 10)
  c1_indices <- which(col.data$genotype_clone == "C1")
  final_indices <- c(c1_indices, selected_c5_indices)
  gene.counts <- gene.counts[, final_indices]
  col.data <- col.data[final_indices, ]
}

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C1_C5_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C1_C5_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C1_C6
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C1")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C6")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C1 and C6

# subset gene.counts for the designated treatment and clone
all(rownames(col.data) == colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C1_C6_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C1_C6_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C2_C3
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C2")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C3")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C2 and C3

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C2_C3_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C2_C3_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C2_C4
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C2")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C4")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C2 and C4

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

c4_indices <- which(col.data$genotype_clone == "C4")
if (length(c4_indices) > 10) {
  set.seed(123)
  selected_c4_indices <- sample(c4_indices, 10)
  c2_indices <- which(col.data$genotype_clone == "C2")
  final_indices <- c(c2_indices, selected_c4_indices)
  gene.counts <- gene.counts[, final_indices]
  col.data <- col.data[final_indices, ]
}

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C2_C4_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C2_C4_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C2_C5
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C2")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C5")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C2 and C5

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

# Subsampling for C5
c5_indices <- which(col.data$genotype_clone == "C5")
if (length(c5_indices) > 10) {
  set.seed(123)
  selected_c5_indices <- sample(c5_indices, 10)
  c2_indices <- which(col.data$genotype_clone == "C2")
  final_indices <- c(c2_indices, selected_c5_indices)
  gene.counts <- gene.counts[, final_indices]
  col.data <- col.data[final_indices, ]
}

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C2_C5_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C2_C5_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C2_C6
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C2")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C6")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C2 and C6

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C2_C6_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C2_C6_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C3_C4
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C3")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C4")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C3 and C4

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

c4_indices <- which(col.data$genotype_clone == "C4")
if (length(c4_indices) > 10) {
  set.seed(123)
  selected_c4_indices <- sample(c4_indices, 10)
  c3_indices <- which(col.data$genotype_clone == "C3")
  final_indices <- c(c3_indices, selected_c4_indices)
  gene.counts <- gene.counts[, final_indices]
  col.data <- col.data[final_indices, ]
}



# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C3_C4_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C3_C4_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C3_C5
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C3")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C5")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C3 and C5

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

# Subsampling for C5
c5_indices <- which(col.data$genotype_clone == "C5")
if (length(c5_indices) > 10) {
  set.seed(123)
  selected_c5_indices <- sample(c5_indices, 10)
  c3_indices <- which(col.data$genotype_clone == "C3")
  final_indices <- c(c3_indices, selected_c5_indices)
  gene.counts <- gene.counts[, final_indices]
  col.data <- col.data[final_indices, ]
}

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C3_C5_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C3_C5_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C3_C6
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C3")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C6")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C3 and C6

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C3_C6_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C3_C6_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C4_C5
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C4")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C5")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C4 and C5

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

#Subsample C4 and C5
c4_indices <- which(col.data$genotype_clone == "C4")
c5_indices <- which(col.data$genotype_clone == "C5")
if (length(c4_indices) > 10) {
  set.seed(123)
  selected_c4_indices <- sample(c4_indices, 10)
  set.seed(123)
  selected_c5_indices <- sample(c5_indices, 10)
  final_indices <- c(selected_c4_indices, selected_c5_indices)
  gene.counts <- gene.counts[, final_indices]
  col.data <- col.data[final_indices, ]
}

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C4_C5_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C4_C5_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C4_C6
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C4")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C6")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C4 and C6

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

c4_indices <- which(col.data$genotype_clone == "C4")
if (length(c4_indices) > 10) {
  set.seed(123)
  selected_c4_indices <- sample(c4_indices, 10)
  c6_indices <- which(col.data$genotype_clone == "C6")
  final_indices <- c(c6_indices, selected_c4_indices)
  gene.counts <- gene.counts[, final_indices]
  col.data <- col.data[final_indices, ]
}


# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C4_C6_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C4_C6_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################
#C5_C6
gene.counts <- read.table("./A7_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
head(gene.counts)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample.info.all.modified.csv", row.names = 4)
head(samples)

m1 <- which(samples[, 8] == "control" & samples[, 12] == "C5")
m2 <- which(samples[, 8] == "control" & samples[, 12] == "C6")

m1 <- c(m1, m2)

samples[m1, ] # should all be control and C5 and C6

# subset gene.counts for the designated treatment and clone
gene.counts <- gene.counts[, m1]

# subset samples for the relevant data
col.data <- samples[m1, ]

# Subsampling for C5
c5_indices <- which(col.data$genotype_clone == "C5")
if (length(c5_indices) > 10) {
  set.seed(123)
  selected_c5_indices <- sample(c5_indices, 10)
  c6_indices <- which(col.data$genotype_clone == "C6")
  final_indices <- c(c6_indices, selected_c5_indices)
  gene.counts <- gene.counts[, final_indices]
  col.data <- col.data[final_indices, ]
}

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
colnames(gene.counts)
rownames(col.data)

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#Construction of DESeqDataSet object and DE analysis
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(clone))
# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")           ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
# isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
seq.name      <- as.character(assembly$seqid[gene_idx])
transcript.id <- assembly$transcript_id[gene_idx]
gene.id       <- assembly$gene_id[gene_idx]
xloc          <- assembly$exon_number[gene_idx]
gene_names  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")

which(duplicated(dds@rowRanges@partitioning@NAMES)) # check that no gene names in dds are duplicates

# ----------------------Normalize the gene counts in dds -----------------------
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

saveRDS(dds, file="./control_results/C5_C6_control.clone.rds")

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)
sink(file = "./control_results/C5_C6_control_DEG_summary.txt")
length(which(res[, 6] < 0.05))
summary(res)
sink(file = NULL)
####################################################################################################################