

setwd("/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/assembled_transcripts_final/")
list.files()
library("DESeq2") 
library("rtracklayer")
library("ggplot2")
library("gplots")
library("seqinr")

#load("steelhead_all.RData")  # steelhead session saved from "deseq2.R"

dds <- readRDS("./DESeq_96_NO_CBSWP_genotype_clone-temp.rds")
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
gtf   <- readGFF("/Users/NJCB/Documents/Purdue/daphnia_RNAseq/reference/new.gtf")             ## read in daphnia GTF
#gtf   <- readGFF("genomic.gtf")             ## read in rainbow trout GTF

m3    <- match(genes$gene_name, gtf$gene_id)  # Many matches simply missing
m4    <- m3[which(is.na(m3) == FALSE) ]       # remove NAs  
gtf[m4, ]  # Of those not missing, almost all are not annotated (argh!) 
gtf.match.back <- gtf[m4, ] # stored as an object to double check position of genes in genome (performed below)


# Creat input file for bedtools
head(genes)
genes.bed <- cbind(as.character(genes$seqid), genes$start, genes$end)

setwd("/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/degs_locations_NO_CBSWP/")
write.table(genes.bed, "96.NO_CBSWP.genotype-clone.temp.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t", append = FALSE)
write.table(genes, "96.NO_CBSWP.genotype-clone.temp.annotated_genes_pre_blast.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t", append = FALSE)



#CHECKS - NOT PART OF ANALYSES ##################################################################################################################################

#---------------------------------------------------------------------------------------------------------#
# Check that first position of merged GTF File contains entire range of sequence (answer: it does!)
OUT <- NULL
m5 <- unique(m2[which(is.na(m2) == FALSE) ])       # remove NAs 
for(i in 1:length(m5)){
  n   <- m5[i]
  n2  <- which(m2 == n) #islolate all positions where transcript id occurs in merged.gtf
  pos <- assembly[n2, ]  
  chromosome <- unique(pos$seqid) # using unique here as check to ensure never matching > 1 chromosome
  positions  <- c(pos$start, pos$end) 
  positions  <- range(positions) # find max and min sequence posistions to isolate from genome and blast
  out        <- cbind(as.character(chromosome), t(positions))
  OUT        <- rbind(OUT, out)
}

check <- cbind(genes, OUT)
check$start - as.numeric(as.character(check[, 15])) # should all = 0
check$end - as.numeric(as.character(check[, 16]))   # should all = 0
#---------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------#
# Check that positions from merged gtf file match those of rainbow trout gtf file (answer: start position always matches; end position sometimes differs a little, should be sufficient for blast)
head(gtf.match.back)
m6 <- match(gtf.match.back$gene, genes$gene_name)
tmp <- genes[m6, ]
tmp$start - gtf.match.back$start
tmp$end - gtf.match.back$end
#---------------------------------------------------------------------------------------------------------#
