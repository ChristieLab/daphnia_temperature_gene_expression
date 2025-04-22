library(EnhancedVolcano)
library(DESeq2)

wd <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/assembled_transcripts_final/")
setwd(wd)


############## 96 ##################
#Read in Deseq matrix
dds <- readRDS("./DESeq_96_NO_CBSWP.genotype_clone.temp.rds")

#calcluate significant results
res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
summary(res)

length(which(res[, 6] < 0.05))
summary(res)


#Volcano_Plot_test
v_96 <- EnhancedVolcano(res,
                lab = NA,
                title = "96 hours",
                subtitle = "Differential expression",
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 5e-2,
                FCcutoff = 3,
                legendPosition = "bottom",
                legendLabSize = 14,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                col = c('gray25', 'gray60', 'deepskyblue3', 'red'),
                colAlpha = 0.8)
v_96

############ 168 #################
#Read in Deseq matrix
dds2 <- readRDS("./DESeq_168_NO_CBSWP.genotype_clone.temp.rds")

#calcluate significant results
res2 <- results(dds2, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
summary(res2)

length(which(res2[, 6] < 0.05))
summary(res2)


#Volcano_Plot_test
v_168 <- EnhancedVolcano(res2,
                        lab = NA,
                        title = "168 hours",
                        subtitle = "Differential expression",
                        x = 'log2FoldChange',
                        y = 'padj',
                        pCutoff = 5e-2,
                        FCcutoff = 3,
                        legendPosition = "bottom",
                        legendLabSize = 14,
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE,
                        col = c('gray25', 'gray60', 'deepskyblue3', 'red'),
                        colAlpha = 0.8)
v_168
######## Save the volcano plots as svg ##############

wdII <- ("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/figures_tables/02_degs/volcano_plots/")
ggsave(paste(wdII, "96_NO_CBSWP_p0.05_log2FC5_volcano_panel.svg", sep = ""), v_96, dpi = 1000, height = 6, width = 8.5)
ggsave(paste(wdII, "168_NO_CBSWP_p0.05_log2FC5_volcano_panel.svg", sep = ""), v_168, dpi = 1000, height = 6, width = 8.5)

