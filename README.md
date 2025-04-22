This repository contains the data and analysis scripts related to the study: "**Delayed transcriptional response of Daphnia pulex to thermal stress**" in review at *G3* 

**Description**

This study investigated the transcriptional response of *Daphnia pulex* to sublethal temperature stress. We aimed to understand how elevated temperatures affect gene expression, providing insights into their resilience and potential impacts on aquatic ecosystems.

## Data

* **Raw sequencing data:** FASTQ files are available on NCBI SRA under BioProject accession number: PRJNA1251057

## Analysis

* **Scripts:**
    * `01_map_and_call_variants`: Scripts for quality control, trimming, alignment, and SNP calling of RNA-Seq reads.
    * `02_DEGs`: Scripts for differential gene expression analysis using DESeq2.
         * `volcano`: Volcano plots from DESeq2 results (Fig2 D & E)
         * `pairwise_comparison_controls/exp`: Scripts for calculating pairwise DEGs between clones in the control and experimental treatment
         * `pairwise_DEGs_FST_corr`: Pairwise control group DEGs correlated with pairwise Weir & Cockerham FST values 
         * `Revigo`: Scripts for Gene Ontology (GO) analysis of DEGs.
    * `experiment`: Script for plotting temperature of tubs during the experiment
    * `FIS`: Scripts for calculating and plotting inbreeding coefficients
    * `fitness`: Script for calculating and plotting mean cumulative offspring throughout the experiment
    * `Fst_clones`: Scripts for calculating pairwise Weir & Cockerham FSt values between clones.
    * `mito_PCA`: Scripts for running and plotting PCA on mitochrondrial SNPs
    * `PCA`: Genome wide SNP PCA



