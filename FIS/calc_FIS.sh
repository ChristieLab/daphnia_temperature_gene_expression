#!/usr/bin/env bash

#SBATCH --job-name=pairwise_clone_WC_FIS
#SBATCH -A beagle
#SBATCH -t 2:00:00
#SBATCH -n 1
#SBATCH --output="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/FIS.out"
#SBATCH --error="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/FIS.err"

DIR="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/FIS"
VCF="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/FIS/all_samples_NO_BQSR_NO_CBSWP_SNPS_filterOut_clean_maf0.05_hwe0_geno0.1_mind0.2_recode_4.2.vcf"
#POPLIST="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/Fst_clone/temp"
#POPREV="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/Fst_clone/tempII"
#CLONES="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/Fst_clone/clone_lists/"

module load biocontainers
module load vcftools/0.1.16

cd ${DIR}

vcftools --vcf ${VCF} --het --out daphnia_FIS

