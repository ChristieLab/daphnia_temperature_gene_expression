#!/usr/bin/env bash

#SBATCH --job-name=daphnia_mito_PCA
#SBATCH -A beagle
#SBATCH -t 2:00:00
#SBATCH -n 1
#SBATCH --output="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/scripts/PCA.out"
#SBATCH --error="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/scripts/PCA.err"

DIR="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/PCA"
VCF="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/PCA/all_samples_NO_BQSR_NO_CBSWP_SNPS_filterOut_clean_maf0.05_hwe0_geno0.1_mind0.2_recode_4.2.vcf"

module load biocontainers
module load vcftools/0.1.16
module load plink2


cd ${DIR}

plink2 --make-bed --vcf ${VCF} --allow-extra-chr --out all_samples_NO_BQSR_NO_CBSWP_SNPS_filterOut_clean_maf0.05_hwe0_geno0.1_mind0.2_recode_4.2
plink2 --make-bpgen --vcf ${VCF} --allow-extra-chr --out all_samples_NO_BQSR_NO_CBSWP_SNPS_filterOut_clean_maf0.05_hwe0_geno0.1_mind0.2_recode_4.2
plink2 --bfile all_samples_NO_BQSR_NO_CBSWP_SNPS_filterOut_clean_maf0.05_hwe0_geno0.1_mind0.2_recode_4.2 --allow-extra-chr --pca