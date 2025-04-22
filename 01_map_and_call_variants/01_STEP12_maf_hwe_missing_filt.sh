#!/usr/bin/env bash

#SBATCH --job-name=maf_hwe_missing
#SBATCH -A beagle
#SBATCH -t 180:00:00
#SBATCH -n 20
#SBATCH --output="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling_NO_BQSR/logs/maf_VCF_logs/slurm-%A_%a.out"
#SBATCH --error="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling_NO_BQSR/logs/maf_VCF_logs/slurm-%A_%a.err"



MASTER_DIR=$1
REFPATH=$2
REFNAME=$3
MAF=$4
HWE=$5
GENO=$6
MIND=$7


#load required modules
module load biocontainers
module load plink2/2.00a2.3
ulimit -s unlimited

cd ${MASTER_DIR}/combined_variants

plink2 --bfile all_samples_NO_BQSR_SNPS_filterOut_clean.vcf --maf ${MAF} --hwe ${HWE} --geno ${GENO} --allow-no-sex --allow-extra-chr --make-bed --out temp
plink2 --bfile temp --mind ${MIND} --allow-no-sex --make-bed --allow-extra-chr --out all_samples_NO_BQSR_SNPS_filterOut_clean_maf${MAF}_hwe${HWE}_geno${GENO}_mind${MIND}

mv *.log /${MASTER_DIR}/logs/maf_VCF_logs/

rm ./temp*

