#!/usr/bin/env bash

#SBATCH --job-name=clean_vcf
#SBATCH -A beagle
#SBATCH -t 180:00:00
#SBATCH -n 20
#SBATCH --output="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling_NO_BQSR/logs/clean_VCF_logs/slurm-%A_%a.out"
#SBATCH --error="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling_NO_BQSR/logs/clean_VCF_logs/slurm-%A_%a.err"



MASTER_DIR=$1
REFPATH=$2
REFNAME=$3



#load required modules
module load biocontainers
module load gatk4/4.3.0.0
module load plink2/2.00a2.3
module load r/4.3.1
ulimit -s unlimited

cd ${MASTER_DIR}/combined_variants

#make a PLINK binary data set
plink2 --vcf ./all_samples_NO_BQSR_SNPS_filterOut.vcf --make-bed --allow-extra-chr --max-alleles 2 -out all_samples_NO_BQSR_SNPS_filterOut.vcf

#copy the log file to logs
mv ./all_samples_NO_BQSR_SNPS_filterOut.vcf.log ${MASTER_DIR}/logs/clean_VCF_logs/

#make a backup copy of the bim file
cp ./all_samples_NO_BQSR_SNPS_filterOut.vcf.bim ./all_samples_NO_BQSR_SNPS_filterOut.vcf.bck.bim

#reformat to run PLINK
awk '{if($2 == ".") {print $1,"\t",$1":"$4,"\t",$3,"\t",$4,"\t",$5,"\t",$6} else {print $0}}' ./all_samples_NO_BQSR_SNPS_filterOut.vcf.bck.bim > all_samples_NO_BQSR_SNPS_filterOut.vcf.bim

#find duplicate snps
#remove variants that are not SNPs
#remove SNPS that are not A, C, T, G
#remove palindromic and strand ambiguous alleles

Rscript /scratch/bell/nbackens/daphnia_RNAseq/scripts/01_map_and_call_variants_NO_BQSR/01_subscripts_NO_BQSR/clean_data.R all_samples_NO_BQSR_SNPS_filterOut.vcf all_samples_NO_BQSR_SNPS_filterOut_clean.vcf

#run plink
plink2 --bfile all_samples_NO_BQSR_SNPS_filterOut.vcf --allow-extra-chr --extract qced_vars.txt --make-bed --out all_samples_NO_BQSR_SNPS_filterOut_clean.vcf
