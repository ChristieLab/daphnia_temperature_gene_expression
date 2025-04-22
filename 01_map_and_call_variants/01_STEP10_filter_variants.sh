#!/usr/bin/env bash

#SBATCH --job-name=filter_variants
#SBATCH -A beagle
#SBATCH -t 180:00:00
#SBATCH -n 20
#SBATCH --output="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling_NO_BQSR/logs/variant_filt_logs/slurm-%A_%a.out"
#SBATCH --error="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling_NO_BQSR/logs/variant_filt_logs/slurm-%A_%a.err"



MASTER_DIR=$1
REFPATH=$2
REFNAME=$3



#load required modules
module load biocontainers
module load gatk4/4.3.0.0
ulimit -s unlimited

#create files for combining and genotyping

cd ${MASTER_DIR}/combined_variants


gatk VariantFiltration \
      -R ${REFPATH}/${REFNAME} \
      -V ${MASTER_DIR}/combined_variants/all_samples_NO_BQSR_genotyped_raw.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QDfilterLessThan2" \
      --filter-expression "MQ < 40.0" --filter-name "MQfilterLessThan40" \
      --filter-expression "FS > 60.0" --filter-name "FSfilterMoreThan60" \
      --genotype-filter-expression "DP < 2" --genotype-filter-name "DP3" \
      --set-filtered-genotype-to-no-call true \
      -O all_samples_NO_BQSR_genotyped_filterKeep.vcf.gz
 
gatk SelectVariants \
      -R ${REFPATH}/${REFNAME} \
      -V ${MASTER_DIR}/combined_variants/all_samples_NO_BQSR_genotyped_filterKeep.vcf.gz \
      --exclude-filtered --exclude-non-variants \
      -O all_samples_NO_BQSR_filterOut.vcf

gatk SelectVariants \
     -R ${REFPATH}/${REFNAME} \
     -V ${MASTER_DIR}/combined_variants/all_samples_NO_BQSR_filterOut.vcf \
     --select-type-to-include SNP \
     -O all_samples_NO_BQSR_SNPS_filterOut.vcf
 
gatk SelectVariants \
     -R ${REFPATH}/${REFNAME} \
     -V ${MASTER_DIR}/combined_variants/all_samples_NO_BQSR_filterOut.vcf \
     --select-type-to-include INDEL \
     -O all_samples_NO_BQSR_INDEL_filterOut.vcf





