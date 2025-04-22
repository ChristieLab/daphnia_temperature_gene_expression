#!/usr/bin/env bash

#SBATCH --job-name=call_variants
#SBATCH -A beagle
#SBATCH -t 180:00:00
#SBATCH -n 5
#SBATCH --output="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling/logs/combine_and_genotype_variant_logs/slurm-%A_%a.out"
#SBATCH --error="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling/logs/combine_and_genotype_variant_logs/slurm-%A_%a.err"



MASTER_DIR=$1
REFPATH=$2
REFNAME=$3



#load required modules
module load biocontainers
module load gatk4/4.3.0.0
ulimit -s unlimited

#create files for combining and genotyping
cd ${MASTER_DIR}
cd ${MASTER_DIR}/variants
readlink -f *.vcf.gz > list.txt
sed 's/^/--variant\t/g' list.txt > arguments.txt
rm ./list.txt
cd ${MASTER_DIR}/reference
cat GCF_021134715_Dpulex.fna | grep '^>' > temp.txt
awk '{print $1}' temp.txt > tempII.txt
sed 's|[>,]||g' tempII.txt > intervals.list
rm ./temp.txt
cd ${MASTER_DIR}

#merge all vcfs
gatk CombineGVCFs -R ${REFPATH}/${REFNAME} --arguments_file ${MASTER_DIR}/variants/arguments.txt -O ${MASTER_DIR}/combined_variants/all_samples_combined.vcf.gz

#genotype VCFs
gatk GenotypeGVCFs -R ${REFPATH}/${REFNAME} --intervals ${MASTER_DIR}/reference/intervals.list -V ${MASTER_DIR}/combined_variants/all_samples_combined.vcf.gz -O ${MASTER_DIR}/combined_variants/all_samples_genotyped_raw.vcf.gz

#validate variants
gatk ValidateVariants -R ${REFPATH}/${REFNAME} -V ${MASTER_DIR}/combined_variants/all_samples_genotyped_raw.vcf.gz--validation-type-to-exclude ALL