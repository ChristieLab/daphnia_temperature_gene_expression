#!/usr/bin/env bash

#SBATCH --job-name=pairwise_clone_WC_Fst
#SBATCH -A beagle
#SBATCH -t 2:00:00
#SBATCH -n 1
#SBATCH --output="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/pairwise.out"
#SBATCH --error="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/pairwise.err"

DIR="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/Fst_clone"
GVCF="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/Fst_clone/all_samples_NO_BQSR_NO_CBSWP_SNPS_filterOut_clean_maf0.05_hwe0_geno0.1_mind0.2_recode_4.2.vcf"
POPLIST="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/Fst_clone/temp"
POPREV="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/Fst_clone/tempII"
CLONES="/scratch/bell/nbackens/daphnia_RNAseq/updated_figures_02-23-25/Fst_clone/clone_lists/"

module load biocontainers
module load vcftools/0.1.16

cd ${DIR}
P1=`cat ${POPLIST}`

for f in $P1
do

P2=`cat $POPREV`

for p in $P2
do
mkdir -p ${DIR}/${f}_${p}_Fst
cd ${DIR}/${f}_${p}_Fst
vcftools --vcf ${GVCF} --weir-fst-pop ${CLONES}/clone_list_${f}.txt --weir-fst-pop ${CLONES}/clone_list_${p}.txt --keep ${CLONES}/clone_list_${f}.txt --keep ${CLONES}/clone_list_${p}.txt --out ${f}_${p}_Fst
cd ${DIR}
done
cd ${DIR}
done





