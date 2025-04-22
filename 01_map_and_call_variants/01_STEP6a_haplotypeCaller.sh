#!/usr/bin/env bash

#SBATCH --job-name=recal_bases
#SBATCH -A beagle
#SBATCH -t 180:00:00
#SBATCH -n 5
#SBATCH --output="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling/logs/variant_call_for_recal_logs/slurm-%A_%a.out"
#SBATCH --error="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling/logs/variant_call_for_recal_logs/slurm-%A_%a.err"
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory = "$SLURM_SUBMIT_DIR


MASTER_DIR=$1
REFPATH=$2
REFNAME=$3
READ_DIR=$4
SUFFIX1=$5
SUFFIX2=$6


#load required modules
module load biocontainers
module load gatk4/4.3.0.0
ulimit -s unlimited

#search all the fastq files from the "data" directory and generate the array
        index=$(( $SLURM_ARRAY_TASK_ID + 1 ))
        read1=$(ls ${READ_DIR}/*${SUFFIX1} | sed -n ${index}p)
        read2=$(ls ${READ_DIR}/*${SUFFIX2} | sed -n ${index}p)
        prefix=${read1%"$SUFFIX1"} # get file name prefix
        ID=${prefix#"${READ_DIR}/"}


echo "sample: ${ID}"

cd ${MASTER_DIR}
mkdir -p ${MASTER_DIR}/variants_for_recal

#The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. 

gatk HaplotypeCaller -R ${REFPATH}/${REFNAME} -I ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq20_markedDup_split_readGroups_split.bam -ERC GVCF -O ${MASTER_DIR}/variants_for_recal/${ID}_first_pass.g.vcf.gz



