#!/bin/sh -l

#SBATCH --job-name=merge_markDups
#SBATCH -A beagle
#SBATCH -t 10:00:00
#SBATCH -n 10
#SBATCH --output="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling/markDups_logs/slurm-%A_%a.out"
#SBATCH --error="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling/markDups_logs/slurm-%A_%a.err"
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
MAPQUAL=$3
READ_DIR=$4
SUFFIX1=$5
SUFFIX2=$6

#load required modules

source /etc/profile.d/modules.sh
module load biocontainers/default
module load sambamba/0.8.2
module load samtools/1.17
module load picard/2.26.10
ulimit -s unlimited


cd ${MASTER_DIR}
mkdir -p test_test


