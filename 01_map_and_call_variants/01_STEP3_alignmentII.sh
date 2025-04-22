#!/bin/sh -l

#SBATCH --job-name=RNAseq_STAR2nd-pass
#SBATCH -A beagle
#SBATCH -t 336:00:00
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH --output="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling/alignment_logsII/slurm-%A_%a.out"
#SBATCH --error="/scratch/bell/nbackens/daphnia_RNAseq/01_snp_calling/alignment_logsII/slurm-%A_%a.err"
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
READ_DIR=$3
GTF=$4
SUFFIX1=$5
SUFFIX2=$6

cd ${MASTER_DIR}

#load required modules
module load biocontainers
module load star/2.7.9a
module load gffread/0.12.7
ulimit -s unlimited

#search all the fastq files from the "data" directory and generate the array
        index=$(( $SLURM_ARRAY_TASK_ID + 1 ))
        read1=$(ls ${READ_DIR}/*${SUFFIX1} | sed -n ${index}p)
        read2=$(ls ${READ_DIR}/*${SUFFIX2} | sed -n ${index}p)
        prefix=${read1%"$SUFFIX1"} # get file name prefix
        ID=${prefix#"${READ_DIR}/"}


echo "sample: ${ID}"
echo "${read1}"
echo "${read2}"


mkdir -p ${MASTER_DIR}/alignment/2-pass
mkdir -p ${MASTER_DIR}/alignment/2-pass/${ID}
cd ${MASTER_DIR}/alignment/2-pass/${ID}

STAR --runThreadN 20 --runMode alignReads --genomeDir ${MASTER_DIR}/alignment/2-pass/genomeDir_2-pass --readFilesIn ${read1} ${read2} --readFilesCommand zcat --sjdbGTFfile ${GTF} --outReadsUnmapped Fastx --outFileNamePrefix ${MASTER_DIR}/alignment/2-pass/${ID}/${ID} 