#!/bin/bash -l

#SBATCH --job-name=degs_locations_search
#SBATCH -A beagle
#SBATCH -t 72:00:00
#SBATCH -n 64
#SBATCH --output="degs_locations_search.out"
#SBATCH --error="degs_locations_search.err"


module load biocontainers/default
module load bedtools/2.31.0
module load eggnog-mapper/2.1.7

cd /scratch/bell/nbackens/daphnia_RNAseq/02_degs/degs_locations_NO_CBSWP
EGG_THREADS="64" #number of threads to use for eggnog-mapper



# .bed file consists of 3 tab delimited columns, chromosome name (excluding carrot and descriptors), start, stop
bedtools getfasta -fi /scratch/bell/nbackens/daphnia_RNAseq/reference/GCF_021134715_Dpulex.fna -bed 96.NO_CBSWP.genotype-clone.temp.bed -fo 96.NO_CBSWP.genotype-clone.temp.fa
bedtools getfasta -fi /scratch/bell/nbackens/daphnia_RNAseq/reference/GCF_021134715_Dpulex.fna -bed 168.NO_CBSWP.genotype-clone.temp.bed -fo 168.NO_CBSWP.genotype-clone.temp.fa


#annotate with eggNog emapper
mkdir -p 96_genotype-clone.temp_NO_CBSWP_emapper_search_default
emapper.py --cpu ${EGG_THREADS} -i 96.NO_CBSWP.genotype-clone.temp.fa --output_dir 96_genotype-clone.temp_NO_CBSWP_emapper_search_default -o 96.NO_CBSWP.genotype-clone.temp.annot --evalue 0.01 --excel --itype genome --genepred search --go_evidence non-electronic 
mkdir -p 168_genotype-clone.temp_NO_CBSWP_emapper_search_default
emapper.py --cpu ${EGG_THREADS} -i 168.NO_CBSWP.genotype-clone.temp.fa --output_dir 168_genotype-clone.temp_NO_CBSWP_emapper_search_default -o 168.NO_CBSWP.genotype-clone.temp.annot --evalue 0.01 --excel --itype genome --genepred search --go_evidence non-electronic 

