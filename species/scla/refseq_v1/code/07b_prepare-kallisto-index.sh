#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_index
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

cd /home/arc78/genomes/species/scla/refseq_v1

mkdir -p results/kallisto

kallisto index -i results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
data/refseq/GCF_013122585.1_ASM1312258v2_rna.fna.gz

