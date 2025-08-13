#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_index
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=50G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load Python 
module load kallisto

cd /home/arc78/genomes/species/bflo/VNyyK

mkdir -p results/kallisto

kallisto index -i results/kallisto/GCF_000003815.2_Bfl_VNyyK_rna.idx \
data/refseq/GCF_000003815.2_Bfl_VNyyK_rna.fna.gz

