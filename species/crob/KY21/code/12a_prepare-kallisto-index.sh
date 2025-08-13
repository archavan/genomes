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

cd /home/arc78/endoderm/species/crob/KY21

mkdir -p results/kallisto

kallisto index -i results/kallisto/KY21_transcripts_with-mito.idx \
results/fasta/KY21_transcripts_with-mito.fasta.gz
