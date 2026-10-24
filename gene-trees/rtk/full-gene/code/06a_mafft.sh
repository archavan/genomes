#!/usr/bin/env bash
#SBATCH --partition=day
#SBATCH --job-name=mafft
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load MAFFT

cd ~/genomes/gene-trees/rtk/full-gene

mkdir -p results/06a_mafft

mafft \
--thread 8 \
--genafpair \
--maxiterate 1000 \
--reorder \
results/05c_seqs-qc/seqs_passed.fa > results/06a_mafft/aligned.fa
