#!/usr/bin/env bash
#SBATCH --partition=day
#SBATCH --job-name=clipkit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

# clipkit gappy mode to remove only columns that are >90% gaps (-g 0.9)

module load miniconda
conda activate iqtree

cd ~/genomes/gene-trees/rtk/full-gene

mkdir -p results/06b_trim

clipkit results/06a_mafft/aligned.fa \
-m gappy \
-g 0.9 \
-o results/06b_trim/trimmed.fa \
-l
