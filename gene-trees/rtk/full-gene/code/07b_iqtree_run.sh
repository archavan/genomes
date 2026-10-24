#!/usr/bin/env bash
#SBATCH --partition=day
#SBATCH --job-name=iqtree
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G
#SBATCH --time=0-06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load miniconda

conda activate iqtree

cd ~/genomes/gene-trees/rtk/full-gene
mkdir -p results/07b_iqtree

iqtree \
-s results/06b_trim/trimmed.fa \
-m MFP \
-B 1000 -bnni \
-T AUTO -ntmax 16 \
--prefix results/07b_iqtree/trimmed
