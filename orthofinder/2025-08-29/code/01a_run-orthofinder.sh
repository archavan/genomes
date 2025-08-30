#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=orthofinder
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36 
#SBATCH --mem-per-cpu=3G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load miniconda
conda activate orthofinder3

cd ~/palmer_scratch/genomes/orthofinder/2025-08-29

orthofinder \
-s ./data/species-tree_abb.txt  \
-o ./results/out \
-f ./results/input-proteomes/
