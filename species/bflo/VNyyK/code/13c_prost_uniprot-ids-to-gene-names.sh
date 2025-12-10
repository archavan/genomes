#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=fetch_uniprot
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8GB
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load R

cd /home/arc78/genomes # must start at project root to use renv

Rscript --no-save species/bflo/VNyyK/code/13c_prost_uniprot-ids-to-gene-names.R
