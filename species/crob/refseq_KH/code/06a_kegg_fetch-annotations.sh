#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --job-name=fetch_kegg
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load R/4.0.3-foss-2020b

cd /home/arc78/endoderm

Rscript --no-save species/crob/genome/refseq_KH/code/06a_kegg_fetch-annotations.R
