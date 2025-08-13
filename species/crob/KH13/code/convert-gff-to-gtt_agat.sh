#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --job-name=convert-gff-to-gtf_agat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load miniconda
conda activate agat

cd /home/arc78/scratch60/endoderm/species/crob/genome/KH13/data/ghost-db

agat_convert_sp_gff2gtf.pl --gff KH.KHGene.2013.gff3 --out KH.KHGene.2013_agat-converted.gtf 2>convert-gff-to-gtf_agat.log

conda deactivate
