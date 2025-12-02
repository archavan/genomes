#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_index
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

cd /home/arc78/genomes/species/mmus/ensembl/v115/results

wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

kallisto index -i Mus_musculus.GRCm39.cdna.all.kallisto.idx Mus_musculus.GRCm39.cdna.all.fa.gz
