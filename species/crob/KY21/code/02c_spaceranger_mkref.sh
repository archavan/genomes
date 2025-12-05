#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=mkref
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu


cd /home/arc78/genomes/species/crob/KY21

gunzip results/fasta/genome_with_mito.fasta.gz
gunzip results/gtf/ky21_fixed_minimal.gtf.gz

cd results

~/applications/spaceranger-4.0.1/spaceranger mkref \
	--genome=spaceranger-4-0-1 \
	--fasta=fasta/genome_with_mito.fasta \
	--genes=gtf/ky21_fixed_minimal.gtf \
	--nthreads=10
	
gzip fasta/genome_with_mito.fasta
gzip gtf/ky21_fixed_minimal.gtf
