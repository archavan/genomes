#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=mkref
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu



cd /home/arc78/genomes/species/bflo/VNyyK/results

~/applications/cellranger-7.1.0/cellranger mkref \
	--genome=cellranger-7.1.0 \
	--fasta=../data/refseq/GCF_000003815.2_Bfl_VNyyK_genomic.fna \
	--genes=gtf/gtf_fixed.gtf \
	--nthreads=10
	
