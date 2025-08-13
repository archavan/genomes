#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=mkref
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu



cd ~/palmer_scratch/endoderm/species/crob/genome/refseq_KH

gunzip data/refseq/GCF_000224145.3_KH_genomic.fna.gz

cd results

~/applications/cellranger-7.1.0/cellranger mkref \
	--genome=cellranger-7.1.0 \
	--fasta=../data/refseq/GCF_000224145.3_KH_genomic.fna \
	--genes=gtf/gtf_fixed.gtf \
	--nthreads=10
	
	
