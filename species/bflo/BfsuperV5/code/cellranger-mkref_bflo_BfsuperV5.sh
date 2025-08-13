#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=mkref_BfsuperV5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

cd /home/arc78/scratch60/macrophages/bflo/genome/BfsuperV5

/home/arc78/applications/cellranger-3.1.0/\
cellranger mkref \
	--genome=cellranger_bflo_BfsuperV5 \
	--fasta=files/BfsuperV5_fix_empty-line-fixed.fa \
	--genes=files/BfsuperV5-4_FIX_spacing-fixed.gtf \
	--nthreads=10
	
