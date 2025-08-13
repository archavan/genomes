#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=mkref_genbank
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

cd /home/arc78/scratch60/macrophages/bflo/genome/genbank

/home/arc78/applications/cellranger-3.1.0/\
cellranger mkref \
	--genome=cellranger_bflo_genbank \
	--fasta=GCA_000003815.1_Version_2/GCA_000003815.1_Version_2_genomic.fna \
	--genes=GCA_000003815.1_Version_2/GCA_000003815.1_Version_2_genomic.gtf \
	--nthreads=10
	