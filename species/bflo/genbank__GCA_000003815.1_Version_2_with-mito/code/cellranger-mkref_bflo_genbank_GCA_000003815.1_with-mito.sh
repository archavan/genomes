#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=mkref_genbank_with_mito
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

cd /home/arc78/macrophages/species/bflo/genome/genbank/GCA_000003815.1_Version_2_with-mito

/home/arc78/applications/cellranger-3.1.0/\
cellranger mkref \
	--genome=cellranger_bflo_genbank_GCA_000003815.1_with-mito \
	--fasta=appended/GCA_000003815.1_Version_2_genomic_with-mito.fasta \
	--genes=appended/GCA_000003815.1_Version_2_genomic_with-mito.gtf \
	--nthreads=10
	