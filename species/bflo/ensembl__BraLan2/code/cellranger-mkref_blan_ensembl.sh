#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=mkref_ensembl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

cd /home/arc78/scratch60/macrophages/bflo/genome/ensembl_BraLan2

/home/arc78/applications/cellranger-3.1.0/\
cellranger mkref \
	--genome=cellranger_blan_ensembl \
	--fasta=Branchiostoma_lanceolatum.BraLan2.dna.toplevel.fa \
	--genes=Branchiostoma_lanceolatum.BraLan2.46.gtf \
	--nthreads=10
	