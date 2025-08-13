#!/usr/bin/env bash
#SBATCH --partition=general
#SBATCH --job-name=mkref
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu



cd /home/arc78/scratch60/endoderm/species/crob/genome/KH13/data/ghost-db

unzip JoinedScaffold.zip
mv JoinedScaffold JoinedScaffold.fasta

~/applications/cellranger-6.1.2/cellranger mkref \
	--genome=KH13.cellranger.6.1.2 \
	--fasta=JoinedScaffold.fasta \
	--genes=KH.KHGene.2013_agat-converted_only-gid-tid.gtf \
	--nthreads=10
	