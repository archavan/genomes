#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=gff_to_gtf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load miniconda
conda activate agat

cd /home/arc78/endoderm/species/crob/KY21/data/ghost-db
gunzip HT.KY21Gene.2.gff3.gz

cd /home/arc78/endoderm/species/crob/KY21

agat_convert_sp_gff2gtf.pl \
--gff data/ghost-db/HT.KY21Gene.2.gff3 \
--out results/gtf/HT.KY21Gene.2_agat.gtf \
1>code/01b_convert-gff-to-gtf_agat.log

gzip results/gtf/HT.KY21Gene.2_agat.gtf

gzip data/ghost-db/HT.KY21Gene.2.gff3

conda deactivate

