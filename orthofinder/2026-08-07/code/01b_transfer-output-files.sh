#!/usr/bin/env bash

mkdir -p ~/Dropbox/research/projects/genomes/orthofinder/2026-08-07/results/out

rsync -avz \
--exclude "Gene_Trees" \
--exclude "MultipleSequenceAlignments" \
--exclude "Orthogroup_Sequences" \
--exclude "WorkingDirectory" \
--exclude "Single_Copy_Orthologue_Sequences" \
arc78@mccleary.ycrc.yale.edu:/home/arc78/palmer_scratch/genomes/orthofinder/2026-08-07/results/out/Results_Aug08/ \
~/Dropbox/research/projects/genomes/orthofinder/2026-08-07/results/out/

# Gene_Trees is excluded in favor of Resolved_Gene_Trees. 
