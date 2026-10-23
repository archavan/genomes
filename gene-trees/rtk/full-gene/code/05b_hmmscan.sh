#!/usr/bin/env bash
set -euo pipefail

BASEDIR="$HOME/Dropbox/research/projects/genomes/gene-trees/rtk/full-gene"
cd "$BASEDIR"
mkdir -p results/05b_hmmscan

hmmdir="$BASEDIR/data/hmm"

cat "$hmmdir"/*.hmm > results/05b_hmmscan/models.hmm
hmmpress -f results/05b_hmmscan/models.hmm

hmmscan \
  --cpu 8 \
  --cut_ga \
  --domtblout results/05b_hmmscan/domains.domtbl \
  results/05b_hmmscan/models.hmm \
  results/05a_seqs-assembled/seqs_assembled.fa > results/05b_hmmscan/hmmscan.out
