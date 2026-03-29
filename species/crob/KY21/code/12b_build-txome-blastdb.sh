#!/usr/bin/env bash
# Build a BLAST nucleotide database from the transcriptome FASTA.

set -euo pipefail

BASEDIR="$(cd "$(dirname "$0")/.." && pwd)"
FASTA_DIR="${BASEDIR}/results/fasta"
FASTA_BASE="KY21_transcripts_with-mito"
FASTA="${FASTA_DIR}/${FASTA_BASE}.fasta.gz"
DB_PATH="${BASEDIR}/results/blastdb/${FASTA_BASE}/db"

mkdir -p "$(dirname "$DB_PATH")"

gunzip -c "$FASTA" |\
  makeblastdb -in - -dbtype nucl \
    -title "$FASTA_BASE" \
    -out "$DB_PATH"
