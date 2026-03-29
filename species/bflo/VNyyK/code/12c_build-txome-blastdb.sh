#!/usr/bin/env bash
# Build a BLAST nucleotide database from the transcriptome FASTA.

set -euo pipefail

BASEDIR="$(cd "$(dirname "$0")/.." && pwd)"
FASTA_DIR="${BASEDIR}/data/refseq"
FASTA_BASE="GCF_000003815.2_Bfl_VNyyK_rna"
FASTA="${FASTA_DIR}/${FASTA_BASE}.fna.gz"
DB_PATH="${BASEDIR}/results/blastdb/${FASTA_BASE}/db"

mkdir -p "$(dirname "$DB_PATH")"

gunzip -c "$FASTA" |\
  makeblastdb -in - -dbtype nucl \
    -title "$FASTA_BASE" \
    -out "$DB_PATH"
