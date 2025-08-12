# Arun Chavan
# 2025-08-11

# background ==================================================================

# Create gene-id and transcript-id / protein-id mapping table, extract a list of
# mitochondrial genes.

# setup =======================================================================
library(tidyverse)
library(here)

genomedir <- here("species/scla/refseq_v1")

read_gxf_as_tibble <- function(.file) {
  rtracklayer::import(.file) |>
    GenomicRanges::as.data.frame() |>
    as_tibble()
}

mt_chr <- "NC_037072.1"

# read data ===================================================================
gtf <- read_gxf_as_tibble(here(genomedir, "results/gtf/gtf_fixed.gtf.gz"))

# gene id and protein id mapping ==============================================
gid_pid <- gtf |>
  select(gene_id, protein_id) |>
  filter(!is.na(protein_id)) |>
  distinct()

write_tsv(gid_pid, here(genomedir, "results/genes/gid-pid_mapping.tsv"))

# gene id and transcript id mapping ===========================================
gid_tid <- gtf %>%
  select(gene_id, transcript_id) |>
  filter(!is.na(transcript_id)) |>
  filter(transcript_id != "") |>
  distinct()

write_tsv(gid_tid, here(genomedir, "results/genes/gid-tid_mapping.tsv"))

# Extract genesets ============================================================
## Mitochondrial genes --------------------------------------------------------
mt_genes <- gtf |>
  filter(seqnames == mt_chr) |>
  select(seqnames, gene_id) |>
  distinct()

write_tsv(mt_genes, here(genomedir, "results/genes/geneset_mitochondrial.tsv"))

# end =========================================================================
