# Arun Chavan
# Started: 2025-12-10

# background ==================================================================

# Calculate lengths of all proteins.

# setup =======================================================================
library(seqinr)
library(tidyverse)
library(here)

genomedir <- here("species/bflo/VNyyK")
resdir <- here(genomedir, "results/genes")

# data ========================================================================
gid_pid <- read_tsv(here(genomedir, "results/genes/gid-pid_mapping.tsv")) |>
  select(gene_id, protein_id)

prot <- read.fasta(
  here(genomedir, "data/refseq/GCF_000003815.2_Bfl_VNyyK_protein.faa"),
  seqtype = "AA"
)

# calculate lengths ===========================================================
protein_length <- enframe(
  lengths(prot),
  name = "protein_id",
  value = "length"
) |>
  left_join(gid_pid) |>
  relocate(gene_id) |>
  arrange(gene_id)

# write =======================================================================
write_tsv(protein_length, here(resdir, "protein-lengths.tsv"))

# end =========================================================================
