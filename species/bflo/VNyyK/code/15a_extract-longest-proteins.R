# Arun Chavan
# Started: 2025-03-06

# background ==================================================================

# Prepare proteome file for orthofinder following https://davidemms.github.io/orthofinder_tutorials/orthofinder-best-practices.html
# 1. get the longest / primary protein for each gene
# 2. rename the fasta entries to gene ids rather than protein ids

# setup =======================================================================
library(tidyverse)
library(here)
library(seqinr)

genomedir <- here("species/bflo/VNyyK")

# data ========================================================================
fa <- read.fasta(
  here(genomedir, "data/refseq/GCF_000003815.2_Bfl_VNyyK_protein.faa"),
  seqtype = "AA"
)
ids <- read_tsv(here(genomedir, "results/genes/gid-pid_mapping.tsv"))

# process =====================================================================
longest <- ids %>%
  mutate(len = map_int(protein_id, ~ length(fa[[.x]]))) %>%
  group_by(gene_id) %>%
  slice_max(order_by = len, n = 1, with_ties = FALSE) %>%
  ungroup()

fa_longest <- fa[longest$protein_id]
names(fa_longest) <- ids$gene_id[match(names(fa_longest), ids$protein_id)]

# orthofinder replaces `:` in gene ids with `_`, so we use use `ncbi.` instead
# of `ncbi:`
names(fa_longest) <- gsub("ncbi:", "ncbi.", names(fa_longest))

# write =======================================================================
write.fasta(
  fa_longest,
  names(fa_longest),
  here(genomedir, "results/genes/longest-proteins.fa")
)

# end =========================================================================
