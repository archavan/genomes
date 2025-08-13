# Arun Chavan
# Started: 2024-12-20

# setup =======================================================================
library(tidyverse)
library(here)

genomedir <- here("species/bflo/VNyyK")
resdir <- here(genomedir, "results/genes")

# data ========================================================================
# tf genes
tf_kegg <- read_tsv(here(genomedir, "results/kegg/kegg_geneset_tf.tsv"))
tf_emapper <- read_tsv(here(
  genomedir,
  "results/emapper/processed/emapper_geneset_tf.tsv"
))

# human transcription factor list from Lambert et al 2018
tf_hs_lambert <- read_csv(
  "https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"
) %>%
  janitor::clean_names() %>%
  filter(is_tf == "Yes")

genes <- read_tsv(here(resdir, "gene-display-names.tsv")) %>%
  select(-display_name)

orth <- read_csv(here(
  genomedir,
  "results/emapper/processed/emapper_orthologs_hsap.csv"
))

# homologs of human transcription factors =====================================
tf_bf_lambert_homologs <- orth %>%
  filter(hs_symbol %in% tf_hs_lambert$hgnc_symbol) %>%
  rename(gene_id = bf_gene_id)

# curate the list =============================================================

## Intersection of all sets
tf_intersect <- list(
  select(tf_bf_lambert_homologs, gene_id),
  select(tf_emapper, gene_id),
  select(tf_kegg, gene_id)
) %>%
  reduce(inner_join, by = "gene_id") %>%
  left_join(genes, by = "gene_id") %>%
  distinct()

## Union of all sets
tf_union <- list(
  select(tf_bf_lambert_homologs, gene_id),
  select(tf_emapper, gene_id),
  select(tf_kegg, gene_id)
) %>%
  reduce(full_join, by = "gene_id") %>%
  left_join(genes, by = "gene_id") %>%
  distinct()

# Write files =================================================================
write_tsv(tf_intersect, here(resdir, "geneset_tf_intersect.tsv"))
write_tsv(tf_union, here(resdir, "geneset_tf_union.tsv"))

# end =========================================================================
