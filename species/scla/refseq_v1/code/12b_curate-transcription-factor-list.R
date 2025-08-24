# Arun Chavan
# Started: 2025-08-24

# setup =======================================================================
library(tidyverse)
library(here)

genomedir <- here("species/scla/refseq_v1")
resdir <- here(genomedir, "results/genes")
orthdir <- here("orthofinder/2025-08-13")

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

orth <- read_tsv(here(orthdir, "results/pairwise-orthologs/hsap_vs_scla.tsv"))
hs_names <- read_tsv(here(orthdir, "results/id-name-mapping/hsap.tsv"))
orth <- orth |>
  left_join(
    y = (hs_names |>
      select(gene_id, symbol) |>
      distinct()),
    by = c("hsap" = "gene_id")
  )

# homologs of human transcription factors =====================================
tf_lambert_homologs <- orth %>%
  filter(symbol %in% tf_hs_lambert$hgnc_symbol) |>
  rename(gene_id = scla)

# curate the list =============================================================

## Intersection of all sets
tf_intersect <- list(
  select(tf_lambert_homologs, gene_id),
  select(tf_emapper, gene_id),
  select(tf_kegg, gene_id)
) %>%
  reduce(inner_join, by = "gene_id") %>%
  left_join(genes, by = "gene_id") %>%
  distinct()

## Union of all sets
tf_union <- list(
  select(tf_lambert_homologs, gene_id),
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
