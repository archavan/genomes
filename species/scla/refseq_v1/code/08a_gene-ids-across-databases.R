#' ---
#' title: "gene ids across databases"
#' author: "Arun Chavan"
#' output: github_document
#' ---
#'
#' Started: 2025-08-12
#'

# background ==================================================================

# Here we put together a table mapping corresponding ids from other databases
# for all ncbi gene ids

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/scla/refseq_v1")

# data ========================================================================
kegg <- read_tsv(here(genomedir, "results/kegg/kegg-genes.tsv"))
ncbi <- read_tsv(here(genomedir, "results/genes/ncbi-genes-table.tsv"))

# consolidate =================================================================
# some databases may have 1:n mapping with ncbi gene id

ids_xdb <- left_join(
  x = select(ncbi, gene_id),
  y = select(kegg, gene_id, kegg_gene_id, ko_id),
  by = "gene_id"
) |>
  distinct()

# write files =================================================================
write_tsv(
  ids_xdb,
  here(genomedir, "results/genes/gene-ids-across-databases.tsv")
)

# end =========================================================================
