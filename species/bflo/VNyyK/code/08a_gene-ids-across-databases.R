#' ---
#' title: "gene ids across databases"
#' author: "Arun Chavan"
#' output: github_document
#' ---
#'
#' Started: 2023-06-08
#'

# background ==================================================================

# At this point we have ncbi gene ids mapped to gene ids in various databases:
# 1. KEGG genes
# 3. Ensembl ids from KEGG and NCBI (we can also get these from biomart)

# Here we put together a table mapping corresponding ids from other databases
# for all ncbi gene ids

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/bflo/VNyyK/")

# data ========================================================================
genes <- read_tsv(here(genomedir, "results/genes/gene-ids.tsv"))
kegg <- read_tsv(here(genomedir, "results/kegg/kegg-genes.tsv"))
ncbi <- read_tsv(here(genomedir, "results/genes/ncbi-genes-table.tsv"))

# consolidate =================================================================
# some databases have 1:n mapping with ncbi gene id
ids_xdb <- genes %>%
  select(gene_id, ncbi_gene_id) %>%
  left_join(
    (kegg %>%
      select(ncbi_gene_id, kegg_gene_id, ko_id)),
    by = "ncbi_gene_id"
  ) %>%
  distinct()

# write files =================================================================
write_tsv(
  ids_xdb,
  here(genomedir, "results/genes/gene-ids-across-databases.tsv")
)

# end =========================================================================
