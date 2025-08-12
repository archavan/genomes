# Arun Chavan
# Started: 2025-08-11

# Background ==================================================================

# We have already downloaded kegg annotations for all genes using REST API and
# saved the results to an .rds file. The KEGG API returns the results for our
# queries in its flat file format described here:
# https://www.kegg.jp/kegg/docs/dbentry.html.

# Here we read in the .rds file and parse the KEGGREST::keggGet() output to
# extract annotations such as Kegg Orthology group, NCBI gene id, etc.

# setup =======================================================================
library(KEGGREST)
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/scla/refseq_v1")

# data ========================================================================
genes <- read_rds(here(genomedir, "results/kegg/kegg_genes.rds"))

# There are 20272 bflo genes in the KEGG genome database. Ensure that all genes
# were queried by our script.
all.equal(nrow(genes), 20272)

# Ensure that KEGGREST did not fail to get info for any genes. If it did, run
# fetch annotations script again but only for the failed genes (it takes several
# hours for it to run for all genes)
all.equal(genes |> filter(map_lgl(errors, ~ !is.null(.x))) |> nrow(), 0)

success <- genes |>
  filter(map_lgl(errors, ~ is.null(.x)))

# fetch failed genes ==========================================================
fetch_gene <- function(gene_id) {
  message(glue("Fetching: {gene_id}"))
  keggGet(dbentries = gene_id)[[1]]
}

insistent_fetch_gene <- purrr::insistently(
  fetch_gene,
  rate = rate_delay(pause = 5, max_times = 5)
)

safe_insistent_fetch_gene <- purrr::safely(insistent_fetch_gene)

failed <- genes |>
  filter(map_lgl(errors, ~ !is.null(.x))) |>
  select(-get, -errors)

failed$out <- map(failed$kegg_gene_id, safe_insistent_fetch_gene)
failed$get <- map(failed$out, "result")
failed$errors <- map(failed$out, "error")
failed$out <- NULL

# bind succeeded and now-succeeded genes genes into one data frame ============
genes <- bind_rows(success, failed)

all.equal(nrow(genes), 20272)
all.equal(genes |> filter(map_lgl(errors, ~ !is.null(.x))) |> nrow(), 0)

# kegg orthology group id =====================================================
genes$ko_id <- map_chr(
  genes$get,
  ~ ifelse(is.null(.x[["ORTHOLOGY"]]), NA, attr(.x[["ORTHOLOGY"]], "names"))
)

# kegg orthology group name ===================================================
genes$ko_name <- map_chr(
  genes$get,
  ~ ifelse(is.null(.x[["ORTHOLOGY"]]), NA, .x[["ORTHOLOGY"]])
)

# refseq_descriptions =========================================================
genes$refseq_description <- map_chr(
  genes$get,
  ~ ifelse(is.null(.x[["NAME"]]), NA, .x[["NAME"]])
)
genes$refseq_description <- str_remove(
  genes$refseq_description,
  "^\\(RefSeq\\) "
)

# db ids ======================================================================
extract_db_id <- function(getres, db) {
  ifelse(
    is.null(getres$DBLINKS),
    NA,
    getres$DBLINKS %>%
      .[str_detect(., db)] %>%
      str_remove(glue("{db}: "))
  )
}

genes <- genes |>
  mutate(
    ncbi_gene_id = map_chr(get, extract_db_id, "NCBI-GeneID")
    # ensembl_gene_id = map_chr(get, extract_db_id, "Ensembl"),
    # uniprot_id = map_chr(get, extract_db_id, "UniProt")
  )

genes <- genes |>
  mutate(gene_id = paste0("ncbi.", ncbi_gene_id)) |>
  select(-ncbi_gene_id)

# kegg pathways ===============================================================
kegg_pathways <- genes |>
  hoist(get, "PATHWAY") |>
  select(kegg_gene_id, PATHWAY) |>
  unnest_longer(col = PATHWAY) |>
  dplyr::rename(kegg_pathway_id = PATHWAY_id, kegg_pathway_name = PATHWAY) |>
  relocate(kegg_pathway_id, .after = kegg_gene_id) |>
  filter(!is.na(kegg_pathway_id) | !is.na(kegg_pathway_name))

# write =======================================================================
genes |>
  dplyr::select(
    gene_id,
    kegg_gene_id,
    ko_id,
    # ensembl_gene_id,
    # uniprot_id,
    refseq_description,
    ko_name
  ) |>
  write_tsv(here(genomedir, "results/kegg/kegg-genes.tsv"))

kegg_pathways |>
  write_tsv(here(genomedir, "results/kegg/kegg-pathways.tsv"))

# end =========================================================================
