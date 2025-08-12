# Arun Chavan
# Started: 2025-08-11

# Background ==================================================================

# Here we will fetch annotations for Styela genes from KEGG.

# Organism id for *Styela clava* on KEGG is T07291. The genome is at
# the following link: https://www.genome.jp/entry/gn:T07291

# KEGG hosts gene annotations for the current refseq assembly of Styela
# GCF_013122585.1, where they annotate the genes using the their
# pipelines KOALA and KoAnn. The main goal of these annotations is to assign a
# Kegg Orthology (KO; K numbers) group to the genes.

# For more on annotation see:
# KEGG genomes database: https://www.genome.jp/kegg/genome.html
# KEGG genes database: https://www.genome.jp/kegg/genes.html
# KEGG annotation stats: https://www.genome.jp/kegg/docs/genes_statistics.html

# For each refseq gene id, KEGG also creates its own gene id. In case of Styela,
# the refseq gene ids are in the format "LOC118403006", starting with "LOC1184"
# followed by 5 digits. The corresponding KEGG ids are the numerical part. For
# each of the KEGG gene ids, we can download the gene names in other databases,
# the Kegg Orthology group that the gene belongs to (the most useful to us),
# etc.

# packages ====================================================================
library(tidyverse)
library(here)
library(glue)
library(KEGGREST)

# setup =======================================================================
genomedir <- here("species/scla/refseq_v1")
genomeid <- "T07291"

# get kegg gene ids ===========================================================
# kegglist returns a named list where names are kegg gene ids and values are
# refseq gene descriptions
genes <- keggList(genomeid) %>%
  enframe() %>%
  rename(kegg_gene_id = name, refseq_description = value)

# fetch gene info =============================================================

fetch_gene <- function(gene_id) {
  message(glue("Fetching: {gene_id}"))
  keggGet(dbentries = gene_id)[[1]]
}

insistent_fetch_gene <- purrr::insistently(
  fetch_gene,
  rate = rate_delay(pause = 5, max_times = 5)
)

safe_insistent_fetch_gene <- purrr::safely(insistent_fetch_gene)

genes$out <- map(genes$kegg_gene_id, safe_insistent_fetch_gene)
genes$get <- map(genes$out, "result")
genes$errors <- map(genes$out, "error")
genes$out <- NULL

final_errors <- purrr::compact(genes$errors)

if (length(final_errors) > 0) {
  message(glue("Errors encountered for {length(final_errors)} genes."))
}

fs::dir_create(here(genomedir, "results/kegg"))
write_rds(genes, here(genomedir, "results/kegg/kegg_genes.rds"))

# end =========================================================================
