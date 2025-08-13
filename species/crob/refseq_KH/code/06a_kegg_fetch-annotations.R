# Arun Chavan
# Started: 2022-08-12

# Background ==================================================================

# Here we will fetch annotations for ciona genes from KEGG.

# details:
#   link to genome: https://www.genome.jp/entry/gn:T01056
#   Ciona intestinalis organism id on kegg gneome: T01056

# KEGG hosts gene annotations for the current refseq assembly of ciona
# CGF_000224145.3_KH, where they annotate the genes using the their
# pipelines KOALA and KoAnn. The main goal of these annotations is to assign a
# Kegg Orthology (KO; K numbers) group to the genes.

# For more on annotation see:
# KEGG genomes database: https://www.genome.jp/kegg/genome.html
# KEGG genes database: https://www.genome.jp/kegg/genes.html
# KEGG annotation stats: https://www.genome.jp/kegg/docs/genes_statistics.html

# For each refseq gene id, KEGG also creates its own gene id. In case of ciona,
# the refseq gene ids are in the format "100176615". The corresponding KEGG ids
# are the numerical part prefixed with the organism name short-hand,
# "cin:100176615". For each of the KEGG gene ids, we can download the gene names
# in other databases, the Kegg Orthology group that the gene belongs to (the
# most useful to us), etc.

# packages ====================================================================
library(tidyverse)
library(here)
library(glue)
library(KEGGREST)

# setup =======================================================================
genomedir <- here("species/crob/refseq_KH/")
genomeid <- "T01056"

# get kegg gene ids ===========================================================
# kegglist returns a named list where names are kegg gene ids and values are
# refseq gene descriptions
genes <- keggList(genomeid) %>%
  enframe() %>%
  rename(kegg_gene_id = name, refseq_description = value)

# fetch gene info =============================================================

fetch_gene_info <- function(gene_id) {
  message(glue("trying {gene_id}"))
  tryCatch(
    keggGet(dbentries = gene_id)[[1]],
    error = function(cond) {
      message(glue("{gene_id} failed: {cond}"))
      return(NA)
    }
  )
}

genes$get <- map(genes$kegg_gene_id, fetch_gene_info)

fs::dir_create(here(genomedir, "results/kegg"))
write_rds(genes, here(genomedir, "results/kegg/kegg_genes.rds"))

# end =========================================================================
