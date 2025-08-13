# Arun Chavan
# Started: 2023-06-07

# Background ==================================================================

# Here we will fetch annotations for amphioxus genes from KEGG.

# Organism id for *Branchiostoma floridae* on KEGG is T01074. The genome is at
# the following link: https://www.genome.jp/dbget-bin/www_bget?gn:T01074

# KEGG hosts gene annotations for the current refseq assembly of amphioxus
# CGF_000003815.2_Bfl_VNyyK, where they annotate the genes using the their
# pipelines KOALA and KoAnn. The main goal of these annotations is to assign a
# Kegg Orthology (KO; K numbers) group to the genes.

# For more on annotation see:
# KEGG genomes database: https://www.genome.jp/kegg/genome.html
# KEGG genes database: https://www.genome.jp/kegg/genes.html
# KEGG annotation stats: https://www.genome.jp/kegg/docs/genes_statistics.html

# For each refseq gene id, KEGG also creates its own gene id. In case of
# amphioxus, the refseq gene ids are in the format "LOC118403006", starting with
# "LOC1184" followed by 5 digits. The corresponding KEGG ids are the numerical
# part prefixed with the organism name short-hand, "bfo:118403006". For each of
# the KEGG gene ids, we can download the gene names in other databases, the Kegg
# Orthology group that the gene belongs to (the most useful to us), etc.

# packages ====================================================================
library(tidyverse)
library(here)
library(glue)
library(KEGGREST)

# setup =======================================================================
genomedir <- here("species/bflo/VNyyK/")
genomeid <- "T01074"

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
