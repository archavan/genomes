# Arun Chavan
# Started: 2025-08-24

# setup =======================================================================
library(tidyverse)
library(here)
library(httr2)
library(jsonlite)
library(glue)

genomedir <- here("species/scla/refseq_v1")

gene_ids <- read_tsv(here(
  genomedir,
  "results/genes/gene-ids-across-databases.tsv"
))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv")) %>%
  select(-display_name)

# data ========================================================================
# fetch_kegg_json <- function(.accession) {
#   url <- glue("https://rest.kegg.jp/get/{.accession}/json")
#   request(url) %>%
#     req_perform() %>%
#     resp_body_json()
# }
# tf_json <- fetch_kegg_json("br:bfo03000")
# nr_json <- fetch_kegg_json("br:bfo03310")

kegg <- read_rds(here(genomedir, "results/kegg/kegg_genes.rds"))

# process =====================================================================
get_genes_with_brite_acccession <- function(.brite) {
  kegg %>%
    mutate(brite = map(get, ~ .[["BRITE"]])) %>%
    unnest(brite) %>%
    filter(str_detect(brite, .brite)) %>%
    select(kegg_gene_id) %>%
    distinct() %>%
    left_join(gene_ids, by = "kegg_gene_id") %>%
    select(gene_id) %>%
    left_join(genes, by = "gene_id")
}


tf_brite <- get_genes_with_brite_acccession("sclv03000")

# write =======================================================================
write_tsv(tf_brite, here(genomedir, "results/kegg/kegg_geneset_tf.tsv"))

# end =========================================================================
