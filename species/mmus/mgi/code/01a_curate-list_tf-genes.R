###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(janitor)

basedir <- here("species/mmus/mgi")

###############################################################################
# data ========================================================================
###############################################################################

mmgenes <- read_csv(here(
  "species/mmus/ensembl/v115/results/mmus_gene-id_gene-name_transcript-id.csv"
)) |>
  select(ensembl_gene_id, external_gene_name) |>
  distinct()

mgi_tf <- read_tsv(here(basedir, "data/queries/mmus_tf-genes_GO0003700.tsv")) |>
  clean_names()

hs_tf <- read_csv(
  "https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"
) |>
  clean_names() |>
  filter(is_tf == "Yes")

###############################################################################
# curate list of TF genes ====================================================
###############################################################################

# get mouse orthologs of human TF genes
mmhs <- read_tsv(
  "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
) |>
  clean_names()

get_mm_orths <- function(x) {
  dbid <- mmhs$db_class_key[which(mmhs$symbol == x)]
  mmhs |>
    filter(db_class_key %in% dbid) |>
    filter(ncbi_taxon_id == 10090) |>
    pull(symbol)
}

tf_hsorth <- hs_tf |>
  select(ensembl_id, hgnc_symbol) |>
  distinct() |>
  mutate(mm_orths = map(hgnc_symbol, get_mm_orths)) |>
  unnest(mm_orths) |>
  select(mm_orths) |>
  distinct()

# combine GO defined TFs and those curated by Lambert

tf_all <- mmgenes |>
  filter(
    external_gene_name %in%
      mgi_tf$symbol |
      external_gene_name %in%
        tf_hsorth$mm_orths
  ) |>
  mutate(
    GO0003700 = if_else(external_gene_name %in% mgi_tf$symbol, "yes", "no"),
    Lambert = if_else(external_gene_name %in% tf_hsorth$mm_orths, "yes", "no")
  ) |>
  arrange(external_gene_name)

# write output ================================================================
write_csv(tf_all, here(basedir, "results/genesets/tf-genes.csv"))

# end =========================================================================
