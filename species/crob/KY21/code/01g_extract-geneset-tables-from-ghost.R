# Arun Chavan
# Started: 2024-12-18

# background ==================================================================

# Extract geneset tables from GHOST: http://ghost.zool.kyoto-u.ac.jp/HT_misc_KY21.html
# These are in turn based on Imai et al 2004: https://doi.org/10.1242/dev.01270

# setup =======================================================================
library(tidyverse)
library(rvest)
library(here)

genomedir <- here("species/crob/KY21")
resdir <- here(genomedir, "results/genes")

ghost_genes <- read_tsv(here(
  genomedir,
  "results/ghost/gene-info_ghostdb_parsed.tsv"
))

clean_up_table <- function(x) {
  x %>%
    rename(
      gene_name = X1,
      clone_used_for_ISH = X2,
      gene_id = X3,
      notes = X4
    ) %>%
    relocate(gene_id, .after = family) %>%
    filter(!is.na(gene_id))
}

# transcription factors =======================================================
tf_non_znf <- read_html("http://ghost.zool.kyoto-u.ac.jp/HT_TF_KY21.html") %>%
  html_table(trim = TRUE, na.strings = c("", "NA"), convert = TRUE) %>%
  .[4:12] %>% # others are header and footer tables
  set_names(c(
    "bHLH",
    "bZIP",
    "Ets",
    "Fox",
    "HMG",
    "Homeobox",
    "Nuclear receptor",
    "T-box",
    "Other TFs"
  )) %>%
  bind_rows(.id = "family") %>%
  clean_up_table()

# They have a list of ZNFs that includes possible TFs and NOT TFs
znf <- read_html("http://ghost.zool.kyoto-u.ac.jp/HT_ZF_KY21.html") %>%
  html_table(trim = TRUE, na.strings = c("", "NA"), convert = TRUE) %>%
  .[4:5] %>%
  set_names(c(
    "possible TFs with ZNF motif",
    "ZNF proteins that are NOT TFs"
  )) %>%
  bind_rows(.id = "family") %>%
  clean_up_table()

tf_znf <- znf %>%
  filter(family == "possible TFs with ZNF motif") %>%
  mutate(family = "ZNF")

tf_clean <- bind_rows(tf_non_znf, tf_znf) %>%
  filter(str_detect(gene_id, "^KY21")) %>% # some don't have a KY21 id
  select(family, gene_id)

# signaling molecules =========================================================
signaling <- read_html("http://ghost.zool.kyoto-u.ac.jp/HT_ST_KY21.html") %>%
  html_table(trim = TRUE, na.strings = c("", "NA"), convert = TRUE) %>%
  .[4:10] %>%
  set_names(c(
    "RTK_MAPK",
    "TGFb",
    "Wnt",
    "Hedgehog",
    "JAK_STAT",
    "NFkB",
    "Notch"
  )) %>%
  bind_rows(.id = "family") %>%
  clean_up_table()

signaling_clean <- signaling %>%
  filter(str_detect(gene_id, "^KY21")) %>%
  select(family, gene_id)

# miscellaneous genes =========================================================

misc <- read_html("http://ghost.zool.kyoto-u.ac.jp/HT_misc_KY21.html") %>%
  html_table(trim = TRUE, na.strings = c("", "NA"), convert = TRUE) %>%
  .[4:7] %>%
  set_names(c(
    "Muscle genes",
    "Protein kinase genes",
    "Protein phosphatase genes",
    "Misc"
  )) %>%
  bind_rows(.id = "family") %>%
  clean_up_table()

misc_clean <- misc %>%
  filter(str_detect(gene_id, "^KY21")) %>%
  select(family, gene_id)


# write files =================================================================
write_tsv(tf_non_znf, here(genomedir, "data/ghost-genesets/TF_non-ZNF.tsv"))
write_tsv(znf, here(genomedir, "data/ghost-genesets/ZNF-proteins.tsv"))
write_tsv(
  signaling,
  here(genomedir, "data/ghost-genesets/signaling-molecules.tsv")
)
write_tsv(misc, here(genomedir, "data/ghost-genesets/misc-genes.tsv"))

write_tsv(tf_clean, here(resdir, "geneset_tf_ghost-curated.tsv"))
write_tsv(signaling_clean, here(resdir, "geneset_signaling_ghost-curated.tsv"))
write_tsv(misc_clean, here(resdir, "geneset_misc_ghost-curated.tsv"))

# end =========================================================================
