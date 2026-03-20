# Arun Chavan
# Started: 2025-12-01

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)

basedir <- here("species/crob/KY21")
orthdir <- here("orthofinder/2025-08-29")

# genes =======================================================================
genes <- read_tsv(here(basedir, "results/genes/gene-display-names.tsv"))

# orthology data ==============================================================
## get mouse gene symbols and descriptions
mgenes <- read_tsv(here(orthdir, "results/id-name-mapping/mmus.tsv")) %>%
  select(gene_id, symbol) %>%
  distinct()

## pairwise orthologs
orth <- read_tsv(here(orthdir, "results/pairwise-orthologs/mmus_vs_crob.tsv"))

orth <- left_join(orth, mgenes, by = c("mmus" = "gene_id"))
orth <- left_join(
  orth,
  select(genes, gene_id, display_name),
  by = c("crob" = "gene_id")
)

hog_labels <- read_tsv(here(orthdir, "results/hog-labels/mmus_vs_crob.tsv"))

# data ========================================================================
mm <- read_tsv(here("species/mmus/gene-modules/functional-modules_gut.tsv"))

# get orthologs ===============================================================

get_orthologs <- function(x) {
  stopifnot(length(x) == 1)
  stopifnot(x %in% orth$symbol)
  orth |>
    filter(symbol == x) |>
    pull(crob) |>
    unique()
}

mm <- mm |>
  select(-comments) |>
  separate_longer_delim(genes, "; ") |>
  mutate(orths = map(genes, possibly(get_orthologs)))

mm_found <- mm |>
  filter(!map_lgl(orths, is.null))

mm_missing <- mm |>
  filter(map_lgl(orths, is.null))

# manual search for genes without an orthofinder ortholog =====================
get_match <- function(x, .col = "display_name") {
  genes |>
    filter(str_detect(.data[[.col]], regex(x, TRUE))) |>
    pull(gene_id) |>
    list()
}

genes |>
  filter(str_detect(display_name, regex("olfm4", TRUE))) |>
  print(n = Inf)

genes |>
  filter(str_detect(symbol, regex("amy2", TRUE))) |>
  print(n = Inf)

get_match("ptf1a", .col = "symbol")

mm_missing <- mm_missing |>
  mutate(
    orths = case_when(
      genes == "Cpe" ~ get_match("carboxypeptidase e"),
      genes == "Scg5" ~ get_match("Scg5"),
      genes == "Cck" ~ get_match("cholecystokinin$"),
      genes == "Dnah11" ~ get_match("dnah11"),
      genes == "Cyp27a1" ~ get_match("Cyp27a1", .col = "symbol"),
      genes == "Ptf1a" ~ get_match("Ptf1a", .col = "symbol"),
      genes == "Amy2a1" ~ get_match("Amy2"),
      genes == "Cpa2" ~ get_match("carboxypeptidase a2"),
      genes == "Pgc" ~ get_match("progastricsin"),
      genes == "Acat1" ~ get_match("Acat1"),
      genes == "Apob" ~ get_match("Apob;"),
      genes == "Slc23a1" ~ get_match("Slc23a1"),
      genes == "Hk3" ~ get_match("hk3"),
      genes == "G6pd" ~ get_match("G6pd"),
      genes == "Ass1" ~ get_match("Ass1"),
      genes == "Hmgcs2" ~ get_match("hmgcs1"),
      genes == "Acat1" ~ get_match("acat1"),
      genes == "Abca2" ~ get_match("abca2"),
      genes == "Ces1" ~ get_match("ces1;"),
      genes == "Ocln" ~ get_match("Ocln"),
      genes == "Muc5ac" ~ get_match("muc5ac"),
      genes == "Muc2" ~ get_match("muc2;"),
      genes == "Lyz1" ~ get_match("lysozyme"),
      genes == "Olfm4" ~ get_match("Olfm4"),
      .default = orths
    )
  )

# combine found and missing ===================================================
mm_all <- bind_rows(mm_found, mm_missing) |>
  unnest(cols = orths, keep_empty = TRUE) |>
  left_join(
    orth |> select(hog, crob) |> distinct(),
    by = c("orths" = "crob")
  ) |>
  left_join(
    select(hog_labels, hog, hog_label = label),
    by = "hog"
  ) |>
  left_join(
    select(genes, gene_id, display_name),
    by = c("orths" = "gene_id")
  ) |>
  arrange(label, module.1, module.2, module.3, module.4) |>
  rename(
    gene_mm = genes,
    ortholog = orths
  )

# write output ================================================================
write_tsv(
  mm_all,
  here(basedir, "results/genes/functional-modules_gut.tsv")
)

# end =========================================================================
