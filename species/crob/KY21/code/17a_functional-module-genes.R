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
get_match <- function(x) {
  genes |>
    filter(str_detect(display_name, regex(x, TRUE))) |>
    pull(gene_id) |>
    list()
}

genes |>
  filter(str_detect(display_name, regex("defensin", TRUE)))

get_match("hk3")

mm_missing <- mm_missing |>
  mutate(
    orths = case_when(
      genes == "Cpe" ~ get_match("carboxypeptidase e"),
      genes == "Scg5" ~ get_match("Scg5"),
      genes == "Cyp27a1" ~ get_match("Cyp27a1"),
      genes == "Ptf1a" ~ get_match("Ptf1a;"),
      genes == "Amy2a1" ~ get_match("Amy2"),
      genes == "Cpa2" ~ get_match("carboxypeptidase a2"),
      genes == "Pgc" ~ get_match("progastricsin"),
      genes == "Hk3" ~ get_match("hk3"),
      genes == "Hmgcs2" ~ get_match("hmgcs1"),
      genes == "Acat1" ~ get_match("acat1"),
      genes == "Abca2" ~ get_match("abca2"),
      genes == "Ces1" ~ get_match("ces1;"),
      genes == "Muc5ac" ~ get_match("muc5ac"),
      genes == "Muc2" ~ get_match("muc2"),
      genes == "Lyz1" ~ get_match("lysozyme"),
      .default = orths
    )
  )

# combine found and missing ===================================================
mm_all <- bind_rows(mm_found, mm_missing)

mm_all <- mm_all |>
  unnest(orths, keep_empty = TRUE) |>
  left_join(
    select(genes, gene_id, display_name),
    by = c("orths" = "gene_id")
  ) |>
  arrange(module.1, module.2, module.3, module.4) |>
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
