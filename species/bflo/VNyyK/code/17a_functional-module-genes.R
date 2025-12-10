# Arun Chavan
# Started: 2025-12-01

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)

basedir <- here("species/bflo/VNyyK")
orthdir <- here("orthofinder/2025-08-29")

# genes =======================================================================
genes <- read_tsv(here(basedir, "results/genes/gene-display-names.tsv"))

# orthology data ==============================================================
## get mouse gene symbols and descriptions
mgenes <- read_tsv(here(orthdir, "results/id-name-mapping/mmus.tsv")) |>
  select(gene_id, symbol) |>
  distinct()

## pairwise orthologs
orth <- read_tsv(here(orthdir, "results/pairwise-orthologs/mmus_vs_bflo.tsv"))
orth <- orth |>
  mutate(bflo = str_replace(bflo, "ncbi\\.", "ncbi\\:")) # orthofinder quirk

orth <- left_join(orth, mgenes, by = c("mmus" = "gene_id"))
orth <- left_join(
  orth,
  select(genes, gene_id, display_name),
  by = c("bflo" = "gene_id")
)

# data ========================================================================
mm <- read_tsv(here("species/mmus/gene-modules/functional-modules_gut.tsv"))

# get orthologs ===============================================================

get_orthologs <- function(x) {
  stopifnot(length(x) == 1)
  stopifnot(x %in% orth$symbol)
  orth |>
    filter(symbol == x) |>
    pull(bflo) |>
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
  filter(str_detect(display_name, regex("leucine-rich", TRUE))) |>
  print(n = Inf)

get_match("hk3")

mm_missing <- mm_missing |>
  mutate(
    orths = case_when(
      genes == "Pcsk1n" ~ get_match("Pcsk1n"),
      genes == "Fabp3" ~ get_match("Fabp3"),
      genes == "Fabp4" ~ get_match("Fabp4"),
      genes == "Fabp5" ~ get_match("Fabp5"),
      genes == "Fabp7" ~ get_match("Fabp7"),
      genes == "Ces1" ~ get_match("carboxylesterase 1-like"),
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
