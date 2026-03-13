# Arun Chavan
# Started: 2025-12-09

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)

basedir <- here("species/scla/refseq_v1")
orthdir <- here("orthofinder/2025-08-29")

# genes =======================================================================
genes <- read_tsv(here(basedir, "results/genes/gene-display-names.tsv"))

# orthology data ==============================================================
## get mouse gene symbols and descriptions
mgenes <- read_tsv(here(orthdir, "results/id-name-mapping/mmus.tsv")) %>%
  select(gene_id, symbol) %>%
  distinct()

## pairwise orthologs
orth <- read_tsv(here(orthdir, "results/pairwise-orthologs/mmus_vs_scla.tsv"))

orth <- left_join(orth, mgenes, by = c("mmus" = "gene_id"))
orth <- left_join(
  orth,
  select(genes, gene_id, display_name),
  by = c("scla" = "gene_id")
)

hog_labels <- read_tsv(here(orthdir, "results/hog-labels/mmus_vs_scla.tsv"))

# data ========================================================================
mm <- read_tsv(here("species/mmus/gene-modules/functional-modules_gut.tsv"))

# get orthologs ===============================================================

get_orthologs <- function(x) {
  stopifnot(length(x) == 1)
  stopifnot(x %in% orth$symbol)
  orth |>
    filter(symbol == x) |>
    pull(scla) |>
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
  filter(str_detect(
    display_name,
    regex("defa5", TRUE)
  )) |>
  print(n = Inf)

genes |>
  filter(str_detect(
    symbol,
    regex("gstp1", TRUE)
  )) |>
  print(n = Inf)

get_match("ncbi.120346634")

mm_missing <- mm_missing |>
  mutate(
    orths = case_when(
      genes == "Cpe" ~ get_match("carboxypeptidase e"),
      genes == "Pcsk1" ~ get_match("Pcsk1"),
      genes == "Gcg" ~ get_match("ncbi.120346634"),
      genes == "Pyy" ~ get_match("Pyy"),
      genes == "Ghrl" ~ get_match("ghrl"),

      genes == "Aqp4" ~ get_match("Aquaporin-4"),
      genes == "Dnah11" ~ get_match("dnah11"),
      genes == "Cela2a" ~ get_match(
        "chymotrypsin-like elastase family member 2A"
      ),
      genes == "Rac1" ~ get_match("ras-related C3 botulinum toxin substrate 1"),
      genes == "Lamp5" ~ get_match("lamp5"),
      genes == "Slc7a9" ~ get_match("slc7a9", .col = "symbol"),
      genes == "Scarb1" ~ get_match("scarb1", .col = "symbol"),
      genes == "Fabp5" ~ get_match("fabp5"),
      genes == "Npc1" ~ get_match("npc1;"),
      genes == "Acat1" ~ get_match("acat1"),
      genes == "Apoa1" ~ get_match("apoa1;"),
      genes == "Abca1" ~ get_match("abca1"),
      genes == "G6pd" ~ get_match("G6pd"),
      genes == "Gstp1" ~ get_match("Gstp1", .col = "symbol"),
      genes == "Ocln" ~ get_match("Ocln"),
      genes == "Muc5ac" ~ get_match(
        "\\[ncbi\\.\\d+\\]; \\[n\\] mucin-5AC-like"
      ),
      genes == "Muc2" ~ get_match("muc2;"),
      genes == "Tff2" ~ get_match("Tff2"),
      genes == "Reg3g" ~ get_match("reg3g"),
      genes == "Lyz1" ~ get_match("lysozyme g-like"),
      genes == "Defa5" ~ get_match("Defa5"),
      .default = orths
    )
  )

# combine found and missing ===================================================
mm_all <- bind_rows(mm_found, mm_missing) |>
  unnest(cols = orths, keep_empty = TRUE) |>
  left_join(
    orth |> select(hog, scla) |> distinct(),
    by = c("orths" = "scla")
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
