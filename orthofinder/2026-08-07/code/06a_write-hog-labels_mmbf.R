# Arun Chavan
# Started: 2026-08-11
#
# Construct informative HOG display labels for mmus vs bflo pairwise orthologs.
# Labels use LCP-based clustering of mouse gene symbols.

# setup =======================================================================
library(tidyverse)
library(here)

orthdir <- here("orthofinder/2026-08-07")
source(here(orthdir, "code/utils_hog-labels.R"))

resdir <- here(orthdir, "results/hog-labels")
fs::dir_create(resdir)

# read data ===================================================================

# pairwise orthologs
pw <- read_tsv(here(orthdir, "results/pairwise-orthologs/mmus_vs_bflo.tsv"))

# mouse gene symbols (deduplicated to gene-level)
symbols <- read_tsv(here(orthdir, "results/id-name-mapping/mmus.tsv")) |>
  select(gene_id, symbol) |>
  distinct() |>
  filter(!duplicated(gene_id))

# join symbols ================================================================
pw <- pw |>
  left_join(symbols, by = c("mmus" = "gene_id"))

# verify one orthology type per HOG ===========================================
type_check <- pw |>
  group_by(hog) |>
  summarize(n_types = n_distinct(type), .groups = "drop") |>
  filter(n_types > 1)

if (nrow(type_check) > 0) {
  stop(
    "The following HOGs have multiple orthology types: ",
    paste(type_check$hog, collapse = ", ")
  )
}

# compute labels ==============================================================
hog_labels <- pw |>
  group_by(hog) |>
  summarize(
    n_mmus = n_distinct(mmus),
    n_bflo = n_distinct(bflo),
    type = unique(type),
    symbols_list = list(unique(symbol)),
    symbols = str_c(sort(unique(symbol)), collapse = ", "),
    .groups = "drop"
  ) |>
  rowwise() |>
  mutate(
    label = format_hog_label(
      syms = symbols_list,
      n_sp1 = n_mmus,
      n_sp2 = n_bflo,
      sp1_abb = "m",
      sp2_abb = "b"
    )
  ) |>
  ungroup() |>
  select(hog, label, n_mmus, n_bflo, type, symbols)

# disambiguate duplicate labels ================================================
# When no gene symbols are available, the label is just "mXbY" — append HOG ID.
# Rare cases (e.g., multiple Gm+Rik pairs) can also produce identical labels.
# Append the HOG ID to disambiguate those too.
dup_labels <- hog_labels$label[duplicated(hog_labels$label)]

hog_labels <- hog_labels |>
  mutate(
    label = if_else(
      label %in% dup_labels | !str_detect(label, "_"),
      paste0(label, "_", hog),
      label
    )
  )

# verify labels are unique ====================================================
stopifnot(!any(duplicated(hog_labels$label)))

# write output ================================================================
write_tsv(hog_labels, here(resdir, "mmus_vs_bflo.tsv"))

# end =========================================================================
