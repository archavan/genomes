# Arun Chavan
# Started: 2025-08-14

# setup =======================================================================
library(ape)
library(tidyverse)
library(here)
library(janitor)

orthdir <- here("orthofinder/2025-08-13")
resdir <- here(orthdir, "results/pairwise-orthologs")
fs::dir_create(resdir)

# species tree ================================================================
tr <- read.tree(here(
  orthdir,
  "results/out/Species_Tree/SpeciesTree_rooted_node_labels.txt"
))

# plot tree
png(
  here(orthdir, "results/species-tree-with-node-labels.png"),
  width = 7,
  height = 7,
  units = "in",
  res = 300
)
plot.phylo(tr, use.edge.length = FALSE, show.node.label = TRUE)
dev.off()

get_node_label_for_pair <- function(.sp1, .sp2) {
  all_labels <- c(tr$tip.label, tr$node.label)
  nodeid <- getMRCA(tr, c(.sp1, .sp2))
  all_labels[nodeid]
}

# read HOGs ====================================================================
# We are using HOGs to assign orthology types. The orthologs files provided by
# orthofinder have problems with one:one orthologs. See:
# https://github.com/davidemms/OrthoFinder/issues/969

# Since v3, orthofinder doesn't seem to create N0.tsv, but for N0 node (the
# root), we can simply use the Orthogroups.tsv file, which contains all
# orthogroups across all species.

df <- expand_grid(sp1 = tr$tip.label, sp2 = tr$tip.label) |>
  filter(sp1 != sp2) |>
  mutate(node = map2_chr(sp1, sp2, ~ get_node_label_for_pair(.x, .y))) |>
  mutate(
    hogfile = case_when(
      node == "N0" ~ here(orthdir, "results/out/Orthogroups/Orthogroups.tsv"),
      TRUE ~
        here(
          orthdir,
          "results/out/Phylogenetic_Hierarchical_Orthogroups",
          paste0(node, ".tsv")
        )
    )
  ) |>
  mutate(hog = map(hogfile, ~ read_tsv(.x) |> clean_names()))

for (i in seq_len(nrow(df))) {
  if (df$node[i] == "N0") {
    df$hog[[i]] <- rename(df$hog[[i]], hog = orthogroup)
  }
}

# functions for extracting orthologs ==========================================
add_orthology_type_info <- function(.orth) {
  x <- .orth |>
    mutate(
      type = case_when(
        !str_detect(sp1, ", ") & !str_detect(orthologs, ", ") ~ "one.to.one",
        !str_detect(sp1, ", ") & str_detect(orthologs, ", ") ~ "one.to.many",
        str_detect(sp1, ", ") & !str_detect(orthologs, ", ") ~ "many.to.one",
        str_detect(sp1, ", ") & str_detect(orthologs, ", ") ~ "many.to.many",
        .default = NA_character_
      )
    )
  stopifnot(!any(is.na(x$type)))
  return(x)
}

unroll_orthologs <- function(.orth) {
  .orth |>
    mutate(
      sp1 = map(sp1, ~ str_split_1(.x, ", ")),
      orthologs = map(orthologs, ~ str_split_1(.x, ", "))
    ) |>
    mutate(
      orth = map2(sp1, orthologs, ~ expand_grid(sp1 = .x, ortholog = .y))
    ) |>
    select(-sp1, -orthologs) |>
    unnest(orth)
}

extract_pairwise_orthologs <- function(.sp1, .sp2, .hog) {
  .hog |>
    select(all_of(c("hog", .sp1, .sp2))) |>
    rename(sp1 = {{ .sp1 }}, orthologs = {{ .sp2 }}) |>
    filter(!is.na(sp1) & !is.na(orthologs)) |>
    add_orthology_type_info() |>
    unroll_orthologs() |>
    relocate(type, .after = last_col())
}

verify_that_one2one_are_not_duplicated <- function(.orth) {
  n_dup <- .orth |>
    filter(type == "one.to.one") |>
    filter(duplicated(sp1) | duplicated(ortholog)) |>
    nrow()
  stopifnot(n_dup == 0)
  return(.orth)
}

rename_columns_to_species_abbr <- function(.orth, .sp1, .sp2) {
  .orth |>
    rename(!!sym(.sp1) := sp1, !!sym(.sp2) := ortholog)
}

extract_and_clean_up_pairwise_orthologs <- function(.sp1, .sp2, .hog) {
  extract_pairwise_orthologs(.sp1, .sp2, .hog) |>
    verify_that_one2one_are_not_duplicated() |>
    rename_columns_to_species_abbr(.sp1, .sp2)
}


# extract pairwise orthologs ==================================================
df$orth <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
  message(paste0(df$sp1[i], " vs ", df$sp2[i]))
  df$orth[[i]] <- extract_and_clean_up_pairwise_orthologs(
    df$sp1[i],
    df$sp2[i],
    df$hog[[i]]
  )
}

# count number of orthologs by type ===========================================
get_n_type <- function(.orth, .type, .focal_spp) {
  .orth |>
    filter(type == .type) |>
    pull(.focal_spp) |>
    unique() |>
    length()
}

df <- df |>
  mutate(
    n_one.to.one = map2_int(orth, sp1, ~ get_n_type(.x, "one.to.one", .y)),
    n_one.to.many = map2_int(orth, sp1, ~ get_n_type(.x, "one.to.many", .y)),
    n_many.to.one = map2_int(orth, sp1, ~ get_n_type(.x, "many.to.one", .y)),
    n_many.to.many = map2_int(orth, sp1, ~ get_n_type(.x, "many.to.many", .y))
  )

# write files =================================================================
for (i in seq_len(nrow(df))) {
  write_tsv(
    df$orth[[i]],
    here(resdir, paste0(df$sp1[i], "_vs_", df$sp2[i], ".tsv"))
  )
}

df |>
  select(
    sp1,
    sp2,
    n_one.to.one,
    n_one.to.many,
    n_many.to.one,
    n_many.to.many
  ) |>
  write_tsv(here(orthdir, "results/number-of-pairwise-orthologs.tsv"))

# end =========================================================================
