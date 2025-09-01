# Arun Chavan
# Started: 2025-08-29

# background ==================================================================

# setup =======================================================================
library(ape)
library(phytools)
library(tidyverse)
library(here)
library(glue)
library(conflicted)
library(janitor)

conflicts_prefer(purrr::map())
conflicts_prefer(dplyr::filter())

basedir <- here("orthofinder/2025-08-29")
resdir <- here(basedir, "results/processed")
fs::dir_create(resdir)

###############################################################################
# plot species tree with node labels ==========================================
###############################################################################

sptr <- read.tree(here(
  basedir,
  "results/out/Species_Tree/SpeciesTree_rooted_node_labels.txt"
))

png(
  here(basedir, "results/species-tree-with-node-labels.png"),
  width = 7,
  height = 6,
  units = "in",
  res = 300
)
plot.phylo(sptr, use.edge.length = FALSE, show.node.label = TRUE)
dev.off()

###############################################################################
# orthogroups =================================================================
###############################################################################

ogs <- read_tsv(here(basedir, "results/out/Orthogroups/Orthogroups.tsv")) |>
  janitor::clean_names()

###############################################################################
# symbol-OG mapping ===========================================================
###############################################################################

# To make it easy to look up gene trees and orthogroups for a gene name, we
# create a mapping between symbols for a well-annotated species and the
# orthogroups they belong to

sym_spp <- "mmus"

symbols <- read_tsv(here(
  basedir,
  "results",
  "id-name-mapping",
  paste0(sym_spp, ".tsv")
)) |>
  select(gene_id, symbol) |>
  distinct() |>
  filter(!duplicated(symbol))

symbol_og <- ogs |>
  select(orthogroup, all_of(sym_spp)) |>
  separate_longer_delim(mmus, ", ") |>
  left_join(symbols, by = set_names("gene_id", sym_spp)) |>
  select(symbol, orthogroup) |>
  distinct() |>
  filter(!is.na(symbol)) |>
  deframe()


###############################################################################
# genetrees ===================================================================
###############################################################################

# save genetrees as a list of phylo objects named by orthogroup ids

genetrees <- read_delim(
  here(basedir, "results/out/Resolved_Gene_Trees/Resolved_Gene_Trees.txt"),
  delim = ": ",
  col_names = FALSE # otherwise, first tree gets omitted
) |>
  deframe() |>
  map(~ read.newick(text = .x))

# When available, replace gene id with symbol in gene tree leaf labels

gids <- list.files(here(basedir, "results", "id-name-mapping")) |>
  (function(x) set_names(x, str_remove(x, ".tsv")))() |>
  map(~ read_tsv(here(basedir, "results", "id-name-mapping", .x))) |>
  bind_rows(.id = "spp") |>
  select(spp, gene_id, symbol) |>
  distinct() |>
  filter(!is.na(symbol)) |>
  mutate(
    tip_id = paste0(spp, "_", gene_id),
    tip_symbol = paste0(spp, "_", symbol)
  )

stopifnot(!any(duplicated(gids$tip_id)))

tip_labels <- gids |>
  select(tip_id, tip_symbol) |>
  deframe()

replace_tip_ids_with_symbols <- function(tr) {
  lbl <- tr$tip.label
  idx <- match(lbl, names(tip_labels))
  hit <- !is.na(idx)
  lbl[hit] <- unname(tip_labels[idx[hit]])
  tr$tip.label <- lbl
  return(tr)
}

genetrees <- map(genetrees, replace_tip_ids_with_symbols, .progress = TRUE)

###############################################################################
# duplication nodes on gene trees =============================================
###############################################################################

nodelabel_to_nodeid <- function(tr, labels) {
  which(tr$node.label %in% labels) + length(tr$tip.label)
}

dups <- read_tsv(here(
  basedir,
  "results/out/Gene_Duplication_Events/Duplications.tsv"
)) |>
  janitor::clean_names()

dup_nodes <- dups |>
  select(orthogroup, gene_tree_node) |>
  group_by(orthogroup) |>
  summarise(nodelabs = list(gene_tree_node)) |>
  ungroup() |>
  mutate(
    nodeids = map2(
      orthogroup,
      nodelabs,
      ~ nodelabel_to_nodeid(genetrees[[.x]], .y)
    )
  ) |>
  # keep nodelabels as names and node ids as values. nodelabels function needs
  # node ids, not node labels
  mutate(nodeids = map2(nodeids, nodelabs, set_names)) |>
  select(orthogroup, nodeids) |>
  deframe()

###############################################################################
# HOGs for all nodes in species tree ==========================================
###############################################################################

# Since v3, orthofinder doesn't seem to create N0.tsv, but for N0 node (the
# root), we can simply use the Orthogroups.tsv file, which contains all
# orthogroups across all species.

hogs <- tibble(node = sort(sptr$node.label)) |>
  mutate(
    hogfile = case_when(
      node == "N0" ~ here(basedir, "results/out/Orthogroups/Orthogroups.tsv"),
      TRUE ~
        here(
          basedir,
          "results/out/Phylogenetic_Hierarchical_Orthogroups",
          paste0(node, ".tsv")
        )
    )
  ) |>
  mutate(hog = map(hogfile, ~ read_tsv(.x) |> clean_names())) |>
  select(-hogfile)

for (i in seq_len(nrow(hogs))) {
  if (hogs$node[i] == "N0") {
    hogs$hog[[i]] <- hogs$hog[[i]] |>
      rename(hog = orthogroup) |>
      mutate(og = hog) |>
      relocate(og, .after = hog) |>
      mutate(gene_tree_parent_clade = "n0") |>
      relocate(gene_tree_parent_clade, .after = og)
  }
}

hogs <- hogs |>
  unnest(cols = hog)

###############################################################################
# HOG node ids for all nodes in respective OG genetree ========================
###############################################################################

hog_ids_are_not_duplicated <- hogs |>
  group_by(og) |>
  summarize(hogs = list(hog)) |>
  ungroup() |>
  mutate(all_unique = map_lgl(hogs, ~ any(duplicated(.x)))) |>
  pull(all_unique) |>
  (function(x) {
    !any(x)
  })()

stopifnot(hog_ids_are_not_duplicated)

hog_nodes <- hogs |>
  filter(og %in% names(genetrees)) |> # some OGs don't seem to have genetrees??
  select(og, node, hog, gene_tree_parent_clade) |>
  distinct() |>
  arrange(og) |>
  mutate(
    genetree_nodeid = map2_int(
      og,
      gene_tree_parent_clade,
      ~ nodelabel_to_nodeid(genetrees[[.x]], .y),
      .progress = TRUE
    )
  )

###############################################################################
# write files =================================================================
###############################################################################

write_rds(symbol_og, here(resdir, "symbol-og-mapping.rds"), compress = "gz")
write_rds(genetrees, here(resdir, "genetrees.rds"), compress = "gz")
write_rds(dup_nodes, here(resdir, "duplication-nodes.rds"), compress = "gz")
write_rds(hogs, here(resdir, "hogs.rds"), compress = "gz")
write_rds(hog_nodes, here(resdir, "hog-node-ids.rds"), compress = "gz")

# end =========================================================================
