# Arun Chavan
# Started: 2025-08-14

# setup =======================================================================
library(phytools)
library(ape)
library(tidyverse)
library(here)
library(glue)
library(conflicted)

conflicts_prefer(purrr::map())
conflicts_prefer(dplyr::filter)

basedir <- here("orthofinder/2025-08-13")
orthdir <- here(basedir, "results/out")
treedir <- here(orthdir, "Resolved_Gene_Trees")

resdir <- here(basedir, "results/gene-tree-plots")
fs::dir_create(resdir)

# data ========================================================================
sptr <- read.tree(here(
  orthdir,
  "Species_Tree/SpeciesTree_rooted_node_labels.txt"
))

genetrees <- read_delim(
  here(treedir, "Resolved_Gene_Trees.txt"),
  delim = ": ",
  col_names = c("og", "tree")
)

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

dups <- read_tsv(here(orthdir, "Gene_Duplication_Events/Duplications.tsv")) |>
  janitor::clean_names()
dup_nodes <- dups |>
  select(orthogroup, gene_tree_node) |>
  group_by(orthogroup) |>
  summarise(nodes = list(gene_tree_node)) |>
  ungroup()

mmgenes <- read_tsv(here(basedir, "results", "id-name-mapping/mmus.tsv")) |>
  select(gene_id, symbol) |>
  distinct()

mmbf <- read_tsv(here(
  basedir,
  "results/pairwise-orthologs/mmus_vs_bflo.tsv"
)) |>
  left_join(mmgenes, by = c("mmus" = "gene_id"))

ogs <- map(
  list.files(here(orthdir, "Phylogenetic_Hierarchical_Orthogroups")),
  ~ read_tsv(here(orthdir, "Phylogenetic_Hierarchical_Orthogroups", .x))
) |>
  bind_rows() |>
  janitor::clean_names() |>
  separate(hog, into = c("sptr_node", NA), sep = "\\.", remove = FALSE)

og_sym <- ogs |>
  select(hog, og, mmus) |>
  separate_longer_delim(mmus, ", ") |>
  left_join(mmgenes, by = c("mmus" = "gene_id")) |>
  distinct()


stop()

# TODO: It's not clear if the HOG node annotations are being plotted at the
# correct node for all trees. Need to investigate this and fix it if necessary.

# functions ===================================================================
get_node_id_for_pair <- function(.sp1, .sp2) {
  all_labels <- c(sptr$tip.label, sptr$node.label)
  nodeid <- getMRCA(sptr, c(.sp1, .sp2))
  all_labels[nodeid]
}

get_og_for_symbol <- function(.symbol) {
  og_sym |>
    filter(symbol == .symbol) |>
    pull(og) |>
    unique()
}

get_hog_nodes_for_og <- function(.og, .sptr_node) {
  ogs |>
    filter(og == .og, sptr_node == .sptr_node) |>
    select(og, hog, gene_tree_parent_clade) |>
    distinct() |>
    pull(gene_tree_parent_clade)
}

get_dup_nodes_for_og <- function(.og) {
  dup_nodes |>
    filter(orthogroup == .og) |>
    pull(nodes) |>
    unlist()
}

replace_tip_ids_with_symbols <- function(.tr) {
  for (i in seq_along(.tr$tip.label)) {
    .tr$tip.label[i] <- ifelse(
      .tr$tip.label[i] %in% gids$tip_id,
      gids$tip_symbol[which(gids$tip_id == .tr$tip.label[i])],
      .tr$tip.label[i]
    )
  }
  return(.tr)
}

nodelabel_to_nodeid <- function(.tr, .label) {
  which(.tr$node.label %in% .label) + length(.tr$tip.label)
}

read_og_tree <- function(.og) {
  read.newick(text = genetrees$tree[genetrees$og == .og])
}

plot_og_tree <- function(
  .og,
  label_hogs = TRUE,
  label_dup_nodes = TRUE,
  node_for_hog = NULL
) {
  x <- read_og_tree(.og)
  x <- replace_tip_ids_with_symbols(x)

  par(mar = c(1, 1, 1, 1))
  plot.phylo(
    x,
    label.offset = 0.01,
    align.tip.label = TRUE,
    x.lim = max(node.depth.edgelength(x)) * 1.3
  )

  if (label_hogs) {
    x_hogs <- get_hog_nodes_for_og(.og, node_for_hog)
    nodelabels(
      node = nodelabel_to_nodeid(x, x_hogs),
      col = "salmon",
      pch = 19,
      cex = 2
    )
  }
  if (label_dup_nodes) {
    x_dup_nodes <- get_dup_nodes_for_og(.og)
    nodelabels(
      node = nodelabel_to_nodeid(x, x_dup_nodes),
      col = "grey",
      pch = 19,
      cex = 2
    )
  }
}

# plot gene trees =============================================================
pdf(file = here(resdir, "Csf1r_HOGs-mmus-bflo.pdf"), 11, 8)
plot_og_tree(
  get_og_for_symbol("Csf1r"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Csf1r_HOGs-mmus-pmar.pdf"), 11, 8)
plot_og_tree(
  get_og_for_symbol("Csf1r"),
  node_for_hog = get_node_id_for_pair("mmus", "pmar")
)
dev.off()

# ===
pdf(file = here(resdir, "Cyp7a1_HOGs-mmus-bflo.pdf"), 11, 8)
plot_og_tree(
  get_og_for_symbol("Cyp7a1"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Cyp7a1_HOGs-mmus-pmar.pdf"), 11, 8)
plot_og_tree(
  get_og_for_symbol("Cyp7a1"),
  node_for_hog = get_node_id_for_pair("mmus", "pmar")
)
dev.off()

# ===
pdf(file = here(resdir, "Ctrb1_HOGs-mmus-bflo.pdf"), 16, 76)
plot_og_tree(
  get_og_for_symbol("Ctrb1"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Ctrb1_HOGs-mmus-pmar.pdf"), 16, 76)
plot_og_tree(
  get_og_for_symbol("Ctrb1"),
  node_for_hog = get_node_id_for_pair("mmus", "pmar")
)
dev.off()

# ===
pdf(file = here(resdir, "Asic5_HOGs-mmus-bflo.pdf"), 11, 8)
plot_og_tree(
  get_og_for_symbol("Asic5"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

# ===
pdf(file = here(resdir, "Slc5a1_HOGs-mmus-bflo.pdf"), 11, 12)
plot_og_tree(
  get_og_for_symbol("Slc5a1"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Slc5a1_HOGs-mmus-pmar.pdf"), 11, 12)
plot_og_tree(
  get_og_for_symbol("Slc5a1"),
  node_for_hog = get_node_id_for_pair("mmus", "pmar")
)
dev.off()

# ===
pdf(file = here(resdir, "Ghr_HOGs-mmus-bflo.pdf"), 11, 14)
plot_og_tree(
  get_og_for_symbol("Ghr"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Ghr_HOGs-mmus-pmar.pdf"), 11, 14)
plot_og_tree(
  get_og_for_symbol("Ghr"),
  node_for_hog = get_node_id_for_pair("mmus", "pmar")
)
dev.off()

# ===
pdf(file = here(resdir, "Stat5a_HOGs-mmus-bflo.pdf"), 11, 6)
plot_og_tree(
  get_og_for_symbol("Stat5a"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Stat5a_HOGs-mmus-hsap.pdf"), 11, 6)
plot_og_tree(
  get_og_for_symbol("Stat5a"),
  node_for_hog = get_node_id_for_pair("mmus", "hsap")
)
dev.off()

# ===
pdf(file = here(resdir, "Wnt5a_HOGs-mmus-bflo.pdf"), 11, 12)
plot_og_tree(
  get_og_for_symbol("Wnt5a"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Wnt5a_HOGs-mmus-hsap.pdf"), 11, 12)
plot_og_tree(
  get_og_for_symbol("Wnt5a"),
  node_for_hog = get_node_id_for_pair("mmus", "hsap")
)
dev.off()

# ===
pdf(file = here(resdir, "Pdgfa_HOGs-mmus-bflo.pdf"), 11, 12)
plot_og_tree(
  get_og_for_symbol("Pdgfa"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Pdgfa_HOGs-mmus-hsap.pdf"), 11, 12)
plot_og_tree(
  get_og_for_symbol("Pdgfa"),
  node_for_hog = get_node_id_for_pair("mmus", "hsap")
)
dev.off()

# ===
pdf(file = here(resdir, "Hgf_HOGs-mmus-bflo.pdf"), 11, 24)
plot_og_tree(
  get_og_for_symbol("Hgf"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Hgf_HOGs-mmus-hsap.pdf"), 11, 24)
plot_og_tree(
  get_og_for_symbol("Hgf"),
  node_for_hog = get_node_id_for_pair("mmus", "hsap")
)
dev.off()

# ===
pdf(file = here(resdir, "Hal_HOGs-mmus-bflo.pdf"), 6, 2)
plot_og_tree(
  get_og_for_symbol("Hal"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

# ===
pdf(file = here(resdir, "Insr_HOGs-mmus-bflo.pdf"), 11, 6)
plot_og_tree(
  get_og_for_symbol("Insr"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Insr_HOGs-mmus-pmar.pdf"), 11, 6)
plot_og_tree(
  get_og_for_symbol("Insr"),
  node_for_hog = get_node_id_for_pair("mmus", "pmar")
)
dev.off()

# ===
pdf(file = here(resdir, "Igf1_HOGs-mmus-bflo.pdf"), 11, 8)
plot_og_tree(
  get_og_for_symbol("Igf1"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

pdf(file = here(resdir, "Igf1_HOGs-mmus-pmar.pdf"), 11, 8)
plot_og_tree(
  get_og_for_symbol("Igf1"),
  node_for_hog = get_node_id_for_pair("mmus", "pmar")
)
dev.off()

# ===
pdf(file = here(resdir, "Amy1_HOGs-mmus-bflo.pdf"), 11, 8)
plot_og_tree(
  get_og_for_symbol("Amy1"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()

# ===
pdf(file = here(resdir, "Lamp1_HOGs-mmus-bflo.pdf"), 11, 8)
plot_og_tree(
  get_og_for_symbol("Lamp1"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)
dev.off()


stop()


# ===
plot_og_tree(
  get_og_for_symbol("Ins2"),
  node_for_hog = get_node_id_for_pair("mmus", "pmar")
)

plot_og_tree(
  get_og_for_symbol("Igfbp2"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)


plot_og_tree(
  get_og_for_symbol("Sox9"),
  node_for_hog = get_node_id_for_pair("mmus", "bflo")
)

# end =========================================================================
