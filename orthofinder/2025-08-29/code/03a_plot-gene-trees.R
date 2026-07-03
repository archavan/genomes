# Arun Chavan
# Started: 2025-08-30

# setup =======================================================================
library(phytools)
library(ape)
library(tidyverse)
library(here)
library(glue)
library(conflicted)

conflicts_prefer(purrr::map())
conflicts_prefer(dplyr::filter())

basedir <- here("orthofinder/2025-08-29")
resdir <- here(basedir, "results/gene-tree-plots")
fs::dir_create(resdir)

# data ========================================================================
sptr <- read.tree(here(
  basedir,
  "results/out/Species_Tree/SpeciesTree_rooted_node_labels.txt"
))

genetrees <- read_rds(here(basedir, "results/processed/genetrees.rds"))
dup_nodes <- read_rds(here(basedir, "results/processed/duplication-nodes.rds"))
sym_og <- read_rds(here(basedir, "results/processed/symbol-og-mapping.rds"))
hogs <- read_rds(here(basedir, "results/processed/hogs.rds"))
hog_nodes <- read_rds(here(basedir, "results/processed/hog-node-ids.rds"))

# functions ===================================================================
get_mrca_node_id <- function(.sp1, .sp2) {
  all_labels <- c(sptr$tip.label, sptr$node.label)
  nodeid <- getMRCA(sptr, c(.sp1, .sp2))
  all_labels[nodeid]
}

get_hog_nodes <- function(.og, .sptr_node) {
  hog_nodes |>
    filter(og == .og, node == .sptr_node) |>
    pull(genetree_nodeid) |>
    unique()
}

plot_gene_tree <- function(
  .symbol = NULL,
  .og = NULL,
  label_dup_nodes = TRUE,
  label_hog_nodes = TRUE,
  node_for_hog = NULL
) {
  if (!is.null(.symbol) & !is.null(.og)) {
    stop("Provide either .symbol or .og, not both")
  }
  if (is.null(.symbol) & !is.null(.og)) {
    og <- .og
  }
  if (!is.null(.symbol) & is.null(.og)) {
    og <- sym_og[[.symbol]]
  }

  x <- genetrees[[og]]

  par(mar = c(1, 1, 1, 1))
  plot.phylo(
    x,
    label.offset = 0.01,
    align.tip.label = TRUE,
    x.lim = max(node.depth.edgelength(x)) * 1.3
  )

  if (label_dup_nodes) {
    nodelabels(node = dup_nodes[[og]], col = "grey", pch = 19, cex = 1.5)
  }

  if (label_hog_nodes) {
    x_hogs <- get_hog_nodes(og, node_for_hog)
    nodelabels(node = x_hogs, col = "salmon", pch = 19, cex = 1.5)
  }
}

# plot gene trees =============================================================
## ===
pdf(file = here(resdir, "Csf1r_HOGs-mmus-bflo.pdf"), 11, 8)
plot_gene_tree("Csf1r", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Csf1r_HOGs-mmus-pmar.pdf"), 11, 8)
plot_gene_tree("Csf1r", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()


## ===
pdf(file = here(resdir, "Cyp7a1_HOGs-mmus-bflo.pdf"), 11, 10)
plot_gene_tree("Cyp7a1", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Cyp7a1_HOGs-mmus-pmar.pdf"), 11, 10)
plot_gene_tree("Cyp7a1", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Cyp7a1_HOGs-crob-scla.pdf"), 11, 10)
plot_gene_tree("Cyp7a1", node_for_hog = get_mrca_node_id("crob", "scla"))
dev.off()


## ===
pdf(file = here(resdir, "Ctrb1_HOGs-mmus-bflo.pdf"), 11, 10)
plot_gene_tree("Ctrb1", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Ctrb1_HOGs-mmus-pmar.pdf"), 11, 10)
plot_gene_tree("Ctrb1", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Ctrb1_HOGs-crob-scla.pdf"), 11, 10)
plot_gene_tree("Ctrb1", node_for_hog = get_mrca_node_id("crob", "scla"))
dev.off()


## ===
pdf(file = here(resdir, "Asic5_HOGs-mmus-bflo.pdf"), 11, 10)
plot_gene_tree("Asic5", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Asic5_HOGs-mmus-pmar.pdf"), 11, 10)
plot_gene_tree("Asic5", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Asic5_HOGs-crob-scla.pdf"), 11, 10)
plot_gene_tree("Asic5", node_for_hog = get_mrca_node_id("crob", "scla"))
dev.off()


## ===
pdf(file = here(resdir, "Abcg1_HOGs-mmus-bflo.pdf"), 11, 7)
plot_gene_tree("Abcg1", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Abcg1_HOGs-mmus-pmar.pdf"), 11, 7)
plot_gene_tree("Abcg1", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Abcg1_HOGs-crob-scla.pdf"), 11, 7)
plot_gene_tree("Abcg1", node_for_hog = get_mrca_node_id("crob", "scla"))
dev.off()


## ===
pdf(file = here(resdir, "Slc5a1_HOGs-mmus-bflo.pdf"), 11, 10)
plot_gene_tree("Slc5a1", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Slc5a1_HOGs-mmus-pmar.pdf"), 11, 10)
plot_gene_tree("Slc5a1", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Slc5a1_HOGs-crob-scla.pdf"), 11, 10)
plot_gene_tree("Slc5a1", node_for_hog = get_mrca_node_id("crob", "scla"))
dev.off()


## ===
pdf(file = here(resdir, "Ghr_HOGs-mmus-bflo.pdf"), 11, 10)
plot_gene_tree("Ghr", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Ghr_HOGs-mmus-pmar.pdf"), 11, 10)
plot_gene_tree("Ghr", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Ghr_HOGs-drer-trub.pdf"), 11, 10)
plot_gene_tree("Ghr", node_for_hog = get_mrca_node_id("drer", "trub"))
dev.off()


## ===
pdf(file = here(resdir, "Igf1_HOGs-mmus-bflo.pdf"), 11, 10)
plot_gene_tree("Igf1", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Igf1_HOGs-mmus-pmar.pdf"), 11, 10)
plot_gene_tree("Igf1", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Igf1_HOGs-crob-scla.pdf"), 11, 10)
plot_gene_tree("Igf1", node_for_hog = get_mrca_node_id("crob", "scla"))
dev.off()


## ===
pdf(file = here(resdir, "Amy1_HOGs-mmus-bflo.pdf"), 11, 10)
plot_gene_tree("Amy1", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Amy1_HOGs-mmus-pmar.pdf"), 11, 10)
plot_gene_tree("Amy1", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Amy1_HOGs-crob-scla.pdf"), 11, 10)
plot_gene_tree("Amy1", node_for_hog = get_mrca_node_id("crob", "scla"))
dev.off()


## ===
pdf(file = here(resdir, "Pparg_HOGs-mmus-bflo.pdf"), 11, 10)
plot_gene_tree("Pparg", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Pparg_HOGs-mmus-pmar.pdf"), 11, 10)
plot_gene_tree("Pparg", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Pparg_HOGs-crob-scla.pdf"), 11, 10)
plot_gene_tree("Pparg", node_for_hog = get_mrca_node_id("crob", "scla"))
dev.off()


## ===
pdf(file = here(resdir, "Hal_HOGs-mmus-bflo.pdf"), 11, 5)
plot_gene_tree("Hal", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Hal_HOGs-mmus-pmar.pdf"), 11, 5)
plot_gene_tree("Hal", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Hal_HOGs-crob-scla.pdf"), 11, 5)
plot_gene_tree("Hal", node_for_hog = get_mrca_node_id("crob", "scla"))
dev.off()


## ===
pdf(file = here(resdir, "Hnf4a_HOGs-mmus-bflo.pdf"), 11, 6)
plot_gene_tree("Hnf4a", node_for_hog = get_mrca_node_id("mmus", "bflo"))
dev.off()

pdf(file = here(resdir, "Hnf4a_HOGs-mmus-pmar.pdf"), 11, 6)
plot_gene_tree("Hnf4a", node_for_hog = get_mrca_node_id("mmus", "pmar"))
dev.off()

pdf(file = here(resdir, "Hnf4a_HOGs-crob-scla.pdf"), 11, 6)
plot_gene_tree("Hnf4a", node_for_hog = get_mrca_node_id("crob", "scla"))
dev.off()

pdf(file = here(resdir, "Cyp3a11.pdf"), 11, 20)
plot_gene_tree("Cyp3a11", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

pdf(file = here(resdir, "Cyp4a10.pdf"), 11, 8)
plot_gene_tree("Cyp4a10", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

pdf(file = here(resdir, "Slc2a1.pdf"), 11, 20)
plot_gene_tree("Slc2a1", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

pdf(file = here(resdir, "Nlrp1b.pdf"), 11, 9)
plot_gene_tree("Nlrp1b", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

pdf(file = here(resdir, "Nlrp1a.pdf"), 11, 11)
plot_gene_tree("Nlrp1a", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

## ===
# RTKs
pdf(file = here(resdir, "Csf1.pdf"), 11, 3)
plot_gene_tree("Csf1", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

pdf(file = here(resdir, "Csf1r.pdf"), 11, 5)
plot_gene_tree("Csf1r", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

pdf(file = here(resdir, "Pdgfa.pdf"), 11, 10)
plot_gene_tree("Pdgfa", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

pdf(file = here(resdir, "Kitl.pdf"), 11, 3)
plot_gene_tree("Kitl", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

pdf(file = here(resdir, "Pdgfra.pdf"), 11, 3)
plot_gene_tree("Pdgfra", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

pdf(file = here(resdir, "Flt1.pdf"), 11, 10)
plot_gene_tree("Flt1", label_dup_nodes = TRUE, label_hog_nodes = FALSE)
dev.off()

# end =========================================================================
