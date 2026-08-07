# Arun Chavan
# Started: 2026-08-07

# background ==================================================================

# setup =======================================================================
library(tidyverse)
library(here)
library(ape)

basedir <- here("orthofinder/2026-08-07")

# data ========================================================================
spp <- read_tsv(here(basedir, "data/species.tsv"))
trr <- read.tree(here(basedir, "data/species.nwk"))
trr <- ladderize(trr)

# generate non-timetree version with spp abbreviations ========================
trr$tip.label <- spp$abb[match(trr$tip.label, spp$binomial)]
trr$edge.length <- NULL
trr$node.label <- NULL

# write tree ==================================================================
write.tree(trr, here(basedir, "data/species-tree_abb.txt"))

# end =========================================================================
