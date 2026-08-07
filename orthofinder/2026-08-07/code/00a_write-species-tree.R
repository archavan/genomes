# Arun Chavan
# Started: 2026-08-07

# background ==================================================================

# setup =======================================================================
library(tidyverse)
library(here)
library(ape)
library(phytools)

basedir <- here("orthofinder/2026-08-07")

# data ========================================================================
spp <- read_tsv(here(basedir, "data/species.tsv"))
trr <- read.tree(here(basedir, "data/species.nwk"))

# generate non-timetree version with spp abbreviations ========================
trr$edge.length <- NULL
trr$node.label <- NULL

# TimeTree could not place Leucoraja erinacea (little skate), so it is missing
# from the downloaded newick file and has to be grafted on manually. It is a
# batoid, so it is placed as sister to the shark Scyliorhinus canicula; together
# they form Elasmobranchii, sister to the holocephalan Callorhinchus milii.
trr <- bind.tip(
  trr,
  tip.label = "Leucoraja_erinacea",
  where = which(trr$tip.label == "Scyliorhinus_canicula")
)

trr <- ladderize(trr)
trr$tip.label <- spp$abb[match(trr$tip.label, spp$binomial)]

stopifnot(
  !anyNA(trr$tip.label),
  setequal(trr$tip.label, spp$abb)
)

# write tree ==================================================================
write.tree(trr, here(basedir, "data/species-tree_abb.txt"))

# end =========================================================================
