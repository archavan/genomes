# Arun Chavan
# Started: 2026-10-22

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(Biostrings)
library(conflicted)

conflicts_prefer(
  dplyr::select,
  dplyr::filter,
  dplyr::slice,
  dplyr::rename
)

basedir <- here("gene-trees/rtk/full-gene")
resdir <- here(basedir, "results/05a_seqs-assembled")
fs::dir_create(resdir)

orthdir <- here("orthofinder/2026-08-07")

###############################################################################
# data ========================================================================
###############################################################################

spp <- read_tsv(here(orthdir, "data/species.tsv")) |>
  select(-ensembl_proteome_fasta)

idmap <- map(
  set_names(spp$abb),
  possibly(function(x) {
    read_tsv(here(
      orthdir,
      "results/id-name-mapping",
      glue("{x}.tsv")
    ))
  })
) |>
  purrr::compact() |>
  bind_rows(.id = "spp") |>
  select(spp, gene_id, symbol) |>
  distinct()

# labels ======================================================================

# Tip label format: <abb>__<symbol>__<gene_id>. Always three fields so labels can
# be parsed by splitting on "__"; gene_id is always present so tips can be linked
# back to the genomes used for single-cell data. Symbols are cleaned to
# characters that are safe in FASTA headers and newick; missing symbols become
# "unnamed".
make_label <- function(.abb, .symbol, .gene_id) {
  symbol <- .symbol |>
    str_replace_all("[^A-Za-z0-9.-]", "-") |>
    replace_na("unnamed")
  paste(.abb, symbol, .gene_id, sep = "__")
}

# seeds =======================================================================
biostring_to_char <- function(x) {
  stopifnot(length(x) == 1)
  as.character(x[[1]])
}

seeds <- read_tsv(here(basedir, "data/seeds/seeds.tsv")) |>
  mutate(
    sequence = map_chr(label, function(.x) {
      bs <- readAAStringSet(here(basedir, glue("data/seeds/{.x}.fasta")))
      biostring_to_char(bs)
    })
  ) |>
  # seed FASTA files are named by the original label (<abb>__<symbol>), so
  # relabel only after reading the sequences. The original label is kept as
  # seed_key because the hits record it in seeded_from.
  mutate(seed_key = label, label = make_label(spp, symbol, gene_id)) |>
  mutate(seqtype = "seed")

## RBH confirmed hits =========================================================

rbh <- read_tsv(here(basedir, "results/blast/hits_rbh-confirmed.tsv")) |>
  mutate(seqtype = "hit") |>
  select(
    gene_id = hit_accession,
    abb = hit_abb,
    seeded_from = seed,
    sequence = hit_sequence,
    seqtype,
    fwd_bitscore,
    seeds_all,
    n_seeds
  ) |>
  left_join(rename(idmap, abb = spp), by = c("abb", "gene_id")) |>
  mutate(label = make_label(abb, symbol, gene_id))

## all sequences together =====================================================
all_seqs <- bind_rows(
  seeds |>
    select(
      label,
      seqtype,
      gene_id,
      abb = spp,
      symbol,
      seed_type,
      rtk_class,
      seed_key,
      sequence
    ),
  rbh |>
    select(
      label,
      seqtype,
      gene_id,
      abb,
      symbol,
      seeded_from,
      fwd_bitscore,
      seeds_all,
      n_seeds,
      sequence
    )
)


###############################################################################
# QC ===============================================================
###############################################################################

stopifnot(!any(duplicated(all_seqs$label)))
stopifnot(!any(duplicated(all_seqs$sequence)))

all_seqs |>
  mutate(has_internal_stop = str_detect(sequence, "\\*")) |>
  filter(has_internal_stop) |>
  nrow()

stopifnot(.Last.value == 0)

###############################################################################
# write =======================================================================
###############################################################################

aa_assembled <- AAStringSet(all_seqs$sequence) |>
  setNames(all_seqs$label)

writeXStringSet(
  aa_assembled,
  here(resdir, "seqs_assembled.fa"),
  width = max(width(aa_assembled))
)

write_tsv(all_seqs, here(resdir, "seqs_assembled.tsv"))

# end #########################################################################
