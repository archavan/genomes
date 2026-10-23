# Arun Chavan
# Started: 2026-10-23

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
  dplyr::first,
  dplyr::rename
)

basedir <- here("gene-trees/rtk/full-gene")
orthdir <- here("orthofinder/2026-08-07")
assembdir <- here(basedir, "results/05a_seqs-assembled")
hmmscandir <- here(basedir, "results/05b_hmmscan")
resdir <- here(basedir, "results/05c_seqs-qc")
fs::dir_create(resdir)

###############################################################################
# data ========================================================================
###############################################################################

seqs <- read_tsv(here(assembdir, "seqs_assembled.tsv"))

# hmmscan --domtblout: target = model, query = sequence. Envelope coordinates
# (env_from/env_to) are positions on the query protein. We read the 21 fixed
# fields and discard the free-text description that follows.
parse_domtbl <- function(.path) {
  cols <- c(
    "model",
    "model_acc",
    "model_len",
    "seq",
    "seq_acc",
    "seq_len",
    "evalue_full",
    "score_full",
    "bias_full",
    "dom_n",
    "dom_of",
    "c_evalue",
    "i_evalue",
    "dom_score",
    "dom_bias",
    "hmm_from",
    "hmm_to",
    "ali_from",
    "ali_to",
    "env_from",
    "env_to"
  )
  read_lines(.path) |>
    discard(function(l) str_starts(l, "#")) |>
    map(function(l) {
      f <- str_split_1(str_squish(l), " ")
      set_names(as.list(f[seq_along(cols)]), cols)
    }) |>
    bind_rows() |>
    type_convert()
}

dom <- parse_domtbl(here(hmmscandir, "domains.domtbl")) |>
  rename(domain = model)

###############################################################################
# QC ==========================================================================
###############################################################################

## check for non-standard characters ==========================================

standard_aa <- "ACDEFGHIKLMNPQRSTVWY"

seqs |>
  filter(str_detect(sequence, glue("[^{standard_aa}]"))) |>
  mutate(
    nonstandard_chars = str_extract_all(sequence, glue("[^{standard_aa}]")),
    nonstandard_chars = map_chr(nonstandard_chars, function(x) {
      paste(sort(unique(x)), collapse = ", ")
    })
  ) |>
  select(seqtype, label, nonstandard_chars)

seqs |>
  filter(str_detect(sequence, glue("[^{standard_aa}]"))) |>
  distinct(label, sequence) |>
  mutate(
    n_X = str_count(sequence, "X"),
    seq_len = nchar(sequence),
    pct_X = round(100 * n_X / seq_len, 2)
  ) |>
  select(label, seq_len, n_X, pct_X)

# Drop low-quality sequences: >1% of residues are X (ambiguous).
seqs <- seqs |>
  filter(str_count(sequence, "X") / nchar(sequence) <= 0.01)

## check lengths ==============================================================

# hits record the seed that retrieved them in seeded_from, which matches the
# seed_key (not the label) of the seed rows
seed_lengths <- seqs |>
  filter(seqtype == "seed") |>
  transmute(seed_key, seed_length = nchar(sequence))

hit_lengths <- seqs |>
  filter(seqtype == "hit") |>
  mutate(hit_length = nchar(sequence)) |>
  left_join(seed_lengths, by = c("seeded_from" = "seed_key")) |>
  mutate(hit_to_seed_length_ratio = hit_length / seed_length) |>
  select(
    label,
    gene_id,
    abb,
    seeded_from,
    hit_length,
    seed_length,
    hit_to_seed_length_ratio
  )

stopifnot(!anyNA(hit_lengths$seed_length))

hit_lengths |>
  ggplot(aes(log2(hit_to_seed_length_ratio))) +
  geom_histogram()

hit_lengths |>
  filter(hit_to_seed_length_ratio < 0.5 | hit_to_seed_length_ratio > 2) |>
  print(n = Inf, width = Inf)

## domain architecture ========================================================

# Collapse all hits per sequence into a left-to-right domain string.
# --cut_ga means every hit is above the curated gathering threshold.
architecture <- dom |>
  arrange(seq, env_from) |>
  mutate(domain = paste0("(", domain, ")")) |>
  summarize(architecture = str_c(domain, collapse = " "), .by = seq)

seqs <- seqs |>
  left_join(architecture, by = c("label" = "seq")) |>
  left_join(
    select(hit_lengths, label, hit_to_seed_length_ratio),
    by = "label"
  ) |>
  mutate(
    has_kinase = coalesce(
      str_detect(architecture, fixed("(PK_Tyr_Ser-Thr)")),
      FALSE
    )
  )

seqs |>
  filter(!has_kinase) |>
  select(label, seqtype, seeded_from, hit_to_seed_length_ratio, architecture) |>
  print(n = Inf, width = Inf)

## remove sequences without a kinase domain ===================================

seqs_passed <- seqs |>
  filter(has_kinase)

###############################################################################
# write =======================================================================
###############################################################################

aa_passed <- AAStringSet(seqs_passed$sequence) |>
  setNames(seqs_passed$label)

writeXStringSet(
  aa_passed,
  here(resdir, "seqs_passed.fa"),
  width = max(width(aa_passed))
)

seqs_passed |>
  select(-has_kinase) |>
  write_tsv(here(resdir, "seqs_passed.tsv"))

architecture |>
  write_tsv(here(resdir, "architecture.tsv"))

# end #########################################################################
