# HCR v3 Probe Design Functions
#
# Functions for designing split-initiator HCR v3 probes from transcript
# sequences. Given a DNA sequence, generates candidate probe pairs: selects
# 52-nt target regions, splits them into two 25-nt halves, reverse complements
# each half, and appends hairpin-specific initiator sequences.
#
# Adapted from Nat Clarke's HCR Probemaker Shiny app:
# https://clarkelab.com/apps/hcr/
#
# Original code by Nat Clarke (Clarke Lab). Modified to extract functions
# for use in a non-interactive pipeline.

###############################################################################
# packages =====================================================================
###############################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(Biostrings)

###############################################################################
# amplifier data ================================================================
###############################################################################

split_init_seqs <- read_tsv(
  here("hcr-probes/data/hcr-v3-initiator-sequences.tsv"),
  show_col_types = FALSE
)

###############################################################################
# functions ====================================================================
###############################################################################

generate_subseqs <- function(.input) {
  seq_length <- nchar(.input)
  if (seq_length < 52) {
    stop(glue(
      "Transcript is too short for probe design ({seq_length} nt). ",
      "Minimum length is 52 nt."
    ))
  }
  data.frame(
    seqname = names(.input),
    start = seq(1, seq_length - 51, 1),
    end = seq(52, seq_length, 1),
    probe_seq = substring(
      .input,
      first = 1:(seq_length - 51),
      last = 52:seq_length
    )
  )
}

# (.)\\1{4,}  = 5+ identical bases in a row (e.g. AAAAA)
# (..)\\2{4,} = 5+ dinucleotide repeats (e.g. ATATATATAT)
# (...)\\3{3,} = 4+ trinucleotide repeats (e.g. ATGATGATGATG)
complexity_filter <- function(.data) {
  .data |>
    filter(!str_detect(probe_seq, "(.)\\1{4,}|(..)\\2{4,}|(...)\\3{3,}"))
}

get_non_overlapping <- function(.df) {
  results <- slice(.df, 0)
  for (row in seq_len(nrow(.df))) {
    data <- .df[row, ]
    compare <- tail(results, 1)
    if (row == 1) {
      results <- rbind(results, data)
    } else if (.df[row, "start"] > compare[1, "end"]) {
      results <- rbind(results, data)
    }
  }
  rows <- nrow(results)
  results$pair_ID <- seq_len(rows)
  results
}

make_fulllength_probes <- function(.input) {
  .input |>
    generate_subseqs() |>
    complexity_filter() |>
    get_non_overlapping()
}

rev_comp_v <- function(.strings) {
  as.character(reverseComplement(DNAStringSet(.strings)))
}

make_split_probes <- function(.fl_probe_df, .amplifier_choice) {
  split_init_seqs_filtered <- split_init_seqs |>
    filter(amplifier == .amplifier_choice)

  if (nrow(split_init_seqs_filtered) == 0) {
    stop("Selected hairpin pair not found in the predefined data.")
  }

  probe_seqs <- .fl_probe_df |>
    mutate(amplifier = .amplifier_choice)

  # Split 52-nt probe into 25-nt left / 2-nt gap / 25-nt right at positions 25
  # and 27. The gap is discarded; each half is reverse-complemented for
  # hybridization.
  left_join(probe_seqs, split_init_seqs_filtered, by = "amplifier") |>
    separate(probe_seq, into = c("p1", "gap", "p2"), sep = c(25, -25)) |>
    mutate(
      probe_1_seq = rev_comp_v(p1),
      probe_2_seq = rev_comp_v(p2)
    ) |>
    mutate(
      odd_ID = paste(pair_ID, "1", sep = "."),
      even_ID = paste(pair_ID, "2", sep = ".")
    ) |>
    select(
      seqname,
      amplifier,
      start,
      end,
      pair_ID,
      odd_ID,
      init_first_half,
      first_spacer,
      probe_1_seq,
      even_ID,
      probe_2_seq,
      second_spacer,
      init_second_half
    )
}

# IDT oPools Base 1 tier is a flat $125 for up to 3300 bases per pool.
# At 45 bases/oligo, 66 oligos (33 probe pairs) = 2970 bases, staying
# within that tier. Beyond 3300 bases, pricing switches to per-base.
# https://www.idtdna.com/pages/products/custom-dna-rna/dna-oligos/custom-dna-oligos/opools-oligo-pools
make_oPools_df <- function(
  .split_probe_df,
  .gene_id,
  .species,
  .max_pairs = 33
) {
  # 33 pairs = 66 oligos x 45 bases = 2970 bases, within IDT oPools Base 1 tier.
  # Uppercase = initiator + probe sequence; lowercase = spacer (IDT convention).
  .split_probe_df |>
    slice_head(n = .max_pairs) |>
    mutate(
      seq_1 = glue("{init_first_half}{tolower(first_spacer)}{probe_1_seq}"),
      seq_2 = glue("{probe_2_seq}{tolower(second_spacer)}{init_second_half}")
    ) |>
    pivot_longer(cols = starts_with("seq_"), values_to = "Sequence") |>
    mutate(
      Pool_name = glue("{.species}_{.gene_id}_{amplifier}_{n() / 2}pp")
    ) |>
    select(Pool_name, Sequence)
}
