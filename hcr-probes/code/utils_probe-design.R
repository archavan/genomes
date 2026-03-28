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
library(seqinr)
library(Biostrings)

###############################################################################
# probeset data ================================================================
###############################################################################

split_init_seqs <- data.frame(
  probeset = c(
    "B1",
    "B2",
    "B3",
    "B4",
    "B5",
    "B7",
    "B9",
    "B10",
    "B11",
    "B13",
    "B14",
    "B15",
    "B17"
  ),
  init_first_half = c(
    "GAGGAGGGCAGCAAACGG",
    "CCTCGTAAATCCTCATCA",
    "GTCCCTGCCTCTATATCT",
    "CCTCAACCTACCTCCAAC",
    "CTCACTCCCAATCTCTAT",
    "CTTCAACCTCCACCTACC",
    "CACGTATCTACTCCACTC",
    "CCTCAAGATACTCCTCTA",
    "CGCTTAGATATCACTCCT",
    "AGGTAACGCCTTCCTGCT",
    "AATGTCAATAGCGAGCGA",
    "CAGATTAACACACCACAA",
    "CGATTGTTTGTTGTGGAC"
  ),
  init_second_half = c(
    "GAAGAGTCTTCCTTTACG",
    "ATCATCCAGTAAACCGCC",
    "CCACTCAACTTTAACCCG",
    "TCTCACCATATTCGCTTC",
    "CTACCCTACAAATCCAAT",
    "TCCAATCCCTACCCTCAC",
    "TCAGCACACTCCCAACCC",
    "CCTACTCGACTACCCTAG",
    "ACGTCGACCACACTCATC",
    "TTATGCTCAACATACAAC",
    "CCCTATATTTCTGCACAG",
    "GGTATCTCGAACACTCTC",
    "GCATGCTAATCGGATGAG"
  ),
  first_spacer = c(
    "aa",
    "aa",
    "tt",
    "aa",
    "aa",
    "aa",
    "ta",
    "tt",
    "at",
    "ta",
    "aa",
    "tt",
    "aa"
  ),
  second_spacer = c(
    "ta",
    "ta",
    "tt",
    "at",
    "aa",
    "at",
    "at",
    "at",
    "ta",
    "tt",
    "at",
    "aa",
    "ta"
  ),
  stringsAsFactors = FALSE
)

###############################################################################
# functions ====================================================================
###############################################################################

generate_subseqs <- function(.input) {
  seq_length <- nchar(.input)
  data.frame(
    seqname = names(.input),
    start = seq(1, seq_length - 51, 1),
    end = seq(52, seq_length, 1),
    probe_seq = substring(
      .input[1:length(.input)],
      first = 1:(seq_length - 51),
      last = 52:seq_length
    )
  )
}

complexity_filter <- function(.data) {
  .data |>
    filter(!str_detect(probe_seq, "(.)\\1{4,}|(..)\\2{4,}|(...)\\3{3,}"))
}

get_non_overlapping <- function(.df) {
  results <- slice(.df, 0)
  for (row in 1:nrow(.df)) {
    data <- .df[row, ]
    compare <- tail(results, 1)
    if (row == 1) {
      results <- rbind(results, data)
    } else if (.df[row, "start"] > compare[1, "end"]) {
      results <- rbind(results, data)
    }
  }
  rows <- nrow(results)
  results$pair_ID <- seq(1, rows, 1)
  results
}

make_fulllength_probes <- function(.input) {
  .input |>
    generate_subseqs() |>
    complexity_filter() |>
    get_non_overlapping()
}

rev_comp <- function(.string) {
  toupper(c2s(rev(comp(s2c(toString(.string))))))
}

rev_comp_v <- Vectorize(rev_comp)

make_split_probes <- function(.fl_probe_df, .probeset_choice) {
  split_init_seqs_filtered <- split_init_seqs |>
    filter(probeset == .probeset_choice)

  if (nrow(split_init_seqs_filtered) == 0) {
    stop("Selected hairpin pair not found in the predefined data.")
  }

  probe_seqs <- .fl_probe_df |>
    mutate(probeset = .probeset_choice)

  full_join(probe_seqs, split_init_seqs_filtered, by = "probeset") |>
    separate(probe_seq, into = c("p1", "gap", "p2"), sep = c(25, -25)) |>
    mutate(across(c("p1", "p2"), rev_comp_v, .names = "{.col}_rc")) |>
    rename(probe_1_seq = p1_rc, probe_2_seq = p2_rc) |>
    mutate(
      odd_ID = paste(pair_ID, "1", sep = "."),
      even_ID = paste(pair_ID, "2", sep = ".")
    ) |>
    select(
      seqname,
      probeset,
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

make_oPools_df <- function(
  .split_probe_df,
  .limit = nrow(.split_probe_df) * 2
) {
  .split_probe_df |>
    mutate(
      seq_1 = glue("{init_first_half}{tolower(first_spacer)}{probe_1_seq}"),
      seq_2 = glue("{probe_2_seq}{tolower(second_spacer)}{init_second_half}")
    ) |>
    pivot_longer(cols = starts_with("seq_"), values_to = "Sequence") |>
    slice_head(n = .limit) |>
    mutate(
      Pool_name = paste(
        seqname,
        min(start),
        max(end),
        probeset,
        glue("{n() / 2}pp"),
        sep = "_"
      )
    ) |>
    select(Pool_name, Sequence)
}
