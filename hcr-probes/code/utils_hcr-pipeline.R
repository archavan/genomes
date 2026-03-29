# HCR Probe Design Pipeline
#
# Main pipeline for designing HCR v3 probes with BLAST-based specificity
# checking. Processes one gene at a time, saving intermediate files for
# manual inspection.
#
# Usage:
#   source(here("hcr-probes/code/utils_hcr-pipeline.R"), local = TRUE)
#   design_probes("ncbi:118414000", species = "bflo", expression = expr_df)
#
# See hcr-probes/docs/how-hcr-v3-probes-work.md for background.

###############################################################################
# packages ====================================================================
###############################################################################

library(tidyverse)
library(Biostrings)
library(glue)
library(writexl)
library(here)
library(conflicted)

conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::slice,
  dplyr::rename,
  dplyr::desc,
  dplyr::collapse,
  dplyr::count
)

source(here("hcr-probes/code/utils_probe-design.R"), local = TRUE)

###############################################################################
# species config ==============================================================
###############################################################################

get_species_config <- function(.species) {
  bflo_dir <- here("species/bflo/VNyyK")
  crob_dir <- here("species/crob/KY21")
  scla_dir <- here("species/scla/refseq_v1")

  # fmt: skip
  configs <- list(
    bflo = list(
      transcriptome = here(bflo_dir, "data/refseq/GCF_000003815.2_Bfl_VNyyK_rna.fna.gz"),
      mapping = here(bflo_dir, "results/genes/gid-tid_mapping.tsv"),
      blastdb = here(bflo_dir, "results/blastdb/GCF_000003815.2_Bfl_VNyyK_rna/db")
    ),
    crob = list(
      transcriptome = here(crob_dir, "results/fasta/KY21_transcripts_with-mito.fasta.gz"),
      mapping = here(crob_dir, "results/genes/gid-tid_mapping.tsv"),
      blastdb = here(crob_dir, "results/blastdb/KY21_transcripts_with-mito/db")
    ),
    scla = list(
      transcriptome = here(scla_dir, "data/refseq/GCF_013122585.1_ASM1312258v2_rna.fna.gz"),
      mapping = here(scla_dir, "results/genes/gid-tid_mapping.tsv"),
      blastdb = here(scla_dir, "results/blastdb/GCF_013122585.1_ASM1312258v2_rna/db")
    )
  )

  if (!.species %in% names(configs)) {
    stop(glue("Species '{.species}' not recognized."))
  }

  configs[[.species]]
}

###############################################################################
# helpers =====================================================================
###############################################################################

sanitize_gene_id <- function(.gene_id) {
  gsub(":", ".", .gene_id)
}

###############################################################################
# list_transcripts ============================================================
###############################################################################

#' List all transcripts for a gene with their lengths
#'
#' @param .gene_id Character. Gene ID (e.g. "ncbi:118414000")
#' @param .species Character. Species abbreviation ("bflo", "crob", "scla")
#' @param .mapping Optional pre-loaded mapping tibble
#' @param .txome Optional pre-loaded DNAStringSet (names already cleaned)
#'
#' @return Tibble with transcript_id, length, and any other columns from mapping
list_transcripts <- function(
  .gene_id,
  .species,
  .mapping = NULL,
  .txome = NULL
) {
  if (is.null(.mapping) || is.null(.txome)) {
    config <- get_species_config(.species)
  }
  if (is.null(.mapping)) {
    .mapping <- read_tsv(config$mapping, show_col_types = FALSE)
  }
  if (is.null(.txome)) {
    .txome <- readDNAStringSet(config$transcriptome)
    names(.txome) <- str_replace(names(.txome), "\\s.*", "")
  }

  gene_txs <- .mapping |>
    filter(gene_id == .env$.gene_id)

  if (nrow(gene_txs) == 0) {
    stop(glue(
      "Gene '{.gene_id}' not found in mapping for species '{.species}'."
    ))
  }

  tx_lengths <- tibble(
    transcript_id = names(.txome),
    length = width(.txome)
  )

  gene_txs |>
    left_join(tx_lengths, by = "transcript_id") |>
    arrange(desc(length))
}

###############################################################################
# step 1: select transcript ===================================================
###############################################################################

select_transcript <- function(
  .gene_id,
  .species,
  .mapping,
  .txome,
  .transcript_id = NULL,
  .expression = NULL
) {
  if (is.null(.transcript_id) && is.null(.expression)) {
    stop("Provide either .expression or .transcript_id.")
  }
  if (!.gene_id %in% .mapping$gene_id) {
    stop(glue("Gene '{.gene_id}' not found in mapping for '{.species}'."))
  }

  gene_txs <- .mapping |>
    filter(gene_id == .env$.gene_id)

  tx_lengths <- tibble(
    transcript_id = names(.txome),
    length = width(.txome)
  )

  selection <- gene_txs |>
    left_join(tx_lengths, by = "transcript_id")

  if (!is.null(.transcript_id)) {
    if (!.transcript_id %in% selection$transcript_id) {
      stop(glue(
        "Transcript '{.transcript_id}' not found for gene '{.gene_id}'."
      ))
    }
    chosen_tid <- .transcript_id
    selection <- selection |>
      mutate(
        tpm = NA_real_,
        selected = transcript_id == chosen_tid,
        selected_by = "manual"
      )
  }

  if (is.null(.transcript_id)) {
    stopifnot(all(c("transcript_id", "tpm") %in% names(.expression)))
    selection <- selection |>
      left_join(
        select(.expression, transcript_id, tpm),
        by = "transcript_id"
      )
    chosen_tid <- selection |>
      filter(!is.na(tpm)) |>
      slice_max(order_by = tpm, n = 1, with_ties = FALSE) |>
      pull(transcript_id)
    if (length(chosen_tid) == 0) {
      stop(glue("No expression data matched transcripts of '{.gene_id}'."))
    }
    selection <- selection |>
      mutate(
        selected = transcript_id == chosen_tid,
        selected_by = "tpm"
      )
  }

  list(
    selection_table = selection,
    transcript_id = chosen_tid,
    sequence = .txome[chosen_tid]
  )
}

###############################################################################
# step 2: design probes =======================================================
###############################################################################

design_probe_candidates <- function(.sequence, .amplifier) {
  make_fulllength_probes(.sequence) |>
    make_split_probes(.amplifier)
}

probes_to_tidy <- function(.split_probes) {
  half1 <- .split_probes |>
    select(id = half1_id, seq = half1_seq)
  half2 <- .split_probes |>
    select(id = half2_id, seq = half2_seq)
  bind_rows(half1, half2) |>
    arrange(id)
}

###############################################################################
# step 3: BLAST ===============================================================
###############################################################################

run_blast <- function(.query_fasta, .blast_db, .evalue = 100) {
  blast_fields <- c(
    "qseqid",
    "sseqid",
    "pident",
    "qcovhsp",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sstrand",
    "evalue",
    "bitscore",
    "qseq",
    "sseq"
  )

  out_file <- tempfile(fileext = ".tsv")

  cmd <- glue(
    "blastn",
    " -query {shQuote(.query_fasta)}",
    " -db {shQuote(.blast_db)}",
    " -outfmt {shQuote(paste('6', paste(blast_fields, collapse = ' ')))}",
    " -task blastn-short",
    " -evalue {.evalue}",
    " -max_hsps 1",
    " -max_target_seqs 100000",
    " -strand both",
    " -out {shQuote(out_file)}"
  )

  exit_code <- system(cmd)
  if (exit_code != 0) {
    stop(glue(
      "BLAST exited with status {exit_code}. ",
      "Check that blastn is installed and the database exists."
    ))
  }

  # If no hits, return empty tibble with correct column names.
  if (file.info(out_file)$size == 0) {
    return(tibble(
      !!!set_names(rep(list(character(0)), length(blast_fields)), blast_fields)
    ))
  }

  read_tsv(out_file, col_names = blast_fields, show_col_types = FALSE)
}

###############################################################################
# step 4: MCH calculation and alignment display ===============================
###############################################################################

# MCH (Max Continuous Homology): longest run of consecutive matching bases
# between query and subject. Used to flag off-target probe hits — a high MCH
# means the probe could cross-hybridize even if overall percent identity is low.
calculate_mch <- function(.qseq, .sseq) {
  q_chars <- strsplit(.qseq, "")[[1]]
  s_chars <- strsplit(.sseq, "")[[1]]
  max_run <- 0L
  current_run <- 0L
  for (i in seq_along(q_chars)) {
    if (toupper(q_chars[i]) == toupper(s_chars[i])) {
      current_run <- current_run + 1L
      if (current_run > max_run) max_run <- current_run
    } else {
      current_run <- 0L
    }
  }
  max_run
}


###############################################################################
# step 5: annotate and filter =================================================
###############################################################################

annotate_blast_results <- function(
  .blast_results,
  .gene_id,
  .species,
  .mapping
) {
  target_tids <- .mapping |>
    filter(gene_id == .env$.gene_id) |>
    pull(transcript_id)

  # build transcript-to-gene lookup
  tid_to_gene <- .mapping |>
    distinct(transcript_id, gene_id) |>
    dplyr::rename(hit_gene = gene_id)

  .blast_results |>
    mutate(sseqid = str_replace(sseqid, "\\s.*", "")) |>
    left_join(tid_to_gene, by = c("sseqid" = "transcript_id")) |>
    mutate(
      hit_gene = replace_na(hit_gene, "unknown"),
      mch = map2_int(qseq, sseq, calculate_mch),
      on_target = sseqid %in% target_tids,
      pair_id = as.integer(str_replace(qseqid, "\\.\\d+$", "")),
      half = as.integer(str_replace(qseqid, "^\\d+\\.", ""))
    )
}

filter_probes <- function(.blast_annotated, .split_probes, .mch_threshold) {
  on_target_ids <- .blast_annotated |>
    filter(on_target) |>
    distinct(qseqid) |>
    pull(qseqid)

  flagged_ids <- .blast_annotated |>
    filter(!on_target, mch >= .mch_threshold) |>
    distinct(qseqid) |>
    pull(qseqid)

  .split_probes |>
    select(pair_id, half1_id, half2_id) |>
    mutate(
      half1_on_target = half1_id %in% on_target_ids,
      half2_on_target = half2_id %in% on_target_ids,
      both_on_target = half1_on_target & half2_on_target,
      half1_flagged = half1_id %in% flagged_ids,
      half2_flagged = half2_id %in% flagged_ids,
      any_flagged = half1_flagged | half2_flagged,
      passed = both_on_target & !any_flagged
    )
}

###############################################################################
# main pipeline ===============================================================
###############################################################################

#' Design HCR v3 probes for a single gene with BLAST specificity checking
#'
#' @param .gene_id Character. Gene ID (e.g. "ncbi:118414000")
#' @param .species Character. Species abbreviation ("bflo", "crob", "scla")
#' @param .amplifier Character. HCR amplifier (e.g. "B1", "B2")
#' @param .transcript_id Character or NULL. Manual transcript override
#' @param .expression Data frame or NULL. Must have transcript_id + tpm columns
#' @param .evalue Numeric. BLAST e-value threshold. Deliberately permissive
#'   (default 100) to catch weak cross-hybridization; MCH filtering handles
#'   specificity downstream.
#' @param .mch_threshold Integer. MCH threshold for off-target flagging
#' @param .outdir Character. Output directory for results
#'
#' @return Invisible list with all intermediate results
design_probes <- function(
  .gene_id,
  .species,
  .amplifier = "B1",
  .transcript_id = NULL,
  .expression = NULL,
  .evalue = 100,
  .mch_threshold = 18,
  .outdir = here("hcr-probes/results")
) {
  config <- get_species_config(.species)
  gene_dir <- here(.outdir, sanitize_gene_id(.gene_id))
  fs::dir_create(gene_dir)

  mapping <- read_tsv(config$mapping, show_col_types = FALSE)
  txome <- readDNAStringSet(config$transcriptome)
  names(txome) <- str_replace(names(txome), "\\s.*", "")

  # step 1: select transcript -------------------------------------------------
  message(glue("Step 1: Selecting transcript for {.gene_id}..."))

  tx <- select_transcript(
    .gene_id = .gene_id,
    .species = .species,
    .transcript_id = .transcript_id,
    .expression = .expression,
    .mapping = mapping,
    .txome = txome
  )

  write_tsv(
    tx$selection_table,
    here(gene_dir, "01_transcript-selection.tsv")
  )
  writeXStringSet(
    tx$sequence,
    here(gene_dir, "02_selected-transcript.fasta")
  )

  message(glue(
    "  Selected transcript: {tx$transcript_id} (length: {width(tx$sequence)} nt)"
  ))

  # step 2: design probes -----------------------------------------------------
  message("Step 2: Designing probe candidates...")

  split_probes <- design_probe_candidates(tx$sequence, .amplifier)
  message(glue("  Generated {nrow(split_probes)} probe pairs"))

  write_tsv(split_probes, here(gene_dir, "03_probe-candidates.tsv"))

  probe_fasta_data <- probes_to_tidy(split_probes)
  probe_seqs <- DNAStringSet(probe_fasta_data$seq)
  names(probe_seqs) <- probe_fasta_data$id
  probe_fasta_path <- here(gene_dir, "03_probe-candidates.fasta")
  writeXStringSet(probe_seqs, probe_fasta_path)

  # step 3: BLAST -------------------------------------------------------------
  message("Step 3: Running BLAST specificity check...")

  blast_results <- run_blast(probe_fasta_path, config$blastdb, .evalue = .evalue)
  message(glue("  BLAST returned {nrow(blast_results)} hits"))

  # step 4: annotate blast results ---------------------------------------------
  message("Step 4: Annotating BLAST results...")

  blast_annotated <- annotate_blast_results(
    blast_results,
    .gene_id,
    .species,
    mapping
  )
  write_tsv(blast_annotated, here(gene_dir, "04_blast-results.tsv"))

  # step 5: filter and output --------------------------------------------------
  message("Step 5: Filtering probes...")

  pair_summary <- filter_probes(blast_annotated, split_probes, .mch_threshold)

  probes_annotated <- split_probes |>
    left_join(
      select(
        pair_summary,
        pair_id,
        half1_on_target,
        half2_on_target,
        both_on_target,
        half1_flagged,
        half2_flagged,
        any_flagged,
        passed
      ),
      by = "pair_id"
    )
  write_tsv(probes_annotated, here(gene_dir, "05_probes-annotated.tsv"))

  passing_probes <- split_probes |>
    filter(pair_id %in% pair_summary$pair_id[pair_summary$passed])

  if (nrow(passing_probes) > 0) {
    opools <- make_oPools_df(passing_probes, sanitize_gene_id(.gene_id), .species)
    write_tsv(opools, here(gene_dir, "06_probes-passing.tsv"))
    opools_xlsx <- opools |> rename(`Pool name` = Pool_name)
    write_xlsx(opools_xlsx, here(gene_dir, "06_probes-passing.xlsx"))
  } else {
    write_tsv(
      tibble(Pool_name = character(), Sequence = character()),
      here(gene_dir, "06_probes-passing.tsv")
    )
  }

  n_total <- nrow(pair_summary)
  n_on_target <- sum(pair_summary$both_on_target)
  n_passed <- sum(pair_summary$passed)

  message(glue(
    "\n  {.gene_id}: {n_passed}/{n_total} probe pairs passed",
    " ({n_on_target} on-target,",
    " {n_total - n_passed} flagged or missing on-target)"
  ))
  message(glue("  Results written to: {gene_dir}"))

  invisible(list(
    gene_id = .gene_id,
    transcript_id = tx$transcript_id,
    selection_table = tx$selection_table,
    split_probes = split_probes,
    blast_results = blast_results,
    blast_annotated = blast_annotated,
    pair_summary = pair_summary,
    probes_annotated = probes_annotated,
    n_passed = n_passed,
    n_total = n_total,
    outdir = gene_dir
  ))
}

###############################################################################
# combine for ordering ========================================================
###############################################################################

#' Combine passing probes from multiple genes into a single IDT oPools xlsx
#'
#' @param .gene_ids Character vector. Gene IDs that have already been run
#'   through design_probes()
#' @param .results_dir Character. Directory containing per-gene results
#' @param .out_path Character. File path for the combined xlsx
#'
#' @return Writes xlsx and returns combined tibble invisibly
combine_opools <- function(
  .gene_ids,
  .results_dir = here("hcr-probes/results"),
  .out_path = here("hcr-probes/results/oPools-order.xlsx")
) {
  combined <- map(.gene_ids, \(gid) {
    gene_dir <- here(.results_dir, sanitize_gene_id(gid))
    tsv_path <- here(gene_dir, "06_probes-passing.tsv")
    df <- read_tsv(tsv_path, show_col_types = FALSE)
    if (nrow(df) == 0) {
      stop(glue("No passing probes for {gid}."))
    }
    df
  }) |>
    bind_rows()

  fs::dir_create(dirname(.out_path))
  combined_xlsx <- combined |> rename(`Pool name` = Pool_name)
  write_xlsx(combined_xlsx, .out_path)

  invisible(combined)
}
