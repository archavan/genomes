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
library(R.utils)
library(glue)
library(here)

source(here("hcr-probes/code/utils_probe-design.R"), local = TRUE)

###############################################################################
# species config ==============================================================
###############################################################################

get_species_config <- function(.species) {
  configs <- list(
    bflo = list(
      transcriptome = here(
        "species/bflo/VNyyK/data/refseq",
        "GCF_000003815.2_Bfl_VNyyK_rna.fna.gz"
      ),
      mapping = here("species/bflo/VNyyK/results/genes/gid-tid_mapping.tsv"),
      blastdb_dir = here("species/bflo/VNyyK/results/blastdb")
    ),
    crob = list(
      transcriptome = here(
        "species/crob/KY21/results/fasta",
        "KY21_transcripts_with-mito.fasta.gz"
      ),
      mapping = here("species/crob/KY21/results/genes/gid-tid_mapping.tsv"),
      blastdb_dir = here("species/crob/KY21/results/blastdb")
    ),
    scla = list(
      transcriptome = here(
        "species/scla/refseq_v1/data/refseq",
        "GCF_013122585.1_ASM1312258v2_rna.fna.gz"
      ),
      mapping = here(
        "species/scla/refseq_v1/results/genes/gid-tid_mapping.tsv"
      ),
      blastdb_dir = here("species/scla/refseq_v1/results/blastdb")
    )
  )

  if (!.species %in% names(configs)) {
    stop(
      "Unknown species: '",
      .species,
      "'. ",
      "Available: ",
      paste(names(configs), collapse = ", ")
    )
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
#'
#' @return Tibble with transcript_id, length, and any other columns from mapping
list_transcripts <- function(.gene_id, .species) {
  config <- get_species_config(.species)
  mapping <- read_tsv(config$mapping, show_col_types = FALSE)
  txome <- readDNAStringSet(config$transcriptome)

  names(txome) <- str_replace(names(txome), "\\s.*", "")

  gene_txs <- mapping |>
    filter(gene_id == .env$.gene_id)

  if (nrow(gene_txs) == 0) {
    stop(glue(
      "Gene '{.gene_id}' not found in mapping for species '{.species}'."
    ))
  }

  tx_lengths <- tibble(
    transcript_id = names(txome),
    length = width(txome)
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
  .transcript_id = NULL,
  .expression = NULL
) {
  config <- get_species_config(.species)
  mapping <- read_tsv(config$mapping, show_col_types = FALSE)
  txome <- readDNAStringSet(config$transcriptome)
  names(txome) <- str_replace(names(txome), "\\s.*", "")

  gene_txs <- mapping |>
    filter(gene_id == .env$.gene_id)

  if (nrow(gene_txs) == 0) {
    stop(glue(
      "Gene '{.gene_id}' not found in mapping for species '{.species}'."
    ))
  }

  tx_lengths <- tibble(
    transcript_id = names(txome),
    length = width(txome)
  )

  selection <- gene_txs |>
    left_join(tx_lengths, by = "transcript_id")

  if (!is.null(.transcript_id)) {
    if (!.transcript_id %in% selection$transcript_id) {
      stop(glue(
        "Transcript '{.transcript_id}' not found for gene '{.gene_id}'.\n",
        "Available: {paste(selection$transcript_id, collapse = ', ')}"
      ))
    }
    selection <- selection |>
      mutate(
        tpm = NA_real_,
        selected = transcript_id == .env$.transcript_id,
        selected_by = "manual"
      )
    chosen_tid <- .transcript_id
  } else if (!is.null(.expression)) {
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
      stop(glue(
        "No expression data matched transcripts of gene '{.gene_id}'.\n",
        "Available transcripts: {paste(selection$transcript_id, collapse = ', ')}"
      ))
    }
    selection <- selection |>
      mutate(
        selected = transcript_id == chosen_tid,
        selected_by = "tpm"
      )
  } else {
    stop(glue(
      "Provide either an expression data frame or an explicit transcript_id.\n",
      "Use list_transcripts('{.gene_id}', '{.species}') to see available transcripts."
    ))
  }

  chosen_seq <- txome[chosen_tid]

  list(
    selection_table = selection,
    transcript_id = chosen_tid,
    sequence = chosen_seq
  )
}

###############################################################################
# step 2: design probes =======================================================
###############################################################################

design_probe_candidates <- function(.sequence, .hairpin) {
  fl_probes <- make_fulllength_probes(.sequence)

  if (nrow(fl_probes) == 0) {
    stop("No valid probe candidates found for this transcript.")
  }

  make_split_probes(fl_probes, .hairpin)
}

probes_to_fasta <- function(.split_probes) {
  half1 <- .split_probes |>
    select(id = odd_ID, seq = probe_1_seq)
  half2 <- .split_probes |>
    select(id = even_ID, seq = probe_2_seq)
  bind_rows(half1, half2) |>
    arrange(id)
}

###############################################################################
# step 3: BLAST ===============================================================
###############################################################################

ensure_blast_db <- function(.config) {
  db_path <- file.path(.config$blastdb_dir, "txome")
  db_file <- glue("{db_path}.ndb")

  if (file.exists(db_file)) {
    message(glue("BLAST database already exists: {db_path}"))
    return(db_path)
  }

  fs::dir_create(.config$blastdb_dir)

  # decompress if gzipped (makeblastdb doesn't read .gz)
  txome_path <- .config$transcriptome
  if (grepl("\\.gz$", txome_path)) {
    decompressed <- file.path(.config$blastdb_dir, "txome_input.fasta")
    if (!file.exists(decompressed)) {
      message("Decompressing transcriptome for BLAST database...")
      gunzip(txome_path, destname = decompressed, remove = FALSE)
    }
    txome_path <- decompressed
  }

  message("Building BLAST database...")
  cmd <- glue(
    "makeblastdb -in {shQuote(txome_path)} -dbtype nucl -out {shQuote(db_path)}"
  )
  result <- system(cmd, intern = TRUE)
  if (!file.exists(db_file)) {
    stop(glue(
      "Failed to build BLAST database. Output:\n{paste(result, collapse = '\n')}"
    ))
  }
  message(glue("BLAST database built: {db_path}"))
  db_path
}

run_blast <- function(.query_fasta, .blast_db, .evalue = 100) {
  outfmt_fields <- paste(
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
    "sseq",
    collapse = " "
  )

  out_file <- tempfile(fileext = ".tsv")

  cmd <- glue(
    "blastn",
    " -query {shQuote(.query_fasta)}",
    " -db {shQuote(.blast_db)}",
    " -outfmt {shQuote(paste('6', outfmt_fields))}",
    " -task blastn-short",
    " -evalue {.evalue}",
    " -max_hsps 1",
    " -max_target_seqs 100000",
    " -strand both",
    " -out {shQuote(out_file)}"
  )

  system(cmd)

  col_names <- c(
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

  if (file.info(out_file)$size == 0) {
    return(tibble(
      !!!set_names(rep(list(character(0)), length(col_names)), col_names)
    ))
  }

  read_tsv(out_file, col_names = col_names, show_col_types = FALSE)
}

###############################################################################
# step 4: MCH calculation and alignment display ===============================
###############################################################################

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

calculate_mch_v <- function(.qseq, .sseq) {
  map2_int(.qseq, .sseq, calculate_mch)
}

format_alignment <- function(
  .qseqid,
  .sseqid,
  .qseq,
  .sseq,
  .qstart,
  .qend,
  .sstart,
  .send,
  .pident,
  .mch
) {
  q_chars <- strsplit(.qseq, "")[[1]]
  s_chars <- strsplit(.sseq, "")[[1]]
  midline <- map2_chr(q_chars, s_chars, \(q, s) {
    if (toupper(q) == toupper(s)) "|" else " "
  })

  q_label <- sprintf("  Query  %-4d ", .qstart)
  s_label <- sprintf("  Sbjct  %-4d ", .sstart)
  m_pad <- strrep(" ", nchar(q_label))

  n_match <- sum(toupper(q_chars) == toupper(s_chars))
  identity_str <- glue("{n_match}/{length(q_chars)} ({round(.pident)}%)")

  glue(
    "{q_label}{.qseq}  {.qend}\n",
    "{m_pad}{paste(midline, collapse = '')}\n",
    "{s_label}{.sseq}  {.send}\n",
    "  MCH: {.mch}   Identity: {identity_str}"
  )
}

generate_alignment_report <- function(
  .blast_annotated,
  .split_probes,
  .gene_id,
  .mapping,
  .mch_threshold
) {
  lines <- character()
  lines <- c(
    lines,
    glue(
      "BLAST Alignment Report for {.gene_id}\n",
      "MCH threshold: {.mch_threshold}\n",
      "{strrep('=', 70)}\n"
    )
  )

  for (pid in sort(unique(.split_probes$pair_ID))) {
    pair <- .split_probes |> filter(pair_ID == pid)
    pair_start <- pair$start
    pair_end <- pair$end

    pair_annot <- .blast_annotated |>
      filter(pair_ID == pid)

    has_offtarget <- any(
      pair_annot$off_target & pair_annot$mch >= .mch_threshold,
      na.rm = TRUE
    )
    status <- if (has_offtarget) "FLAGGED" else "PASSED"

    lines <- c(
      lines,
      glue(
        "\n=== Probe pair {pid} (position {pair_start}-{pair_end}) -- {status} ===\n"
      )
    )

    for (half in c(1, 2)) {
      half_id <- if (half == 1) pair$odd_ID else pair$even_ID
      half_hits <- .blast_annotated |>
        filter(qseqid == half_id)

      lines <- c(lines, glue("\n  Probe half {half} ({half_id}):\n"))

      if (nrow(half_hits) == 0) {
        lines <- c(lines, "    No BLAST hits.\n")
        next
      }

      on_target_hits <- half_hits |> filter(!off_target)
      off_target_hits <- half_hits |> filter(off_target)

      if (nrow(on_target_hits) > 0) {
        tids <- paste(unique(on_target_hits$sseqid), collapse = ", ")
        lines <- c(lines, glue("    On-target hits: {tids}\n"))
      }

      if (nrow(off_target_hits) > 0) {
        for (j in seq_len(nrow(off_target_hits))) {
          hit <- off_target_hits[j, ]
          flag <- if (hit$mch >= .mch_threshold) " *** FLAGGED ***" else ""
          hit_gene <- .mapping |>
            filter(transcript_id == hit$sseqid) |>
            pull(gene_id) |>
            unique() |>
            paste(collapse = ", ")
          if (length(hit_gene) == 0 || hit_gene == "") {
            hit_gene <- "unknown"
          }

          lines <- c(
            lines,
            glue(
              "\n    Off-target: {hit$sseqid} (gene: {hit_gene}){flag}\n"
            )
          )

          aln <- format_alignment(
            hit$qseqid,
            hit$sseqid,
            hit$qseq,
            hit$sseq,
            hit$qstart,
            hit$qend,
            hit$sstart,
            hit$send,
            hit$pident,
            hit$mch
          )
          lines <- c(lines, glue("{aln}\n"))
        }
      } else {
        lines <- c(lines, "    No off-target hits.\n")
      }
    }
  }

  paste(lines, collapse = "\n")
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

  .blast_results |>
    mutate(sseqid = str_replace(sseqid, "\\s.*", "")) |>
    mutate(
      mch = calculate_mch_v(qseq, sseq),
      on_target = sseqid %in% target_tids,
      off_target = !on_target,
      pair_ID = as.integer(str_replace(qseqid, "\\.\\d+$", "")),
      half = as.integer(str_replace(qseqid, "^\\d+\\.", ""))
    )
}

filter_probes <- function(.blast_annotated, .split_probes, .mch_threshold) {
  .split_probes |>
    select(pair_ID, odd_ID, even_ID) |>
    mutate(
      half1_on_target = map_lgl(odd_ID, \(id) {
        any(.blast_annotated$qseqid == id & .blast_annotated$on_target)
      }),
      half2_on_target = map_lgl(even_ID, \(id) {
        any(.blast_annotated$qseqid == id & .blast_annotated$on_target)
      }),
      both_on_target = half1_on_target & half2_on_target,
      half1_flagged = map_lgl(odd_ID, \(id) {
        any(
          .blast_annotated$qseqid == id &
            .blast_annotated$off_target &
            .blast_annotated$mch >= .mch_threshold
        )
      }),
      half2_flagged = map_lgl(even_ID, \(id) {
        any(
          .blast_annotated$qseqid == id &
            .blast_annotated$off_target &
            .blast_annotated$mch >= .mch_threshold
        )
      }),
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
#' @param .hairpin Character. HCR hairpin channel (e.g. "B1", "B2")
#' @param .transcript_id Character or NULL. Manual transcript override
#' @param .expression Data frame or NULL. Must have transcript_id + tpm columns
#' @param .blast_db Character or NULL. Path to pre-built BLAST db
#' @param .evalue Numeric. BLAST e-value threshold (default 100, permissive)
#' @param .mch_threshold Integer. MCH threshold for off-target flagging
#' @param .outdir Character. Output directory for results
#'
#' @return Invisible list with all intermediate results
design_probes <- function(
  .gene_id,
  .species,
  .hairpin = "B1",
  .transcript_id = NULL,
  .expression = NULL,
  .blast_db = NULL,
  .evalue = 100,
  .mch_threshold = 14,
  .outdir = here("hcr-probes/results")
) {
  config <- get_species_config(.species)
  gene_dir <- file.path(.outdir, sanitize_gene_id(.gene_id))
  fs::dir_create(gene_dir)

  mapping <- read_tsv(config$mapping, show_col_types = FALSE)

  # step 1: select transcript -------------------------------------------------
  message(glue("Step 1: Selecting transcript for {.gene_id}..."))

  tx <- select_transcript(
    .gene_id = .gene_id,
    .species = .species,
    .transcript_id = .transcript_id,
    .expression = .expression
  )

  write_tsv(
    tx$selection_table,
    file.path(gene_dir, "01_transcript-selection.tsv")
  )
  writeXStringSet(
    tx$sequence,
    file.path(gene_dir, "02_selected-transcript.fasta")
  )

  message(glue(
    "  Selected transcript: {tx$transcript_id} (length: {width(tx$sequence)} nt)"
  ))

  # step 2: design probes -----------------------------------------------------
  message("Step 2: Designing probe candidates...")

  split_probes <- design_probe_candidates(tx$sequence, .hairpin)
  message(glue("  Generated {nrow(split_probes)} probe pairs"))

  write_tsv(split_probes, file.path(gene_dir, "03_probe-candidates.tsv"))

  probe_fasta_data <- probes_to_fasta(split_probes)
  probe_seqs <- DNAStringSet(probe_fasta_data$seq)
  names(probe_seqs) <- probe_fasta_data$id
  probe_fasta_path <- file.path(gene_dir, "03_probe-candidates.fasta")
  writeXStringSet(probe_seqs, probe_fasta_path)

  # step 3: BLAST -------------------------------------------------------------
  message("Step 3: Running BLAST specificity check...")

  if (is.null(.blast_db)) {
    .blast_db <- ensure_blast_db(config)
  }

  blast_results <- run_blast(probe_fasta_path, .blast_db, .evalue = .evalue)
  write_tsv(blast_results, file.path(gene_dir, "04_blast-results.tsv"))
  message(glue("  BLAST returned {nrow(blast_results)} hits"))

  # step 4: annotate and generate alignment report -----------------------------
  message("Step 4: Calculating MCH and generating alignment report...")

  blast_annotated <- annotate_blast_results(
    blast_results,
    .gene_id,
    .species,
    mapping
  )

  report <- generate_alignment_report(
    blast_annotated,
    split_probes,
    .gene_id,
    mapping,
    .mch_threshold
  )
  writeLines(report, file.path(gene_dir, "04_blast-alignments.txt"))

  # step 5: filter and output --------------------------------------------------
  message("Step 5: Filtering probes...")

  pair_summary <- filter_probes(blast_annotated, split_probes, .mch_threshold)

  probes_annotated <- split_probes |>
    left_join(
      select(
        pair_summary,
        pair_ID,
        half1_on_target,
        half2_on_target,
        both_on_target,
        half1_flagged,
        half2_flagged,
        any_flagged,
        passed
      ),
      by = "pair_ID"
    )
  write_tsv(probes_annotated, file.path(gene_dir, "05_probes-annotated.tsv"))

  passing_probes <- split_probes |>
    filter(pair_ID %in% pair_summary$pair_ID[pair_summary$passed])

  if (nrow(passing_probes) > 0) {
    opools <- make_oPools_df(passing_probes)
    write_tsv(opools, file.path(gene_dir, "06_probes-passing.tsv"))
  } else {
    write_tsv(
      tibble(Pool_name = character(), Sequence = character()),
      file.path(gene_dir, "06_probes-passing.tsv")
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
