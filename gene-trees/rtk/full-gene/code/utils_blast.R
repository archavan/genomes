# Arun Chavan

#' Submit a blastp query to the NCBI BLAST API and return parsed JSON results
#'
#' @param .query Amino acid sequence string to use as the query.
#' @param .hitlist_size Maximum number of hits to return.
#' @param .expect E-value threshold for reporting hits.
#' @param .entrez_query Entrez search term to restrict the database, e.g.
#'   `"Chordata[organism]"`.
#' @param .poll_interval Seconds between status checks. NCBI asks for no more
#'   than one poll per minute per RID.
#' @param .blast_url Base URL for the NCBI BLAST CGI endpoint.
#' @param .word_size Minimum word size for seeding alignments. Default 3
#'   (blastp default); increase to 6 for faster, less sensitive searches.
#' @param .blast_email Email address sent with each request, as required by
#'   NCBI's usage policy.
#'
#' @return A list parsed from the BLAST JSON2_S response.
run_blast <- function(
  .query,
  .hitlist_size = 100,
  .expect = 0.001,
  .entrez_query,
  .word_size = 3,
  .poll_interval = 60,
  .blast_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
  .blast_email = "arun.chavan@yale.edu"
) {
  # stage 1: submit
  resp <- request(.blast_url) |>
    req_body_form(
      CMD = "Put",
      PROGRAM = "blastp",
      DATABASE = "refseq_protein",
      QUERY = .query,
      ENTREZ_QUERY = .entrez_query,
      HITLIST_SIZE = as.character(.hitlist_size),
      EXPECT = as.character(.expect),
      WORD_SIZE = as.character(.word_size),
      EMAIL = .blast_email
    ) |>
    req_timeout(60) |>
    req_options(http_version = 2) |>
    req_perform()

  body <- resp_body_string(resp)
  rid <- str_extract(body, "(?<=RID = )\\S+")
  stopifnot(!is.na(rid))
  rtoe <- str_extract(body, "(?<=RTOE = )\\d+") |> as.integer()

  message(glue("submitted: RID = {rid}, estimated wait = {rtoe}s"))

  # stage 2: poll
  start_time <- Sys.time()
  Sys.sleep(max(rtoe, .poll_interval))

  repeat {
    status_resp <- request(.blast_url) |>
      req_url_query(CMD = "Get", RID = rid, FORMAT_OBJECT = "SearchInfo") |>
      req_timeout(60) |>
      req_options(http_version = 2) |>
      req_perform()

    status <- resp_body_string(status_resp) |> str_extract("(?<=Status=)\\w+")
    elapsed <- round(as.numeric(difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )))
    message(glue("  [{elapsed}s elapsed] status: {status}"))

    if (is.na(status)) {
      stop(glue("BLAST job {rid}: unreadable status response"))
    }
    if (status == "READY") {
      break
    }
    if (status %in% c("FAILED", "UNKNOWN")) {
      stop(glue("BLAST job {rid}: {status}"))
    }

    Sys.sleep(.poll_interval)
  }

  # stage 3: retrieve
  message("retrieving results...")
  result_resp <- request(.blast_url) |>
    req_url_query(
      CMD = "Get",
      RID = rid,
      FORMAT_TYPE = "JSON2_S",
      HITLIST_SIZE = as.character(.hitlist_size)
    ) |>
    req_timeout(120) |>
    req_options(http_version = 2) |>
    req_perform()

  result <- resp_body_json(result_resp)
  message(glue("done ({elapsed}s total)"))
  result
}

#' Download the reference assembly protein sequences for a taxon via NCBI datasets
#'
#' Looks up the reference genome assembly for the taxon, downloads its protein
#' FASTA via the datasets CLI, and writes it to .out_fasta.
#'
#' @param .taxid NCBI taxonomy ID.
#' @param .out_fasta Path to output FASTA file.
download_refseq_proteome <- function(.taxid, .out_fasta) {
  summary_out <- system2(
    path.expand("~/Dropbox/research/apps/datasets"),
    args = c(
      "summary",
      "genome",
      "taxon",
      .taxid,
      "--reference",
      "--as-json-lines"
    ),
    stdout = TRUE
  )
  gcf <- jsonlite::fromJSON(summary_out[1])$accession
  stopifnot(length(gcf) == 1)
  message(glue("using assembly {gcf} for txid{.taxid}"))

  zip_file <- tempfile(fileext = ".zip")
  on.exit(unlink(zip_file))

  exit_code <- system2(
    path.expand("~/Dropbox/research/apps/datasets"),
    args = c(
      "download",
      "genome",
      "accession",
      gcf,
      "--include",
      "protein",
      "--filename",
      zip_file
    )
  )
  if (exit_code != 0) {
    stop(glue("datasets download exited with code {exit_code}"))
  }

  tmp_dir <- tempfile()
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  unzip(zip_file, exdir = tmp_dir)

  faa <- list.files(
    tmp_dir,
    pattern = "protein\\.faa$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(faa) != 1) {
    all_files <- list.files(tmp_dir, recursive = TRUE)
    stop(glue(
      "{length(faa)} *_protein.faa found in datasets zip for txid{.taxid}. ",
      "Files present: {paste(all_files, collapse = ', ')}"
    ))
  }
  file.copy(faa, .out_fasta)
  invisible(.out_fasta)
}

#' Run local blastp for a forward search, returning columns matching parse_hits_to_table()
#'
#' @param .query Amino acid sequence string.
#' @param .db Path to a BLAST+ database (without file extension).
#' @param .hitlist_size Maximum number of hits to return.
#' @param .word_size Word size for seeding alignments.
#' @param .evalue E-value threshold.
#' @param .nthreads Number of CPU threads.
#'
#' @return A tibble with columns: accession, title, bitscore, evalue, length,
#'   query_from, query_to, hit_from, hit_to, align_len, gaps, identity.
run_blast_local_forward <- function(
  .query,
  .db,
  .hitlist_size = 100,
  .word_size = 3,
  .evalue = 0.001,
  .nthreads = 4
) {
  query_file <- tempfile(fileext = ".fasta")
  out_file <- tempfile(fileext = ".tsv")
  on.exit({
    unlink(query_file)
    unlink(out_file)
  })

  writeLines(paste0(">query\n", .query), query_file)

  exit_code <- system2(
    path.expand("~/Dropbox/research/apps/ncbi-blast-2.17.0+/bin/blastp"),
    args = c(
      "-query", query_file,
      "-db", .db,
      "-out", out_file,
      "-outfmt", shQuote("6 sacc stitle bitscore evalue slen qstart qend sstart send length gaps nident"),
      "-max_target_seqs", as.character(.hitlist_size),
      "-word_size", as.character(.word_size),
      "-evalue", as.character(.evalue),
      "-num_threads", as.character(.nthreads)
    )
  )
  if (exit_code != 0) {
    stop(glue("blastp exited with code {exit_code}"))
  }

  empty <- tibble(
    accession  = character(),
    title      = character(),
    bitscore   = numeric(),
    evalue     = numeric(),
    length     = integer(),
    query_from = integer(),
    query_to   = integer(),
    hit_from   = integer(),
    hit_to     = integer(),
    align_len  = integer(),
    gaps       = integer(),
    identity   = integer()
  )

  if (!file.exists(out_file) || file.size(out_file) == 0) {
    return(empty)
  }

  read_tsv(
    out_file,
    col_names = c(
      "accession", "title", "bitscore", "evalue", "length",
      "query_from", "query_to", "hit_from", "hit_to",
      "align_len", "gaps", "identity"
    ),
    col_types = "ccddiiiiiiii",
    show_col_types = FALSE
  )
}

#' Run local blastp for a reciprocal best-hit check
#'
#' @param .query Amino acid sequence string.
#' @param .db Path to a BLAST+ database (without file extension).
#' @param .hitlist_size Maximum number of hits to return.
#' @param .word_size Word size for seeding alignments.
#' @param .evalue E-value threshold.
#' @param .nthreads Number of CPU threads.
#'
#' @return A one-row tibble with rbh_accession, rbh_title, rbh_bitscore, rbh_evalue.
run_blast_local_rbh <- function(
  .query,
  .db,
  .hitlist_size = 5,
  .word_size = 6,
  .evalue = 0.001,
  .nthreads = 4
) {
  query_file <- tempfile(fileext = ".fasta")
  out_file <- tempfile(fileext = ".tsv")
  on.exit({
    unlink(query_file)
    unlink(out_file)
  })

  writeLines(paste0(">query\n", .query), query_file)

  exit_code <- system2(
    path.expand("~/Dropbox/research/apps/ncbi-blast-2.17.0+/bin/blastp"),
    args = c(
      "-query", query_file,
      "-db", .db,
      "-out", out_file,
      "-outfmt", shQuote("6 sacc stitle bitscore evalue"),
      "-max_target_seqs", as.character(.hitlist_size),
      "-word_size", as.character(.word_size),
      "-evalue", as.character(.evalue),
      "-num_threads", as.character(.nthreads)
    )
  )
  if (exit_code != 0) {
    stop(glue("blastp exited with code {exit_code}"))
  }

  if (!file.exists(out_file) || file.size(out_file) == 0) {
    return(tibble(
      rbh_accession = NA_character_,
      rbh_title     = NA_character_,
      rbh_bitscore  = NA_real_,
      rbh_evalue    = NA_real_
    ))
  }

  read_tsv(
    out_file,
    col_names = c("rbh_accession", "rbh_title", "rbh_bitscore", "rbh_evalue"),
    col_types = "ccdd",
    show_col_types = FALSE
  ) |>
    slice(1)
}
