# Arun Chavan
# Started: 2026-10-22

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(Biostrings)

basedir <- here("gene-trees/rtk/full-gene")
resdir <- here(basedir, "results/blast/hits")
orthdir <- here("orthofinder/2026-08-07")

source(here(basedir, "code/utils_blast.R"), local = TRUE)

###############################################################################
# data ========================================================================
###############################################################################

spp <- read_tsv(here(orthdir, "data/species.tsv")) |>
  select(-ensembl_proteome_fasta)

seeds <- fs::dir_ls(here(basedir, "data/seeds"), glob = "*.fasta") |>
  map(~ str_remove(.x, here(basedir, "data/seeds/"))) |>
  (function(x) set_names(x, str_remove(x, ".fasta")))() |>
  map(~ readAAStringSet(here(basedir, "data/seeds", .x))) |>
  map(~ as.character(.x[[1]]))

###############################################################################
# set up local BLAST databases ================================================
###############################################################################

db_dir <- path.expand("~/Downloads/blastdbs")
fs::dir_create(db_dir)

for (i in unique(spp$abb)) {
  db_path <- file.path(db_dir, as.character(i))
  if (
    fs::file_exists(glue("{db_path}.psq")) ||
      fs::file_exists(glue("{db_path}.pdb"))
  ) {
    message(glue("BLAST db for txid{i} already exists, skipping"))
    rm(i, db_path)
    next
  }
  fasta_file <- here(orthdir, "results/input-proteomes", glue("{i}.fa"))
  stopifnot(fs::file_exists(fasta_file))
  message(glue("indexing BLAST database for {i}..."))
  exit_code <- system2(
    path.expand("~/Dropbox/research/apps/ncbi-blast-2.17.0+/bin/makeblastdb"),
    args = c("-in", fasta_file, "-dbtype", "prot", "-out", db_path)
  )
  stopifnot(exit_code == 0)
  rm(i, db_path, fasta_file, exit_code)
}

###############################################################################
# run blast ===================================================================
###############################################################################

res_df <- expand_grid(seed = names(seeds), abb = spp$abb) |>
  left_join(spp, by = "abb") |>
  mutate(resfile = glue("{seed}/{abb}.rds"), hits = vector("list", n()))

for (i in names(seeds)) {
  fs::dir_create(here(resdir, i))
  rm(i)
}

for (i in seq_len(nrow(res_df))) {
  rds_path <- here(resdir, res_df$resfile[i])
  if (fs::file_exists(rds_path)) {
    message(glue("loading cached result: {res_df$resfile[i]}"))
    res_df$hits[[i]] <- read_rds(rds_path)
    rm(i, rds_path)
    next
  }
  message(glue(
    "running BLAST for seed {res_df$seed[i]} against taxid {res_df$abb[i]}"
  ))
  db_path <- file.path(db_dir, as.character(res_df$abb[i]))
  res_df$hits[[i]] <- run_blast_local_forward(
    .query = seeds[[res_df$seed[i]]],
    .db = db_path
  )
  write_rds(res_df$hits[[i]], rds_path, compress = "gz")
  rm(i, rds_path, db_path)
}

###############################################################################
# parse hits ==================================================================
###############################################################################

hits_all <- res_df |>
  unnest(col = hits) |>
  mutate(
    pct_identity = (identity / align_len) * 100,
    query_len = map_int(seed, ~ seeds[[.x]] |> nchar()),
    query_coverage = 100 * (query_to - query_from + 1) / query_len
  ) |>
  select(-resfile)

write_tsv(hits_all, here(basedir, "results/blast/hits_all.tsv"))

# end #########################################################################
