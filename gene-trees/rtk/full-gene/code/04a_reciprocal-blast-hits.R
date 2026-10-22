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
  dplyr::filter,
  dplyr::select,
  dplyr::rename,
  dplyr::slice,
  base::setdiff
)

basedir <- here("gene-trees/rtk/full-gene")
orthdir <- here("orthofinder/2026-08-07")

source(here(basedir, "code/utils_blast.R"), local = TRUE)

###############################################################################
# data ========================================================================
###############################################################################

hits_all <- read_tsv(here(basedir, "results/blast/hits_all.tsv"))
spp <- read_tsv(here(orthdir, "data/species.tsv")) |>
  select(-ensembl_proteome_fasta)
seeds_meta <- read_tsv(here(basedir, "data/seeds/seeds.tsv"))

###############################################################################
# deduplicate =================================================================
###############################################################################
hits_all_dedup <- hits_all |>
  slice_max(
    order_by = bitscore,
    n = 1,
    by = c(seed, abb, accession),
    with_ties = FALSE
  )

###############################################################################
# limit to top hits ===========================================================
###############################################################################

hits_top <- hits_all_dedup |>
  slice_max(order_by = bitscore, n = 25, by = c(seed, abb)) |>
  filter(bitscore > 100)

###############################################################################
# fetch sequences =============================================================
###############################################################################

hits_seqs <- spp |>
  mutate(
    faa = map(
      abb,
      ~ readAAStringSet(here(
        orthdir,
        "results/input-proteomes",
        glue("{.x}.fa")
      ))
    )
  )

hits_top <- hits_top |>
  mutate(
    sequence = map2_chr(
      abb,
      accession,
      function(x, y) {
        prot <- hits_seqs$faa[[which(hits_seqs$abb == x)]]
        as.character(prot[[y]])
      },
      .progress = TRUE
    )
  )

###############################################################################
# setup for reciprocal blast ==================================================
###############################################################################

rbh <- hits_top |>
  select(
    seed,
    hit_accession = accession,
    hit_sequence = sequence
  ) |>
  separate(
    seed,
    into = c("seed_abb", "seed_symbol"),
    sep = "__",
    remove = FALSE
  )

###############################################################################
# setup local blast databases for reciprocal blast ============================
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
# reciprocal blast ============================================================
###############################################################################

rbh_cache_file <- here(basedir, "results/blast/rbh_cache.tsv")

if (!file.exists(rbh_cache_file)) {
  write_tsv(
    tibble(
      hit_accession = character(),
      seed_abb = integer(),
      rbh_accession = character(),
      rbh_title = character(),
      rbh_bitscore = numeric(),
      rbh_evalue = numeric()
    ),
    rbh_cache_file
  )
}

rbh_cache <- read_tsv(rbh_cache_file)

rbh_pending <- rbh |>
  distinct(seed_abb, hit_accession, hit_sequence) |>
  anti_join(rbh_cache, by = c("hit_accession", "seed_abb"))

for (i in seq_len(nrow(rbh_pending))) {
  if (i == 1 || i %% 5 == 0) {
    message(glue("processing {i} of {nrow(rbh_pending)}"))
  }
  db_path <- here(db_dir, as.character(rbh_pending$seed_abb[i]))
  top_hit <- run_blast_local_rbh(
    .query = rbh_pending$hit_sequence[i],
    .db = db_path
  )
  write_tsv(
    bind_cols(
      tibble(
        hit_accession = rbh_pending$hit_accession[i],
        seed_abb = rbh_pending$seed_abb[i]
      ),
      top_hit
    ),
    rbh_cache_file,
    append = TRUE
  )
}

rbh_cache <- read_tsv(rbh_cache_file)

###############################################################################
# confirm reciprocal best hits ================================================
###############################################################################

family_gids <- as.character(seeds_meta$gene_id)

hits_final <- rbh |>
  left_join(
    select(
      hits_top,
      seed,
      hit_accession = accession,
      hit_abb = abb,
      fwd_bitscore = bitscore
    ),
    by = c("seed", "hit_accession")
  ) |>
  left_join(rbh_cache, by = c("hit_accession", "seed_abb")) |>
  mutate(
    is_rbh = !is.na(rbh_accession) &
      map_lgl(rbh_accession, function(x) x %in% family_gids)
  )

rbh_confirmed <- hits_final |>
  filter(is_rbh) |>
  filter(hit_accession != rbh_accession) |>
  mutate(
    seeds_all = str_c(seed, collapse = ";"),
    n_seeds = n(),
    .by = c(hit_abb, hit_accession)
  ) |>
  arrange(desc(fwd_bitscore)) |>
  slice_head(n = 1, by = c(hit_abb, hit_accession))


write_tsv(rbh_confirmed, here(basedir, "results/blast/hits_rbh-confirmed.tsv"))

# end #########################################################################
