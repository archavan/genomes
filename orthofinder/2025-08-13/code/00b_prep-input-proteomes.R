# Arun Chavan
# Started: 2025-08-13

# background ==================================================================
# Prep input proteomes for orthofinder: use only the longest protein per gene.

# setup =======================================================================
library(Biostrings)
library(tidyverse)
library(here)
library(glue)

basedir <- here("orthofinder/2025-08-13")

# data ========================================================================
spp <- read_tsv(here(basedir, "data/species.tsv"))

###############################################################################
# prepare proteomes from ensembl ==============================================
###############################################################################

spp_ens <- filter(spp, !is.na(ensembl_proteome_fasta))

spp_ens <- spp_ens |>
  mutate(fa = map(ensembl_proteome_fasta, readAAStringSet))

# get gene_id and protein_id mapping ==========================================
pull_id_mapping <- function(.fa) {
  ann <- tibble(annot = names(.fa), protein_length = width(.fa))

  idmap <- ann |>
    mutate(
      protein_id = str_extract(annot, "^(.+?) ", group = 1),
      gene_id = str_extract(annot, "gene:(.+?) ", group = 1),
      symbol = str_extract(annot, "gene_symbol:(.+?) ", group = 1)
    ) |>
    select(gene_id, symbol, protein_id, protein_length, annot) |>
    arrange(gene_id, desc(protein_length)) |>
    distinct()

  return(idmap)
}

spp_ens <- spp_ens |>
  mutate(ids = map(fa, pull_id_mapping))

fs::dir_create(here(basedir, "results/id-name-mapping"))
walk2(
  spp_ens$ids,
  spp_ens$abb,
  ~ write_tsv(
    select(.x, gene_id, symbol, protein_id, protein_length),
    here(basedir, "results/id-name-mapping", paste0(.y, ".tsv"))
  )
)

# subset to longest proteins ==================================================
spp_ens <- spp_ens |>
  mutate(
    longest_proteins = map(
      ids,
      ~ .x |>
        group_by(gene_id) |>
        slice_max(protein_length, n = 1, with_ties = FALSE) |>
        ungroup() |>
        pull(annot)
    )
  )

subset_fa <- function(.fa, .annot) {
  .fa[.annot]
}

spp_ens <- spp_ens |>
  mutate(fa_longest = map2(fa, longest_proteins, subset_fa))

# use gene_id as header =======================================================

replace_header_with_gene_id <- function(.fa) {
  names(.fa) <- str_extract(names(.fa), "gene:(.+?) ", group = 1)
  return(.fa)
}

spp_ens <- spp_ens |>
  mutate(fa_longest_geneid = map(fa_longest, replace_header_with_gene_id))

# sanity checks ===============================================================
headers_are_unique <- function(.fa) {
  length(unique(names(.fa))) == length(names(.fa))
}

stopifnot(
  all(map_lgl(spp_ens$fa_longest_geneid, headers_are_unique))
)

stopifnot(
  all(map_lgl(
    spp_ens$fa_longest,
    ~ headers_are_unique(replace_header_with_gene_id(.x))
  ))
)

stopifnot(
  !all(map_lgl(
    spp_ens$fa,
    ~ headers_are_unique(replace_header_with_gene_id(.x))
  ))
)

###############################################################################
# non ensembl species =========================================================
###############################################################################

spp_other <- filter(spp, is.na(ensembl_proteome_fasta))

spp_other <- spp_other |>
  mutate(
    version = case_when(
      abb == "crob" ~ "KY21",
      abb == "scla" ~ "refseq_v1",
      abb == "bflo" ~ "VNyyK"
    )
  )

spp_other <- spp_other |>
  mutate(
    fa_longest_geneid = map2(
      abb,
      version,
      ~ readAAStringSet(here(
        "species",
        .x,
        .y,
        "results/genes/longest-proteins.fa"
      ))
    )
  )

###############################################################################
# combine spp_ens and spp_other ===============================================
###############################################################################

spp_all <- bind_rows(
  select(spp_ens, abb, fa_longest_geneid),
  select(spp_other, abb, fa_longest_geneid)
)

fs::dir_create(here(basedir, "results/input-proteomes"))
walk2(
  spp_all$fa_longest_geneid,
  spp_all$abb,
  ~ writeXStringSet(
    .x,
    filepath = here(basedir, "results/input-proteomes", paste0(.y, ".fa")),
    format = "fasta"
  )
)

# end =========================================================================
