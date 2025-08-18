# Arun Chavan
# Started: 2025-08-15

# background ==================================================================

# Get gene names for uniprot ids:

# The protein names in the PROST output are the names of the proteins on
# Uniprot, not the genes. For example, the protein PERT_PIG is thryoid
# peroxidase, for which protein symbol on uniprot is PERT, but the gene name is
# TPO. We need to get gene names for all best hits.

# setup =======================================================================
library(tidyverse)
library(here)
library(jsonlite)

genomedir <- here("species/scla/refseq_v1")
resdir <- here(genomedir, "results/prost/processed")

# data ========================================================================
best_hits_gene_level <- read_tsv(here(resdir, "best-hits.tsv"))

# get gene names ==============================================================
read_uniprot_page <- function(.uniprot_id) {
  url <- sprintf("https://rest.uniprot.org/uniprotkb/%s.json", .uniprot_id)
  tryCatch(
    read_json(url),
    error = function(e) {
      message(sprintf(
        "Failed to read UniProt ID %s: %s",
        .uniprot_id,
        e$message
      ))
      NULL
    }
  )
}

insistent_read_uniprot_page <- purrr::insistently(
  read_uniprot_page,
  rate = rate_delay(pause = 5, max_times = 5)
)


get_gene_name <- function(.uniprot_id) {
  pg <- insistent_read_uniprot_page(.uniprot_id)

  if (is.null(pg)) {
    return(NA_character_)
  }

  if (is.null(pg$genes) || length(pg$genes) == 0) {
    return(NA_character_)
  }

  # extract the primary gene name
  primary <- pg$genes[[1]]$geneName$value
  if (is.null(primary)) {
    return(NA_character_)
  }

  return(primary)
}

best_hits_with_gene_names <- best_hits_gene_level |>
  mutate(gene_name = map_chr(uniprot_id, get_gene_name, .progress = TRUE)) |>
  relocate(gene_name, .after = protein_name)


# write files =================================================================
write_tsv(
  best_hits_with_gene_names,
  here(resdir, "best-hits_with-gene-names.tsv")
)

# end =========================================================================
