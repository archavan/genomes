# Arun Chavan
# Started: 2024-11-07

# background ==================================================================

# Get gene names for uniprot ids:

# The protein names in the PROST output are the names of the proteins on
# Uniprot, not the genes. For example, the protein PERT_PIG is thryoid
# peroxidase, for which protein symbol on uniprot is PERT, but the gene name is
# TPO. We need to get gene names for all best hits.

# setup =======================================================================
library(tidyverse)
library(here)
library(httr)
library(parallel)
library(parallelly)

genomedir <- here("species/crob/KY21")
resdir <- here(genomedir, "results/prost/processed")

# data ========================================================================
best_hits_gene_level <- read_tsv(here(resdir, "best-hits_gene-level.tsv"))

# get gene names ==============================================================

get_gene_name <- function(uniprot_id) {
  url <- paste0("https://www.uniprot.org/uniprot/", uniprot_id, ".txt")
  response <- GET(url)

  if (response$status_code == 200) {
    x <- content(response, as = "text", encoding = "UTF-8")
    lines <- strsplit(x, "\n")[[1]]
    gene_name_line <- grep("^GN\\s+Name", lines, value = TRUE)
    if (length(gene_name_line) > 0) {
      gene_name <- gsub("(^GN\\s+Name=)(.+?)[ ;].*", "\\2", gene_name_line)
      return(gene_name)
    } else {
      return(NA_character_)
    }
  } else {
    return(NA_character_)
  }
}

best_hits_with_gene_names <- best_hits_gene_level %>%
  mutate(
    gene_name = mclapply(
      uniprot_id,
      FUN = get_gene_name,
      mc.cores = availableCores() - 1
    )
  ) %>%
  mutate(gene_name = map_chr(gene_name, `[[`, 1)) %>%
  relocate(gene_name, .after = protein_name)

# write files =================================================================
write_tsv(
  best_hits_with_gene_names,
  here(resdir, "best-hits_gene-level_with-gene-names.tsv")
)

# end =========================================================================
