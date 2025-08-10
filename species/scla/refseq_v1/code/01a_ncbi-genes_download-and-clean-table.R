# Arun Chavan
# Started: 2025-08-07

# background ==================================================================
# Download gene info from ncbi genes and clean up the file for further use.

# setup =======================================================================
library(tidyverse)
library(here)
library(janitor)

genomedir <- here("species/scla/refseq_v1")
resdir <- here(genomedir, "results/genes")
fs::dir_create(resdir)

# download data ===============================================================
# I downloaded gene info table from NCBI Genes database from the following link manually:
# https://www.ncbi.nlm.nih.gov/datasets/gene/taxon/7725/
# select columns: all columns
# on date: 2025-08-07

# read data ===================================================================
genes <- read_tsv(
  here(genomedir, "data/ncbi-genes/gene-info_raw.tsv"),
  col_types = list(
    `NCBI GeneID` = col_character(),
    `Annotation Genomic Range Accession` = col_character(),
    `Annotation Genomic Range Start` = col_integer(),
    `Annotation Genomic Range Stop` = col_integer()
  )
) %>%
  clean_names()

# clean up ====================================================================
## get rid of columns for which all values are NA
genes <- genes[, apply(genes, 2, function(x) !all(is.na(x)))]

## get rid of columns we don't need
genes <- select(genes, -taxonomic_name)

genes <- genes %>%
  rename(
    seqname = annotation_genomic_range_accession,
    begin = annotation_genomic_range_start,
    end = annotation_genomic_range_stop
  )

genes <- genes %>%
  arrange(seqname, ncbi_gene_id)

## add ncbi prefix to gene ids
genes <- genes %>%
  mutate(gene_id = paste0("ncbi.", ncbi_gene_id)) %>%
  relocate(gene_id) %>%
  select(-ncbi_gene_id)

# write cleaned up file =======================================================
write_tsv(genes, here(genomedir, "results/genes/ncbi-genes-table.tsv"))

# end =========================================================================
