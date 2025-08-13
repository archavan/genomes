# Arun Chavan
# Started: 2022-06-22

# background ==================================================================
# Download gene info from ncbi genes and clean up the file for further use. 

# setup =======================================================================
# packages 
library(tidyverse)
library(here)

genomedir <- here("species/crob/genome/refseq_KH")

# download data ===============================================================
# I downloaded gene info table from NCBI Genes database from the following link manually:
# https://www.ncbi.nlm.nih.gov/data-hub/gene/table/taxon/7719/
# select columns: all columns
# on date: 2022-06-22

# read data ===================================================================
genes <- read_tsv(here(genomedir, "data/ncbi-genes/gene-info_raw.tsv"))

# clean up ====================================================================
## clean up column names
clean_names <- function(x) {
  x <- str_replace_all(x, "\\s", "_")
  x <- tolower(x)
  x <- str_replace_all(x, "geneid", "gene_id")
  x <- str_replace_all(x, "ids", "id")
  return(x)
}

names(genes) <- clean_names(names(genes))

## get rid of columns for which all values are NA
genes <- genes[, apply(genes, 2, function(x) !all(is.na(x)))]

## since these are all genes for ciona, we don't really need the columns that
## have information abot which species these genes come from
genes <- genes %>% 
  select(-taxonomic_name,
         -common_name,
         -taxonomic_id) 

# write cleaned up file =======================================================
write_tsv(genes, here(genomedir, "results/genes/ncbi-genes-table.tsv"))

# end =========================================================================
