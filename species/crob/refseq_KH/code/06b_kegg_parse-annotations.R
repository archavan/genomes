# Arun Chavan
# Started: 2022-08-14  

# Background ==================================================================

# We have already downloaded kegg annotations for all genes using REST API and
# saved the results to an .rds file. The KEGG API returns the results for our
# queries in its flat file format described here:
# https://www.kegg.jp/kegg/docs/dbentry.html.

# Here we read in the .rds file and parse the KEGGREST::keggGet() output to
# extract annotations such as Kegg Orthology group, NCBI gene id, etc.

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/crob/genome/refseq_KH/")

# data ========================================================================
genes <- read_rds(here(genomedir, "results/kegg/kegg_genes.rds"))
logfile <- read_lines(here(genomedir, "results/kegg/06a_kegg_fetch-annotations.out"))

# There are 14141 ciona genes in the KEGG genome database. Ensure that all genes
# were queried by our script.
all.equal(nrow(genes), 14141) 


# Ensure that KEGGREST did not fail to get info for any genes. If it did, run
# fetch annotations script again but only for the failed genes (it takes several
# hours for it to run for all genes)
all.equal(sum(grepl("failed", logfile)), 0)


# kegg orthology group id =====================================================
genes$ko_id <- map_chr(genes$get, 
                       ~ ifelse(is.null(.x[["ORTHOLOGY"]]),
                                NA,
                                attr(.x[["ORTHOLOGY"]], "names")))

# kegg orthology group name ===================================================
genes$ko_name <- map_chr(genes$get, 
                         ~ ifelse(is.null(.x[["ORTHOLOGY"]]),
                                  NA,
                                  .x[["ORTHOLOGY"]]))

# db ids ======================================================================
extract_db_id <- function(getres, db) {
  ifelse(is.null(getres$DBLINKS),
         NA,
         getres$DBLINKS %>% 
           .[str_detect(., db)] %>% 
           str_remove(glue("{db}: ")))
}

genes <- genes %>% 
  mutate(ncbi_gene_id = map_chr(get, extract_db_id, "NCBI-GeneID"),
         ensembl_gene_id = map_chr(get, extract_db_id, "Ensembl"),
         uniprot_id = map_chr(get, extract_db_id, "UniProt")) 

# kegg pathways ===============================================================
kegg_pathways <- genes %>% 
  hoist(get, "PATHWAY") %>% 
  select(kegg_gene_id, PATHWAY) %>% 
  unnest_longer(col = PATHWAY) %>% 
  dplyr::rename(kegg_pathway_id = PATHWAY_id, kegg_pathway_name = PATHWAY) %>% 
  relocate(kegg_pathway_id, .after = kegg_gene_id) %>% 
  filter(!is.na(kegg_pathway_id) | !is.na(kegg_pathway_name))

# write =======================================================================
genes %>% 
  dplyr::select(ncbi_gene_id, 
                kegg_gene_id, 
                ko_id, 
                ensembl_gene_id,
                uniprot_id,
                refseq_description, 
                ko_name) %>% 
  write_tsv(here(genomedir, "results/kegg/kegg-genes.tsv"))

kegg_pathways %>% 
  write_tsv(here(genomedir, "results/kegg/kegg-pathways.tsv"))

# end =========================================================================

