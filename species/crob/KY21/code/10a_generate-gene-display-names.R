# Arun Chavan
# Started: 2024-11-05

# background ==================================================================

# Generate gene display names for plotting combining symbols and descriptions
# from KY21, emapper, and PROST.

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/crob/genome/KY21")

###############################################################################
# data ========================================================================
###############################################################################

## ghost ======================================================================
ghost <- read_tsv(here(genomedir, "results/ghost/gene-info_ghostdb_parsed.tsv"))

ghost <- ghost %>% 
  rename(ghost_symbol = symbol,
         ghost_description = description)

## emapper ====================================================================
# with emapper, we kept all symbols and descriptions associated with each gene
# id. We need to pick one, so we will just keep the first one.

retain_first_entry <- function(x) {
  x %>% 
    group_by(gene_id) %>% 
    slice_head(n = 1) %>% 
    ungroup()
}

emap_symbol <- read_tsv(here(genomedir, "results/emapper/processed/emapper_predicted-protein-name.tsv"))
emap_desc <- read_tsv(here(genomedir, "results/emapper/processed/emapper_description.tsv"))

emap_symbol <- retain_first_entry(emap_symbol)
emap_desc <- retain_first_entry(emap_desc)

emap <- full_join(emap_symbol, emap_desc, by = "gene_id") %>% 
  rename(emap_symbol = predicted_protein_name,
         emap_description = description)

## prost ======================================================================
prost <- read_tsv(here(genomedir, "results/prost/processed/best-hits_gene-level_with-gene-names.tsv"))

prost <- prost %>% 
  select(gene_id, 
         prost_symbol = gene_name,
         prost_description = protein_type)

###############################################################################
# generate gene display names =================================================
###############################################################################

genes <- reduce(list(ghost, emap, prost), left_join, by = "gene_id")

# add source information to columns ===========================================
add_source_prefix <- function(x, prefix) {
  ifelse(!is.na(x), glue("{prefix}{x}"), NA)
}

genes <- genes %>% 
  mutate(ghost_symbol = add_source_prefix(ghost_symbol, "[g] "),
         ghost_description = add_source_prefix(ghost_description, "[g] "),
         hsap_homolog_symbol = add_source_prefix(hsap_homolog_symbol, "[h] "),
         hsap_homolog_description = add_source_prefix(hsap_homolog_description, "[h] "),
         emap_symbol = add_source_prefix(emap_symbol, "[e] "),
         emap_description = add_source_prefix(emap_description, "[e] "),
         prost_symbol = add_source_prefix(prost_symbol, "[p] "),
         prost_description = add_source_prefix(prost_description, "[p] "))


# coalesce ====================================================================
# we use the order of columns in coalesce() to set preference: prefer
# ghost_symbol if exists, if not use hsap_homolog_symbol, then emapper, and prost
display_names <- genes %>% 
  mutate(symbol = coalesce(ghost_symbol, 
                           hsap_homolog_symbol,
                           emap_symbol,
                           prost_symbol),
         description = coalesce(ghost_description, 
                                hsap_homolog_description,
                                emap_description,
                                prost_description)) %>% 
  select(gene_id, symbol, description)

# do all genes have some symbol or description?
display_names %>% 
  filter(is.na(symbol) | is.na(description))

display_names <- display_names %>% 
  mutate(description = ifelse(is.na(description),
                              "uncharacterized",
                              description)) 

# merge symbol and description into a single name =============================
display_names <- display_names %>% 
  mutate(display_name = paste0(
    "[", gene_id, "]", "; ",
    ifelse(is.na(symbol), "", paste0(symbol, "; ")),
    ifelse(is.na(description), "", paste0(description))
  )) 

# write files =================================================================
write_tsv(display_names, here(genomedir, "results/genes/gene-display-names.tsv"))

# end =========================================================================
