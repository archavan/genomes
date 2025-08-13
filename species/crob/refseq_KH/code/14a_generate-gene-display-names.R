# Arun Chavan
# Started: 2022-08-15

# background ==================================================================
# We have gene annotations from 4 sources: 
#
# 1. RefSeq
# 2. KEGG 
# 3. emapper
# 4. PROST 
#
# Here we try to merge them into a single table. Generally with the following 
# order of preference but see code for details:
# 
# 1. refseq
# 2. kegg
# 3. emapper
# 4. prost

# set up ======================================================================
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/crob/genome/refseq_KH/")

# data ========================================================================
gene_ids <- read_tsv(here(genomedir, "results/genes/gene-ids.tsv"))
ncbi <- read_tsv(here(genomedir, "results/genes/ncbi-genes-table.tsv"))
kegg <- read_tsv(here(genomedir, "results/kegg/kegg-genes.tsv"))
emap.ko <- read_tsv(here(genomedir, "results/emapper/processed/emapper_kegg_ko.tsv"))
emap.name <- read_tsv(here(genomedir, "results/emapper/processed/emapper_predicted-protein-name.tsv"))
prost <- read_tsv(here(genomedir, "results/prost/processed/best-hits_gene-level_with-gene-names.tsv"))

###############################################################################
# symbols =====================================================================
###############################################################################

symbols <- gene_ids %>% 
  select(gene_id, ncbi_gene_id) 

## add ncbi symbols -----------------------------------------------------------
ncbi.symbol <- ncbi %>% 
  select(ncbi_gene_id, symbol) %>% 
  filter(!is.na(symbol)) %>% 
  filter(str_detect(symbol, "^LOC", negate = TRUE)) %>% 
  distinct() %>% 
  dplyr::rename(symbol_ncbi = symbol)

symbols <- left_join(symbols, ncbi.symbol, by = "ncbi_gene_id")

## add emapper predicted protein name -----------------------------------------

# 26 ncbi gene ids have multiple predicted names by emapper. For these genes, I
# will manually pick the best one that makes the most sense given the ncbi
# description or kegg name for those genes

# workflow:
# 1. If one of the names matches the refseq description, use it. If not,
# 2. If one of the names matches the kegg name, use it. if not,
# 3. look up gene id on www.orthodb.org; if one of the genes matches the name 
#    of the orthology group, use it. If not, 
# 4. free web search and use whatever information available to make in informed 
#    decision. If non found,
# 5. use the first name.

emap.name.unique <- emap.name %>% 
  group_by(gene_id) %>% 
  summarise(names = list(predicted_protein_name)) %>% 
  mutate(n_names = map_int(names, length)) %>% 
  filter(n_names == 1) %>% 
  unnest_longer(col = names) %>% 
  select(gene_id, symbol_emapper = names, -n_names)

emap.name.multi <- emap.name %>% 
  group_by(gene_id) %>% 
  summarise(names = list(predicted_protein_name)) %>% 
  mutate(n_names = map_int(names, length)) %>% 
  arrange(-n_names) %>% 
  mutate(pasted = map_chr(names, ~ paste0(.x, collapse = ","))) %>% 
  filter(n_names > 1) %>% 
  left_join(select(gene_ids, ncbi_gene_id, gene_id), by = "gene_id") %>% 
  left_join(select(kegg, ncbi_gene_id, refseq_description, ko_name), by = "ncbi_gene_id") %>% 
  unique()

emap.name.multi <- emap.name.multi %>% 
  mutate(picked_name = case_when(
    refseq_description == "syntenin-1 isoform X3" ~ "SDCBP",
    refseq_description == "guanine nucleotide-binding protein alpha-4 subunit-like isoform X2" ~ "GNAI3",
    refseq_description == "glycerol-3-phosphate phosphatase-like isoform X2" ~ "PGP",
    refseq_description == "small integral membrane protein 8 isoform X2" ~ "Smim8",
    refseq_description == "major facilitator superfamily domain-containing protein 3 isoform X1" ~ "MFSD3",
    refseq_description == "zinc finger matrin-type protein 2 isoform X2" ~ "ZMAT2",
    refseq_description == "uncharacterized protein LOC100177386" ~ "nvd",
    refseq_description == "trehalase isoform X1" ~ "TREH",
    refseq_description == "uncharacterized protein LOC100178097 isoform X1" ~ "tnxb",
    refseq_description == "uncharacterized protein LOC100179778 isoform X1" ~ "SNAPC1",
    refseq_description == "protein odr-4 homolog isoform X1" ~ "ODR4",
    refseq_description == "uncharacterized protein LOC100181469 isoform X2" ~ "MTG2",
    refseq_description == "B9 domain-containing protein 1-like isoform X3" ~ "B9D1",
    refseq_description == "uncharacterized protein LOC100181989 isoform X2" ~ "prr14",
    refseq_description == "protein-glucosylgalactosylhydroxylysine glucosidase-like isoform X4" ~ "ATHL1",
    refseq_description == "uncharacterized protein LOC100184499 isoform X1" ~ "STAC2",
    refseq_description == "uncharacterized protein C18orf19 homolog A isoform X1" ~ "FAM210A",
    refseq_description == "UDP-N-acetylglucosamine transporter isoform X3" ~ "SLC35A3",
    refseq_description == "probable E3 ubiquitin-protein ligase RNF144A isoform X2" ~ "RNF144A",
    refseq_description == "N-acylethanolamine-hydrolyzing acid amidase-like isoform X1" ~ "NAAA",
    refseq_description == "integrin alpha-4-like isoform X1" ~ "ITGAX",
    refseq_description == "cyclin-K-like isoform X1" ~ "CCNK",
    refseq_description == "uncharacterized protein LOC104266258 isoform X1" ~ "PTPRQ",
    refseq_description == "mediator of RNA polymerase II transcription subunit 24-like isoform X1" ~ "MED24",
    refseq_description == "solute carrier family 22 member 6-A-like isoform X1" ~ "SLC22A6",
    refseq_description == "ci-zf(c3h/u1like)-1; zinc finger protein Ci-ZF(C3H/U1like)-1" ~ "ci-zf(c3h/u1like)-1"
  ))

emap.name.multi <- emap.name.multi %>% 
  select(gene_id, symbol_emapper = picked_name) 

symbols <- left_join(symbols, 
                     bind_rows(emap.name.unique,
                               emap.name.multi),
                     by = "gene_id")

## add symbols extracted from emapper ko name ----------------------------------
# Some gene ids have multiple KO ids associated with them. I don't think it's
# strictly true with KO ids, but for the sake of making a decision, let's assume
# that the larger KO ids are more specific and pick the largest one out of
# multiple KO ids.
emap.ko.name <- emap.ko %>% 
  group_by(gene_id) %>% 
  arrange(desc(kegg_ko_id), .by_group = TRUE) %>% 
  slice_head(n = 1) %>% 
  ungroup()

# only deal with gene ids for which we haven't found a symbol yet.
emap.ko.name <- emap.ko.name %>% 
  filter(gene_id %in% (symbols %>% 
                         filter(is.na(symbol_ncbi) & is.na(symbol_emapper)) %>% 
                         pull(gene_id)))

# separate out gene ids that have multiple symbols listed and unique ones
emap.ko.name.unique <- emap.ko.name %>% 
  separate(kegg_ko_name, into = c("symbol", "description"), ";\\s+") %>% 
  mutate(symbol = map(symbol, ~ str_split(.x, ",\\s+")[[1]])) %>% 
  mutate(n_names = map_int(symbol, length)) %>% 
  filter(n_names == 1) %>% 
  unnest_longer(col = symbol) %>% 
  select(gene_id, symbol_emap.ko = symbol) %>% 
  distinct()

emap.ko.name.multi <- emap.ko.name %>% 
  separate(kegg_ko_name, into = c("symbol", "description"), ";\\s+") %>% 
  mutate(symbol = map(symbol, ~ str_split(.x, ",\\s+")[[1]])) %>% 
  mutate(n_names = map_int(symbol, length)) %>% 
  filter(n_names > 1)

# for the multi ones, we will pick the first symbol. For most rows that I
# spot-checked, the first one seems to match the description.
emap.ko.name.multi <- emap.ko.name.multi %>% 
  mutate(symbol_emap.ko = map_chr(symbol, `[[`, 1)) %>% 
  select(gene_id, symbol_emap.ko)

# add to symbols
symbols <- left_join(symbols,
                     bind_rows(emap.ko.name.unique,
                               emap.ko.name.multi),
                     by = "gene_id")


## add symbols from PROST -----------------------------------------------------
prost.symbols <- prost %>% 
  select(gene_id, symbol_prost = gene_name) %>% 
  filter(!is.na(symbol_prost)) %>% 
  distinct() 

# add to symbols
symbols <- left_join(symbols, prost.symbols, by = "gene_id")

## add source prefixes --------------------------------------------------------
add_source_prefix <- function(x, prefix) {
  ifelse(is.na(x), NA_character_, paste0(prefix, x))
}

symbols <- symbols %>% 
  mutate(symbol_ncbi = add_source_prefix(symbol_ncbi, "[n] "),
         symbol_emapper = add_source_prefix(symbol_emapper, "[e] "),
         symbol_emap.ko = add_source_prefix(symbol_emap.ko, "[k] "),
         symbol_prost = add_source_prefix(symbol_prost, "[p] "))

## put it all togehter --------------------------------------------------------
symbols <- symbols %>% 
  mutate(symbol = coalesce(symbol_ncbi, symbol_emapper, symbol_emap.ko, symbol_prost)) %>% 
  select(gene_id, symbol) 

## add a source column --------------------------------------------------------
symbols <- symbols %>% 
  mutate(source_symbol = case_when(
    str_detect(symbol, "^\\[n\\]") ~ "ncbi",
    str_detect(symbol, "^\\[e\\]") ~ "emapper",
    str_detect(symbol, "^\\[k\\]") ~ "kegg",
    str_detect(symbol, "^\\[p\\]") ~ "prost",
    .default = NA_character_
  ))

## summary ---------------------------------------------------------------------
nrow(symbols)
#> [1] 16297

sum(is.na(symbols$symbol))
#> 3245

# We have symbols for 13052 gene ids now and NAs for 3245 genes. 

###############################################################################
# descriptions ================================================================
###############################################################################

# make one table with descriptions from the NCBI and KEGG
# 1. NCBI has gene descriptions
# 2. KEGG:
#       - KEGG KO name is a description, but
#       - KEGG also provides a column called refseq description.
#           - For the most part, these are identical to the NCBI descriptions
#           - Some are different. Some contain suffixes like "isoform X1" etc, 
#             indicating that these were pulled from protein descriptions. We 
#             can get rid of those suffixes to make them essentially the same as 
#             NCBI description. 
#           - These are all prefixed with the gene symbol, which we have to str_remove
#       - There several transcription factors for which NCBI description and 
#         refseq description only says "transcription factor", but KEGG name is 
#         actually more informative, e.g, Insulin promoting factor-1, and so on. 
#         In these cases, we want to prioritize using KO name instead of NCBI 
#         description even if it exists. 
# 3. PROST

descriptions <- gene_ids %>% 
  select(gene_id, ncbi_gene_id) %>% 
  left_join((ncbi %>% 
               select(ncbi_gene_id, description_refseq = description) %>% 
               distinct()),
            by = "ncbi_gene_id") %>% 
  left_join((kegg %>% 
               group_by(ncbi_gene_id) %>% 
               arrange(desc(ko_id), .by_group = TRUE) %>% 
               slice_head(n = 1) %>% 
               ungroup() %>% 
               select(ncbi_gene_id, 
                      description_ko.refseq = refseq_description, 
                      description_ko = ko_name) %>% 
               mutate(description_ko.refseq = str_remove(description_ko.refseq, "^.+?;\\s+"))),
            by = "ncbi_gene_id") %>% 
  left_join((prost %>% 
               select(gene_id, description_prost = protein_type) %>% 
               distinct()),
            by = "gene_id")
  
# remove extraneous prefixes and suffixes from the refseq description column from kegg
descriptions <- descriptions %>% 
  mutate(description_ko.refseq = str_remove(description_ko.refseq, "\\s+isoform.+$")) %>% 
  mutate(description_ko.refseq = str_remove(description_ko.refseq, "^LOW QUALITY PROTEIN:\\s+")) %>% 
  mutate(description_ko.refseq = str_replace(description_ko.refseq, "^uncharacterized protein", "uncharacterized")) %>% 
  mutate(description_ko.refseq = str_remove(description_ko.refseq, "\\s+precursor$"))

# The descriptions from refseq and the refseq descriptions from KEGG are
# identical except in some prefixes etc. We will firth make a unified
# description column from these two columns that we will use as the refseq
# description

descriptions <- descriptions %>% 
  mutate(description_refseq = ifelse(str_detect(description_refseq, "uncharacterized"),
                                     NA_character_, 
                                     description_refseq)) %>% 
  mutate(description_ko.refseq = ifelse(str_detect(description_ko.refseq, "uncharacterized"),
                                        NA_character_, 
                                        description_ko.refseq))

descriptions <- descriptions %>% 
  mutate(description_refseq_combined = coalesce(description_refseq, 
                                                description_ko.refseq)) %>% 
  relocate(description_refseq_combined, .after = description_ko.refseq) %>% 
  select(gene_id, description_refseq_combined, description_ko, description_prost)

## put them all together ------------------------------------------------------
descriptions <- descriptions %>% 
  mutate(description_refseq_combined = add_source_prefix(description_refseq_combined, "[n] "),
         description_ko = add_source_prefix(description_ko, "[k] "),
         description_prost = add_source_prefix(description_prost, "[p] ")) %>% 
  mutate(description = coalesce(description_refseq_combined, description_ko, description_prost))


## better transcription factor annotations ------------------------------------

# There are some genes that have refseq description as "transcription factor
# protein", but their KO description is more specific. We will use the KO
# description for these genes.

descriptions <- descriptions %>% 
  mutate(
    description = case_when(
      description == "[n] transcription factor protein" & 
        !is.na(description_ko) ~ description_ko,
      description == "[n] transcription factor protein" &
        is.na(description_ko) &
        !is.na(description_prost) ~ description_prost,
      .default = description
    )
  )

# we still have 2 genes with "transcription factor protein" description and NA
# ko_name as well as NA prost name. We will manually look them up on orthodb and
# update their description.

descriptions <- descriptions %>% 
  mutate(description = case_when(
    gene_id == "ncbi:778715" ~ "[m] Zinc finger C2H2-type", # m for manual
    gene_id == "ncbi:778929" ~ "[n] transcription factor protein",
    .default = description
  )) 


sum(is.na(descriptions$description))
#> [1] 2421

## replace NAs with uncharacterized -------------------------------------------
descriptions <- descriptions %>% 
  mutate(description = ifelse(is.na(description),
                              "uncharacterized",
                              description))

## add a source column --------------------------------------------------------
descriptions <- descriptions %>% 
  mutate(source_description = case_when(
    str_detect(description, "^\\[n\\]") ~ "ncbi",
    str_detect(description, "^\\[k\\]") ~ "kegg",
    str_detect(description, "^\\[p\\]") ~ "prost",
    str_detect(description, "^\\[m\\]") ~ "manual",
    .default = NA_character_
  )) 


## finalize descriptions ------------------------------------------------------
sum(descriptions$description == "uncharacterized")
# now we have a description for all genes, with 2421 genes labelled as
# "uncharacterized". Now we get rid of all columns except the final description
# and gene id

descriptions <- descriptions %>% 
  select(gene_id, description, source_description)


###############################################################################
# final merged table ==========================================================
###############################################################################

merged <- left_join(symbols, descriptions, by = "gene_id")

merged <- merged %>% 
  mutate(display_name = paste0(
    "[", gene_id, "]", "; ",
    ifelse(is.na(symbol), "", paste0(symbol, "; ")),
    ifelse(is.na(description), "", paste0(description))
  )) %>% 
  select(gene_id, symbol, description, display_name) 

###############################################################################
# write files =================================================================
###############################################################################
write_tsv(merged, here(genomedir, "results/genes/gene-display-names.tsv"))
write_tsv(symbols, here(genomedir, "results/genes/gene-symbols.tsv"))
write_tsv(descriptions, here(genomedir, "results/genes/gene-descriptions.tsv"))

# end =========================================================================
