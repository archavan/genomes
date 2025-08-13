# Arun Chavan
# started: 2024-11-05

# background ==================================================================

# Process eggnog mapper output to extract various annotations by for all
# gene_ids, e.g. GO terms, KEGG ids, predicted protein names etc. The emapper
# output is at the protein id level, so we have to summarize it to the gene
# level.

# info about output files for this version of emapper (2.1.7)
# https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.8#output-format

# setup =======================================================================
library(GO.db)
library(KEGGREST)
library(rtracklayer)
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/crob/genome/KY21/")
resdir <- here(genomedir, "results/emapper/processed")
fs::dir_create(resdir)

read_gxf_as_tibble <- function(.file) {
  rtracklayer::import(.file) %>% 
    GenomicRanges::as.data.frame() %>% 
    as_tibble()
}

# genes data ==================================================================
gid.pid <- read_tsv(here(genomedir, "results/genes/gid-pid_mapping.tsv"))
genes.gtf <- read_gxf_as_tibble(here(genomedir, "results/gtf/ky21_fixed.gtf.gz"))

# emapper output ==============================================================
emap <- read_tsv(here(genomedir, "results/emapper/out/out.emapper.annotations"),
                 col_names = FALSE,
                 comment = "#")

# from: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.8#annotations-file
emap_headers <- c(
  # serch hit fields
  "query_name",
  "seed_eggNOG_ortholog",
  "seed_ortholog_evalue",
  "seed_ortholog_score",
  # orthologous group fields
  "eggNOG_OGs",
  "max_annot_lvl",
  "COG_cat",
  "Description",
  # transferred annotation fields
  "Preferred_name",
  "GOs",
  "EC",
  "KEGG_ko",
  "KEGG_Pathway",
  "KEGG_Module",
  "KEGG_Reaction",
  "KEGG_rclass",
  "BRITE",
  "KEGG_TC",
  "CAZy",
  "BiGG_Reaction",
  "PFAMs"
)

colnames(emap) <- emap_headers

# the emapper output is at the protein id level -- we need to add gene ids to
# the emapper output so that we summarize at the level of genes
emap <- gid.pid %>% 
  select(gene_id, protein_id) %>% 
  left_join(emap, by = c("protein_id" = "query_name"))



# functions for processing the output =========================================

#' Extract specified column from emapper output for all genes
#' If the column contains multiple values for the same genes (GO terms, e.g.)
#'  split the values, and unnest the df on them to produce a tidy format output
#'
#' @param df emapper output that has been joined with gene_ids
#' @param eggnog_col column to extract (quoted)
#' @param rename.to rename the column to (quoted)
#'
#' @return data frame in tidy format with two columns, gene_id and values from the extracted column

unnest_to_gene_level <- function(df = emap, eggnog_col, rename.to = NULL) {
  
  eggnog_col <- sym(eggnog_col)

  # GO IDs/KEGG ids are present as a character string for each protein id. Here, we
  # combine all GO/KEGG IDs associated with all protein products of a gene, str_split
  # them, and find unique ones. Then, we unnest the df on the list column with GO/KEGG
  # IDs so that we end up with a df with two columns in tidy format: gene_id and
  # go/KEGG_id
  
  df <- df %>%
    dplyr::select(gene_id, {{ eggnog_col }}) %>%
    group_by(gene_id) %>%
    summarize({{ eggnog_col }} := list({{ eggnog_col }})) %>%
    ungroup() %>%
    mutate(
      {{ eggnog_col }} := {{ eggnog_col }} %>%
        map_chr(str_flatten, collapse = ",") %>%
        str_split(",") %>%
        map(unique)
    ) %>%
    dplyr::filter(!is.na({{ eggnog_col }})) %>%
    unnest(col = {{ eggnog_col }})
  
  # emapper output contains dashes when it doesn't have any information. There
  # are NAs too, but it looks like NAs are used when no seed was found in the
  # eggnog database and therefore there is not information at all on the query.
  # For cases where only some information is missing, you have dashes instead.
  # We want to get rid of these rows too because they are essentially NAs.
  
  df <- df %>%
    filter(!{{ eggnog_col }} == "-")
  
  
  if (!is.null(rename.to)) {
    rename.to  <- sym(rename.to)
    df <- dplyr::rename(df, {{ rename.to }} := {{ eggnog_col }})
  }
  
  return(df)
  
}

#' Fetch name for a given kegg_id
#' (Actually, keggList() works better for what we need to do here, so we don't need this function)
#' @param kegg_id 
#'
#' @return named character string
 
fetch_kegg_name <- function(kegg_id) { 
  # we have to use tryCatch because REST errors out if it doesn't find a kegg_id
  out <- tryCatch( 
    keggGet(dbentries = kegg_id)[[1]]$NAME,
    error = function(cond) {
      message(glue("{kegg_id} failed: {cond}"))
      return(NA)
    }
  )    
  return(out)
}


#' Fetch ids and names for all entries from a given kegg database
#'
#' @param database one of KEGGREST::listDatabases()
#' @param prefix ids are returned with a prefix e.g. "path:" for pathway ids, 
#'   which needs to be removed to match with emap output
#'
#' @return tibble columns kegg id and kegg name

fetch_kegg_list <- function(database, prefix = NULL) {
  
  df <- keggList(database) %>% 
    enframe() # kegglist returns of named list
  
  if (!is.null(prefix)) {
    df <- mutate(df, name = str_remove(name, prefix))
  }  
  
  names(df)[names(df) == "name"] <- glue("kegg_{database}_id")
  names(df)[names(df) == "value"] <- glue("kegg_{database}_name")
  
  return(df)
  
}

# GO ids and terms ============================================================
emap.go <- unnest_to_gene_level(df = emap, eggnog_col = "GOs", rename.to = "go_id")

## Now we add GO term to this df. We get the GO terms associated with GO ids using {GO.db}. 
### First create a dataframe with GOID, TERM, and ONTOLOGY mapping
godb <- AnnotationDbi::select(GO.db, 
                              keys = unique(emap.go$go_id), 
                              columns = c("ONTOLOGY", "TERM"), 
                              keytype = "GOID")
godb <- godb %>% 
  dplyr::filter(!is.na(TERM)) %>% # emapper also returns some obsolete GO IDs for which GO.db returns NA for term
  dplyr::mutate(TERM = glue("[{ONTOLOGY}] {TERM}")) %>% 
  dplyr::rename(go_term = TERM) %>% 
  dplyr::select(-ONTOLOGY)

  
### join with emap.go
emap.go <- left_join(emap.go, godb, by = c("go_id" = "GOID")) %>% 
  dplyr::filter(!is.na(go_term)) # get rid of rows with obsolete go_ids

write_tsv(emap.go, here(resdir, "emapper_go.tsv"))

# KEGG database lists =========================================================
kegg.ko      <- fetch_kegg_list("ko")
kegg.pathway <- fetch_kegg_list("pathway", "path:")
kegg.module  <- fetch_kegg_list("module", "md:")
kegg.brite   <- fetch_kegg_list("brite", "br:")

# Almost all ids from these databases that are output by emapper are found in
# the lists we fetched from KEGG REST API. For kegg.pathway, however, we only
# fetch patway ids that start with "map", and about 388 of the kegg pathway ids
# output by empapper are not not found in the list we fetched. These turn out to
# be pathway ids that start with "ko:". For pathways, the KEGGREST::keggGet()
# method ends up finding all pathway ids. However, the ones that start with
# "ko:" are essentially the same as those that start with "map:", and as far as
# I could figure out, the only difference between them is that map: ids are
# manually curated (https://www.genome.jp/kegg/kegg3a.html).  So, for pathways,
# we will still use KEGGREST::keggList() method and only retain mpa: ids and
# ignore ko: ids.

# KEGG pathway ================================================================
emap.kegg.pathway <- unnest_to_gene_level(df = emap, 
                                          eggnog_col = "KEGG_Pathway", 
                                          rename.to = "kegg_pathway_id")

emap.kegg.pathway <- left_join(emap.kegg.pathway, 
                               kegg.pathway, 
                               by = "kegg_pathway_id") %>% 
  dplyr::filter(!is.na(kegg_pathway_name))

write_tsv(emap.kegg.pathway, 
          here(resdir, "emapper_kegg_pathway.tsv"))

# KEGG KO =====================================================================
emap.kegg.ko <- unnest_to_gene_level(df = emap, 
                                     eggnog_col = "KEGG_ko", 
                                     rename.to = "kegg_ko_id")
emap.kegg.ko <- emap.kegg.ko %>% 
  mutate(kegg_ko_id = str_remove(kegg_ko_id, "ko:"))

emap.kegg.ko <- left_join(emap.kegg.ko,
                          kegg.ko,
                          by = "kegg_ko_id") %>% 
  dplyr::filter(!is.na(kegg_ko_name))

write_tsv(emap.kegg.ko, here(resdir, "emapper_kegg_ko.tsv"))

# KEGG module =================================================================
emap.kegg.module <- unnest_to_gene_level(df = emap, 
                                         eggnog_col = "KEGG_Module", 
                                         rename.to = "kegg_module_id")

emap.kegg.module <- left_join(emap.kegg.module, 
                              kegg.module,
                              by = "kegg_module_id") %>% 
  dplyr::filter(!is.na(kegg_module_name))

write_tsv(emap.kegg.module, here(resdir, "emapper_kegg_module.tsv"))

# BRITE =======================================================================
emap.kegg.brite <- unnest_to_gene_level(df = emap, 
                                        eggnog_col = "BRITE", 
                                        rename.to = "kegg_brite_id")

emap.kegg.brite <- left_join(emap.kegg.brite, 
                             kegg.brite,
                             by = "kegg_brite_id") %>% 
  dplyr::filter(!is.na(kegg_brite_name))

write_tsv(emap.kegg.brite, here(resdir, "emapper_kegg_brite.tsv"))

# Predicted protein name ======================================================
emap.protein.name <- unnest_to_gene_level(df = emap,
                                          eggnog_col = "Preferred_name",
                                          rename.to = "predicted_protein_name")

write_tsv(emap.protein.name, 
          here(resdir, "emapper_predicted-protein-name.tsv"))

# PFAM domains ================================================================
emap.pfam <- unnest_to_gene_level(df = emap, 
                                  eggnog_col = "PFAMs",
                                  rename.to = "pfam_domain")

write_tsv(emap.pfam, here(resdir, "emapper_pfam-domains.tsv"))

# description =================================================================
# 69 genes have more than one descriptions associated with them because
# different proteins asociated with those genes have different descriptions. We
# will write them all to file -- they are in the tidy format.
emap.descriptions <- emap %>% 
  select(gene_id, Description) %>% 
  rename(description = Description) %>% 
  distinct() %>% 
  filter(!is.na(description)) %>% 
  filter(!description == "-")

write_tsv(emap.descriptions, here(resdir, "emapper_description.tsv"))

# enzyme EC code ==============================================================
emap.ec <- unnest_to_gene_level(df = emap, 
                                eggnog_col = "EC", 
                                rename.to = "enzyme_ec")

write_tsv(emap.ec, here(resdir, "emapper_enzyme-ec.tsv"))


###############################################################################
# extract genesets ############################################################
###############################################################################

# TF genes ====================================================================
# We will define TF genes as those belonging to either of the following:
# 1) GO:0003700 (DNA-binding transcription factor activity) https://www.ebi.ac.uk/QuickGO/term/GO:0003700 
# 2) BRITE ko03000 (Transcription factors): https://www.genome.jp/kegg-bin/get_htext?ko03000.keg

tf.go <- dplyr::filter(emap.go, go_id == "GO:0003700") %>% 
  select(gene_id) %>% 
  distinct() # 539 genes

tf.brite <- dplyr::filter(emap.kegg.brite, kegg_brite_id == "ko03000") %>% 
  select(gene_id) %>% 
  distinct() # 365 genes

# STATS:
# intersect = 288 genes
# union = 616 genes.

# Take union of the unique genes and extract those rows from emapper output
tf <- emap %>% 
  dplyr::filter((gene_id %in% tf.go$gene_id) | (gene_id %in% tf.brite$gene_id)) %>% 
  dplyr::select(gene_id, Preferred_name) %>% 
  dplyr::rename(predicted_protein_name = Preferred_name) %>% 
  dplyr::distinct() %>% 
  mutate(predicted_protein_name = ifelse(predicted_protein_name == "-", 
                                         NA_character_, 
                                         predicted_protein_name))

write_tsv(tf, here(resdir, "emapper_geneset_tf.tsv"))

# end =========================================================================

