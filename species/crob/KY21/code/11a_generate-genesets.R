# Arun Chavan
# Started: 2024-11-08

# background ==================================================================

# Generate lists of genes belonging to a specific geneset.

# setup =======================================================================
library(tidyverse)
library(here)

genomedir <- here("species/crob/genome/KY21")
resdir <- here(genomedir, "results/genes")

# data ========================================================================
genes <- read_tsv(here(resdir, "gene-display-names.tsv"))

go <- read_tsv(here(resdir, "go_coalesced.tsv.gz"))
keggbrite <- read_tsv(here(genomedir, "results/emapper/processed/emapper_kegg_brite.tsv"))

# TF genes ====================================================================
# We will define TF genes as those belonging to either of the following:
# 1) GO:0003700 (DNA-binding transcription factor activity) https://www.ebi.ac.uk/QuickGO/term/GO:0003700 
# 2) BRITE ko03000 (Transcription factors): https://www.genome.jp/kegg-bin/get_htext?ko03000.keg

tf_go <- filter(go, go_id == "GO:0003700") %>% 
  select(gene_id) %>% 
  distinct() # 668 genes

tf_brite <- filter(keggbrite, kegg_brite_id == "ko03000") %>% 
  select(gene_id) %>% 
  distinct() # 355 genes

# STATS:
# intersect =  genes
# union = 723 genes

# Take union of the unique genes and extract those rows from emapper output
tf <- genes %>% 
  filter((gene_id %in% tf_go$gene_id) | (gene_id %in% tf_brite$gene_id)) %>% 
  select(-display_name)

write_tsv(tf, here(resdir, "geneset_tf.tsv"))

# Neuropeptides ===============================================================
# GO:0160041 - [MF] neuropeptide activity
# GO:0005184 -[MF] neuropeptide hormone activity

np_go <- go %>% 
  filter(go_id == "GO:0005184") %>% 
  select(gene_id) %>%
  distinct() 

np <- genes %>%
  filter(gene_id %in% np_go$gene_id) %>%
  select(-display_name)

write_tsv(np, here(resdir, "geneset_neuropeptides.tsv"))

# All hormones ================================================================
# GO:0005179 [MF] hormone activity

hormones_go <- go %>% 
  filter(go_id == "GO:0005179") %>% 
  select(gene_id) %>%
  distinct()

hormones <- genes %>%
  filter(gene_id %in% hormones_go$gene_id) %>%
  select(-display_name)

write_tsv(hormones, here(resdir, "geneset_hormones.tsv"))

# GPCRs =======================================================================
# GO:0004930 - [MF] G protein-coupled receptor activity

gpcr_go <- go %>% 
  filter(go_id == "GO:0004930") %>% 
  select(gene_id) %>%
  distinct()

gpcr <- genes %>%
  filter(gene_id %in% gpcr_go$gene_id) %>%
  select(-display_name)

write_tsv(gpcr, here(resdir, "geneset_gpcr.tsv"))

# end =========================================================================