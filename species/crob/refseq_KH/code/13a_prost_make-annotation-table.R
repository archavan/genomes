# Arun Chavan
# Started: 2023-12-19

# background ==================================================================

# Here, we generate a table of gene id, protein id, symbol, and description. The
# purpose of this is to see which genes are still uncharacterized after emapper
# etc, so that Kejue can run them through his PLM pipeline.

# setup =======================================================================
library(tidyverse)
library(here)
library(seqinr)

genomedir <- here("species/crob/refseq_KH/")
resdir <- here(genomedir, "results/prost")
fs::dir_create(resdir)

# data ========================================================================
genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
pids <- read_tsv(here(genomedir, "results/genes/gid-pid_mapping.tsv"))
proteome <- read.fasta(
  here(genomedir, "data/refseq/GCF_000224145.3_KH_protein.faa.gz"),
  seqtype = "AA"
)

# create table ================================================================
genes <- select(genes, -display_name)

# add protein ids
genes <- genes %>%
  left_join(select(pids, -symbol), by = "gene_id") %>%
  relocate(protein_id)

# filter out ncRNA genes
genes <- genes %>%
  filter(!is.na(protein_id))

uncharacterized_coding <- genes %>%
  filter(str_detect(description, "uncharacterized"))

# make sure that all proteins in this list are in the fasta file ==============

stopifnot(sum(!genes$protein_id %in% names(proteome)) == 0)

# write tables ================================================================

write_tsv(genes, here(resdir, "genes_coding_all.tsv"))
write_tsv(
  uncharacterized_coding,
  here(resdir, "genes_coding_uncharacterized.tsv")
)

# end =========================================================================
