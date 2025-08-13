# Arun Chavan
# 2023-06-06

# background ==================================================================

# Create gene-id and transcript-id / protein-id mapping table, extract a list of
# mitochondrial genes.

# setup =======================================================================
library(tidyverse)
library(here)

genomedir <- here("species/bflo/genome/VNyyK/")

# read data ===================================================================
gtf <- read_tsv(here(genomedir, "results/gtf/gtf_fixed_parsed.tsv"), 
                col_types = cols(
                  .default = col_character(),
                  start = col_double(),
                  end = col_double(),
                  width = col_double(),
                  exon_number = col_double()
                ))

# gene id and protein id mapping ==============================================
gid_pid <- gtf %>% 
  select(gene_id, symbol, protein_id) %>% 
  filter(!is.na(protein_id)) %>% 
  distinct()

write_tsv(gid_pid, here(genomedir, "results/genes/gid-pid_mapping.tsv"))

# gene id and transcript id mapping ===========================================
gid_tid <- gtf %>% 
  select(gene_id, symbol, transcript_id) %>% 
  filter(!is.na(transcript_id)) %>% 
  distinct()

write_tsv(gid_tid, here(genomedir, "results/genes/gid-tid_mapping.tsv"))

# Extract genesets ============================================================
## Mitochondrial genes --------------------------------------------------------
# Mitochondrial genome sequence name: https://www.ncbi.nlm.nih.gov/nuccore/NC_000834.1
mt_chr <- "NC_000834.1" 

mt_genes <- gtf %>% 
  filter(seqnames == mt_chr) %>% 
  select(seqnames, gene_id, symbol) %>% 
  distinct()

write_tsv(mt_genes, here(genomedir, "results/genes/geneset_mitochondrial.tsv"))

# end =========================================================================


