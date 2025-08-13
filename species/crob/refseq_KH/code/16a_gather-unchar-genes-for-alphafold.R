# Arun Chavan
# Started: 2024-04-03

# background ==================================================================

# There are many genes that are still uncharacterized, even after emapper and
# prost. Kejue will run alphafold on them.

# setup =======================================================================
library(seqinr)
library(tidyverse)
library(here)

genomedir <- here("species/crob/genome/refseq_KH/")
resdir <- here(genomedir, "results/alphafold")
fs::dir_create(here(resdir, "input"))

# data ========================================================================
genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
plen <- read_tsv(here(genomedir, "results/genes/protein-lengths.tsv"))
pfa <- read.fasta(here(genomedir, 
                       "data/refseq/GCF_000224145.3_KH_protein.faa.gz"), 
                  seqtype = "AA")

# get sequences for uncharacterized genes =====================================

# there are 1,688 uncharacterized non-coding RNA genes too but with alphafolf we
# can only use the coding ones

unchar_coding_genes_longest_protein <- genes %>% 
  filter(description == "uncharacterized") %>%
  left_join(plen) %>% 
  filter(!is.na(protein_id)) %>%  # non-coding RNAs
  group_by(gene_id) %>% 
  slice_max(length, n = 1, with_ties = FALSE) %>% 
  ungroup()

# there are 733 uncharacterized coding genes.

pfa_unchar <- pfa[unchar_coding_genes_longest_protein$protein_id]

# write files =================================================================
write_tsv(unchar_coding_genes_longest_protein, 
          here(resdir, "input/unchar_coding-genes_longest-protein.tsv"))

write.fasta(pfa_unchar, 
            names(pfa_unchar),
            here(resdir, "input/unchar_coding-genes_longest-protein.fa"))

# end =========================================================================