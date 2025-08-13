# From the gtf file:
# 1. Calculate transcript lengths
# 2. Calculate median length of all transcripts by gene_id

# this is needed for running {goseq}.

# packages ====================================================================
library(tidyverse)
library(here)

genomedir <- here("species/bflo/genome/VNyyK/")

# gtf =========================================================================
gtf <- read_tsv(here(genomedir, "results/gtf/gtf_fixed_parsed.tsv"), 
                col_types = cols(
                  .default = col_character(),
                  start = col_double(),
                  end = col_double(),
                  width = col_double(),
                  exon_number = col_double()
                ))

# calculate median gene length ================================================
glen <- gtf %>% 
  select(gene_id, transcript_id, type, strand, start, end) %>% 
  filter(type == "exon") %>% 
  mutate(length = end - start + 1) %>% 
  group_by(gene_id, transcript_id) %>% 
  summarize(length = sum(length)) %>% 
  ungroup() %>% 
  group_by(gene_id) %>% 
  summarize(length = round(median(length))) %>% 
  ungroup()

write_tsv(glen, here(genomedir, "results/genes/gene-length-median.tsv"))

# end =========================================================================
