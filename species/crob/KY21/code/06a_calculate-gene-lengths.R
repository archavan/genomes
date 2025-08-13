# From the gtf file:
# 1. Calculate transcript lengths
# 2. Calculate median length of all transcripts by gene_id

# this is needed for running {goseq}.

# packages ====================================================================
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(here)

genomedir <- here("species/crob/genome/KY21")
                  
# gtf =========================================================================
read_gxf_as_tibble <- function(.file) {
  rtracklayer::import(.file) %>% 
    GenomicRanges::as.data.frame() %>% 
    as_tibble()
} 

gtf <- read_gxf_as_tibble(here(genomedir, "results/gtf/ky21_fixed_minimal.gtf.gz"))

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


