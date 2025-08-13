# parse gtf file to create a dataframe with gene information.
# started: 2022-06-20

# packages ====================================================================
library(rtracklayer)
library(tidyverse)
library(here)
library(bioseq)

# setup =======================================================================
genomedir <- here("species/crob/genome/KH13/")

# parse gtf ===================================================================
gtf <- rtracklayer::import(
  here(genomedir, "data/ghost-db/KH.KHGene.2013_agat-converted.gtf")
) %>% 
  GenomicRanges::as.data.frame() %>% 
  as_tibble()

# write to file
write_tsv(gtf, here(genomedir, "results/gtf/gtf_parsed.tsv"))

# gene-id and protein-id mapping ==============================================
# It's hard to tell from the names of the fields in the gtf attributes column
# which one is really, but the one labelled `target` is the one that has ids
# that look like protein ids from the headers of the fasta file.

pid <- gtf %>% 
  select(gene_id, Target) %>% 
  filter(!is.na(Target)) %>% 
  separate(Target, into = c("protein_id", NA, NA), sep = "\\s+") %>% 
  distinct() %>% 
  rename(ghost_gene_id = gene_id, ghost_protein_id = protein_id)


# 1. There are still about 4000 proteins that are in the fasta but not in the gtf
#    file. 
# 2. don't know why and I don't know what to do with them. For now I am
#    ignoring them until I know better.
# 3. It appeara that most of the protein ids in the fasta that are absent in the
#    gtf have "ND" in their 5th field, which stands for "Not Determined". These
#    are supposed to be proteins that have been computationally predicted. That
#    is perhaps why they are missing from the gtf. 

pfasta <- read_fasta(here(genomedir, "data/ghost-db/KH.KHGene.2013.Longest.protein.fasta"), 
                     type = "AA")
sum(names(pfasta) %in% pid$ghost_protein_id)

# write gid pid to files
fs::dir_create(here(genomedir, "results/genes"))
write_tsv(pid, here(genomedir, "results/genes/gid-pid-mapping.tsv"))

### end =======================================================================
