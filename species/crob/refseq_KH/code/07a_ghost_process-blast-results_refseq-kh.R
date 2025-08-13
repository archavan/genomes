#' ---
#' title: "process refseq-ghosh blast results"
#' author: "Arun Chavan"
#' output: github_document
#' ---
#'
#' Started: 2022-08-15

# background ==================================================================
# Here we process blast for KH annotations between refseq and ghost-db.
# Shweta Ramdas ran the blast and sent me the output file.

# NCBI : 21,096 proteins
# Ghost: 55,017 proteins

# Notes about the run from Shweta:
# 1. Most of the NCBI sequences seem to have a match in KH, so that's good. But
#    there are 587 that are missing.
# 2. I took one of the four proteins that hadn't been found in BLAST and BLASTed
#    this protein only. It shows an e-value of 1e-04, marginally significant.
#    And its original e value had been much lower. I think this reflect the
#    properties of e-values, and how they depend on what is being queried against
#    what, since it is trying to simulate some kind of null from the database we
#    are using. I think I trust my original results a little more.
# 3. I tried the inverse of what I did--use the NCBI proteins as the query and
#    the KH proteins as the database. This might change the properties of the
#    e-values. So the two ways of doing blast (switching up the query and the
#    database) both gave the same number of significant NCBI protein hits. I have
#    a little more confidence in these results. I would use them, and for the
#    four extra proteins for which we don't have a result, we can use the top
#    result from the blast on the ghost webpage and justify it saying: they have
#    a reasonably high e-value, it just got missed genome-wide.

# setup =======================================================================
# packages
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/crob/refseq_KH/")

# blast output ================================================================
blast <- read_tsv(
  here(genomedir, "results/ghost/blast_ghost-vs-refseq_kh_results_1e-06.out"),
  col_names = c(
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
  )
)
# blastn output format: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

# genomic data ================================================================
gene_ids <- read_tsv(here(genomedir, "results/genes/gene-ids.tsv"))
gid_pid <- read_tsv(here(genomedir, "results/genes/gid-pid_mapping.tsv"))
ghost_pid <- read_tsv(here(
  "species/crob/KH13/results/genes/gid-pid-mapping.tsv"
))

# exploratory analysis of the blast output ====================================
# number of significant refseq hits for each ghost protein id
blast %>%
  group_by(qseqid) %>%
  summarise(sseqid = list(sseqid)) %>%
  ungroup() %>%
  mutate(n_hits = map_int(sseqid, length)) %>%
  select(-sseqid) %>%
  summary()

#> qseqid             n_hits
#> Length:48672       Min.   :   1.00
#> Class :character   1st Qu.:   1.00
#> Mode  :character   Median :   3.00
#>                    Mean   :  39.92
#>                    3rd Qu.:  17.00
#>                    Max.   :4025.00

# filter  =====================================================================
# for each ghost protein, keep only one ncbi protein, one with the highest
# percent identity
top_hit <- blast %>%
  filter(evalue < 10^-6) %>% # already filtered, but just to make sure
  group_by(qseqid) %>%
  arrange(desc(pident), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

# distribution of percent identities for top hits vs all hits
(left_join(
  x = blast,
  y = (top_hit %>%
    mutate(hit = "top") %>%
    select(qseqid, sseqid, hit)),
  by = c("qseqid", "sseqid")
) %>%
  ggplot(aes(pident)) +
  geom_histogram(aes(fill = hit), color = "white", size = 0.25) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())) %>%
  ggsave(
    filename = here(
      genomedir,
      "results/ghost/blast_ghost-vs-refseq_kh_percent-identity-distribution.pdf"
    ),
    plot = .,
    width = 6,
    height = 4
  )

# create gene id mapping between refseq and kh based on protein id mapping ====
rs_gh <- left_join(
  gid_pid,
  select(top_hit, qseqid, sseqid),
  by = c("protein_id" = "sseqid")
) %>%
  rename(ghost_protein_id = qseqid, ncbi_gene_id = gene_id) %>%
  left_join(ghost_pid, by = "ghost_protein_id") %>%
  select(ncbi_gene_id, ghost_gene_id) %>%
  filter(!is.na(ghost_gene_id)) %>%
  distinct()

# write files =================================================================
write_tsv(
  top_hit,
  here(genomedir, "results/ghost/blast_ghost-vs-refseq_kh_top-hits.tsv")
)
write_tsv(
  rs_gh,
  here(genomedir, "results/ghost/refseq-ghost-gene-id-mapping.tsv")
)

# end =========================================================================
