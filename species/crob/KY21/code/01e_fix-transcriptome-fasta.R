# Arun Chavan
# Started: 2024-11-07

# background ==================================================================

# Here we combine the fasta files from KY21 and from Ensembl (mitochondrial genes)

# setup =======================================================================
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/crob/genome/KY21")

read_gxf_as_tibble <- function(.file) {
  rtracklayer::import(.file) %>% 
    GenomicRanges::as.data.frame() %>% 
    as_tibble()
} 

# data ========================================================================
gtf <- read_gxf_as_tibble(here(genomedir, "results/gtf/ky21_fixed.gtf.gz"))

tx_ky21 <- readDNAStringSet(here(genomedir, "data/ghost-db/HT.KY21Gene.2.fasta.gz"))
cdna_ens <- readDNAStringSet("https://ftp.ensembl.org/pub/release-113/fasta/ciona_intestinalis/cdna/Ciona_intestinalis.KH.cdna.all.fa.gz")
ncrna_ens <- readDNAStringSet("https://ftp.ensembl.org/pub/release-113/fasta/ciona_intestinalis/ncrna/Ciona_intestinalis.KH.ncrna.fa.gz")

# get mitochondrial transcripts from ensembl fasta ============================
strip_fasta_description <- function(x) {
  str_replace(x, "(.+?)\\s.*", "\\1")
}

remove_id_version <- function(x) {
  str_replace(x, "(ENSCINT\\d+)\\.\\d+", "\\1")
}

names(cdna_ens) <- names(cdna_ens) %>% 
  strip_fasta_description() %>% 
  remove_id_version()

names(ncrna_ens) <- names(ncrna_ens) %>% 
  strip_fasta_description() %>% 
  remove_id_version()

mt_tx_ids <- gtf %>% 
  filter(seqnames == "MT") %>% 
  filter(!is.na(transcript_id)) %>% 
  pull(transcript_id) %>% 
  unique()

# Not all mitochondrial trascript ids are actually present in the cdna fasta
# from ensembl. Non-coding genes (tRNA, rRNA) genes are in the ncrna fasta file.

cdna_mt <- cdna_ens[mt_tx_ids[mt_tx_ids %in% names(cdna_ens)]]
ncrna_mt <- ncrna_ens[mt_tx_ids[mt_tx_ids %in% names(ncrna_ens)]]

tx_mt <- c(cdna_mt, ncrna_mt)

# remove transcripts missing from GTF =========================================

# There are 304 transcripts that are present in the KY21 transcriptome fasta
# file but are not to be found in the GTF file. We will remove them from the
# transcriptome fasta file consistency with the GTF.

names(tx_ky21)[!names(tx_ky21) %in% gtf$transcript_id[!is.na(gtf$transcript_id)]] # 304
gtf$transcript_id[!is.na(gtf$transcript_id)][!gtf$transcript_id[!is.na(gtf$transcript_id)] %in% names(tx_ky21)] # 0, except for ENSEMBL transcripts from MT chromosome

tx_ky21 <- tx_ky21[names(tx_ky21) %in% gtf$transcript_id[!is.na(gtf$transcript_id)]]

# append ======================================================================
tx_all <- c(tx_mt, tx_ky21)

stopifnot(length(tx_all) == length(unique(gtf$transcript_id[!is.na(gtf$transcript_id)])))

# write file ==================================================================
writeXStringSet(tx_all, 
                here(genomedir, 
                     "results/fasta/KY21_transcripts_with-mito.fasta.gz"), 
                format = "fasta", 
                compress = "gzip")

# end =========================================================================
