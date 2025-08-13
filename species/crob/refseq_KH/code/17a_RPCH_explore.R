# Arun Chavan
# Started: 2024-05-09

# background ==================================================================

# Investigate the ciona gene that prost aligns with red pigment concentrating
# hormone (RPCH) and happens to be the most highly expressed gene in the
# dataset.

# setup =======================================================================
library(seqinr)
library(tidyverse)
library(here)
library(Seurat)
library(httr2)

genomedir <- here("species/crob/genome/refseq_KH")
datadir <- here("species/crob/10x/run-02/refseq_KH")

resdir <- here(genomedir, "results/neuropeptides/rpch")
fs::dir_create(resdir)

# genomic data ================================================================
genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
plen <- read_tsv(here(genomedir, "results/genes/protein-lengths.tsv"))
fa <- read.fasta(here(genomedir, 
                      "data/refseq/GCF_000224145.3_KH_protein.faa.gz"), 
                 seqtype = "AA")

prost <- read_tsv(here(genomedir, "results/prost/processed/best-hits.tsv"))

# single cell data ============================================================
cr <- read_rds(here(datadir, "data/processed/03_clustered/clustered_v1.rds"))

aggregate_counts <- function(x) {
  aggcts <- GetAssayData(x, slot = "counts") %*% model.matrix(~ 0 + Idents(x))
  colnames(aggcts) <- gsub("Idents\\(x\\)", "g", colnames(aggcts))
  aggcts <- as.matrix(aggcts)
  return(aggcts)
}

agg_cts_by_gene <- rowSums(aggregate_counts(cr))

genes$agg_cts <- agg_cts_by_gene[genes$gene_id]
genes <- arrange(genes, -agg_cts)

# get sequences of RPCH isoforms ==============================================
# The highest expressed gene in the entire dataset is ncbi:100178806, which is
# annotated by PROST as Red pigment-concentrating prohormone. 

# Here we get sequences for all protein isoforms of this gene. There are 8.
rpch_id <- names(which.max(agg_cts_by_gene))
rpch_isoforms <- filter(plen, gene_id == rpch_id)
rpch_prot_seq <- fa[rpch_isoforms$protein_id]

# get sequence of the PROST hit RPCH ==========================================
rpch_prost <- filter(prost, gene_id == rpch_id)

# UNIPROT best hit for the longest protein isoform
rpch_uniprot_hit <- rpch_prost %>% 
  left_join(plen, by = c("gene_id", "protein_id")) %>% 
  slice_max(length, n = 1, with_ties = FALSE) %>% 
  pull(uniprot_id)

uniprot_req <- request(paste0(
  "https://rest.uniprot.org/uniprotkb/", rpch_uniprot_hit
))

resp <- req_perform(uniprot_req)
resp_json <- resp_body_json(resp)

rpch_carma_seq <- as.SeqFastaAA(s2c(resp_json$sequence$value), 
                                name = resp_json$uniProtkbId)

# add both ciona and C. maema sequences to a single file to run BLAST =========
rpch_prot_seq$RPCH_CARMA <- rpch_carma_seq

# write fasta file ============================================================
write.fasta(rpch_prot_seq,
            names = names(rpch_prot_seq),
            file.out = here(resdir, "RPCH_seq_crob-and-cmae.fa"))


# end =========================================================================

uniprot_to_fasta <- function(.uniprot_id) {
  req <- request(paste0(
    "https://rest.uniprot.org/uniprotkb/", .uniprot_id
  ))
  
  resp <- req_perform(req)
  resp_json <- resp_body_json(resp)
  
  fa <- as.SeqFastaAA(s2c(resp_json$sequence$value), 
                      name = resp_json$uniProtkbId)
  
  return(fa)
}

uniprot_to_fasta(rpch_uniprot_hit)

new_al <- vector("list")

new_al$RPCH_CARMA <- uniprot_to_fasta(rpch_uniprot_hit)
new_al$XP_002131878.2 <- rpch_prot_seq$XP_002131878.2
write.fasta(new_al,
            names = names(new_al),
            file.out = here(resdir, "RPCH_seq_crob-and-cmae_longest.fa"))
