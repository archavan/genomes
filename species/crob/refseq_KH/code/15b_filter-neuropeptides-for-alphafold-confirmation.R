# Arun Chavan
# Started: 2024-03-28

# background ==================================================================

# Here we subset the neuropeptide geneset to a smaller number (~50, ideally) of
# genes that Kejue can test against AlphaFold to confirm the annotations.

# setup =======================================================================
library(seqinr)
library(tidyverse)
library(here)
library(Seurat)

genomedir <- here("species/crob/genome/refseq_KH")
datadir <- here("species/crob/10x/run-02/refseq_KH")

# data ========================================================================
np <- read_tsv(here(genomedir, "results/genes/geneset_neuropeptides.tsv"))
cr <- read_rds(here(datadir, "data/processed/03_clustered/clustered_v1.rds"))
plen <- read_tsv(here(genomedir, "results/genes/protein-lengths.tsv"))
fa <- read.fasta(here(genomedir, 
                      "data/refseq/GCF_000224145.3_KH_protein.faa.gz"), 
                 seqtype = "AA")

aggregate_counts <- function(x) {
  aggcts <- GetAssayData(x, slot = "counts") %*% model.matrix(~ 0 + Idents(x))
  colnames(aggcts) <- gsub("Idents\\(x\\)", "g", colnames(aggcts))
  aggcts <- as.matrix(aggcts)
  return(aggcts)
}

agg_cts_by_gene <- rowSums(aggregate_counts(cr))
exp_genes <- names(agg_cts_by_gene[which(agg_cts_by_gene > 0)])

# filter geneset for AlphaFold confirmation ===================================

# add aggregated counts information to help make decision about which genes to
# prioritize. i deally we want only 50 genes or so to test against AlphaFold.

np <- np %>% 
  mutate(agg_cts = map_int(gene_id, ~agg_cts_by_gene[.x])) %>% 
  arrange(-agg_cts)

np_1 <- np %>% 
  filter(agg_cts > 0) %>%
  filter(str_detect(symbol, "^\\[p\\]") | 
           str_detect(description, "^\\[p\\]")) %>%  # annotated with prost
  filter(!str_detect(description, "[Uu]ncharacterized")) %>%  
  filter(!str_detect(description, "[Ii]nterleukin")) %>%
  filter(!str_detect(description, "[Ii]nterferon")) %>%
  filter(!str_detect(description, "[Mm]ucin")) %>%
  filter(!str_detect(description, "[Ff]ragment")) %>%
  filter(!str_detect(description, "[Aa]polipoprotein")) %>%
  filter(!str_detect(description, "epidermal growth factor")) %>%
  filter(!str_detect(description, "slit homolog 2 protein")) %>%
  filter(!str_detect(description, "[Tt]ransmembrane")) %>%
  filter(!str_detect(description, "[Ss]erine protease")) %>%
  filter(!str_detect(description, "[Mm]eiosis-specific")) %>%
  filter(!str_detect(description, "bypass of stop codon protein")) %>%
  filter(!str_detect(description, "von Willebrand factor")) %>%
  filter(!str_detect(description, "hedgehog-interacting protein-like")) %>%
  filter(!str_detect(description, "Phospholipase A2 inhibitor")) %>%
  filter(!str_detect(description, "ALKAL")) %>%
  filter(!str_detect(description, "UL131A")) %>% 
  filter(!str_detect(description, "^\\[n\\]")) %>%
  filter(!str_detect(description, "UPF0375 protein C08F11.12")) %>%
  filter(!str_detect(description, "Meiotically up-regulated gene 117")) %>%
  filter(!str_detect(description, "Tan_10cys"))

nrow(np_1)
#> [1] 71

# We still have 71 genes. we need to filter further. We will simply pick top 50
# most expressed genes.

np_top <- slice_max(np_1, order_by = agg_cts, n = 50)

# generate fasta file for AlphaFold ===========================================
# We will use the longest protein for each gene. 

np_top <- left_join(np_top, plen, by = c("gene_id" = "gene_id")) %>% 
  group_by(gene_id) %>% 
  slice_max(order_by = length, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  arrange(-agg_cts) %>% 
  print(n = Inf)

np_top_fa <- fa[np_top$protein_id] 

# write files =================================================================
fs::dir_create(here(genomedir, "results/prost/alphafold"))
write.fasta(np_top_fa, 
            names(np_top_fa),
            here(genomedir, 
                 "results/prost/alphafold/neuropeptides.fa"))
write_tsv(np_top, here(genomedir, "results/prost/alphafold/neuropeptides.tsv"))

# end =========================================================================
