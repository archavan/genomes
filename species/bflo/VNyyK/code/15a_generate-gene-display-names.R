# Arun Chavan
# Started: 2023-06-08

# background ==================================================================
# We have gene annotations from 4 sources: RefSeq (just the gene descriptions),
# KEGG, emapper, and PROST. Here we try to merge them into a single table. The
# primary purpose of doing this is so that we have a single gene description for
# each gene id, which we can append to the gene ids when making gene expression
# dotplots.

# We will prefer the Refseq description; if that's missing, we will use KEGG
# description (based KO group mapped by KEGG annotation, not from emapper). From
# emapper, we have GO id mapping etc but not descriptions. But we will also use
# predicted protein names from emapper.

# set up ======================================================================
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/bflo/VNyyK/")

###############################################################################
# data ========================================================================
###############################################################################

gene_ids <- read_tsv(here(genomedir, "results/genes/gene-ids.tsv"))

ncbi <- read_tsv(
  here(genomedir, "results/genes/ncbi-genes-table.tsv"),
  col_types = cols(ncbi_gene_id = col_character())
) |>
  mutate(ncbi_gene_id = paste0("ncbi:", ncbi_gene_id))

kegg <- read_tsv(here(genomedir, "results/kegg/kegg-genes.tsv"))

emap_symbol <- read_tsv(here(
  genomedir,
  "results/emapper/processed/emapper_predicted-protein-name.tsv"
))

emap_desc <- read_tsv(here(
  genomedir,
  "results/emapper/processed/emapper_description.tsv"
))

prost <- read_tsv(here(
  genomedir,
  "results/prost/processed/best-hits_with-gene-names.tsv"
))

ncbi <- ncbi |>
  mutate(
    gene_type = case_when(
      gene_type == "PROTEIN_CODING" ~ "protein_coding",
      gene_type == "PSEUDO" ~ "pseudo",
      .default = gene_type
    )
  )

# emap.ko <- read_tsv(here(
#   genomedir,
#   "results/emapper/processed/emapper_kegg_ko.tsv"
# ))
# emap.name <- read_tsv(here(
#   genomedir,
#   "results/emapper/processed/emapper_predicted-protein-name.tsv"
# ))

## clean up emapper symbols ===================================================
emap_symbol_unique <- emap_symbol |>
  group_by(gene_id) |>
  summarize(predicted_protein_name = list(predicted_protein_name)) |>
  ungroup() |>
  mutate(n_symbols = map_int(predicted_protein_name, length)) |>
  filter(n_symbols == 1) |>
  select(-n_symbols) |>
  unnest_longer(col = predicted_protein_name)

emap_symbol_multi <- emap_symbol |>
  group_by(gene_id) |>
  summarise(predicted_protein_name = list(predicted_protein_name)) |>
  ungroup() |>
  mutate(n_symbols = map_int(predicted_protein_name, length)) |>
  filter(n_symbols > 1) |>
  select(-n_symbols) |>
  left_join(ncbi, by = c("gene_id" = "ncbi_gene_id")) |>
  mutate(
    emapper_symbol = map_chr(predicted_protein_name, paste0, collapse = ";")
  ) |>
  select(
    gene_id,
    emapper_symbol,
    ncbi_description = description
  )

# generate case_when statements and copy-paste into the mutate query below
emap_symbol_multi |>
  mutate(
    cw = glue(
      "gene_id == \"{gene_id}\" & emapper_symbol == \"{emapper_symbol}\" ~ TRUE, # {ncbi_description}"
    )
  ) |>
  pull(cw) |>
  writeLines()

#fmt: skip
emap_symbol_multi_filtered <- emap_symbol_multi |> 
  mutate(emapper_symbol = str_split(emapper_symbol, ";")) |>
  unnest_longer(col = emapper_symbol) |> 
  mutate(
    retain = case_when(
      gene_id == "ncbi:118403812" & emapper_symbol == "CLASP1"  ~ TRUE, # CLIP-associating protein 1-B-like
      gene_id == "ncbi:118404105" & emapper_symbol == "TPD52"   ~ TRUE, # tumor protein D54-like
      gene_id == "ncbi:118404778" & emapper_symbol == "YWHAB"   ~ TRUE, # 14-3-3 protein zeta-like
      gene_id == "ncbi:118406164" & emapper_symbol == "SPATA7"  ~ TRUE, # spermatogenesis-associated protein 7 homolog
      gene_id == "ncbi:118406541" & emapper_symbol == "NAV2"    ~ TRUE, # neuron navigator 2-like
      gene_id == "ncbi:118406715" & emapper_symbol == "SCN1A"   ~ TRUE, # sodium channel protein type 4 subunit alpha-like
      gene_id == "ncbi:118407675" & emapper_symbol == "NUDT17"  ~ TRUE, # nucleoside diphosphate-linked moiety X motif 17-like
      gene_id == "ncbi:118407854" & emapper_symbol == "RBFOX3"  ~ TRUE, # RNA binding protein fox-1 homolog 3-like
      gene_id == "ncbi:118409294" & emapper_symbol == "BNC2"    ~ TRUE, # zinc finger protein basonuclin-2-like
      gene_id == "ncbi:118409440" & emapper_symbol == "TMEM234" ~ TRUE, # transmembrane protein 234-like
      gene_id == "ncbi:118410506" & emapper_symbol == "nrx-1"   ~ TRUE, # contactin-associated protein-like 2
      gene_id == "ncbi:118410823" & emapper_symbol == "SPOCK2"  ~ TRUE, # testican-2-like
      gene_id == "ncbi:118411124" & emapper_symbol == "SYNDIG1" ~ TRUE, # uncharacterized LOC118411124
      gene_id == "ncbi:118411158" & emapper_symbol == "MEGF10"  ~ TRUE, # multiple epidermal growth factor-like domains protein 11
      gene_id == "ncbi:118411163" & emapper_symbol == "ANKRD17" ~ TRUE, # ankyrin repeat domain-containing protein 17-like
      gene_id == "ncbi:118411573" & emapper_symbol == "CAMK2D"  ~ TRUE, # calcium/calmodulin-dependent protein kinase type II subunit delta-like
      gene_id == "ncbi:118411981" & emapper_symbol == "HINT1"   ~ TRUE, # histidine triad nucleotide-binding protein 1-like
      gene_id == "ncbi:118412562" & emapper_symbol == "NRG2"    ~ TRUE, # pro-neuregulin-2, membrane-bound isoform-like
      gene_id == "ncbi:118412595" & emapper_symbol == "ERLIN2"  ~ TRUE, # erlin-2-B-like
      gene_id == "ncbi:118413292" & emapper_symbol == "NFYA"    ~ TRUE, # nuclear transcription factor Y subunit alpha-like
      gene_id == "ncbi:118413358" & emapper_symbol == "SRSF10"  ~ TRUE, # serine/arginine-rich splicing factor 10-like
      gene_id == "ncbi:118413730" & emapper_symbol == "PPP2R1A" ~ TRUE, # serine/threonine-protein phosphatase 2A 65 kDa regulatory subunit A
      gene_id == "ncbi:118413855" & emapper_symbol == "FGFR3"   ~ TRUE, # fibroblast growth factor receptor 2-like
      gene_id == "ncbi:118414456" & emapper_symbol == "KLHL6"   ~ TRUE, # kelch-like protein 24
      gene_id == "ncbi:118414743" & emapper_symbol == "REEP2"   ~ TRUE, # receptor expression-enhancing protein 2-like
      gene_id == "ncbi:118415340" & emapper_symbol == "IMPDH1"  ~ TRUE, # inosine-5'-monophosphate dehydrogenase 1b-like
      gene_id == "ncbi:118416933" & emapper_symbol == "RBM5"    ~ TRUE, # RNA-binding protein 5-like
      gene_id == "ncbi:118418184" & emapper_symbol == "ROBO1"   ~ TRUE, # roundabout homolog 1-like
      gene_id == "ncbi:118418954" & emapper_symbol == "CDX1"    ~ TRUE, # homeobox protein CHOX-CAD-like
      gene_id == "ncbi:118419074" & emapper_symbol == "ASB11"   ~ TRUE, # ankyrin repeat and SOCS box protein 13-like
      gene_id == "ncbi:118421334" & emapper_symbol == "BANF2"   ~ TRUE, # barrier-to-autointegration factor B-like
      gene_id == "ncbi:118421376" & emapper_symbol == "PIGC"    ~ TRUE, # phosphatidylinositol N-acetylglucosaminyltransferase subunit C-like
      gene_id == "ncbi:118422770" & emapper_symbol == "PAAF1"   ~ TRUE, # proteasomal ATPase-associated factor 1-like
      gene_id == "ncbi:118425250" & emapper_symbol == "TMEM138" ~ TRUE, # transmembrane protein 138-like
      gene_id == "ncbi:118425437" & emapper_symbol == "SPRYD4"  ~ TRUE, # E3 ubiquitin-protein ligase TRIM41-like
      gene_id == "ncbi:118425463" & emapper_symbol == "MSI2"    ~ TRUE, # RNA-binding protein Musashi homolog 2-like
      gene_id == "ncbi:118425543" & emapper_symbol == "MLXIP"   ~ TRUE, # MLX-interacting protein-like
      gene_id == "ncbi:118425676" & emapper_symbol == "AUTS2"   ~ TRUE, # autism susceptibility gene 2 protein-like
      gene_id == "ncbi:118426824" & emapper_symbol == "MINK1"   ~ TRUE, # mitogen-activated protein kinase kinase kinase kinase 4-like
      gene_id == "ncbi:118427201" & emapper_symbol == "NOVA1"   ~ TRUE, # RNA-binding protein Nova-1-like
      gene_id == "ncbi:118427250" & emapper_symbol == "LRCH1"   ~ TRUE, # leucine-rich repeat and calponin homology domain-containing protein 1-like
      gene_id == "ncbi:118427510" & emapper_symbol == "CPZ"     ~ TRUE, # mucin-5AC-like
      gene_id == "ncbi:118427601" & emapper_symbol == "COLEC11" ~ TRUE, # EMI domain-containing protein 1-like
      gene_id == "ncbi:118427794" & emapper_symbol == "DLG1"    ~ TRUE, # disks large homolog 2-like
      gene_id == "ncbi:118427949" & emapper_symbol == "FAM89A"  ~ TRUE, # protein FAM89A-like
      gene_id == "ncbi:118428320" & emapper_symbol == "MBNL1"   ~ TRUE, # muscleblind-like protein 3
      gene_id == "ncbi:118429098" & emapper_symbol == "PITPNC1" ~ TRUE, # cytoplasmic phosphatidylinositol transfer protein 1-like
      gene_id == "ncbi:118429842" & emapper_symbol == "ACTG2"   ~ TRUE, # actin, alpha cardiac muscle 2-like
      gene_id == "ncbi:118430371" & emapper_symbol == "PGLYRP1" ~ TRUE, # peptidoglycan-recognition protein SC2-like
      gene_id == "ncbi:118431774" & emapper_symbol == "dock7"   ~ TRUE # dedicator of cytokinesis protein 7-like
    )
  ) |> 
  filter(retain) |> 
  select(gene_id, predicted_protein_name = emapper_symbol)

emap_symbol_clean <- bind_rows(emap_symbol_unique, emap_symbol_multi_filtered)

## clean up emapper descriptions ==============================================
# we will just use the first description for each gene id for simplicity

retain_first_entry <- function(x) {
  x |>
    group_by(gene_id) |>
    slice_head(n = 1) |>
    ungroup()
}

emap_desc_clean <- retain_first_entry(emap_desc)

## put together emapper symbols and descriptions ==============================
emap <- full_join(emap_symbol_clean, emap_desc_clean, by = "gene_id") |>
  rename(emap_symbol = predicted_protein_name, emap_description = description)

stopifnot(!any(duplicated(emap$gene_id)))

## prost ======================================================================
prost <- prost |>
  select(gene_id, prost_symbol = gene_name, prost_description = protein_type) |>
  distinct()

stopifnot(!any(duplicated(prost$gene_id)))

###############################################################################
# generate gene display names =================================================
###############################################################################
add_source_prefix <- function(x, prefix) {
  ifelse(!is.na(x), glue("{prefix}{x}"), NA)
}

genes <- reduce(
  list(
    select(
      ncbi,
      gene_id = ncbi_gene_id,
      ncbi_symbol = symbol,
      ncbi_description = description,
      ncbi_gene_type = gene_type
    ),
    emap,
    prost
  ),
  left_join,
  by = "gene_id"
) |>
  mutate(
    ncbi_symbol = add_source_prefix(ncbi_symbol, "[n] "),
    ncbi_description = add_source_prefix(ncbi_description, "[n] "),
    emap_symbol = add_source_prefix(emap_symbol, "[e] "),
    emap_description = add_source_prefix(emap_description, "[e] "),
    prost_symbol = add_source_prefix(prost_symbol, "[p] "),
    prost_description = add_source_prefix(prost_description, "[p] ")
  ) |>
  mutate(
    ncbi_symbol = case_when(
      str_detect(ncbi_symbol, "LOC\\d+") ~ NA_character_,
      .default = ncbi_symbol
    )
  ) |>
  mutate(
    ncbi_description = case_when(
      str_detect(ncbi_description, "uncharacterized") ~ NA_character_,
      .default = ncbi_description
    )
  )

# coalesce ====================================================================
# we use the order of columns in coalesce() to set preference: prefer
# ncbi if exists, if not use emapper, then emapper, and prost
display_names <- genes |>
  mutate(
    symbol = coalesce(
      ncbi_symbol,
      emap_symbol,
      prost_symbol
    ),
    description = coalesce(
      ncbi_description,
      prost_description, # prost descriptions are better than emapper's
      emap_description
    )
  ) |>
  select(gene_id, symbol, description, ncbi_gene_type)

# do all genes have some symbol or description?
display_names |>
  filter(is.na(symbol) | is.na(description))

display_names <- display_names |>
  mutate(
    description = ifelse(
      is.na(description),
      glue("uncharacterized ({ncbi_gene_type})"),
      description
    )
  )

display_names <- select(display_names, -ncbi_gene_type)

# merge symbol and description into a single name =============================
display_names <- display_names |>
  mutate(
    display_name = paste0(
      "[",
      gene_id,
      "]",
      "; ",
      ifelse(is.na(symbol), "", paste0(symbol, "; ")),
      ifelse(is.na(description), "", paste0(description))
    )
  )

print(display_names, n = 500)

# write files =================================================================
write_tsv(
  display_names,
  here(genomedir, "results/genes/gene-display-names.tsv")
)

# end =========================================================================
