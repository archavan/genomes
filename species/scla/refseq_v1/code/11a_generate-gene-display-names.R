# Arun Chavan
# Started: 2025-08-18

# background ==================================================================

# Generate gene display names for plotting combining symbols and descriptions
# from ncbi, emapper, and PROST.

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)

genomedir <- here("species/scla/refseq_v1")
orthdir <- here("orthofinder/2025-08-29")

###############################################################################
# data ========================================================================
###############################################################################

## orthofinder ================================================================
orth <- read_tsv(here(orthdir, "results/pairwise-orthologs/mmus_vs_scla.tsv"))
hog_labels <- read_tsv(here(orthdir, "results/hog-labels/mmus_vs_scla.tsv"))

orth <- orth |>
  left_join(select(hog_labels, hog, label), by = "hog") |>
  select(gene_id = scla, hog_label = label) |>
  distinct()

## ncbi =======================================================================
ncbi <- read_tsv(here(genomedir, "results/genes/ncbi-genes-table.tsv"))

ncbi <- ncbi |>
  mutate(
    gene_type = case_when(
      gene_type == "PROTEIN_CODING" ~ "protein-coding",
      gene_type == "PSEUDO" ~ "pseudo",
      .default = gene_type
    )
  )

## kegg =======================================================================
kegg <- read_tsv(here(genomedir, "results/kegg/kegg-genes.tsv"))

## emapper ====================================================================
emap_symbol <- read_tsv(here(
  genomedir,
  "results/emapper/processed/emapper_predicted-protein-name.tsv"
))
emap_desc <- read_tsv(here(
  genomedir,
  "results/emapper/processed/emapper_description.tsv"
))

## prost ======================================================================
prost <- read_tsv(here(
  genomedir,
  "results/prost/processed/best-hits_with-gene-names.tsv"
))

prost <- prost |>
  select(gene_id, prost_symbol = gene_name, prost_description = protein_type)

stopifnot(nrow(filter(prost, duplicated(gene_id))) == 0)

###############################################################################
# clean up emapper output =====================================================
###############################################################################

## clean up emapper symbols ===================================================
# For emapper we kept all symbols and descriptions for all genes. We will pick
# one depending on what makes most sense based on ncbi description

emap_symbol_unique <- emap_symbol |>
  group_by(gene_id) |>
  summarise(predicted_protein_name = list(predicted_protein_name)) |>
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
  left_join(ncbi, by = "gene_id") |>
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
      gene_id == "ncbi.120325390" & emapper_symbol == "CYP2C8"   ~ TRUE, # cytochrome P450 2C42-like
      gene_id == "ncbi.120325915" & emapper_symbol == "DIS3"     ~ TRUE, # exosome complex exonuclease RRP44-like
      gene_id == "ncbi.120325923" & emapper_symbol == "HHIPL1"   ~ TRUE, # HHIP-like protein 1
      gene_id == "ncbi.120327742" & emapper_symbol == "GYG1"     ~ TRUE, # glycogenin-1-like
      gene_id == "ncbi.120328013" & emapper_symbol == "ACTR5"    ~ TRUE, # actin-related protein 5-like
      gene_id == "ncbi.120328025" & emapper_symbol == "WNK1"     ~ TRUE, # uncharacterized LOC120328025
      gene_id == "ncbi.120328537" & emapper_symbol == "RFT1"     ~ TRUE, # protein RFT1 homolog
      gene_id == "ncbi.120328564" & emapper_symbol == "HS3ST1"   ~ TRUE, # heparan sulfate glucosamine 3-O-sulfotransferase 4-like
      gene_id == "ncbi.120329022" & emapper_symbol == "WAS"      ~ TRUE, # neural Wiskott-Aldrich syndrome protein-like
      gene_id == "ncbi.120329047" & emapper_symbol == "ACER3"    ~ TRUE, # alkaline ceramidase 3-like
      gene_id == "ncbi.120329435" & emapper_symbol == "MYO18A"   ~ TRUE, # unconventional myosin-XVIIIa-like
      gene_id == "ncbi.120330165" & emapper_symbol == "PRODH"    ~ TRUE, # hydroxyproline dehydrogenase-like
      gene_id == "ncbi.120330337" & emapper_symbol == "TRAK2"    ~ TRUE, # trafficking kinesin-binding protein 2-like
      gene_id == "ncbi.120330351" & emapper_symbol == "PIAS1"    ~ TRUE, # E3 SUMO-protein ligase PIAS1-like
      gene_id == "ncbi.120330690" & emapper_symbol == "SULT1B1"  ~ TRUE, # sulfotransferase family cytosolic 1B member 1-like
      gene_id == "ncbi.120330785" & emapper_symbol == "PROM1"    ~ TRUE, # prominin-1-A-like
      gene_id == "ncbi.120331594" & emapper_symbol == "ZDHHC5"   ~ TRUE, # palmitoyltransferase ZDHHC5-A-like
      gene_id == "ncbi.120331656" & emapper_symbol == "WSCD1"    ~ TRUE, # WSC domain-containing protein 1-like
      gene_id == "ncbi.120331905" & emapper_symbol == "DTX3"     ~ TRUE, # serine/threonine protein kinase OSK1-like
      gene_id == "ncbi.120331979" & emapper_symbol == "ATL2"     ~ TRUE, # uncharacterized LOC120331979
      gene_id == "ncbi.120331988" & emapper_symbol == "ZDHHC5"   ~ TRUE, # palmitoyltransferase ZDHHC5-A-like
      gene_id == "ncbi.120332134" & emapper_symbol == "RAB11B"   ~ TRUE, # ras-related protein Rab-11B-like
      gene_id == "ncbi.120332346" & emapper_symbol == "TBX21"    ~ TRUE, # T-box transcription factor TBX6-like
      gene_id == "ncbi.120332681" & emapper_symbol == "CERS4"    ~ TRUE, # ceramide synthase 5-like
      gene_id == "ncbi.120332819" & emapper_symbol == "TCF12"    ~ TRUE, # transcription factor 12-like
      gene_id == "ncbi.120333353" & emapper_symbol == "ATL1"     ~ TRUE, # atlastin-like
      gene_id == "ncbi.120333359" & emapper_symbol == "ODC1"     ~ TRUE, # ornithine decarboxylase-like
      gene_id == "ncbi.120333420" & emapper_symbol == "FMO1"     ~ TRUE, # dimethylaniline monooxygenase [N-oxide-forming] 1-like
      gene_id == "ncbi.120334652" & emapper_symbol == "TSPAN8"   ~ TRUE, # tetraspanin-8-like
      gene_id == "ncbi.120335643" & emapper_symbol == "NRXN3"    ~ TRUE, # neurexin-3a-like
      gene_id == "ncbi.120337510" & emapper_symbol == "ZNFX1"    ~ TRUE, # NFX1-type zinc finger-containing protein 1-like
      gene_id == "ncbi.120337595" & emapper_symbol == "HS3ST4"   ~ TRUE, # heparan sulfate glucosamine 3-O-sulfotransferase 1-like
      gene_id == "ncbi.120337928" & emapper_symbol == "ODC1"     ~ TRUE, # antizyme inhibitor 2-like
      gene_id == "ncbi.120337945" & emapper_symbol == "FAM134A"  ~ TRUE, # reticulophagy regulator 3-like
      gene_id == "ncbi.120337953" & emapper_symbol == "SLC6A11"  ~ TRUE, # sodium- and chloride-dependent GABA transporter 3-like
      gene_id == "ncbi.120338197" & emapper_symbol == "JAK2"     ~ TRUE, # tyrosine-protein kinase JAK2-like
      gene_id == "ncbi.120338234" & emapper_symbol == "PARD3"    ~ TRUE, # partitioning defective 3 homolog
      gene_id == "ncbi.120338418" & emapper_symbol == "MFGE8"    ~ TRUE, # lactadherin-like
      gene_id == "ncbi.120338655" & emapper_symbol == "CNKSR2"   ~ TRUE, # CNK3/IPCEF1 fusion protein-like
      gene_id == "ncbi.120338717" & emapper_symbol == "GLIS3"    ~ TRUE, # zinc finger protein GLIS3-like
      gene_id == "ncbi.120339696" & emapper_symbol == "CDK17"    ~ TRUE, # cyclin-dependent kinase 17-like
      gene_id == "ncbi.120339859" & emapper_symbol == "BHLHE40"  ~ TRUE, # hairy/enhancer-of-split related with YRPW motif-like protein
      gene_id == "ncbi.120340200" & emapper_symbol == "WNT5A"    ~ TRUE, # protein Wnt-5-like
      gene_id == "ncbi.120340441" & emapper_symbol == "SMARCC2"  ~ TRUE, # SWI/SNF complex subunit SMARCC2-like
      gene_id == "ncbi.120340663" & emapper_symbol == "TFPI"     ~ TRUE, # papilin-like
      gene_id == "ncbi.120341158" & emapper_symbol == "FLNC"     ~ TRUE, # filamin-C-like
      gene_id == "ncbi.120342233" & emapper_symbol == "HSPA8"    ~ TRUE, # heat shock 70 kDa protein II-like
      gene_id == "ncbi.120342265" & emapper_symbol == "MYL3"     ~ TRUE, # myosin catalytic light chain, smooth muscle-like
      gene_id == "ncbi.120342610" & emapper_symbol == "HNRNPM"   ~ TRUE, # heterogeneous nuclear ribonucleoprotein M-like
      gene_id == "ncbi.120342629" & emapper_symbol == "HOXA4"    ~ TRUE, # homeobox protein Hox-A3-like
      gene_id == "ncbi.120342695" & emapper_symbol == "PPP1R2"   ~ TRUE, # protein phosphatase inhibitor 2-like
      gene_id == "ncbi.120343458" & emapper_symbol == "HSPA2"    ~ TRUE, # heat shock cognate 71 kDa protein-like
      gene_id == "ncbi.120344451" & emapper_symbol == "CHAT"     ~ TRUE, # probable vesicular acetylcholine transporter-B
      gene_id == "ncbi.120345149" & emapper_symbol == "FMO5"     ~ TRUE, # flavin-containing monooxygenase 5-like
      gene_id == "ncbi.120345170" & emapper_symbol == "FGFR4"    ~ TRUE, # fibroblast growth factor receptor 1-like
      gene_id == "ncbi.120345183" & emapper_symbol == "TNC"      ~ TRUE, # ryncolin-1-like
      gene_id == "ncbi.120346013" & emapper_symbol == "FAM204A"  ~ TRUE, # protein FAM204A-like
      gene_id == "ncbi.120346135" & emapper_symbol == "CPB2"     ~ TRUE, # carboxypeptidase B-like
      gene_id == "ncbi.120346362" & emapper_symbol == "ADAMTS12" ~ TRUE, # A disintegrin and metalloproteinase with thrombospondin motifs 12-like
      gene_id == "ncbi.120346370" & emapper_symbol == "TIAL1"    ~ TRUE, # nucleolysin TIAR-like
      gene_id == "ncbi.120346645" & emapper_symbol == "KCNB1"    ~ TRUE, # potassium voltage-gated channel subfamily F member 1-like
      gene_id == "ncbi.120346677" & emapper_symbol == "GCNT1"    ~ TRUE, # beta-1,3-galactosyl-O-glycosyl-glycoprotein beta-1,6-N-acetylglucosaminyltransferase 3-like
      gene_id == "ncbi.120347172" & emapper_symbol == "PPP3CA"   ~ TRUE, # serine/threonine-protein phosphatase 2B catalytic subunit
      gene_id == "ncbi.120348134" & emapper_symbol == "RUFY2"    ~ TRUE, # RUN and FYVE domain-containing protein 2-like
      gene_id == "ncbi.120348254" & emapper_symbol == "PAPLN"    ~ TRUE, # papilin-like
      gene_id == "ncbi.120348402" & emapper_symbol == "FGR"      ~ TRUE, # tyrosine-protein kinase Yes-like
      gene_id == "ncbi.120348683" & emapper_symbol == "HTR1A"    ~ TRUE, # 5-hydroxytryptamine receptor 2B-like
      .default = FALSE
    )
  ) |> 
  filter(retain) |> 
  select(gene_id, predicted_protein_name = emapper_symbol)


emap_symbol_clean <- bind_rows(emap_symbol_unique, emap_symbol_multi_filtered)

## clean up emapper descriptions ==============================================
# we will just reatain the first one for duplicates for simplicity

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

stopifnot(nrow(filter(emap, duplicated(gene_id))) == 0)

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
      gene_id,
      ncbi_symbol = symbol,
      ncbi_description = description,
      ncbi_gene_type = gene_type
    ),
    orth,
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
# ghost_symbol if exists, if not use hsap_homolog_symbol, then emapper, and prost
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
  select(gene_id, hog_label, symbol, description, ncbi_gene_type)

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
    display_hog = ifelse(
      is.na(hog_label),
      NA,
      paste0(
        "[",
        gene_id,
        "]; ",
        hog_label,
        "; ",
        ifelse(is.na(description), "", description)
      )
    )
  ) |>
  mutate(
    display_non_hog = paste0(
      "[",
      gene_id,
      "]",
      "; ",
      ifelse(is.na(symbol), "", paste0(symbol, "; ")),
      ifelse(is.na(description), "", paste0(description))
    )
  ) |>
  mutate(display_name = coalesce(display_hog, display_non_hog))

print(display_names, n = 500)

# write files =================================================================
write_tsv(
  display_names |>
    select(gene_id, hog_label, symbol, description, display_name),
  here(genomedir, "results/genes/gene-display-names.tsv")
)

# end =========================================================================
