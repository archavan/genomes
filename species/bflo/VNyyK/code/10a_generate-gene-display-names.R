# Arun Chavan
# Started: 2023-06-08

# background ==================================================================
# We have gene annotations from 3 sources: RefSeq (just the gene descriptions),
# KEGG, and emapper. Here we try to merge them into a single table. The primary
# purpose of doing this is so that we have a single gene description for each
# gene id, which we can append to the gene ids when making gene expression
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

# data ========================================================================
gene_ids <- read_tsv(here(genomedir, "results/genes/gene-ids.tsv"))
ncbi <- read_tsv(here(genomedir, "results/genes/ncbi-genes-table.tsv"))
kegg <- read_tsv(here(genomedir, "results/kegg/kegg-genes.tsv"))
emap.ko <- read_tsv(here(
  genomedir,
  "results/emapper/processed/emapper_kegg_ko.tsv"
))
emap.name <- read_tsv(here(
  genomedir,
  "results/emapper/processed/emapper_predicted-protein-name.tsv"
))

###############################################################################
# symbols =====================================================================
###############################################################################

symbols <- gene_ids %>%
  select(gene_id, ncbi_gene_id)

## add ncbi symbols -----------------------------------------------------------
ncbi.symbol <- ncbi %>%
  select(ncbi_gene_id, symbol) %>%
  filter(!is.na(symbol)) %>%
  filter(str_detect(symbol, "^LOC", negate = TRUE)) %>%
  distinct() %>%
  dplyr::rename(ncbi_symbol = symbol)

symbols <- left_join(symbols, ncbi.symbol, by = "ncbi_gene_id")

## add emapper predicted protein name -----------------------------------------

# 50 ncbi gene ids have multiple predicted names by emapper. For these genes, I
# will manually pick the best one that makes the most sense given the ncbi
# description or kegg name for those genes

# workflow:
# 1. If one of the names matches the refseq description, use it. If not,
# 2. If one of the names matches the kegg name, use it. if not,
# 3. look up gene id on www.orthodb.org; if one of the genes matches the name
#    of the orthology group, use it. If not,
# 4. free web search and use whatever information available to make in informed
#    decision. If non found,
# 5. use the first name.

emap.name.unique <- emap.name %>%
  group_by(gene_id) %>%
  summarise(names = list(predicted_protein_name)) %>%
  mutate(n_names = map_int(names, length)) %>%
  filter(n_names == 1) %>%
  unnest_longer(col = names) %>%
  select(gene_id, emapper_symbol = names, -n_names)

emap.name.multi <- emap.name %>%
  group_by(gene_id) %>%
  summarise(names = list(predicted_protein_name)) %>%
  mutate(n_names = map_int(names, length)) %>%
  arrange(-n_names) %>%
  mutate(pasted = map_chr(names, ~ paste0(.x, collapse = ","))) %>%
  filter(n_names > 1) %>%
  left_join(select(gene_ids, ncbi_gene_id, gene_id), by = "gene_id") %>%
  left_join(
    select(kegg, ncbi_gene_id, refseq_description, ko_name),
    by = "ncbi_gene_id"
  ) %>%
  unique()

emap.name.multi <- emap.name.multi %>%
  mutate(
    picked_name = case_when(
      refseq_description == "tumor protein D52-like isoform X1" ~ "TPD52",
      refseq_description == "CLIP-associating protein 1-like isoform X1" ~
        "CLASP1",
      refseq_description == "14-3-3 protein zeta-like isoform X1" ~
        "14-3-3zeta",
      refseq_description ==
        "spermatogenesis-associated protein 7-like isoform X1" ~
        "SPATA7",
      refseq_description == "neuron navigator 3-like isoform X1" ~ "NAV2",
      refseq_description ==
        "sodium channel protein type 3 subunit alpha-like isoform X1" ~
        "SCN1A",
      refseq_description ==
        "protein tyrosine phosphatase domain-containing protein 1-like isoform X1" ~
        "PTPDC1",
      refseq_description ==
        "RNA binding protein fox-1 homolog 3-like isoform X1" ~
        "RBFOX2",
      refseq_description == "zinc finger protein basonuclin-2-like isoform X1" ~
        "BNC2",
      refseq_description ==
        "protein phosphatase 1K, mitochondrial-like isoform X1" ~
        "PPM1K",
      refseq_description == "contactin-associated protein-like 2 isoform X1" ~
        "CNTNAP2",
      refseq_description == "uncharacterized protein LOC118410823 isoform X1" ~
        "SPOCK2",
      refseq_description == "uncharacterized protein LOC118411124 isoform X1" ~
        "SYNDIG1",
      refseq_description ==
        "multiple epidermal growth factor-like domains protein 10 isoform X1" ~
        "MEGF10",
      refseq_description ==
        "ankyrin repeat domain-containing protein 17-like isoform X1" ~
        "ANKRD17",
      refseq_description ==
        "calcium/calmodulin-dependent protein kinase type II delta chain-like isoform X1" ~
        "CAMK2D",
      refseq_description ==
        "histidine triad nucleotide-binding protein 1-like isoform X1" ~
        "HINT1",
      refseq_description == "uncharacterized protein LOC118412562 isoform X1" ~
        "NRG2",
      refseq_description == "erlin-2-B-like isoform X1" ~ "ERLIN2",
      refseq_description ==
        "nuclear transcription factor Y subunit alpha-like isoform X1" ~
        "NFYA",
      refseq_description ==
        "serine/arginine-rich splicing factor 10-like isoform X1" ~
        "SRSF10",
      refseq_description ==
        "serine/threonine-protein phosphatase 2A 65 kDa regulatory subunit A alpha isoform-like isoform X1" ~
        "PPP2R1A",
      refseq_description ==
        "fibroblast growth factor receptor 2-like isoform X1" ~
        "FGFR3",
      refseq_description == "kelch-like protein 24 isoform X1" ~ "KLHL6",
      refseq_description == "uncharacterized protein LOC118414743 isoform X1" ~
        "REEP2",
      refseq_description ==
        "inosine-5'-monophosphate dehydrogenase 1b-like isoform X1" ~
        "IMPDH1",
      refseq_description == "RNA-binding protein 10-like isoform X1" ~ "RBM10",
      refseq_description == "roundabout homolog 1-like isoform X1" ~ "ROBO1",
      refseq_description == "homeobox protein CHOX-CAD-like isoform X1" ~
        "CDX1",
      refseq_description ==
        "ankyrin repeat and SOCS box protein 13-like isoform X1" ~
        "ASB5",
      refseq_description ==
        "guanine nucleotide-binding protein subunit beta-like protein 1 isoform X1" ~
        "GNB1L",
      refseq_description ==
        "phosphatidylinositol N-acetylglucosaminyltransferase subunit C-like isoform X1" ~
        "PIGC",
      refseq_description ==
        "nuclear receptor subfamily 0 group B member 2-like isoform X1" ~
        "NR0B2",
      refseq_description ==
        "cancer-related nucleoside-triphosphatase homolog isoform X1" ~
        "NTPCR",
      refseq_description ==
        "E3 ubiquitin-protein ligase TRIM41-like isoform X1" ~
        "FADD",
      refseq_description ==
        "RNA-binding protein Musashi homolog 2-like isoform X1" ~
        "MSI2",
      refseq_description == "MLX-interacting protein-like isoform X1" ~ "MLXIP",
      refseq_description ==
        "autism susceptibility gene 2 protein homolog isoform X1" ~
        "AUTS2",
      refseq_description ==
        "mitogen-activated protein kinase kinase kinase kinase 4-like isoform X1" ~
        "MINK1",
      refseq_description == "RNA-binding protein Nova-1-like isoform X1" ~
        "NOVA1",
      refseq_description ==
        "leucine-rich repeat and calponin homology domain-containing protein 1-like isoform X1" ~
        "LRCH1",
      refseq_description == "mucin-17-like isoform X1" ~ "CPZ",
      refseq_description == "uncharacterized protein LOC118427601 isoform X1" ~
        "UGT2A1",
      refseq_description == "disks large homolog 4-like isoform X1" ~ "DLG1",
      refseq_description == "protein FAM89A-like isoform X1" ~ "FAM89A",
      refseq_description == "muscleblind-like protein 3 isoform X1" ~ "MBNL1",
      refseq_description ==
        "cytoplasmic phosphatidylinositol transfer protein 1-like isoform X1" ~
        "PITPNC1",
      refseq_description == "actin-5-like isoform X1" ~ "ACTG1",
      refseq_description ==
        "trans-aconitate 2-methyltransferase-like isoform X1" ~
        "jhamt",
      refseq_description ==
        "dedicator of cytokinesis protein 7-like isoform X1" ~
        "dock7"
    )
  )

emap.name.multi <- emap.name.multi %>%
  select(gene_id, emapper_symbol = picked_name)

symbols <- left_join(
  symbols,
  bind_rows(emap.name.unique, emap.name.multi),
  by = "gene_id"
)

# add symbols extracted from emapper ko name ----------------------------------
# Some gene ids have multiple KO ids associated with them. I don't think it's
# strictly true with KO ids, but for the sake of making a decision, let's assume
# that the larger KO ids are more specific and pick the largest one out of
# multiple KO ids.
emap.ko.name <- emap.ko %>%
  group_by(gene_id) %>%
  arrange(desc(kegg_ko_id), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

# only deal with gene ids for which we haven't found a symbol yet.
emap.ko.name <- emap.ko.name %>%
  filter(
    gene_id %in%
      (symbols %>%
        filter(is.na(ncbi_symbol) & is.na(emapper_symbol)) %>%
        pull(gene_id))
  )

# separate out gene ids that have multiple symbols listed and unique ones
emap.ko.name.unique <- emap.ko.name %>%
  separate(kegg_ko_name, into = c("symbol", "description"), ";\\s+") %>%
  mutate(symbol = map(symbol, ~ str_split(.x, ",\\s+")[[1]])) %>%
  mutate(n_names = map_int(symbol, length)) %>%
  filter(n_names == 1) %>%
  unnest_longer(col = symbol) %>%
  select(gene_id, emap_ko_symbol = symbol) %>%
  distinct()

emap.ko.name.multi <- emap.ko.name %>%
  separate(kegg_ko_name, into = c("symbol", "description"), ";\\s+") %>%
  mutate(symbol = map(symbol, ~ str_split(.x, ",\\s+")[[1]])) %>%
  mutate(n_names = map_int(symbol, length)) %>%
  filter(n_names > 1)

# for the multi ones, we will pick the first symbol. For most rows that I
# spot-checked, the first one seems to match the description.
emap.ko.name.multi <- emap.ko.name.multi %>%
  mutate(emap_ko_symbol = map_chr(symbol, `[[`, 1)) %>%
  select(gene_id, emap_ko_symbol)

# add to symbols
symbols <- left_join(
  symbols,
  bind_rows(emap.ko.name.unique, emap.ko.name.multi),
  by = "gene_id"
)

## put it all togehter --------------------------------------------------------
symbols <- symbols %>%
  mutate(
    symbol = case_when(
      !is.na(ncbi_symbol) ~ ncbi_symbol,
      !is.na(emapper_symbol) ~ emapper_symbol,
      !is.na(emap_ko_symbol) ~ emap_ko_symbol,
      TRUE ~ NA_character_
    )
  )

symbols <- symbols %>%
  select(gene_id, symbol)

## summary ---------------------------------------------------------------------
nrow(symbols)
#> [1] 29755

sum(is.na(symbols$symbol))
#> [1] 11504

# We have symbols for 18251 gene ids now and NAs for 11504 genes.

###############################################################################
# descriptions ================================================================
###############################################################################

# make one table with descriptions from the NCBI and KEGG
# 1. NCBI has gene descriptions
# 2. KEGG:
#       - KEGG KO name is a description, but
#       - KEGG also provides a column called refseq description.
#           - For the most part, these are identical to the NCBI descriptions
#           - Some are different. Some contain suffixes like "isoform X1" etc,
#             indicating that these were pulled from protein descriptions. We
#             can get rid of those suffixes to make them essentially the same as
#             NCBI description.
#       - There several transcription factors for which NCBI description and
#         refseq description only says "transcription factor", but KEGG name is
#         actually more informative, e.g, Insulin promoting factor-1, and so on.
#         In these cases, we want to prioritize using KO name instead of NCBI
#         description even if it exists.

descriptions <- gene_ids %>%
  select(gene_id, ncbi_gene_id) %>%
  left_join(
    (ncbi %>%
      select(ncbi_gene_id, description) %>%
      distinct()),
    by = "ncbi_gene_id"
  ) %>%
  left_join(
    (kegg %>%
      group_by(ncbi_gene_id) %>%
      arrange(desc(ko_id), .by_group = TRUE) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(ncbi_gene_id, refseq_description, ko_name)),
    by = "ncbi_gene_id"
  )

# remove extraneous prefixes and suffixes from the refseq description column
descriptions <- descriptions %>%
  mutate(
    refseq_description = str_remove(refseq_description, "\\s+isoform.+$")
  ) %>%
  mutate(
    refseq_description = str_remove(
      refseq_description,
      "^LOW QUALITY PROTEIN:\\s+"
    )
  ) %>%
  mutate(
    refseq_description = str_replace(
      refseq_description,
      "^uncharacterized protein",
      "uncharacterized"
    )
  ) %>%
  mutate(
    refseq_description = str_replace(
      refseq_description,
      "^uncharacterized LOC\\d+",
      "uncharacterized"
    )
  ) %>%
  mutate(refseq_description = str_remove(refseq_description, "\\s+precursor$"))

## start filling into a column for final descriptions -------------------------
# if ncbi == refseq, and ncbi != uncharacterized, use ncbi description
descriptions <- descriptions %>%
  mutate(
    final = ifelse(
      (description == refseq_description) &
        (str_detect(description, "uncharacterized", negate = TRUE)),
      description,
      NA_character_
    )
  )

sum(is.na(descriptions$final))
#> [1] 9741

# if ncbi and refseq both uncharacterized and ko_name not NA, use KO name
descriptions <- descriptions %>%
  mutate(
    final = ifelse(
      is.na(final),
      ifelse(
        str_detect(description, "uncharacterized") &
          str_detect(refseq_description, "uncharacterized") &
          !is.na(ko_name),
        ko_name,
        NA_character_
      ),
      final
    )
  )

sum(is.na(descriptions$final))
#> [1] 9258

# if ncbi == refseq == transcription factor protein and KO is not NA (and often a
# more specific description of the TF), use KO name
# this was true for ciona but not for amphioxus
# descriptions <- descriptions %>%
#   mutate(final = ifelse(
#     description == "transcription factor protein" &
#       refseq_description == "transcription factor protein" &
#       !is.na(ko_name),
#     ko_name,
#     final
#   ))

# descriptions <- descriptions %>%
#   mutate(final = ifelse(
#     description == "transcription factor" &
#       refseq_description == "transcription factor" &
#       !is.na(ko_name),
#     ko_name,
#     final
#   ))

# if ncbi != refseq, and ncbi != uncharacterized, use NCBI
descriptions <- descriptions %>%
  mutate(
    final = ifelse(
      is.na(final),
      ifelse(
        description != refseq_description &
          str_detect(description, "uncharacterized", negate = TRUE),
        description,
        NA_character_
      ),
      final
    )
  )

sum(is.na(descriptions$final))
#> [1] 8710

# if ncbi != refseq, and ncbi == uncharacterized, refseq != uncharacterized, use refseq
descriptions <- descriptions %>%
  mutate(
    final = ifelse(
      is.na(final),
      ifelse(
        description != refseq_description &
          str_detect(description, "uncharacterized") &
          str_detect(refseq_description, "uncharacterized", negate = TRUE),
        refseq_description,
        NA_character_
      ),
      final
    )
  )

sum(is.na(descriptions$final))
#> [1] 8706

# if ncbi == uncharacterized, refseq == uncharacterized, and ko_name == NA, use "uncharacterized"
descriptions <- descriptions %>%
  mutate(
    final = ifelse(
      is.na(final),
      ifelse(
        str_detect(description, "uncharacterized") &
          str_detect(refseq_description, "uncharacterized") &
          is.na(ko_name),
        "uncharacterized",
        NA_character_
      ),
      final
    )
  )

sum(is.na(descriptions$final))
#> [1] 2894

# for the rest, if ncbi != uncharacterized, use ncbi; else uncharacterized
descriptions <- descriptions %>%
  mutate(
    final = ifelse(
      is.na(final),
      ifelse(
        str_detect(description, "uncharacterized", negate = TRUE),
        description,
        "uncharacterized"
      ),
      final
    )
  )

sum(is.na(descriptions$final))
#> [1] 0

## finalize descriptions ------------------------------------------------------
sum(descriptions$final == "uncharacterized")
# now we have a description for all genes, with 8265 genes labelled as
# "uncharacterized". Now we get rid of all columns except the final description
# and gene id

descriptions <- descriptions %>%
  select(gene_id, final) %>%
  rename(description = final)


###############################################################################
# final merged table ==========================================================
###############################################################################

merged <- left_join(symbols, descriptions, by = "gene_id")

merged <- merged %>%
  mutate(
    display_name = paste0(
      paste0("[", gene_id, "]"),
      ifelse(is.na(symbol), "", paste0(" ", symbol, " |")),
      ifelse(is.na(description), "", paste0(" ", description))
    )
  )


# save ========================================================================
write_tsv(merged, here(genomedir, "results/genes/gene-display-names.tsv"))

# end =========================================================================
