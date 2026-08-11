# Arun Chavan
# 2026-04-01

###############################################################################
# setup =======================================================================
###############################################################################

source("hcr-probes/code/utils_hcr-pipeline.R")

genomedir <- here("species/crob/KY21")
bulkdir <- here(
  "~/Dropbox/research/projects/endoderm/species/crob/bulk/ARCA123/KY21"
)

resdir <- here("hcr-probes/results/crob")

###############################################################################
# data ========================================================================
###############################################################################

tpm <- read_csv(here(bulkdir, "data/03_compiled/tpm_by-transcript-id.csv")) |>
  column_to_rownames("transcript_id") |>
  rowMeans() |>
  enframe(name = "transcript_id", "tpm")

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

###############################################################################
# probes ======================================================================
###############################################################################

## stromal clusters ===========================================================

# cluster_v1: 6, 19
design_probes(
  .gene_id = "KY21.Chr7.743",
  .species = "crob",
  .amplifier = "B1",
  .expression = tpm,
  .evalue = 20,
  .outdir = here(resdir, "likely-stromal_cr6.19")
)

# cluster_v1: 35
design_probes(
  .gene_id = "KY21.Chr1.567",
  .species = "crob",
  .amplifier = "B1",
  .expression = tpm,
  .evalue = 20,
  .outdir = here(resdir, "likely-stromal_cr35")
)

# cluster_v1: 28, 36
walk(
  c("KY21.Chr12.747"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B1",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "likely-stromal_cr28.36")
  )
)

## likely enteric phagocytes ==================================================

# cluster_v1: 2, 4, 31
walk(
  c("KY21.Chr14.387", "KY21.Chr12.262", "KY21.Chr10.498", "KY21.Chr6.173"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B1",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "likely-enteric-phagocyte_cr2.4.31")
  )
)

## clusters 10 and 17 =========================================================
walk(
  c("KY21.Chr13.184", "KY21.Chr10.33", "KY21.Chr8.932"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B1",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "cr10.17")
  )
)

## clusters 5 and 29 ==========================================================
walk(
  c("KY21.Chr11.710", "KY21.Chr1.2158", "KY21.Chr1.2159", "KY21.Chr1.270"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B3",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "cr5.29")
  )
)

## clusters 8 and 21 (hepatocyte-like) ========================================
walk(
  c("KY21.Chr8.628"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B1",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "cr8.21")
  )
)

## clusters 12 and 13 (absorptive) ============================================

walk(
  c(
    "KY21.Chr2.1138",
    "KY21.Chr8.1221"
  ),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B2",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "cr12.13")
  )
)

walk(
  c(
    "KY21.Chr3.750"
  ),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B4",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "cr12.13")
  )
)

## cluster 23 =================================================================
walk(
  c("KY21.Chr1.1392", "KY21.Chr7.852", "KY21.Chr7.426"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B3",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "cr23")
  )
)

## cluster 30 =================================================================
walk(
  c("KY21.Chr8.1350", "KY21.Chr13.294"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B3",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "cr30")
  )
)

## immune clusters ============================================================

walk(
  c("KY21.Chr4.869", "KY21.Chr2.1024"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B4",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "immune_cr11.25")
  )
)

walk(
  c("KY21.Chr1.1090"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B2",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "immune_cr7")
  )
)

walk(
  c("KY21.Chr1.1264"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B3",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "immune_cr18")
  )
)

walk(
  c("KY21.Chr11.714", "KY21.Chr7.236"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B4",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "immune_cr22")
  )
)

walk(
  c("KY21.Chr7.840", "KY21.Chr11.668"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B1",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "immune_cr26")
  )
)

walk(
  c("KY21.Chr13.78"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B1",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "immune_cr34")
  )
)

walk(
  c("KY21.Chr12.538", "KY21.Chr11.222"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B2",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "immune_cr11.15.20")
  )
)

walk(
  c("KY21.Chr6.542"),
  ~ design_probes(
    .gene_id = .x,
    .species = "crob",
    .amplifier = "B1",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "immune_m3c1@Spi_1.b.c")
  )
)

###############################################################################
# combined order sheets =======================================================
###############################################################################

combine_opools(
  .dirs = paste0(
    here(resdir),
    "/",
    c(
      "likely-enteric-phagocyte_cr2.4.31/KY21.Chr12.262",
      "cr8.21/KY21.Chr8.628",
      "cr12.13/KY21.Chr2.1138",
      "cr5.29/KY21.Chr1.2158",
      "cr23/KY21.Chr7.852",
      "immune_cr11.25/KY21.Chr2.1024",
      "cr10.17/KY21.Chr8.932",
      "immune_cr7/KY21.Chr1.1090",
      "immune_cr18/KY21.Chr1.1264",
      "immune_cr22/KY21.Chr11.714",
      "cr30/KY21.Chr8.1350"
    )
  ),
  .out_path = here("hcr-probes/docs/orders/2026-06-24/combined-opools.xlsx")
)


combine_opools(
  .dirs = paste0(
    here(resdir),
    "/",
    c(
      "cr12.13/KY21.Chr3.750",
      "cr12.13/KY21.Chr8.1221"
    )
  ),
  .out_path = here(
    "hcr-probes/docs/orders/2026-08-11/combined-opools-crob.xlsx"
  )
)
