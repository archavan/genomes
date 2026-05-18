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
  c("KY21.Chr13.184", "KY21.Chr10.33"),
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
    .amplifier = "B1",
    .expression = tpm,
    .evalue = 20,
    .outdir = here(resdir, "cr5.29")
  )
)

## clusters 8 and 21 ===========================================================
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

## immune clusters ============================================================

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
