# Arun Chavan
# 2026-04-01

###############################################################################
# setup =======================================================================
###############################################################################

source("hcr-probes/code/utils_hcr-pipeline.R")

genomedir <- here("species/bflo/VNyyK")
bulkdir <- here(
  "~/Dropbox/research/projects/endoderm/species/bflo/bulk/ARCA118"
)

resdir <- here("hcr-probes/results/bflo")

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

## stromal cell clusters: gut =================================================

# cluster_v1: 5, 21, 26, 29
design_probes(
  .gene_id = "ncbi:118410441",
  .species = "bflo",
  .amplifier = "B1",
  .expression = tpm,
  .outdir = here(resdir, "likely-stromal-gut")
)

# cluster_v1: 5, 21, 29
design_probes(
  .gene_id = "ncbi:118414160",
  .species = "bflo",
  .amplifier = "B2",
  .expression = tpm,
  .outdir = here(resdir, "likely-stromal-gut")
)

## stromal cell clusters: mdiv ================================================

# cluster_v1: 17
design_probes(
  .gene_id = "ncbi:118423395",
  .species = "bflo",
  .amplifier = "B1",
  .expression = tpm,
  .outdir = here(resdir, "likely-stromal-mdiv")
)

# cluster_v1: 13
design_probes(
  .gene_id = "ncbi:118418336",
  .species = "bflo",
  .amplifier = "B2",
  .expression = tpm,
  .outdir = here(resdir, "likely-stromal-mdiv")
)

## coelomocyte ================================================================

# cluster_v1: 33
design_probes(
  .gene_id = "ncbi:118426472",
  .species = "bflo",
  .amplifier = "B2",
  .expression = tpm,
  .outdir = here(resdir, "likely-coelomocytes_bf33")
)

# cluster_v1: 33
design_probes(
  .gene_id = "ncbi:118414102",
  .species = "bflo",
  .amplifier = "B1",
  .expression = tpm,
  .outdir = here(resdir, "likely-coelomocytes_bf33")
)

## Spic =======================================================================

walk(
  c("ncbi:118414344", "ncbi:118412785"),
  ~ design_probes(
    .gene_id = .x,
    .species = "bflo",
    .amplifier = "B1",
    .expression = tpm,
    .outdir = here(resdir, "m3b2@Spi_1.b.c")
  )
)

# ncbi:118414344	m3b2_Spi(1,b,c)
# ncbi:118412785	m3b2_Spi(1,b,c)

## atrial =====================================================================

# cluster_v1: 3, 22
design_probes(
  .gene_id = "ncbi:118414386",
  .species = "bflo",
  .amplifier = "B1",
  .expression = tpm,
  .outdir = here(resdir, "atrial")
)

## Kupffer-like ===============================================================

# cluster_v1: 18
design_probes(
  .gene_id = "ncbi:118412380",
  .species = "bflo",
  .amplifier = "B1",
  .expression = tpm,
  .outdir = here(resdir, "kupffer-like")
)
