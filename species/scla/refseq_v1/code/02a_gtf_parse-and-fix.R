# Arun Chavan
# Started: 2025-08-07

# background ==================================================================
# Here we read in the GTF file, parse it (which happens automatically when we
# import using {rtracklayer}), fix several gene id issues, and filter out
# pseudogenes.

# setup =======================================================================
library(rtracklayer)
library(tidyverse)
library(here)

genomedir <- here("species/scla/refseq_v1")
mt_chr <- "NC_037072.1"

read_gxf_as_tibble <- function(.file) {
  rtracklayer::import(.file) |>
    GenomicRanges::as.data.frame() |>
    as_tibble()
}

# data ========================================================================
gtf <- read_gxf_as_tibble(
  here(genomedir, "data/refseq/GCF_013122585.1_ASM1312258v2_genomic.gtf.gz")
)

genes <- read_tsv(here(genomedir, "results/genes/ncbi-genes-table.tsv"))

###############################################################################
# fix gene ids ================================================================
###############################################################################

gtf <- gtf |>
  mutate(gene_id_new = NA_character_) |>
  relocate(gene_id_new, .after = gene_id)

# This is a typical problem with refseq GTF files, where the gene_id attribute
# doesn't actually have the gene ids, but often symbols or something else. The
# gene ids are actually in the db_xref attribute. This is in the format,
# datatase:numeric_id.

# some observations ===========================================================

# 1. Most rows have a GeneID in the db_xref attribute, but some have an ID from
#    a different db.

gtf |>
  separate(db_xref, into = c("db", "db_id"), sep = ":") |>
  count(db)
#> # A tibble: 4 × 2
#>   db                    n
#>   <chr>             <int>
#> 1 GeneID           722751
#> 2 RFAM                560
#> 3 UniProtKB/TrEMBL     39
#> 4 <NA>                 78

# 2. There are 78 rows that don't have a db_xref at all. These rows are all from
#    the mitochondrial chromosome.

gtf |>
  filter(is.na(db_xref)) |>
  count(seqnames)
#> # A tibble: 1 × 2
#>   seqnames        n
#>   <fct>       <int>
#> 1 NC_037072.1    78

# fix mitochondrial gene ids ==================================================
# The gene_id column for mitochondrial gene entries (not codon, CDS, etc) holds
# the locus tag and db_xref holds the gene id is GeneID:numeric_id format. They
# both have 1:1 mapping, so we can fill in GeneID:numeric_id for all rows
# associated with a given gene based on the relationship between locus tag and
# dbxref.

mt_idmap <- gtf |>
  filter(seqnames == mt_chr) |>
  filter(type == "gene") |>
  select(gene_id, db_xref) |>
  distinct()

get_gene_id_for_mt_locus_tag <- function(.locus_tag) {
  mt_idmap$db_xref[which(mt_idmap$gene_id == .locus_tag)]
}

for (i in mt_idmap$gene_id) {
  gtf$gene_id_new[which(gtf$gene_id == i)] <- get_gene_id_for_mt_locus_tag(i)
}


# add gene ids for genes that have gene ids in db_xref ========================

gtf <- gtf |>
  mutate(
    gene_id_new = map2_chr(
      gene_id_new,
      db_xref,
      function(.x, .y) {
        ifelse(is.na(.x) & grepl("^GeneID", .y), .y, .x)
      },
      .progress = TRUE
    )
  )

# non GeneID db_xref entries ==================================================

# At this point we have gene ids for all gene entries. There are 23364 genes in
# the ncbi genes table and also in the gtf file, all of which have a gene_id_new
# entry now. The only rows that are still NA for gene_id_new are non-gene
# entries (transcripts, exons, etc) that have a separate db_xref entry, e.g.
# from RFAM. For these rows we can simply inherit the gene_id_new entry from the
# gene row. The other option is that, for these rows the gene_id column actually
# holds the symbols starting with LOC1203XXXXXX where the number after LOC is
# the same as the numeric part of GeneID:numeric_id in db_xref, so we can just
# use the symbol and then replace LOC with GeneID:.

gtf |>
  filter(type == "gene", !is.na(gene_id_new)) |>
  select(gene_id, gene_id_new, db_xref) |>
  distinct()
#> # A tibble: 23,364 × 3
#>    gene_id      gene_id_new      db_xref
#>    <chr>        <chr>            <chr>
#>  1 LOC120330027 GeneID:120330027 GeneID:120330027
#>  2 LOC120326291 GeneID:120326291 GeneID:120326291
#>  3 LOC120343894 GeneID:120343894 GeneID:120343894
#>  4 LOC120329136 GeneID:120329136 GeneID:120329136
#>  5 LOC120328199 GeneID:120328199 GeneID:120328199
#>  6 LOC120325368 GeneID:120325368 GeneID:120325368
#>  7 LOC120327205 GeneID:120327205 GeneID:120327205
#>  8 LOC120330924 GeneID:120330924 GeneID:120330924
#>  9 LOC120331548 GeneID:120331548 GeneID:120331548
#> 10 LOC120334032 GeneID:120334032 GeneID:120334032
#> # ℹ 23,354 more rows

gtf |>
  filter(is.na(gene_id_new)) |>
  count(type)
#> # A tibble: 2 × 2
#>   type           n
#>   <fct>      <int>
#> 1 transcript   280
#> 2 exon         280

gtf |>
  filter(gene_id %in% gtf$gene_id[which(is.na(gene_id_new))]) |>
  select(seqnames, start, end, gene_id, gene_id_new, db_xref) |>
  filter(!is.na(gene_id_new)) |>
  mutate(
    gene_id = str_remove(gene_id, "^LOC"),
    gene_id_new = str_remove(gene_id_new, "^GeneID:")
  ) |>
  filter(gene_id != gene_id_new)
#> A tibble: 0 × 6

gtf <- gtf |>
  mutate(
    gene_id_new = map2_chr(
      gene_id_new,
      gene_id,
      function(.x, .y) {
        ifelse(is.na(.x) & grepl("^LOC", .y), .y, .x)
      }
    )
  )

# make gene_id prefixes consistent ============================================

# Now we have all gene_id_new entries populated. They have either
# GeneID:numeric_id or LOC1203XXXXXX format. We just replace the LOC or GeneID:
# with ncbi. now.

stopifnot(!any(is.na(gtf$gene_id_new)))

table(str_remove(gtf$gene_id_new, "\\d+"))
#>
#> GeneID:     LOC
#>  722868     560

gtf <- gtf |>
  mutate(
    gene_id_new = str_replace(gene_id_new, "^LOC", "ncbi."),
    gene_id_new = str_replace(gene_id_new, "^GeneID:", "ncbi.")
  )

# replace gene_id with gene_id_new ============================================
gtf <- gtf |>
  select(-gene_id) |>
  rename(gene_id = gene_id_new)

###############################################################################
# write files =================================================================
# #############################################################################

fs::dir_create(here(genomedir, "results/gtf"))

gtf |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
  rtracklayer::export(here(genomedir, "results/gtf/gtf_fixed.gtf"))

system(paste0("gzip ", here(genomedir, "results/gtf/gtf_fixed.gtf")))

# end =========================================================================
