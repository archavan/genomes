# Arun Chavan
# Started: 2026-10-21

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(Biostrings)

basedir <- here("gene-trees/rtk/full-gene")

# use species and proteomes from orthofinder run
# in the proteome file, already:
#   1. longest protein chosen
#   2. using gene_id as header instead of protein_id
orthdir <- here("orthofinder/2026-08-07")
idmap_mmus <- read_tsv(here(orthdir, "results/id-name-mapping/mmus.tsv"))
proteome_mmus <- readAAStringSet(here(
  orthdir,
  "results/input-proteomes/mmus.fa"
))

###############################################################################
# define seeds and outgroups ==================================================
###############################################################################

# class info: https://en.wikipedia.org/wiki/Receptor_tyrosine_kinase
seeds <- tribble(
  ~symbol  , ~seed_type , ~rtk_class ,
  "Csf1r"  , "ingroup"  , "III"      ,
  "Pdgfra" , "ingroup"  , "III"      ,
  "Pdgfrb" , "ingroup"  , "III"      ,
  "Kit"    , "ingroup"  , "III"      ,
  "Kdr"    , "ingroup"  , "IV"       ,
  "Flt1"   , "ingroup"  , "IV"       ,
  "Flt3"   , "ingroup"  , "III"      ,
  "Flt4"   , "ingroup"  , "IV"       ,
  "Fgfr1"  , "outgroup" , "V"        ,
  "Fgfr2"  , "outgroup" , "V"        ,
  "Fgfr3"  , "outgroup" , "V"        ,
  "Fgfr4"  , "outgroup" , "V"        ,
  "Tek"    , "outgroup" , "XI"       ,
  "Tie1"   , "outgroup" , "XI"
) |>
  left_join(
    y = (idmap_mmus |>
      select(gene_id, symbol) |>
      distinct()),
    by = "symbol"
  ) |>
  mutate(spp = "mmus") |>
  mutate(label = paste0(spp, "__", symbol))

stopifnot(!any(duplicated(seeds$gene_id)), !any(is.na(seeds$gene_id)))

###############################################################################
# fetch sequences =============================================================
###############################################################################

seeds <- seeds |>
  mutate(sequence = map(gene_id, ~ proteome_mmus[.x])) |>
  mutate(
    sequence = map2(sequence, label, function(x, y) {
      names(x) <- y
      return(x)
    })
  )

###############################################################################
# write sequences and metadata ================================================
###############################################################################

fs::dir_create(here(basedir, "data/seeds"))
seeds_maxwidth <- max(map_int(seeds$sequence, nchar))
for (i in seq_len(nrow(seeds))) {
  writeXStringSet(
    seeds$sequence[[i]],
    filepath = here(basedir, "data/seeds", glue("{seeds$label[i]}.fasta")),
    width = seeds_maxwidth
  )
}
write_tsv(select(seeds, -sequence), here(basedir, "data/seeds/seeds.tsv"))

# end #########################################################################
